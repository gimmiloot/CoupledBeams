#!/usr/bin/env python3
"""Warn about selected Russian scientific-prose style risks.

The linter is deliberately heuristic. It never rewrites input and does not
assess scientific correctness. Exit codes: 0 = no warnings, 1 = warnings,
2 = usage, input, encoding, or allowlist error.
"""

from __future__ import annotations

import argparse
import bisect
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import re
import sys
from typing import Iterable


MAX_SENTENCE_WORDS = 35
MAX_FIGURE_PARAGRAPH_NUMBERS = 4

INFLATED_PHRASES = (
    "в рамках настоящего комплексного исследования",
    "следует особо подчеркнуть",
    "представляется целесообразным отметить",
    "всесторонне исследовано",
    "уникальный подход",
    "инновационный метод",
    "значительный научный вклад",
    "открывает широкие перспективы",
    "обеспечивает глубокое понимание",
    "результаты однозначно доказывают",
)

COINED_PHRASES = (
    "комплексная спектральная картина",
    "многофакторная параметрическая обусловленность",
    "спектральная реорганизация",
    "модально-спектральная перестройка",
    "ориентационно-индуцированная связанность",
    "геометрически обусловленная модальная трансформация",
    "параметрически индуцированный эффект",
    "мера продольной геометрической асимметрии",
)

SUSPECT_TERM_VARIANTS = (
    "спектральные показатели",
    "кластеризация ветвей",
    "угол ориентации материала",
    "поворот материальных осей",
    "ориентационный параметр",
    "направление анизотропии",
    "параметр асимметрии длин",
)

INTRO_CLICHES = (
    "таким образом",
    "следует отметить",
    "необходимо отметить",
    "важно отметить",
    "следует подчеркнуть",
    "в данной работе",
    "полученные результаты показывают",
)

NOMINAL_HEADS = {
    "анализ",
    "исследование",
    "определение",
    "изменение",
    "описание",
    "сравнение",
    "влияние",
    "формирование",
    "построение",
    "оценка",
    "характер",
    "перестройка",
}

WORD_RE = re.compile(r"[^\W\d_]+(?:-[^\W\d_]+)*", re.UNICODE)
NUMBER_RE = re.compile(r"(?<![\w])[-+]?\d+(?:[.,]\d+)?\s*%?", re.UNICODE)
FIGURE_WORD_RE = re.compile(
    r"\b(?:рис(?:\.|унок|унке|унка)?|график(?:е|а)?|табл(?:\.|ица|ице|ицы)?)\b",
    re.IGNORECASE,
)
RANGE_RE = re.compile(
    r"\bс\b[^.!?;\n]{1,60}\bдо\b[^.!?;\n]{0,50},\s*а\b"
    r"[^.!?;\n]{0,60}\bс\b[^.!?;\n]{1,60}\bдо\b",
    re.IGNORECASE,
)
FINITE_VERB_RE = re.compile(
    r"(?:ется|ются|ится|атся|ятся|ет|ют|ит|ат|ят|ем|им|ли|лся|лась|лось|лись)$",
    re.IGNORECASE,
)
NOUNISH_RE = re.compile(
    r"(?:ние|ния|нию|нием|ция|ции|цию|цией|ость|ости|остью|ство|ства|ству|ством|"
    r"ого|его|ой|ей|ых|их|а|я|ы|и|ов|ев)$",
    re.IGNORECASE,
)
RUSSIAN_LETTER_RE = re.compile(r"[А-Яа-яЁё]")

MASKED_LATEX_ENVIRONMENTS = (
    "equation",
    "equation*",
    "align",
    "align*",
    "alignat",
    "alignat*",
    "gather",
    "gather*",
    "multline",
    "multline*",
    "flalign",
    "flalign*",
    "displaymath",
    "math",
    "split",
    "verbatim",
    "verbatim*",
    "lstlisting",
    "minted",
)
_LATEX_ENVIRONMENT_ALTERNATION = "|".join(
    re.escape(name) for name in sorted(MASKED_LATEX_ENVIRONMENTS, key=len, reverse=True)
)
LATEX_ENVIRONMENT_RE = re.compile(
    r"\\begin\s*\{\s*(?P<environment>"
    + _LATEX_ENVIRONMENT_ALTERNATION
    + r")\s*\}.*?\\end\s*\{\s*(?P=environment)\s*\}",
    re.DOTALL,
)


@dataclass(frozen=True)
class WarningItem:
    code: str
    start: int
    end: int
    message: str
    subject: str


def normalise(text: str) -> str:
    return " ".join(text.casefold().split())


def is_russian_text(text: str) -> bool:
    """Return whether a span is predominantly Russian prose."""

    words = WORD_RE.findall(text)
    cyrillic_words = sum(bool(RUSSIAN_LETTER_RE.search(word)) for word in words)
    return cyrillic_words >= 3 and cyrillic_words * 2 >= len(words)


def mask_span(chars: list[str], start: int, end: int) -> None:
    for index in range(start, end):
        if chars[index] not in "\r\n":
            chars[index] = " "


def mask_pattern(chars: list[str], pattern: re.Pattern[str]) -> None:
    for match in pattern.finditer("".join(chars)):
        mask_span(chars, match.start(), match.end())


def mask_latex_comments(chars: list[str]) -> None:
    """Mask TeX comments, respecting the parity of preceding backslashes."""

    source = "".join(chars)
    for index, char in enumerate(source):
        if char != "%":
            continue
        backslashes = 0
        cursor = index - 1
        while cursor >= 0 and source[cursor] == chr(92):
            backslashes += 1
            cursor -= 1
        if backslashes % 2:
            continue
        end = index
        while end < len(source) and source[end] not in "\r\n":
            end += 1
        mask_span(chars, index, end)


def mask_latex(text: str) -> str:
    """Mask comments, selected environments, maths, and commands in place."""

    chars = list(text)
    mask_latex_comments(chars)
    mask_pattern(chars, LATEX_ENVIRONMENT_RE)
    patterns = (
        re.compile(r"\$\$.*?\$\$", re.DOTALL),
        re.compile(r"(?<!\\)\$(?:\\.|[^$])*?(?<!\\)\$", re.DOTALL),
        re.compile(r"\\\[.*?\\\]", re.DOTALL),
        re.compile(r"\\\(.*?\\\)", re.DOTALL),
        re.compile(r"\\(?:cite|citep|citet|ref|eqref|label|url|href)\*?(?:\[[^]]*\])?\{[^}]*\}"),
    )
    for pattern in patterns:
        mask_pattern(chars, pattern)
    mask_pattern(chars, re.compile(r"\\[A-Za-zА-Яа-яЁё@]+\*?"))
    return "".join(chars)


def iter_paragraphs(text: str) -> Iterable[tuple[int, int, str]]:
    for match in re.finditer(r"(?s)(?:^|(?<=\n))\s*\S.*?(?=\n\s*\n|\Z)", text):
        start, end = match.span()
        while start < end and text[start].isspace():
            start += 1
        while end > start and text[end - 1].isspace():
            end -= 1
        if start < end:
            yield start, end, text[start:end]


def protect_sentence_periods(text: str) -> str:
    chars = list(text)
    for match in re.finditer(r"(?<=\d)\.(?=\d)", text):
        chars[match.start()] = "¤"
    abbreviations = (
        r"\bрис\.",
        r"\bтабл\.",
        r"\bстр\.",
        r"\bгл\.",
        r"\bт\.\s*е\.",
        r"\bт\.\s*д\.",
        r"\bт\.\s*п\.",
        r"\bим\.",
    )
    for pattern in abbreviations:
        for match in re.finditer(pattern, text, re.IGNORECASE):
            for index in range(match.start(), match.end()):
                if chars[index] == ".":
                    chars[index] = "¤"
    for match in re.finditer(r"\b[А-ЯЁA-Z]\.\s*[А-ЯЁA-Z]\.", text):
        for index in range(match.start(), match.end()):
            if chars[index] == ".":
                chars[index] = "¤"
    return "".join(chars)


def iter_sentences(text: str) -> Iterable[tuple[int, int, str]]:
    protected = protect_sentence_periods(text)
    boundaries = [0]
    boundaries.extend(
        match.end()
        for match in re.finditer(r"[.!?](?=\s+(?:[А-ЯЁA-Z]|\)))", protected)
    )
    boundaries.append(len(text))
    for left, right in zip(boundaries, boundaries[1:]):
        while left < right and text[left].isspace():
            left += 1
        while right > left and text[right - 1].isspace():
            right -= 1
        if left < right and WORD_RE.search(text[left:right]):
            yield left, right, text[left:right]


def phrase_warnings(text: str, phrases: Iterable[str], code: str, message: str) -> list[WarningItem]:
    warnings: list[WarningItem] = []
    for phrase in phrases:
        for match in re.finditer(re.escape(phrase), text, re.IGNORECASE):
            warnings.append(
                WarningItem(code, match.start(), match.end(), message, normalise(match.group()))
            )
    return warnings


def nominalisation_warnings(text: str, base_offset: int = 0) -> list[WarningItem]:
    warnings: list[WarningItem] = []
    tokens = list(WORD_RE.finditer(text))
    used_until = -1
    for index, token in enumerate(tokens):
        if token.start() < used_until or token.group().casefold() not in NOMINAL_HEADS:
            continue
        window = tokens[index : index + 8]
        if len(window) < 5:
            continue
        selected = [window[0]]
        for candidate in window[1:]:
            gap = text[selected[-1].end() : candidate.start()]
            if re.search(r"[.!?;:\n]", gap) or FINITE_VERB_RE.search(candidate.group()):
                break
            selected.append(candidate)
        nounish = sum(
            item.group().casefold() in NOMINAL_HEADS or bool(NOUNISH_RE.search(item.group()))
            for item in selected
        )
        if len(selected) >= 5 and nounish >= 4:
            local_start, local_end = selected[0].start(), selected[-1].end()
            subject = normalise(text[local_start:local_end])
            warnings.append(
                WarningItem(
                    "NOM001",
                    base_offset + local_start,
                    base_offset + local_end,
                    "возможная цепочка отглагольных или отвлечённых существительных",
                    subject,
                )
            )
            used_until = local_end
    return warnings


def repeated_cliche_warnings(text: str, paragraphs: list[tuple[int, int, str]]) -> list[WarningItem]:
    occurrences: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for phrase in INTRO_CLICHES:
        for match in re.finditer(re.escape(phrase), text, re.IGNORECASE):
            paragraph_index = next(
                (i for i, (start, end, _) in enumerate(paragraphs) if start <= match.start() < end),
                None,
            )
            if paragraph_index is None:
                continue
            occurrences[phrase].append((match.start(), match.end(), paragraph_index))

    warnings: list[WarningItem] = []
    for phrase, items in occurrences.items():
        paragraph_seen: dict[int, int] = defaultdict(int)
        for document_index, (start, end, paragraph_index) in enumerate(items):
            paragraph_seen[paragraph_index] += 1
            repeated_in_paragraph = paragraph_index >= 0 and paragraph_seen[paragraph_index] >= 2
            repeated_in_document = document_index >= 2
            if repeated_in_paragraph or repeated_in_document:
                warnings.append(
                    WarningItem(
                        "REP001",
                        start,
                        end,
                        "повторяющееся вводное клише; проверьте необходимость повтора",
                        normalise(phrase),
                    )
                )
    return warnings


def service_marker_warnings(text: str) -> list[WarningItem]:
    patterns = (
        re.compile(r"\b(?:TODO|FIXME|TBD|NEEDS_CHECK|proved_here|needs_caution|not_reached_yet|branch-ready)\b", re.IGNORECASE),
        re.compile(r"\[(?:вставить|уточнить|проверить)[^\]]*\]", re.IGNORECASE),
        re.compile(r"\bздесь\s+будет\b", re.IGNORECASE),
        re.compile(r"\?{3,}"),
        re.compile(r"\{\{[^{}]+\}\}"),
        re.compile(r"<<[^<>]+>>"),
    )
    warnings: list[WarningItem] = []
    for pattern in patterns:
        for match in pattern.finditer(text):
            warnings.append(
                WarningItem(
                    "META001",
                    match.start(),
                    match.end(),
                    "возможная служебная метка в научном тексте",
                    normalise(match.group()),
                )
            )
    return warnings


def collect_warnings(raw_text: str) -> list[WarningItem]:
    text = mask_latex(raw_text)
    paragraphs = list(iter_paragraphs(text))
    russian_paragraphs = [
        (start, end, paragraph)
        for start, end, paragraph in paragraphs
        if is_russian_text(paragraph)
    ]
    warnings: list[WarningItem] = []

    for paragraph_start, _, paragraph in paragraphs:
        for local_start, local_end, sentence in iter_sentences(paragraph):
            if not is_russian_text(sentence):
                continue
            start = paragraph_start + local_start
            end = paragraph_start + local_end
            word_count = len(WORD_RE.findall(sentence))
            if word_count > MAX_SENTENCE_WORDS:
                warnings.append(
                    WarningItem(
                        "LEN001",
                        start,
                        end,
                        f"предложение содержит {word_count} слов; проверьте возможность разделения",
                        normalise(raw_text[start:end]),
                    )
                )
            warnings.extend(nominalisation_warnings(sentence, start))

    warnings.extend(
        phrase_warnings(text, INFLATED_PHRASES, "LEX001", "нежелательное оценочное или усложняющее выражение")
    )
    warnings.extend(
        phrase_warnings(text, COINED_PHRASES, "LEX002", "возможный придуманный или неоправданно усложнённый термин")
    )
    warnings.extend(repeated_cliche_warnings(text, russian_paragraphs))

    for paragraph_start, _, paragraph in russian_paragraphs:
        for match in RANGE_RE.finditer(paragraph):
            start = paragraph_start + match.start()
            end = paragraph_start + match.end()
            warnings.append(
                WarningItem(
                    "RNG001",
                    start,
                    end,
                    "двойной пересказ диапазонов; сформулируйте основную зависимость",
                    normalise(raw_text[start:end]),
                )
            )

    for start, end, paragraph in russian_paragraphs:
        numbers = list(NUMBER_RE.finditer(paragraph))
        if FIGURE_WORD_RE.search(paragraph) and len(numbers) >= MAX_FIGURE_PARAGRAPH_NUMBERS:
            warnings.append(
                WarningItem(
                    "NUM001",
                    start,
                    end,
                    f"абзац о рисунке или таблице содержит {len(numbers)} числовых значений; проверьте, не заменяют ли они анализ",
                    normalise(raw_text[start:end]),
                )
            )

    warnings.extend(service_marker_warnings(text))

    coined_ranges = [
        (warning.start, warning.end) for warning in warnings if warning.code == "LEX002"
    ]
    for variant in SUSPECT_TERM_VARIANTS:
        for match in re.finditer(re.escape(variant), text, re.IGNORECASE):
            if any(left <= match.start() and match.end() <= right for left, right in coined_ranges):
                continue
            warnings.append(
                WarningItem(
                    "TERM001",
                    match.start(),
                    match.end(),
                    "возможный новый вариант термина; проверьте формулировку по текущей утверждённой рукописи и документации проекта",
                    normalise(match.group()),
                )
            )

    unique: dict[tuple[str, int, int], WarningItem] = {}
    for warning in warnings:
        unique[(warning.code, warning.start, warning.end)] = warning
    return sorted(unique.values(), key=lambda item: (item.start, item.code, item.end))


def read_utf8(path: str) -> tuple[str, str]:
    if path == "-":
        data = sys.stdin.buffer.read()
        source = "<stdin>"
    else:
        data = Path(path).read_bytes()
        source = path
    return data.decode("utf-8-sig"), source


def read_allowlist(paths: Iterable[str]) -> set[tuple[str, str]]:
    allowed: set[tuple[str, str]] = set()
    for path_string in paths:
        path = Path(path_string)
        text = path.read_bytes().decode("utf-8-sig")
        for line_number, raw_line in enumerate(text.splitlines(), 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if "\t" not in raw_line:
                raise ValueError(f"{path}:{line_number}: expected CODE<TAB>exact subject")
            code, subject = raw_line.split("\t", 1)
            code = code.strip().upper()
            subject = normalise(subject)
            if not re.fullmatch(r"(?:[A-Z]{3,5}\d{3}|\*)", code) or not subject:
                raise ValueError(f"{path}:{line_number}: invalid allowlist entry")
            allowed.add((code, subject))
    return allowed


def line_column(line_starts: list[int], offset: int) -> tuple[int, int]:
    line_index = bisect.bisect_right(line_starts, offset) - 1
    return line_index + 1, offset - line_starts[line_index] + 1


def short_fragment(text: str, start: int, end: int, limit: int = 120) -> str:
    fragment = " ".join(text[start:end].split())
    if len(fragment) > limit:
        fragment = fragment[: limit - 1].rstrip() + "…"
    return fragment


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Warn about selected Russian scientific-writing style risks without changing the input."
    )
    parser.add_argument("path", nargs="?", default="-", help="UTF-8 file, or -/omitted for stdin")
    parser.add_argument(
        "--allowlist",
        action="append",
        default=[],
        metavar="FILE",
        help="UTF-8 file with CODE<TAB>exact-subject entries; may be repeated",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    args = build_parser().parse_args(argv)
    try:
        text, source = read_utf8(args.path)
        allowed = read_allowlist(args.allowlist)
    except (OSError, UnicodeDecodeError, ValueError) as error:
        print(f"style_lint: {error}", file=sys.stderr)
        return 2

    warnings = [
        warning
        for warning in collect_warnings(text)
        if (warning.code, warning.subject) not in allowed
        and ("*", warning.subject) not in allowed
    ]
    line_starts = [0]
    line_starts.extend(match.end() for match in re.finditer("\n", text))
    for warning in warnings:
        line, column = line_column(line_starts, warning.start)
        fragment = short_fragment(text, warning.start, warning.end)
        print(f"{source}:{line}:{column}: {warning.code} {warning.message}: «{fragment}»")
    return 1 if warnings else 0


if __name__ == "__main__":
    raise SystemExit(main())
