from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPT = SCRIPT_DIR / "style_lint.py"
sys.path.insert(0, str(SCRIPT_DIR))
import style_lint as lint  # noqa: E402


def warnings_for(text: str) -> list[lint.WarningItem]:
    return lint.collect_warnings(text)


def codes_for(text: str) -> list[str]:
    return [warning.code for warning in warnings_for(text)]


def run_cli(*args: str, stdin: str | bytes = b"") -> subprocess.CompletedProcess[bytes]:
    data = stdin.encode("utf-8") if isinstance(stdin, str) else stdin
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        input=data,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        timeout=10,
    )


class StyleLintTests(unittest.TestCase):
    def test_long_russian_sentence_warns_len001(self) -> None:
        text = " ".join(["Исследование"] * 36) + "."
        self.assertEqual(codes_for(text).count("LEN001"), 1)

    def test_long_english_sentence_has_no_len001(self) -> None:
        text = " ".join(["analysis"] * 40) + "."
        self.assertNotIn("LEN001", codes_for(text))

    def test_mixed_tex_skips_english_abstract(self) -> None:
        english = " ".join(["analysis"] * 40) + "."
        russian = "Рассматривается " + " ".join(["зависимость"] * 36) + "."
        text = "\\begin{abstract}\n" + english + "\n\\end{abstract}\n\n" + russian
        warnings = warnings_for(text)
        length_warnings = [warning for warning in warnings if warning.code == "LEN001"]
        self.assertEqual(len(length_warnings), 1)
        self.assertEqual(length_warnings[0].start, text.index("Рассматривается"))
        self.assertNotIn("NOM001", [warning.code for warning in warnings])

    def test_russian_with_latin_terms_is_still_checked(self) -> None:
        text = (
            "В расчёте MAC и FEM по модели Euler--Bernoulli "
            + " ".join(["последовательно"] * 32)
            + " сопоставляются полученные русские результаты."
        )
        self.assertTrue(lint.is_russian_text(text))
        self.assertIn("LEN001", codes_for(text))

    def test_all_requested_environments_are_masked(self) -> None:
        body = (
            " ".join(["исследование"] * 40)
            + " 1 2 3 4 5 TODO спектральные показатели"
        )
        for environment in lint.MASKED_LATEX_ENVIRONMENTS:
            with self.subTest(environment=environment):
                if environment == "align*":
                    opener = "\\begin\r\n  { align* }"
                    closer = "\\end\r\n  { align* }"
                elif environment == "minted":
                    opener = "\\begin{minted}[linenos]{python}"
                    closer = "\\end{minted}"
                else:
                    opener = f"\\begin{{{environment}}}"
                    closer = f"\\end{{{environment}}}"
                text = opener + "\r\n" + body + "\r\n" + closer + "\r\nКраткий русский текст."
                masked = lint.mask_latex(text)
                self.assertEqual(len(masked), len(text))
                self.assertEqual(
                    [index for index, char in enumerate(masked) if char in "\r\n"],
                    [index for index, char in enumerate(text) if char in "\r\n"],
                )
                self.assertEqual(codes_for(text), [])

    def test_commented_environment_opener_does_not_hide_prose(self) -> None:
        text = (
            "% \\begin{equation}\n"
            "Следует особо подчеркнуть полученный результат.\n"
            "\\end{equation}"
        )
        self.assertIn("LEX001", codes_for(text))

    def test_comment_escape_uses_backslash_parity(self) -> None:
        even = (
            "\\\\% \\begin{equation}\n"
            "Следует особо подчеркнуть полученный результат.\n"
            "\\end{equation}"
        )
        odd = (
            "\\% \\begin{equation}\n"
            "Следует особо подчеркнуть полученный результат.\n"
            "\\end{equation}"
        )
        self.assertIn("LEX001", codes_for(even))
        self.assertNotIn("LEX001", codes_for(odd))

    def test_equation_does_not_trigger_len_or_num(self) -> None:
        body = " ".join(["исследование"] * 40) + " 1 2 3 4 5 6"
        text = (
            "На рис. 1 показана зависимость.\n"
            "\\begin { equation }\n"
            + body
            + "\n\\end { equation }"
        )
        codes = codes_for(text)
        self.assertNotIn("LEN001", codes)
        self.assertNotIn("NUM001", codes)

    def test_align_star_does_not_trigger_len_or_num(self) -> None:
        body = " ".join(["исследование"] * 40) + " 1 2 3 4 5 6"
        text = (
            "На рис. 1 показана зависимость.\n"
            "\\begin{align*}\n"
            + body
            + "\n\\end{align*}"
        )
        codes = codes_for(text)
        self.assertNotIn("LEN001", codes)
        self.assertNotIn("NUM001", codes)

    def test_masking_preserves_following_warning_position(self) -> None:
        text = (
            "\\begin{equation}\n"
            "x_1=1+2+3\n"
            "x_2=4+5+6\n"
            "\\end{equation}\n\n"
            "    Следует особо подчеркнуть полученный результат."
        )
        masked = lint.mask_latex(text)
        self.assertEqual(len(masked), len(text))
        self.assertEqual(
            [index for index, char in enumerate(masked) if char == "\n"],
            [index for index, char in enumerate(text) if char == "\n"],
        )
        warning = next(item for item in warnings_for(text) if item.code == "LEX001")
        self.assertEqual(warning.start, text.index("Следует"))
        line_starts = [0]
        line_starts.extend(index + 1 for index, char in enumerate(text) if char == "\n")
        self.assertEqual(lint.line_column(line_starts, warning.start), (6, 5))

    def test_figure_value_catalogue_warns_num001(self) -> None:
        text = "На рис. 2 показаны значения 1,0, 2,0, 3,0 и 4,0 для выбранных случаев."
        self.assertEqual(codes_for(text).count("NUM001"), 1)

    def test_masked_formula_numbers_do_not_count_for_num001(self) -> None:
        text = (
            "На рис. 2 показана рассчитанная зависимость.\n"
            "\\begin{equation}\n1+2+3+4+5+6=21\n\\end{equation}"
        )
        self.assertNotIn("NUM001", codes_for(text))

    def test_figure_caption_remains_visible(self) -> None:
        text = (
            "\\begin{figure}\n"
            "\\caption{На рис. 4 показаны значения 1, 2, 3 и 4 для выбранных случаев.}\n"
            "\\end{figure}"
        )
        self.assertIn("NUM001", codes_for(text))

    def test_term001_message_is_neutral(self) -> None:
        warning = next(
            item
            for item in warnings_for("Введены спектральные показатели для сравнения.")
            if item.code == "TERM001"
        )
        self.assertIn(
            "проверьте формулировку по текущей утверждённой рукописи и документации проекта",
            warning.message,
        )
        self.assertNotIn("собственные частоты", warning.message)
        self.assertNotIn("отношение длин стержней", warning.message)

    def test_exact_and_wildcard_allowlist(self) -> None:
        text = "Введены спектральные показатели для сравнения."
        warning = next(item for item in warnings_for(text) if item.code == "TERM001")
        with tempfile.TemporaryDirectory() as directory:
            allowlist = Path(directory) / "allowlist.txt"
            for code in ("TERM001", "*"):
                with self.subTest(code=code):
                    allowlist.write_bytes(f"{code}\t{warning.subject}\n".encode("utf-8"))
                    result = run_cli("-", "--allowlist", str(allowlist), stdin=text)
                    self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8"))
                    self.assertEqual(result.stdout, b"")

    def test_cli_exit_codes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            malformed = Path(directory) / "malformed.txt"
            malformed.write_bytes("TERM001 без табуляции\n".encode("utf-8"))
            cases = (
                ("clean", run_cli("-", stdin="Краткий русский текст."), 0),
                ("warning", run_cli("-", stdin="TODO"), 1),
                (
                    "malformed_allowlist",
                    run_cli("-", "--allowlist", str(malformed), stdin="Краткий русский текст."),
                    2,
                ),
                ("invalid_utf8", run_cli("-", stdin=b"\xff"), 2),
            )
            for name, result, expected in cases:
                with self.subTest(name=name):
                    self.assertEqual(result.returncode, expected)

    def test_meta_marker_warns(self) -> None:
        self.assertEqual(codes_for("Здесь оставлена метка TODO.").count("META001"), 1)

    def test_language_gate_skips_english_dominated_paragraph(self) -> None:
        text = (
            " ".join(["analysis"] * 30)
            + " следует отметить затем следует отметить изменение с 1 до 2, а значение с 3 до 4 "
            + "на рис. 5 приведены 1, 2, 3 и 4."
        )
        codes = codes_for(text)
        for code in ("REP001", "RNG001", "NUM001"):
            self.assertNotIn(code, codes)

    def test_russian_paragraph_heuristics_still_warn(self) -> None:
        text = (
            "Следует отметить, что частота изменилась с 1 до 2, а амплитуда с 3 до 4. "
            "Следует отметить, что сравнение относится к рассмотренной сетке."
        )
        codes = codes_for(text)
        self.assertIn("REP001", codes)
        self.assertIn("RNG001", codes)


if __name__ == "__main__":
    unittest.main()
