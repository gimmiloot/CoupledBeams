# Corpus map

The four sources below were used with equal authority. Equality means that no paper is treated as the preferred authorial model. Evidence is nevertheless limited by what each file actually contains: a source cannot confirm Russian body-prose syntax when its body is in English.

The original files are provenance only. Future use of the skill does not require access to them.

## A1. Published coupled-rods article

- File at analysis time: `Статья-Дорофеев-2025.pdf`.
- Title: «Собственные колебания сопряжённых стержней».
- Russian material used: abstract; introduction; problem statement; analytical solution; angle study; small-angle analysis; conclusion.
- Elements analysed: narrowing of the introduction, order of the eigenvalue formulation, transitions between equations, discussion of frequency curves and modes, limited conclusions.
- Excluded: author and publisher data, bibliography formatting, journal boilerplate, mathematical glyphs damaged during extraction, typographical errors.
- Technical issue: the PDF was produced with legacy font encoding. Extracted Cyrillic is usable after checking, but Greek letters, punctuation, and formulas are not reliable enough to reproduce from the text layer.

## A2. Geometric-parameters manuscript

- File at analysis time: `Дорофеев_текст.docx`.
- Visible Russian title: «Частоты собственных колебаний сопряжённых стержней. Влияние геометрических параметров».
- Corpus-mapping note: the supplied prompt associated this file with the length-ratio title. Inspection of the document shows that the two titles were interchanged in the prompt; the file is identified here by its actual contents.
- Russian material used: title, abstract, keywords, and only ordinary Russian prose in the front matter.
- Other material used: full-paper composition and paragraph functions in the English body, without treating English constructions as Russian syntax.
- Elements analysed: abstract order, introduction-to-purpose transition, problem-statement sequence, results organisation, conclusion sequence.
- Excluded: English wording as a Russian phrase model, journal placeholders, author data, funding, bibliography, dates and DOI/EDN blanks, internal publisher notes.
- Technical issue: OOXML extraction is clean for prose, but the document contains about 160 OfficeMath objects; extracted equations must not be used as authoritative mathematical text.

## A3. Length-ratio article

- File at analysis time: `Dorofeev_beams_2.pdf`.
- Visible Russian title: «Частоты собственных колебаний стержневой конструкции. Влияние отношения длин стержней».
- Corpus-mapping note: this is the length-ratio paper even though the supplied prompt assigned the geometric-parameters title to this filename.
- Russian material used: pages 1--10, including abstract, introduction, problem statement, analytical solution, results, mode discussion, and conclusion.
- Elements analysed: parameter introduction, formula-to-plot-to-interpretation sequence, comparison with FEM, mode-shape discussion, cautious future work.
- Excluded: English duplicate on pages 11--12, bibliography formatting, captions as body-prose templates, numerical values and notation as universal style.
- Technical issue: the text layer is readable, but discretionary hyphens, page breaks, and several draft punctuation or grammar errors were filtered out.

## A4. Composite-beam-theory article

- File at analysis time: `Dorofeev.pdf`.
- Title: «Спектры композитной стержневой конструкции. Влияние выбора теории стержня на собственные частоты».
- Russian material used: pages 1--9, including abstract, introduction, model hierarchy, problem statement, comparison criterion, results, and conclusion.
- Elements analysed: purpose-led literature synthesis, sequential definition of models, criterion-led comparisons, separation of observed discrepancy from its supported mechanism.
- Excluded: English duplicate on pages 10--11, author/publisher matter, bibliography formatting, table-layout defects, model-specific symbols as reusable style.
- Technical issue: extracted overbars, indices, derivatives, and large brackets can collapse. Formulas were treated only as structural landmarks and never copied from extraction.

## Corpus-wide exclusions

- OCR, extraction, punctuation, grammar, and encoding errors;
- publisher instructions, placeholders, headers, author biographies, funding, and submission dates;
- journal-specific layout, bibliography, transliteration, and English-translation syntax;
- unverified terminology differences between drafts;
- scientific results, formulas, parameter values, and notation as stylistic templates.
