# Reddy 2004: симметрично слоистые балки

- Источник: J. N. Reddy, *Mechanics of Laminated Composite Plates and Shells: Theory and Analysis*, 2nd ed., 2004.
- PDF: `docs/literature/pdf/EB__Mechanics_of_Laminated_Composite_Plates_and_Shells_-JN_Reddy.pdf`.
- SHA-256: `A8D4D0FA67C7073D6EB48903B868BFA157875DC60F54484376A5B510547B37EA`.
- Число PDF-страниц: 855.
- Соответствие страниц: `pdf_index_zero_based = printed_page + 22`.
- Роль в проекте: основной источник для lamina/laminate stiffness, FSDT-редукции одного симметрично слоистого стержня и Table 4.3.3.

## Используемые положения

Для плоского напряжённого состояния книга задаёт reduced lamina stiffness \(Q\) и transformed stiffness \(\bar Q\). Мембранные компоненты упорядочены как \((xx,yy,xy)\). Transverse-shear block упорядочен как \((yz,xz)\), поэтому в балочной редукции \(xz\)-жёсткость связана с \(A_{55}\).

Для FSDT принята кинематическая связь

\[
\gamma_{xz}=w_{0,x}+\phi_x.
\]

Матрицы ламината и массовые моменты определены интегрированием по толщине:

\[
(A_{ij},B_{ij},D_{ij})
=\int_{-h/2}^{h/2}\bar Q_{ij}(1,z,z^2)\,dz,
\qquad
(I_0,I_1,I_2)
=\int_{-h/2}^{h/2}(1,z,z^2)\rho\,dz.
\]

Для рассматриваемых симметричных stacking sequences должны выполняться \(B\approx0\) и \(I_1\approx0\). Эти условия проверяются численно и не подменяют исходные определения.

## Карта страниц

| Содержание | Печатная страница | PDF index, zero-based |
|---|---:|---:|
| Преобразованные упругие характеристики, (2.3.17)–(2.3.18) | 96–97 | 118–119 |
| \(Q\), \(\bar Q\), membrane/shear ordering, (2.4.4)–(2.4.8) | 101 | 123 |
| Координата по толщине и нумерация слоёв | 113 | 135 |
| \(I_0,I_1,I_2\), (3.3.20c) | 122 | 144 |
| \(A,B,D\), (3.3.38) | 128 | 150 |
| FSDT shear strain, (3.4.3) | 134 | 156 |
| Shear correction и shear resultants, (3.4.10), (3.4.19), (3.4.22) | 135, 138, 139 | 157, 160, 161 |
| Вывод \(K=5/6\) для однородного прямоугольного сечения | 136 | 158 |
| Stacking convention и условие симметрии | 142, 148 | 164, 170 |
| Effective beam properties, (4.2.6), (4.2.8) | 168 | 190 |
| Beam boundary conditions, (4.2.11c) | 169 | 191 |
| Материальные отношения, (4.2.25) | 176 | 198 |
| Первые classical characteristic roots, Table 4.2.3 | 184 | 206 |
| Table 4.2.4 | 186 | 208 |
| Beam FSDT reduction, (4.3.1)–(4.3.9) | 187–188 | 209–210 |
| Free-vibration equations, (4.3.40)–(4.3.44) | 197 | 219 |
| Characteristic relations, (4.3.45)–(4.3.51) | 198 | 220 |
| Нормировка, (4.3.52) | 199 | 221 |
| Table 4.3.3 | 200 | 222 |
| Figures 4.3.3–4.3.4 | 201 | 223 |

## Source qualifications для Table 4.3.3

У Reddy координата \(z_R\) положительна вниз: верхняя поверхность
соответствует \(-h/2\), нижняя — \(+h/2\), а слои нумеруются сверху вниз.
Проект хранит стопку снизу вверх в координате \(z_P=-z_R\). Это соглашение
Python API, а не буквальная координатная конвенция книги. Для текущих
симметричных палиндромных стопок разворот не меняет benchmark; при будущем
рассмотрении \(B\ne0\) преобразование потребуется учитывать явно.

В cross-ply block напечатано `(90/0)_s`. Table 4.2.4, Figures 4.3.3–4.3.4 и численные значения Table 4.3.3 согласуются с `(0/90)_s = [0/90/90/0]`. В проекте используются

```text
printed_label = "(90/0)_s"
used_label = "(0/90)_s"
correction_status = "CORRECTED_BY_INTERNAL_SOURCE_CROSSCHECK"
```

Для воспроизведения таблицы принято \(K=5/6\) с provenance `INFERRED_FROM_INTERNAL_NUMERICAL_CONSISTENCY`. Книга явно выводит это значение для однородного прямоугольного сечения, но не задаёт его повторно в caption таблицы. Поэтому \(K\) остаётся явным входом и не переносится на произвольный ламинат как универсальная константа.

Для \(0^\circ\) и \(90^\circ\) таблица содержит buckling row и четыре frequency rows. Две classical-characteristic rows обозначаются как `source_classical_characteristic_with_RI` и `source_classical_characteristic_without_RI`; они не называются pure CLPT frequencies. Для `(0/90)_s` и `(45/-45)_s` напечатаны только buckling row, `fsdt_frequency_with_RI` и `source_classical_characteristic_with_RI`. Отсутствующие no-RI rows не реконструируются.

В первой колонке Table 4.3.3 однонаправленные случаи буквально подписаны
`0` и `90`. Machine contract сохраняет эти строки в `printed_label`; знаки
градуса используются только в `used_label` и в научном тексте.

При расчёте source-classical-characteristic rows используются именно
округлённые константы Table 4.2.3: \(\pi\), \(4{,}730\) и \(1{,}875\) для HH,
CC и CF соответственно. Они относятся только к воспроизведению процедуры
таблицы; аналитические проверки используют независимо найденные полноточные
корни.

Промежуточная clamped-boundary строка Eq. (4.3.50a) имеет статус `PRINTED_INTERMEDIATE_FORMULA_INCONSISTENCY`. Она не служит программным источником. Clamped boundary matrix выводится из физического first-order state, а Eq. (4.3.51) используется для вторичной проверки.

Полная числовая транскрипция и машинные статусы находятся в `tests/data/reddy_ch4_table_4_3_3.json`. Принятые решения подробно записаны в `docs/laminated_beams/reddy_ch4_source_contract.md`.

## Границы применимости

Источник используется только для одного прямого симметрично слоистого стержня. Он не задаёт coupled-rods joint, beta sweep, torsion, damping или FEM-модель проекта. Несимметричный случай \(B\ne0\) в текущую реализацию не входит.
