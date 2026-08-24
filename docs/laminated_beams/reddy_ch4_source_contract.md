# Source contract: Reddy Chapter 4, symmetric laminated beam

## Назначение и статус

Документ фиксирует источник для модели одного прямого симметрично слоистого стержня. Он отделяет напечатанные данные Reddy от принятых при воспроизведении уточнений. Текущий статус источника:

```text
RLB-0-SOURCE: PASS_WITH_DOCUMENTED_SOURCE_RECONSTRUCTION
```

Статус разрешает реализацию RLB-0A/B/C. Он не означает, что все численные проверки решателя уже выполнены.

Машиночитаемая транскрипция находится в `tests/data/reddy_ch4_table_4_3_3.json`.

## Идентификация источников

Основной источник:

- J. N. Reddy, *Mechanics of Laminated Composite Plates and Shells: Theory and Analysis*, 2nd ed., 2004;
- файл: `docs/literature/pdf/EB__Mechanics_of_Laminated_Composite_Plates_and_Shells_-JN_Reddy.pdf`;
- SHA-256: `A8D4D0FA67C7073D6EB48903B868BFA157875DC60F54484376A5B510547B37EA`;
- число PDF-страниц: 855;
- соответствие страниц: `pdf_index_zero_based = printed_page + 22`.

Независимый источник для проверки одномерных балансов и энергетически сопряжённых величин:

- В. В. Елисеев, *Механика упругих тел*, 1999;
- файл: `docs/literature/pdf/ELISEEV V.V._ MEXANIKA UPRUGIX TEL, 1999, 341s.pdf`;
- SHA-256: `615DE3610AEFC7967FEB1DF28A607AC09ADE7BDA9F756F0575D8BFC8A0606EF2`;
- число PDF-страниц: 341;
- соответствие страниц: `pdf_index_zero_based = printed_page - 1`.

Источник Елисеева не заменяет балочную редукцию Reddy и не задаёт значения Table 4.3.3.

## Карта формул Reddy

| Содержание | Печатная страница | PDF index, zero-based |
|---|---:|---:|
| Преобразованные упругие характеристики, (2.3.17)–(2.3.18) | 96–97 | 118–119 |
| Матрицы \(Q\), \(\bar Q\) и ordering компонент, (2.4.4)–(2.4.8) | 101 | 123 |
| Координата по толщине и нумерация слоёв | 113 | 135 |
| Массовые моменты \(I_0,I_1,I_2\), (3.3.20c) | 122 | 144 |
| Матрицы \(A,B,D\), (3.3.38) | 128 | 150 |
| Знак \(\gamma_{xz}=w_{0,x}+\phi_x\), (3.4.3) | 134 | 156 |
| Shear correction и transverse-shear resultants, (3.4.10), (3.4.19), (3.4.22) | 135, 138, 139 | 157, 160, 161 |
| Вывод \(K=5/6\) для однородного прямоугольного сечения и оговорка для общего ламината | 136 | 158 |
| Правило записи stacking sequence | 142 | 164 |
| Условие \(B=0\) для симметричного ламината | 148 | 170 |
| Эффективная изгибная жёсткость и массовые моменты балки, (4.2.6), (4.2.8) | 168 | 190 |
| Граничные условия, (4.2.11c) | 169 | 191 |
| Отношения свойств материала, (4.2.25) | 176 | 198 |
| Первые classical characteristic roots, Table 4.2.3 | 184 | 206 |
| Балочная редукция FSDT, (4.3.1)–(4.3.9) | 187–188 | 209–210 |
| Уравнения свободных колебаний, (4.3.40)–(4.3.44) | 197 | 219 |
| Частотные ветви и характеристические уравнения, (4.3.45)–(4.3.51) | 198 | 220 |
| Нормировка, (4.3.52) | 199 | 221 |
| Table 4.3.3 | 200 | 222 |

Для проверки cross-ply identity используются также Table 4.2.4 на печатной странице 186 (PDF index 208) и Figures 4.3.3–4.3.4 на печатной странице 201 (PDF index 223).

## Координаты, порядок компонент и формулы

У Reddy ось \(x\) направлена вдоль балки, \(y\) — по ширине, а координата
\(z_R\) по толщине положительна вниз. Верхняя поверхность имеет координату
\(z_R=-h/2\), нижняя — \(z_R=h/2\); слои нумеруются сверху вниз при
возрастании \(z_R\).

Python API использует отдельную project storage convention:
\(z_P=-z_R\), стопка перечисляется снизу вверх при возрастании \(z_P\),
\(z_{P,0}=-h/2\) и \(z_{P,N}=h/2\). Для текущих палиндромных симметричных
стопок это преобразование не меняет \(A\), \(D\) и benchmark frequencies.
Оно должно учитываться явно в будущей модели с \(B\ne0\). Длина, ширина,
полная толщина и угол слоя называются `length`, `width`, `thickness` и
`theta_deg`.

Membrane- и bending-компоненты имеют порядок \((xx,yy,xy)\). Порядок
transverse-shear-компонент равен \((yz,xz)\). Поэтому используемая балкой
\(xz\)-компонента является второй компонентой shear matrix.

Для ортотропного слоя

\[
\nu_{21}=\nu_{12}\frac{E_2}{E_1},\qquad
\Delta=1-\nu_{12}\nu_{21},
\]

\[
Q=
\begin{bmatrix}
E_1/\Delta & \nu_{12}E_2/\Delta & 0\\
\nu_{12}E_2/\Delta & E_2/\Delta & 0\\
0&0&G_{12}
\end{bmatrix}.
\]

При \(m=\cos\theta\), \(n=\sin\theta\) ненулевые элементы преобразованной
матрицы имеют вид

\[
\begin{aligned}
\bar Q_{11}&=Q_{11}m^4+2(Q_{12}+2Q_{66})m^2n^2+Q_{22}n^4,\\
\bar Q_{22}&=Q_{11}n^4+2(Q_{12}+2Q_{66})m^2n^2+Q_{22}m^4,\\
\bar Q_{12}&=(Q_{11}+Q_{22}-4Q_{66})m^2n^2+Q_{12}(m^4+n^4),\\
\bar Q_{16}&=(Q_{11}-Q_{12}-2Q_{66})m^3n
 -(Q_{22}-Q_{12}-2Q_{66})mn^3,\\
\bar Q_{26}&=(Q_{11}-Q_{12}-2Q_{66})mn^3
 -(Q_{22}-Q_{12}-2Q_{66})m^3n,\\
\bar Q_{66}&=(Q_{11}+Q_{22}-2Q_{12}-2Q_{66})m^2n^2
 +Q_{66}(m^4+n^4).
\end{aligned}
\]

В порядке \((yz,xz)\) transverse-shear block определяется равенствами

\[
\bar Q_{44}=G_{23}m^2+G_{13}n^2,\qquad
\bar Q_{55}=G_{23}n^2+G_{13}m^2,\qquad
\bar Q_{45}=(G_{13}-G_{23})mn.
\]

Для слоёв \(k=1,\ldots,N\)

\[
\begin{aligned}
A&=\sum_k \bar Q^{(k)}(z_k-z_{k-1}),\\
B&=\frac12\sum_k \bar Q^{(k)}(z_k^2-z_{k-1}^2),\\
D&=\frac13\sum_k \bar Q^{(k)}(z_k^3-z_{k-1}^3),\\
A_s&=\sum_k \bar Q_s^{(k)}(z_k-z_{k-1}),\\
I_j&=\int_{-h/2}^{h/2}\rho(z)z^j\,dz,\qquad j=0,1,2.
\end{aligned}
\]

Коэффициент \(K\) не входит в \(A_s\). Он применяется один раз при
одномерной редукции. После исключения свободных поперечных результирующих

\[
A_{\mathrm{axial}}=
\frac{\text{width}}{(A^{-1})_{11}},\qquad
D_{\mathrm{bending}}=
\frac{\text{width}}{(D^{-1})_{11}},\qquad
S_{\mathrm{shear}}=
K\,\frac{\text{width}}
{\boldsymbol e_{xz}^{\mathsf T}A_s^{-1}\boldsymbol e_{xz}},
\]

\[
m=\text{width}\,I_0,\qquad J=\text{width}\,I_2.
\]

Каждая эффективная жёсткость независимо вычисляется также через Schur
complement. Для допускаемой модели требуется \(B\simeq0\) и \(I_1\simeq0\).
Общий случай \(B\ne0\) не реализуется.

Знак деформации сдвига подтверждён по Eq. (3.4.3):

\[
\gamma_{xz}=w_{0,x}+\phi_x.
\]

При внутреннем обозначении \(\psi_b=\phi_x\) physical bending state равен

\[
\boldsymbol y_b=[w,\psi_b,Q,M]^{\mathsf T},
\]

а его уравнения имеют вид

\[
w'=\frac{Q}{S_{\mathrm{shear}}}-\psi_b,\qquad
\psi_b'=\frac{M}{D_{\mathrm{bending}}},\qquad
Q'=-m\omega^2w,\qquad
M'=Q-J\omega^2\psi_b.
\]

Используются граничные условия

- HH: \(w(0)=M(0)=w(L)=M(L)=0\);
- CC: \(w(0)=\psi_b(0)=w(L)=\psi_b(L)=0\);
- CF: \(w(0)=\psi_b(0)=Q(L)=M(L)=0\).

Кинетическая нормировка изгибной формы задаётся равенством

\[
\int_0^L\left(mw^2+J\psi_b^2\right)\,dx=1,
\]

а энергетическая проверка — равенством

\[
\int_0^L\left[D_{\mathrm{bending}}(\psi_b')^2
+S_{\mathrm{shear}}(w'+\psi_b)^2\right]\,dx
=\omega^2\int_0^L\left(mw^2+J\psi_b^2\right)\,dx.
\]

## Утверждённая реконструкция источника

### Cross-ply label

В Table 4.3.3 напечатано `(90/0)_s`. Значения таблицы, Table 4.2.4 и Figures 4.3.3–4.3.4 согласуются с последовательностью `(0/90)_s = [0/90/90/0]`. Поэтому контракт хранит обе метки:

```text
printed_label = "(90/0)_s"
used_label = "(0/90)_s"
correction_status = "CORRECTED_BY_INTERNAL_SOURCE_CROSSCHECK"
```

Текст и изображение книги не изменяются.

### Shear correction factor

Для benchmark Table 4.3.3 принимается

\[
K=\frac{5}{6},
\]

с provenance `INFERRED_FROM_INTERNAL_NUMERICAL_CONSISTENCY`. Reddy явно получает \(5/6\) для однородного прямоугольного сечения, но не повторяет это значение в caption Table 4.3.3. Поэтому \(K\) остаётся явным входным параметром и не считается универсальным свойством произвольного ламината. Для всех benchmark cases используется одно значение; отдельный подбор \(K\) запрещён.

Контроль source audit при \(a/h=100\) дал для buckling rows: \(20.460706\) вместо \(20.461\) для \(0^\circ\) и \(0.822061\) вместо \(0.822\) для \(90^\circ\). Основной предварительный gate использует напечатанные HH frequencies with rotary inertia для всех четырёх laminates и \(a/h=100,20,10\).

### Row roles

Для \(0^\circ\) и \(90^\circ\) напечатаны пять строк:

1. `buckling_load`;
2. `fsdt_frequency_with_RI`;
3. `source_classical_characteristic_with_RI`;
4. `fsdt_frequency_without_RI`;
5. `source_classical_characteristic_without_RI`.

Для `(0/90)_s` и `(45/-45)_s` напечатаны только строки 1–3. Строки 4–5 имеют `source_value = null`, `source_status = "NOT_PUBLISHED_IN_TABLE"` и не входят в source PASS/FAIL. Их численные значения не реконструируются. Собственный no-RI расчёт допускается только как project diagnostic.

Buckling rows сохраняются как source transcription. Buckling solver в рассматриваемой реализации не вводится.

Для строк `source_classical_characteristic_*` используются напечатанные в
Table 4.2.3 первые характеристические корни: \(\pi\) для HH, \(4{,}730\) для
CC и \(1{,}875\) для CF. Три знака после запятой сохраняются намеренно: это
часть воспроизводимой процедуры источника. Полноточные корни применяются в
отдельных аналитических gates и не подменяются табличными округлениями.

### Printed Eq. (4.3.50a)

Промежуточная clamped-boundary строка регистрируется со статусом

```text
PRINTED_INTERMEDIATE_FORMULA_INCONSISTENCY
```

Она не используется для сборки решателя. Boundary matrix выводится непосредственно из physical first-order state и условий \(w=\psi_b=0\). Eq. (4.3.51) служит только вторичной проверкой transfer-matrix roots.

### Однонаправленные ламинаты

\(0^\circ\) и \(90^\circ\) представляются одним слоем полной толщины. Обязательная partition-invariance check сопоставляет такой слой с произвольным числом одинаково ориентированных слоёв той же суммарной толщины. Должны совпадать \(A,B,D\), shear matrix, \(I_0,I_1,I_2\), reduced beam properties и frequencies.

## Полная транскрипция Table 4.3.3

В каждой ячейке значения следуют порядку \(a/h=100,20,10\). `HH`, `CC` и `CF` обозначают hinged–hinged, clamped–clamped и clamped–free соответственно. Все числа напечатаны с тремя десятичными знаками; абсолютная printed tolerance равна \(0.0005\).

Однонаправленные блоки в первой колонке книги буквально подписаны `0` и
`90`, без знака градуса. Поэтому machine contract хранит эти строки в
`printed_label`, а обозначения `0°` и `90°` — в `used_label`.

### \(0^\circ\), Tier A

| Row role | HH | CC | CF |
|---|---|---|---|
| `buckling_load` | 20.461 / 18.304 / 13.768 | 80.655 / 55.070 / 27.656 | 5.134 / 4.987 / 4.576 |
| `fsdt_frequency_with_RI` | 14.210 / 13.430 / 11.635 | 31.899 / 25.327 / 17.212 | 5.070 / 4.930 / 4.528 |
| `source_classical_characteristic_with_RI` | 14.210 / 13.430 / 11.635 | 32.110 / 28.506 / 22.140 | 5.070 / 4.965 / 4.675 |
| `fsdt_frequency_without_RI` | 14.211 / 13.441 / 11.657 | 31.824 / 24.636 / 16.680 | 5.063 / 4.813 / 4.229 |
| `source_classical_characteristic_without_RI` | 14.211 / 13.441 / 11.657 | 32.113 / 28.547 / 22.186 | 5.070 / 4.966 / 4.680 |

### \(90^\circ\), Tier A

| Row role | HH | CC | CF |
|---|---|---|---|
| `buckling_load` | 0.822 / 0.812 / 0.784 | 3.283 / 3.135 / 2.747 | 0.205 / 0.205 / 0.203 |
| `fsdt_frequency_with_RI` | 2.848 / 2.829 / 2.771 | 6.450 / 6.260 / 5.761 | 1.015 / 1.012 / 1.004 |
| `source_classical_characteristic_with_RI` | 2.848 / 2.829 / 2.771 | 6.454 / 6.356 / 6.079 | 1.015 / 1.012 / 1.005 |
| `fsdt_frequency_without_RI` | 2.848 / 2.832 / 2.781 | 6.449 / 6.232 / 5.681 | 1.015 / 1.009 / 0.993 |
| `source_classical_characteristic_without_RI` | 2.848 / 2.832 / 2.781 | 6.455 / 6.370 / 6.125 | 1.015 / 1.013 / 1.006 |

### `(0/90)_s`, Tier B

В первой колонке исходной книги этот блок подписан `(90/0)_s`.

| Row role | HH | CC | CF |
|---|---|---|---|
| `buckling_load` | 18.015 / 15.689 / 11.179 | 70.748 / 44.716 / 20.800 | 4.525 / 4.362 / 3.922 |
| `fsdt_frequency_with_RI` | 13.334 / 12.434 / 10.488 | 29.857 / 22.672 / 14.837 | 4.758 / 4.594 / 4.132 |
| `source_classical_characteristic_with_RI` | 13.334 / 12.434 / 10.488 | 30.106 / 26.041 / 19.504 | 4.759 / 4.636 / 4.307 |
| `fsdt_frequency_without_RI` | NOT_PUBLISHED | NOT_PUBLISHED | NOT_PUBLISHED |
| `source_classical_characteristic_without_RI` | NOT_PUBLISHED | NOT_PUBLISHED | NOT_PUBLISHED |

### `(45/-45)_s`, Tier B

| Row role | HH | CC | CF |
|---|---|---|---|
| `buckling_load` | 1.436 / 1.419 / 1.369 | 5.737 / 5.478 / 4.802 | 0.359 / 0.358 / 0.355 |
| `fsdt_frequency_with_RI` | 3.765 / 3.739 / 3.663 | 8.526 / 8.275 / 7.616 | 1.341 / 1.338 / 1.326 |
| `source_classical_characteristic_with_RI` | 3.765 / 3.739 / 3.663 | 8.531 / 8.402 / 8.036 | 1.341 / 1.338 / 1.328 |
| `fsdt_frequency_without_RI` | NOT_PUBLISHED | NOT_PUBLISHED | NOT_PUBLISHED |
| `source_classical_characteristic_without_RI` | NOT_PUBLISHED | NOT_PUBLISHED | NOT_PUBLISHED |

## Benchmark tiers

- Tier A содержит \(0^\circ\), \(90^\circ\), все напечатанные frequency rows, три значения \(a/h\) и три типа граничных условий.
- Tier B содержит исправленный `(0/90)_s`, `(45/-45)_s`, принятое \(K=5/6\) и только фактически напечатанные frequency rows.

Единое принятое значение \(K=5/6\) с указанным inferred provenance
используется в обоих tiers; отдельного подбора для Tier A нет.

Tier A и Tier B сравниваются и отчётно агрегируются раздельно. Для каждой comparison row сохраняются `printed_label`, `used_label`, `row_role`, `K`, `K_provenance`, `source_value`, `computed_value`, `absolute_error`, `relative_error`, `printed_tolerance` и `status`.

## Нормировка и ограничения

Table 4.3.3 использует

\[
\bar N=N_{xx}^{0}\frac{a^2}{E_2h^3},
\qquad
\bar\omega=\omega_1a^2\sqrt{\frac{I_0}{E_2h^3}}.
\]

Материальные отношения имеют вид

\[
\frac{E_1}{E_2}=25,
\qquad G_{12}=G_{13}=0.5E_2,
\qquad G_{23}=0.2E_2,
\qquad \nu_{12}=0.25.
\]

Контракт не распространяется на coupled rods, angular joint, torsion, damping, FEM и несимметричные ламинаты. Заданные source qualifications не должны скрываться или трактоваться как ошибка решателя.
