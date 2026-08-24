# Проверка одного симметрично слоистого стержня Reddy

## Результат

Проверка выполнена 24 августа 2026 г. для замороженного source contract и
единого benchmark-входа \(K=5/6\).

| Этап | Статус | Основание |
|---|---|---|
| RLB-0-SOURCE | `PASS_WITH_DOCUMENTED_SOURCE_RECONSTRUCTION` | хеши и страницы записаны; 144 опубликованных значения проверены визуально; 36 отсутствующих значений оставлены `null` |
| RLB-0A-BENDING | `PARTIAL_PASS` | аналитические gates прошли; 86 из 108 source comparisons находятся в printed tolerance |
| RLB-0B-AXIAL | `PASS` | первые шесть FF и FC частот и форм совпали с аналитическими решениями |
| RLB-0C-COMBINED-UNION | `PASS` | для двух физических случаев проверено по 6 axial и 6 bending членов, а также точное и близкое вырождения |
| OVERALL | `PARTIAL_PASS` | RLB-0A не имеет полного source-table PASS |

Полные построчные данные находятся в
`results/laminated_beams/reddy_symmetric_single_beam/`. Итоговый статус не
повышался за счёт изменения \(K\), плотности, \(E_2\), нормировки, граничных
условий или допуска.

## Источники и соглашения

Основной PDF Reddy имеет SHA-256
`A8D4D0FA67C7073D6EB48903B868BFA157875DC60F54484376A5B510547B37EA`.
Использованы печатные страницы 96–97, 101, 113, 122, 128, 134–136, 138–139,
142, 148, 168–169, 176, 184, 186–188, 197–201; соответствующие zero-based
PDF indices равны 118–119, 123, 135, 144, 150, 156–158, 160–161, 164, 170,
190–191, 198, 206, 208–210, 219–223.

Независимый PDF Елисеева имеет SHA-256
`615DE3610AEFC7967FEB1DF28A607AC09ADE7BDA9F756F0575D8BFC8A0606EF2`.
Использованы печатные страницы 141–143, 151–155, 162, 164–166, 240,
242–243; zero-based indices 140–142, 150–154, 161, 163–165, 239, 241–242.

У Reddy координата \(z_R\) положительна вниз. Проектная стопка хранится
снизу вверх при \(z_P=-z_R\). В Table 4.3.3 однонаправленные labels
напечатаны как `0` и `90`; научные обозначения `0°` и `90°` хранятся в
`used_label`.

Cross-ply block напечатан как `(90/0)_s`, но Table 4.2.4 и Figures
4.3.3–4.3.4 соответствуют `(0/90)_s=[0/90/90/0]`. В расчёте сохранены

```text
printed_label = "(90/0)_s"
used_label = "(0/90)_s"
correction_status = "CORRECTED_BY_INTERNAL_SOURCE_CROSSCHECK"
```

## Модель

Порядок membrane- и bending-компонент равен \((xx,yy,xy)\), порядок
transverse shear — \((yz,xz)\). Матрицы \(A\), \(B\), \(D\), \(A_s\) и
массовые моменты \(I_0,I_1,I_2\) интегрируются по слоям. Для допускаемой
редукции численно проверяются \(B\simeq0\) и \(I_1\simeq0\).

Эффективные свойства вычисляются через compliance element и независимо через
Schur complement:

\[
A_{\mathrm{axial}}=\frac{\mathrm{width}}{(A^{-1})_{11}},\qquad
D_{\mathrm{bending}}=\frac{\mathrm{width}}{(D^{-1})_{11}},\qquad
S_{\mathrm{shear}}=K\frac{\mathrm{width}}
{\boldsymbol e_{xz}^{\mathsf T}A_s^{-1}\boldsymbol e_{xz}}.
\]

Изгибный state и его уравнения:

\[
\boldsymbol y_b=[w,\psi_b,Q,M]^{\mathsf T},
\]

\[
w'=Q/S_{\mathrm{shear}}-\psi_b,\quad
\psi_b'=M/D_{\mathrm{bending}},\quad
Q'=-m\omega^2w,\quad
M'=Q-J\omega^2\psi_b.
\]

Граничные условия имеют вид HH: \(w=M=0\), CC: \(w=\psi_b=0\), CF:
\(w(0)=\psi_b(0)=Q(L)=M(L)=0\). Printed Eq. (4.3.50a) не используется
при сборке. Eq. (4.3.51) является только вторичной проверкой.

Продольный state

\[
\boldsymbol y_a=[u,N]^{\mathsf T},\qquad
u'=N/A_{\mathrm{axial}},\qquad N'=-m\omega^2u
\]

получен в проекте из symmetric CLT и стандартной динамики стержня. Combined
state \([u,w,\psi_b,N,Q,M]^{\mathsf T}\) является точным блочно-диагональным
объединением этих подсистем.

## Численное основание для \(K=5/6\)

Все 12 предварительно утверждённых HH/RI controls прошли допуск
\(0{,}0005\).

| Laminate | \(a/h\) | Source | Calculated | \(|\Delta|\) | Статус |
|---|---:|---:|---:|---:|---|
| \(0^\circ\) | 100 | 14.210 | 14.209948 | 0.000052 | PASS |
| \(0^\circ\) | 20 | 13.430 | 13.429632 | 0.000368 | PASS |
| \(0^\circ\) | 10 | 11.635 | 11.635330 | 0.000330 | PASS |
| \(90^\circ\) | 100 | 2.848 | 2.848290 | 0.000290 | PASS |
| \(90^\circ\) | 20 | 2.829 | 2.828859 | 0.000141 | PASS |
| \(90^\circ\) | 10 | 2.771 | 2.770977 | 0.000023 | PASS |
| \((0/90)_s\) | 100 | 13.334 | 13.333579 | 0.000421 | PASS |
| \((0/90)_s\) | 20 | 12.434 | 12.434104 | 0.000104 | PASS |
| \((0/90)_s\) | 10 | 10.488 | 10.487532 | 0.000468 | PASS |
| \((45/-45)_s\) | 100 | 3.765 | 3.765021 | 0.000021 | PASS |
| \((45/-45)_s\) | 20 | 3.739 | 3.739370 | 0.000370 | PASS |
| \((45/-45)_s\) | 10 | 3.663 | 3.662954 | 0.000046 | PASS |

Это одно глобальное значение \(K\), а не отдельный fit по laminate или case.
Оно является явным параметром API и сохраняется вместе с provenance во всех
generated CSV/JSON; те же сведения записаны в metadata двух PNG.

## Использованные и пропущенные строки Table 4.3.3

| Набор | Buckling transcription | Frequency comparisons | Missing rows |
|---|---:|---:|---:|
| \(0^\circ\), \(90^\circ\), Tier A | 18 | 72 | 0 |
| \((0/90)_s\), \((45/-45)_s\), Tier B | 18 | 36 | 36 |
| Итого | 36 | 108 | 36 |

Для Tier A использованы четыре frequency roles:
`fsdt_frequency_with_RI`, `source_classical_characteristic_with_RI`,
`fsdt_frequency_without_RI` и
`source_classical_characteristic_without_RI`. Для Tier B использованы только
первые две напечатанные frequency roles с rotary inertia. Отсутствующие
`fsdt_frequency_without_RI` и
`source_classical_characteristic_without_RI` имеют `source_value=null`, статус
`NOT_PUBLISHED_IN_TABLE` и исключены из PASS/FAIL. Buckling solver не
реализован; 36 buckling values сохранены только как транскрипция.

Строки `source_classical_characteristic_*` рассчитаны с напечатанными в Table
4.2.3 корнями \(\pi\), \(4{,}730\) и \(1{,}875\). Полноточные аналитические
gates решаются независимо.

## Сверка Table 4.3.3

| Tier | PASS | FAIL | Всего | Максимальное \(|\Delta|\) |
|---|---:|---:|---:|---:|
| A | 57 | 15 | 72 | 0.703491 |
| B | 29 | 7 | 36 | 0.046374 |
| Итого | 86 | 22 | 108 | 0.703491 |

Структура 22 расхождений:

| Tier | Row role | BC | Число | Максимальное \(|\Delta|\) |
|---|---|---|---:|---:|
| A | `fsdt_frequency_with_RI` | CF | 4 | 0.031813 |
| A | `fsdt_frequency_without_RI` | CC | 6 | 0.703491 |
| A | `fsdt_frequency_without_RI` | CF | 5 | 0.336560 |
| B | `fsdt_frequency_with_RI` | CC | 2 | 0.000629 |
| B | `fsdt_frequency_with_RI` | CF | 4 | 0.046374 |
| B | `source_classical_characteristic_with_RI` | CC | 1 | 0.000628 |

Таким образом, 21 расхождение относится к direct FSDT roots, одно — к
source-classical diagnostic. Все 72 direct transfer roots, включая допустимые
project no-RI diagnostics для отсутствующих multilayer rows, приняты SVD-gate.
Матрицы CC и CF собраны непосредственно из физического state. Изменение
printed Eq. (4.3.50a) или подгонка входов не выполнялись.

Полная построчная таблица, включая `printed_label`, `used_label`, `row_role`,
\(K\), provenance, обе частоты, ошибки, допуск и статус, записана в
`bending_source_comparison.csv`.

## Аналитические и численные gates

| Gate | Максимальная ошибка или residual | Порог | Статус |
|---|---:|---:|---|
| \(B\), \(I_1\), scaled symmetry | \(1.5051\cdot10^{-18}\) | \(10^{-12}\) | PASS |
| Compliance / Schur | \(5.9690\cdot10^{-16}\) | \(10^{-11}\) | PASS |
| Partition invariance, \(0^\circ\) и \(90^\circ\) | \(2.1316\cdot10^{-16}\) | \(10^{-12}\) | PASS |
| HH, \(n=1,2,3\), обе ветви: analytic / transfer | 0 | \(10^{-10}\) | PASS |
| HH boundary residual | \(1.1523\cdot10^{-14}\) | \(10^{-9}\) | PASS |
| HH energy identity | \(1.5772\cdot10^{-15}\) | \(10^{-8}\) | PASS |
| HH \(1-\mathrm{MAC}\) | \(2.2204\cdot10^{-16}\) | \(10^{-8}\) | PASS |
| Eq. (4.3.51) / independent CC transfer root | 0 | \(10^{-10}\) | PASS |
| Axial analytic / transfer, первые шесть FF и FC | 0 | \(10^{-10}\) | PASS |
| Axial boundary residual | \(1.8931\cdot10^{-15}\) | \(10^{-9}\) | PASS |
| Combined union, isolated roots | \(2.8333\cdot10^{-14}\) | \(10^{-9}\) | PASS |
| Independent combined inventories | \(6.1735\cdot10^{-14}\) | \(10^{-9}\) | PASS |
| Combined energy identity | \(1.0018\cdot10^{-9}\) | \(10^{-8}\) | PASS |
| Euler–Bernoulli large-\(S\), \(J=0\) limit | \(3.1064\cdot10^{-13}\) | \(10^{-10}\) | PASS |
| \(E_2,\rho\) scale invariance of \(\bar\omega\) | 0 | \(10^{-12}\) | PASS |

Оба combined cases используют full 6×6 endpoint matrix: cross-ply FF+CC и
angle-ply FC+CF. В каждом случае подтверждены шесть axial и шесть bending
членов. Для двух построенных near-degenerate cases cluster multiplicity равна
2, каждый корень имеет nullity 1, а leakage противоположного семейства равен
нулю. В точном вырождении nullity равна 2; residual разделённого axial/bending
basis равен \(2.1204\cdot10^{-16}\), его mass inner product равен нулю.
Произвольная смесь SVD-векторов не интерпретируется как физическая связь.

## Качество корней и обусловленность

| Набор | Принято | Максимальный singular ratio | Максимальный конечный condition number | Оценки \(\infty\) |
|---|---:|---:|---:|---:|
| Bending, source diagnostics и HH gates | 78/78 | \(7.0417\cdot10^{-10}\) | \(1.8014\cdot10^{16}\) | 9 |
| Combined full 6×6 | 28/28 | \(8.6771\cdot10^{-17}\) | \(2.8210\cdot10^{18}\) | 4 |

Большие и бесконечные оценки condition number сохранены в CSV и не скрыты.
Значение \(\infty\) возникает, когда вычисленный \(\sigma_{\min}\) равен нулю
в машинной арифметике у искомой сингулярной boundary matrix. Приём корня
требует одновременно конечной положительной частоты, SVD-gate, допустимого
boundary residual и отсутствия rejected status.

## Формы и generated outputs

Изгибные и продольные формы нормированы по массе. Для изгиба проверены сетки
401 и 801 точка, энергетическое равенство и детерминированный знак первой
существенной амплитуды \(w(x)\). Combined kinetic shares вычислены интегрально,
а не назначены по имени семейства.

Созданы ровно две диагностические фигуры:

- `bending_modes.png` — первые три mass-normalized lower-branch HH modes;
- `combined_spectrum.png` — axial и bending families для двух случаев.

`--plot-only` построил их из сохранённых таблиц без импорта solver. Каталог
results подтверждён как Git-ignored.

## Ограничения

Результат не распространяется на два сопряжённых стержня, угловой узел,
beta sweep, несимметричный ламинат, \(B\ne0\), torsion, damping, FEM или
complex roots. Arbitrary vectors в вырожденном подпространстве не являются
свидетельством продольно-изгибной связанности. Глобальные проектные индексы и
журнал изменений не обновлялись, поскольку они содержат защищённые локальные
изменения пользователя.
