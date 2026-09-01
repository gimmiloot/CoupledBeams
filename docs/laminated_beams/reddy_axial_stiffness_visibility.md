# RLB-2H: visibility of axial stiffness under geometric coupling

## Scope

Этап RLB-2H строит ordinary frequency map по политике
[`frequency-map-v1`](../numerics/frequency_map_computation_policy.md) для
двух симметрично слоистых идентичных плеч Reddy. Внешняя геометрия и массовые
параметры фиксированы. Меняется только reduced axial stiffness
\(A_{\mathrm{beam}}\) через перераспределение in-plane жёсткости между
наружными и внутренними слоями.

## Frozen contract

\[
\mu=\tau=0,\qquad L_1=L_2=l=1,\qquad L_{\mathrm{total}}=2,
\]

\[
b_1=b_2=0.20,\qquad h_1=h_2=0.05,\qquad K=5/6.
\]

Каждое плечо содержит четыре равных слоя толщины
\(h_{\mathrm{ply}}=0.0125\). Во всех слоях используется угол
\(0^\circ\). Рассматриваются только два значения угла узла:
\(\beta=0^\circ\) и \(\beta=30^\circ\).

Базовый материал \(M_0\):

\[
E_{1,0}=1.1,\qquad E_{2,0}=0.9,\qquad \nu_{12,0}=0.3,
\]

\[
G_{12,0}=G_{13,0}=G_{23,0}=1/2.6,\qquad \rho_0=1.
\]

## A-only redistribution

Главный sweep parameter:

\[
\alpha_A=\frac{A_{\mathrm{beam}}}{A_0},\qquad
\alpha_A=0.70,0.72,\ldots,1.30.
\]

Стопка обоих плеч:

```text
OUTER / INNER / INNER / OUTER
```

Скалярные множители:

\[
s_{\mathrm{outer}}=\frac{4-\alpha_A}{3},\qquad
s_{\mathrm{inner}}=\frac{7\alpha_A-4}{3}.
\]

Масштабируются только in-plane параметры слоя:

\[
E_1=sE_{1,0},\qquad E_2=sE_{2,0},\qquad G_{12}=sG_{12,0},
\]

при неизменных
\(\nu_{12}=\nu_{12,0}\),
\(G_{13}=G_{13,0}\),
\(G_{23}=G_{23,0}\) и
\(\rho=\rho_0\).

Для четырёх равных 0°-слоёв аналитически фиксируются

\[
\frac{s_{\mathrm{outer}}+s_{\mathrm{inner}}}{2}=\alpha_A,
\qquad
\frac{7s_{\mathrm{outer}}+s_{\mathrm{inner}}}{8}=1.
\]

Следовательно ожидается:

\[
A(\alpha_A)=\alpha_A A_0,\qquad D(\alpha_A)=D_0,
\]

\[
S(\alpha_A)=S_0,\qquad m(\alpha_A)=m_0,\qquad J(\alpha_A)=J_0,
\]

а также \(B\simeq0\) и \(I_1\simeq0\).

## Constitutive gate

Constitutive gate получил статус `PASS` на всех 31 значениях
\(\alpha_A\). Получены диапазоны
\(A/A_0=[0.7,1.3]\),
\(D/D_0=[0.9999999999999997,1]\) и
\(S/S_0=m/m_0=J/J_0=1\).

Максимальные относительные residuals равны:

| Проверка | Максимальный residual |
| --- | ---: |
| scaling матрицы \(A\) | \(2.093\cdot10^{-16}\) |
| инвариантность матрицы \(D\) | \(2.386\cdot10^{-16}\) |
| формула \(A/A_0=\alpha_A\) | \(3.828\cdot10^{-16}\) |
| совокупная инвариантность \(D,S,m,J,K,b\) | \(3.696\cdot10^{-16}\) |
| \(B\) | \(1.508\cdot10^{-17}\) |
| \(I_1\) | \(0\) |

Таким образом, в reduced beam properties на принятой сетке меняется только
\(A_{\mathrm{beam}}\). Геометрия, transverse shear stiffness, погонная масса
и вращательная инерция сохраняются.

## Local frequency-map instance

```yaml
frequency_map_policy: frequency-map-v1
calculation_mode: fast_plot
spectrum_semantics: sorted_positions
sweep_parameter: alpha_A
parameter_grid: 0.70:0.02:1.30
continuation_anchor: 1.00
continuation_paths:
  - 1.00:0.98:0.70
  - 1.00:1.02:1.30
K_plot: 8
K_guard: 9
guard_root_role: completeness_only
neighbour_audit: enabled
local_repair_policy: triggered_only
strict_audit_default: false
branch_tracking: false
mac: false
mode_shapes: false
energy_analysis: false
```

В каждой точке сохраняются independently sorted positions 1--8 и root 9 как
completeness guard. Отдельный intentional search или audit spectral tail выше
root 9 не выполняется, а строки выше guard не экспортируются. Текущая версия
entry point отклоняет distinct pre-trim slot выше root 9; исторические уже
завершённые группы после добавления этого gate не пересчитывались. Predictor
используется только для локализации search windows; exported frequencies всегда
получаются из characteristic matrix.

## beta = 0 mechanism and reference

При \(\beta=0\) и одинаковых свойствах двух плеч конструкция должна
совпадать с прямой fixed-fixed Reddy-балкой общей длины
\(L_{\mathrm{total}}=2\), имеющей ту же reduced section и однородной вдоль
длины. Продольная и bending-подсистемы остаются разделёнными:

- isolated bending spectrum должен быть invariant по \(\alpha_A\);
- axial frequencies должны следовать
  \(\omega_{\mathrm{axial}}\propto\sqrt{\alpha_A}\);
- при используемой нормировке
  \(\Lambda_{\mathrm{axial}}\propto\alpha_A^{1/4}\).

Изменение independently sorted positions при \(\beta=0\) допускается только
из-за движения axial family относительно постоянных bending levels.

Для этой проверки не выполнялся новый независимый поиск bending roots.
Завершённая coupled characteristic group при \(\alpha_A=1\), содержащая
позиции 1--9, была повторно использована без пересчёта. Её частоты ранее
получены determinant/SVD refiner из coupled characteristic matrix. Одна
продольная частота, попавшая в первые девять позиций, сопоставлена с exact
axial family. Остальные восемь частот образуют единственный bending reference.

Для каждого \(\alpha_A\) exact axial family объединяется с этим bending
reference. При объединении учитываются multiplicity, total nullity и cluster
centre. Получено 279 reference rows, то есть 31 группа по девять позиций.
Новых вызовов root solver для этой reference-конструкции не было.

В точках \(\alpha_A=0.70,1.00,1.30\) direct fixed-fixed boundary matrix
общей длины 2 независимо проверялась на сингулярность на частотах, заданных
exact axial и единственным bending reference. Direct boundary matrix не
локализовала эти частоты независимо. Coupled roots также не пересчитывались.
Поэтому эта проверка подтверждает boundary singularity, но не является вторым
независимым частотным решением. Direct diagnostic отдельно сохраняет detected
nullity каждой поданной частоты; cluster grouping при этом задаётся subsystem
union и не считается независимо локализованным direct inventory.

Обе beta=0 проверки получили статус `PASS`. Ошибка инвариантности bending
matrix равна \(2.125\cdot10^{-14}\), а ошибка закона
\(\Lambda_{\mathrm{axial}}\propto\alpha_A^{1/4}\) составляет
\(1.891\cdot10^{-16}\). Максимальная cluster-aware относительная разность
между coupled spectrum и exact axial+bending union равна
\(1.330\cdot10^{-11}\). На тех же заданных union frequencies все три direct
boundary matrices прошли singularity quality gate. Отдельная разность между
независимо локализованными direct и coupled frequencies не вычислялась,
поскольку direct localization не выполнялась.

## beta = 30 finite-grid result

При \(\beta=30^\circ\) используется полный coupled determinant без
artificial decoupling. Научный вопрос ограничен конечной grid и формулируется
осторожно: может ли изменение только \(A_{\mathrm{beam}}\), не влияющее на
isolated bending subsystem, менять independently sorted coupled spectrum
через geometric rigid joint.

Для описания результата введены величины

\[
\delta_k=
\frac{\Lambda_k(1.30)-\Lambda_k(0.70)}
{\max\left(|\Lambda_k(0.70)|,|\Lambda_k(1.30)|\right)},
\]

\[
d_k=\max_{\alpha_A}
\frac{|\Lambda_k(\alpha_A)-\Lambda_k(1)|}{|\Lambda_k(1)|}.
\]

| Sorted position \(k\) | \(\delta_k\) | \(d_k\) | Поведение на grid |
| ---: | ---: | ---: | --- |
| 1 | \(7.260932\cdot10^{-5}\) | \(4.737544\cdot10^{-5}\) | `NONDECREASING` |
| 2 | \(4.636169\cdot10^{-2}\) | \(2.882694\cdot10^{-2}\) | `NONDECREASING` |
| 3 | \(7.544155\cdot10^{-2}\) | \(4.021077\cdot10^{-2}\) | `NONDECREASING` |
| 4 | \(3.550492\cdot10^{-4}\) | \(2.423897\cdot10^{-4}\) | `NONDECREASING` |
| 5 | \(2.219910\cdot10^{-2}\) | \(1.114106\cdot10^{-2}\) | `NONDECREASING` |
| 6 | \(4.577981\cdot10^{-2}\) | \(4.194196\cdot10^{-2}\) | `NONDECREASING` |
| 7 | \(1.154930\cdot10^{-1}\) | \(6.538451\cdot10^{-2}\) | `NONDECREASING` |
| 8 | \(8.643971\cdot10^{-3}\) | \(4.874540\cdot10^{-3}\) | `NONDECREASING` |

На принятой сетке все восемь independently sorted positions не убывают, но
их относительное изменение существенно различается. Наибольшее изменение
между концами сетки получено для позиции 7, а наименьшее -- для позиции 1.
Следовательно, в рассматриваемой модели изменение только
\(A_{\mathrm{beam}}\) видно в coupled sorted spectrum при
\(\beta=30^\circ\). Этот вывод не задаёт modal identity.

## Normalization

Используется прежняя reference normalization:

\[
\Omega=\omega l^2\sqrt{\frac{\rho_0A_{\mathrm{ref}}}{E_0I_{y0}}},
\qquad
\Lambda=\sqrt{\Omega},
\]

\[
A_{\mathrm{ref}}=b_0h_0,\qquad
I_{y0}=\frac{b_0h_0^3}{12},
\qquad
E_0=\rho_0=l=1.
\]

Здесь `reference_area = b0*h0` и `axial_stiffness = A_beam` разделены
именами, чтобы нормировка не скрывала исследуемый эффект.

## Frequency-map quality and performance

Обе карты завершены: получены 31/31 группы при \(\beta=0^\circ\) и 31/31
группы при \(\beta=30^\circ\). Inventory содержит 558 строк `BASE` и 62
root-9 guards. Корни выше root 9 не вычислялись.

Neighbour audit отметил 14 точек. Во всех случаях локальный пересчёт получил
статус `REPRODUCED_AFTER_LOCAL_REPAIR`; сохранено 126 строк
`LOCAL_REFINEMENT`. Неразрешённых точек нет. Сглаживание и использование
predictor или interpolation вместо частоты не применялись.

Минимальные расстояния между соседними sorted positions равны:

- \(\Delta\Lambda_{\min}=0.0100513\) при \(\beta=0^\circ\) и
  \(\alpha_A=0.82\) между позициями 6 и 7;
- \(\Delta\Lambda_{\min}=0.244122\) при \(\beta=30^\circ\) и
  \(\alpha_A=1.30\) между позициями 7 и 8.

Эти интервалы являются только диагностикой близости соседних упорядоченных
частот. Они не устанавливают crossing, veering или обмен модальным характером.

Максимальные значения \(\sigma_{\min}/\sigma_{\max}\) и boundary residual
среди строк `BASE` равны соответственно \(4.673\cdot10^{-17}\) и
\(3.797\cdot10^{-16}\). Минимальный right margin root 9 по \(\Omega\) равен
2. Root-quality gate получил статус `PASS`.

Шесть production anchors дали conservative ETA 206.7 s для 56 оставшихся
точек при лимите 2700 s. Сохранённая сумма measured times равна 66.4 s и
является нижней границей. Времена и evaluation counters для 14 локальных
пересчётов не были сохранены до прерванной первой финализации, поэтому эти
точки не пересчитывались только ради метрик. Известная нижняя граница числа
determinant и SVD evaluations равна 248193 для каждого типа вычислений.
Максимальный RSS всего workflow составил 139624448 bytes, то есть 133.2 MiB.

Frozen production-physics hashes и result-tree hashes этапов RLB-2A, RLB-2D,
RLB-2E, RLB-2F и RLB-2G сохранились.

## Reproducibility and outputs

Основной или missing-only запуск:

```powershell
python scripts/analysis/laminated_beams/sweep_reddy_axial_stiffness_visibility.py --missing-only
```

Повторная отрисовка готовых CSV без поиска корней:

```powershell
python scripts/analysis/laminated_beams/sweep_reddy_axial_stiffness_visibility.py --plot-only
```

Вывод contract manifest в stdout без спектрального расчёта и изменения
готовых outputs:

```powershell
python scripts/analysis/laminated_beams/sweep_reddy_axial_stiffness_visibility.py --manifest-only
```

Targeted regressions:

```powershell
python -m pytest -q -p no:cacheprovider tests/test_reddy_axial_stiffness_visibility.py
```

Результаты находятся в
`results/laminated_beams/reddy_axial_stiffness_visibility/`:

- `spectrum_roots.csv`;
- `section_properties.csv`;
- `neighbour_audit.csv`;
- `beta0_subsystem_reference.csv`;
- `benchmark.json` и `checkpoint.json`;
- `lambda_vs_A_ratio_beta0_beta30.png`;
- `beta0_decoupled_reference.png`;
- `report.md` и `run_manifest.json`.

## Status and limitations

Все пять status gates получили `PASS`:

```text
RLB-2H-CONSTITUTIVE-A-ONLY: PASS
RLB-2H-BETA0-DIRECT-SUBSYSTEM-REFERENCE: PASS
RLB-2H-FREQUENCY-MAP: PASS
RLB-2H-ROOT-INVENTORY: PASS
OVERALL: PASS
```

Даже при `PASS` результат относится только к:

- четырём равным 0°-слоям;
- двум identical arms;
- \(\beta\in\{0^\circ,30^\circ\}\);
- конечной grid \(\alpha_A=0.70:0.02:1.30\);
- independently sorted positions 1--8 и root 9 guard.

Этап не устанавливает modal identity, localization, veering, branch exchange
или обобщение на другие геометрии, углы, stacking sequences и несимметричные
ламинаты. Branch tracking, mode shapes, MAC, energy partition, Ritz, FEM и
certified audit не выполнялись.
