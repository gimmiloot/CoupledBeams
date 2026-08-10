# AP-0: screening применимости спектрального приближения при theta = 2 deg

> Этот диагностический этап относится только к прямоугольным
> анизотропным/ортотропным стержням Chapter 2 Ярцева. Он не использует
> круглые изотропные стержни, isotropic-steel defaults, article
> safe-prefix/certification workflows или FEM.

## Назначение и границы этапа

AP-0 проверяет на фиксированной малой сетке, насколько прямоугольная
ортотропная модель Эйлера--Бернулли с кручением Сен-Венана воспроизводит
первые шесть independently sorted значений `Lambda` полной Chapter-2 модели
при `theta1=theta2=2 deg`. Computational workflow имеет статус `PASS`;
научная классификация отдельных геометрий записывается отдельно и не является
computational PASS/FAIL.

Этап не определяет точную допустимую границу по `theta`, не уточняет сетку по
`theta` или `beta`, не отслеживает modal descendants и не применяет MAC,
mode shapes, energy classification либо FEM. T2 используется как
одномерный reference для сравнения моделей, а не как утверждение об абсолютной
точности реальной трёхмерной конструкции.

## Три модели

- `T2`: HMS/DX-209, elastic, `theta1=theta2=2 deg`, Chapter-2
  `state_corrected` Timoshenko bending, generalized rectangular torsion,
  rotated material properties и `Sbar16` coupling.
- `T0`: ортотропный предел той же Chapter-2 модели при `theta=0`.
- `EB0`: существующий validated rectangular orthotropic
  Euler--Bernoulli/Saint-Venant comparator при `theta=0`.

Все модели используют `book_slope_clamp` и неизменённый `J_book(beta)`.
Сравниваются sorted spectral positions 1--6; root 7 служит completeness guard
и входит в расчёт gaps 1--2, ..., 6--7. Совпадающий индекс не означает общий
modal descendant.

## Объёмно-сохраняющая геометрия

Зафиксированы `l=0.4 m`, `a0=0.005 m`, `b0=0.020 m` и `kappa=4`, то есть
`b_i=4 a_i`. Параметры определены как

```text
L1 = l(1-mu)
L2 = l(1+mu)
tau = (a2-a1)/(a1+a2)
s = a0/sqrt(1+tau^2+2*mu*tau)
a1 = s(1-tau)
a2 = s(1+tau)
b1 = 4*a1
b2 = 4*a2
```

Используется только знаковый параметр `tau`; старый thickness parameter в
этом workflow отсутствует. Hard gates требуют положительных размеров,
`L1+L2=0.8 m`, `b_i/a_i=4` с tolerance `1e-14` и
`abs(V-V0)/V0 <= 1e-12`, где `V0=8.0e-5 m^3`. Hard slenderness threshold не
вводился; `L_i/a_i` и `L_i/b_i` сохраняются только как diagnostics.

Сетка содержит ровно `mu={0,0.25,0.5}` и `tau={-0.2,0,0.2}`. Полученные
размеры и volume residuals:

| geometry | mu | tau | L1, m | L2, m | a1, mm | a2, mm | b1, mm | b2, mm | volume residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `mu_p0p00_tau_m0p20` | 0 | -0.2 | 0.4 | 0.4 | 5.883484 | 3.922323 | 23.533936 | 15.689291 | 1.694e-16 |
| `mu_p0p00_tau_p0p00` | 0 | 0 | 0.4 | 0.4 | 5.000000 | 5.000000 | 20.000000 | 20.000000 | 0 |
| `mu_p0p00_tau_p0p20` | 0 | 0.2 | 0.4 | 0.4 | 3.922323 | 5.883484 | 15.689291 | 23.533936 | 1.694e-16 |
| `mu_p0p25_tau_m0p20` | 0.25 | -0.2 | 0.3 | 0.5 | 6.188527 | 4.125685 | 24.754110 | 16.502740 | 1.694e-16 |
| `mu_p0p25_tau_p0p00` | 0.25 | 0 | 0.3 | 0.5 | 5.000000 | 5.000000 | 20.000000 | 20.000000 | 0 |
| `mu_p0p25_tau_p0p20` | 0.25 | 0.2 | 0.3 | 0.5 | 3.746343 | 5.619515 | 14.985373 | 22.478059 | 0 |
| `mu_p0p50_tau_m0p20` | 0.5 | -0.2 | 0.2 | 0.6 | 6.546537 | 4.364358 | 26.186147 | 17.457431 | 1.694e-16 |
| `mu_p0p50_tau_p0p00` | 0.5 | 0 | 0.2 | 0.6 | 5.000000 | 5.000000 | 20.000000 | 20.000000 | 1.694e-16 |
| `mu_p0p50_tau_p0p20` | 0.5 | 0.2 | 0.2 | 0.6 | 3.592106 | 5.388159 | 14.368424 | 21.552636 | 1.694e-16 |

Во всех случаях volume равен `8.0e-5 m^3`, а масса при
`rho=1580 kg/m^3` равна `0.1264 kg`.

## Сетка beta и fixed normalization

Screening grid: `beta=0,5,...,90 deg`, всего 19 точек. Для всех моделей и
геометрий применяется одна reference normalization

```text
Lambda_ref = (rho*A0*omega^2*l^4/(E_x0*I_y0))^(1/4),
A0 = a0*b0,
I_y0 = a0^3*b0/12,
E_x0 = 191 GPa,
rho = 1580 kg/m^3.
```

Физическая T2 model при этом использует фактические rotated properties при
2 deg. Reference normalization не подменяет их и не зависит от текущей
геометрии.

## Метрики

Для `k=1,...,6`:

```text
delta_beam,k   = abs(Lambda_k(T0)-Lambda_k(EB0))/Lambda_k(T0)
delta_orient,k = abs(Lambda_k(T2)-Lambda_k(T0))/Lambda_k(T2)
delta_simpl,k  = abs(Lambda_k(T2)-Lambda_k(EB0))/Lambda_k(T2)
```

`Delta_beam`, `Delta_orient` и `Delta_simpl` являются максимумами по шести
позициям. Научный screening criterion применяется только к
`Delta_simpl_max=max_beta Delta_simpl`: значение не более `0.10` получает
метку `WITHIN_10_PERCENT_ON_SCREENING_GRID`, иначе
`EXCEEDS_10_PERCENT_ON_SCREENING_GRID`.

Для T0 и T2 по roots 1--7 вычисляются

```text
g_k = 2*(Lambda[k+1]-Lambda[k])/(Lambda[k+1]+Lambda[k])
g_min_0 = min_k g_k(T0)
G_nearest = g_kstar(T2)-g_kstar(T0), kstar=argmin_k g_k(T0)
G_open = max_k (g_k(T2)-g_k(T0))
G_close = min_k (g_k(T2)-g_k(T0)).
```

Новый эмпирический критерий по gaps не вводился.

## Root quality, fast solver и reuse

Fast-solver validation evidence имеет `FAST_SOLVER_PASS` и ошибку относительно
legacy oracle менее `1e-8`. AP-0 переиспользует exact transfer cache,
predictor/local windows, cluster handling, global inventories, обязательные
anchors, atomic checkpoints и существующий root-quality finalizer. Если
local/predictor path не проходит, screening-local coordinator сохраняет уже
проверенный global inventory; generic fast-solver checkpoint contract не
изменён.

Каждый spectrum содержит семь конечных положительных строго возрастающих
roots. Gate принимает scaled quality или normalized physical-raw quality при
determinant и singular residual не более `1e-8`; rejected roots запрещены.
Итог: 3591 accepted roots, `0` rejected, максимум determinant residual
`1.6594e-14`, максимум relative singular residual `1.7552e-12`.

Точно переиспользованы пять families:

- `mu=0,tau=0`: T2 из Figure 10, T0 и EB0 из Figure 3;
- `mu=0.25,tau=0`: T0 и EB0 из Figure 5.

Для всех reuse rows beta/frequency/Lambda arrays совпадают с source arrays
через `array_equal`; максимальный residual повторной fixed normalization
равен `2.90e-16`. Figure 6 не использован, поскольку его сечения не
удовлетворяют текущему similar-section contract `b_i=4a_i`.

В полном запуске рассчитаны 22 новые и переиспользованы 5 families. Recorded
scientific family runtime — `244.919 s`; выполнены 154 global anchors,
219 global inventories, 46 automatic global fallbacks и 411248 transfer
exponential evaluations.

## Результаты девяти геометрий

| geometry | Delta_beam,max | Delta_orient,max | Delta_simpl,max | beta/mode simpl | min baseline gap | beta/pair | max G_open | beta/pair | classification |
|---|---:|---:|---:|---|---:|---|---:|---|---|
| `mu_p0p00_tau_m0p20` | 0.010444 | 0.032138 | 0.032177 | 25/1 | 0.002747 | 70/5--6 | 0.036060 | 50/4--5 | WITHIN_10_PERCENT_ON_SCREENING_GRID |
| `mu_p0p00_tau_p0p00` | 0.010631 | 0.021827 | 0.025442 | 0/4 | 0.000225 | 75/4--5 | 0.033906 | 0/2--3 | WITHIN_10_PERCENT_ON_SCREENING_GRID |
| `mu_p0p00_tau_p0p20` | 0.010444 | 0.026582 | 0.035693 | 15/5 | 0.002747 | 70/5--6 | 0.042676 | 20/5--6 | WITHIN_10_PERCENT_ON_SCREENING_GRID |
| `mu_p0p25_tau_m0p20` | 0.009029 | 0.024951 | 0.025153 | 20/1 | 0.011274 | 0/3--4 | 0.026696 | 90/3--4 | WITHIN_10_PERCENT_ON_SCREENING_GRID |
| `mu_p0p25_tau_p0p00` | 0.010631 | 0.027554 | 0.029521 | 35/4 | 0.003489 | 0/5--6 | 0.043738 | 90/4--5 | WITHIN_10_PERCENT_ON_SCREENING_GRID |
| `mu_p0p25_tau_p0p20` | 0.009253 | 0.026109 | 0.029974 | 0/5 | 0.011566 | 85/6--7 | 0.038254 | 90/6--7 | WITHIN_10_PERCENT_ON_SCREENING_GRID |
| `mu_p0p50_tau_m0p20` | 0.010477 | 0.028643 | 0.031590 | 35/3 | 0.014149 | 0/5--6 | 0.037648 | 90/3--4 | WITHIN_10_PERCENT_ON_SCREENING_GRID |
| `mu_p0p50_tau_p0p00` | 0.010631 | 0.026183 | 0.030360 | 35/3 | 0.003489 | 0/5--6 | 0.037195 | 90/3--4 | WITHIN_10_PERCENT_ON_SCREENING_GRID |
| `mu_p0p50_tau_p0p20` | 0.010703 | 0.028513 | 0.033125 | 40/3 | 0.000378 | 0/5--6 | 0.036321 | 55/3--4 | WITHIN_10_PERCENT_ON_SCREENING_GRID |

На finite grid все `9/9` геометрий находятся в пределах 10%; превышений нет.
Максимальное наблюдённое `Delta_simpl_max=0.035693` существенно ниже 0.10.
Это не устанавливает точную boundary по `theta` и ничего не утверждает о
значениях между beta-grid points.

## Проверка гипотезы о плотности исходного спектра

По 171 точке Spearman correlation `g_min_0` с `Delta_orient` равна `+0.3770`,
а с `G_open` — `+0.05335`. В lowest-gap quartile median значения равны
`Delta_orient=0.01965` и `G_open=0.02951`; на остальных точках — `0.02479` и
`0.03196`. Поэтому предварительная тенденция «меньший исходный gap — больший
orientation effect/opening» на этой finite grid не подтверждается как
монотонная. Это описательная проверка без regression law и без утверждения
причинности.

## Candidates

- stable: `mu=0.25,tau=-0.2`, `Delta_simpl_max=0.025153`;
- sensitive: `mu=0,tau=0.2`, `Delta_simpl_max=0.035693`;
- borderline: `NO_CLOSE_BORDERLINE_CASE_FOUND`; ближайший случай находится
  от порога 0.10 на `0.064307`;
- smallest baseline gap: `mu=0,tau=0`, `0.00022477`;
- largest gap opening: `mu=0.25,tau=0`, `0.043738`.

Эти случаи только записаны в `candidate_cases.json`; следующий gate не
запускался.

## Запуск и outputs

```text
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --smoke
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --resume
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --reuse-data
```

Smoke использует три геометрии и `beta={0,45,90 deg}`, пишет в отдельный
`theta2_smoke/` и не является scientific baseline. Полные CSV, JSON, report,
atomic family checkpoints и три diagnostic figures находятся в
`results/anisotropic_rods/yartsev_ch2_spectral_applicability_screening/theta2_small_grid/`.
Figures S1--S3 — это heatmap `Delta_simpl_max`, heatmap maximum `G_open` и
scatter `g_min_0` versus `Delta_orient`; PDF являются vector, PNG записаны с
300 dpi. Отдельные 27 `Lambda(beta)` figures не создаются.

Supervisor Figures 1--12, их `plot_manifest.json` и `report.md` не
перезаписываются. Точный theta boundary, local beta refinement, MAC, shapes,
energy analysis, FEM и 3D FEM остаются вне AP-0.

## AP-1 — theta=5 deg same-grid screening

AP-1 повторяет AP-0 на точно тех же девяти геометриях и на той же сетке
`beta=0,5,...,90 deg`, заменяя основную модель T2 на T5:

- `T5`: Chapter-2 `state_corrected` Timoshenko/generalized torsion,
  `theta1=theta2=5 deg`, actual rotated HMS/DX-209 properties и `Sbar16`;
- `T0`: неизменённый theta-zero Chapter-2 reference;
- `EB0`: неизменённый rectangular orthotropic EB/Saint-Venant comparator.

Clamp, joint, seven-root inventory, six compared sorted positions, fixed
`Lambda_ref` и quality thresholds полностью совпадают с AP-0. T5 является
одномерным reference для сравнения моделей, а не абсолютной 3D truth model.

### Geometry и reuse

`theta5_small_grid/geometry_manifest.csv` является byte-identical копией
AP-0 manifest. Независимая проверка подтверждает девять geometry IDs,
`L1+L2=0.8 m`, `b_i/a_i=4`, массу `0.1264 kg` и maximum volume residual
`1.694e-16`.

AP-1 не пересчитывает reference spectra:

- все девять T0 и девять EB0 families точно переиспользованы из
  `theta2_small_grid/screening_spectra.csv`;
- baseline T5 `mu=0,tau=0` точно переиспользован из supervisor Figure 7;
- только остальные восемь T5 families рассчитаны заново.

Reuse audit содержит 19 строк. Frequency, `Lambda`, root status, quality basis
и quality residual arrays совпадают с источниками точно. Pointwise
`Delta_beam` также array-equal AP-0 во всех 171 точках. AP-1 checkpoint
fingerprint включает target theta и source hashes, поэтому theta2 checkpoint
не может быть принят как theta5 family.

### Theta=2 versus theta=5

| geometry | Delta beam | Delta orient, 2 | Delta orient, 5 | increment | Delta simpl, 2 | Delta simpl, 5 | increment | beta/mode at 5 | transition |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| `mu_p0p00_tau_m0p20` | 0.010444 | 0.032138 | 0.113165 | +0.081027 | 0.032177 | 0.114324 | +0.082147 | 0/5 | WITHIN_AT_2_EXCEEDS_AT_5 |
| `mu_p0p00_tau_p0p00` | 0.010631 | 0.021827 | 0.101097 | +0.079270 | 0.025442 | 0.104277 | +0.078835 | 0/2 | WITHIN_AT_2_EXCEEDS_AT_5 |
| `mu_p0p00_tau_p0p20` | 0.010444 | 0.026582 | 0.108316 | +0.081734 | 0.035693 | 0.119683 | +0.083990 | 5/5 | WITHIN_AT_2_EXCEEDS_AT_5 |
| `mu_p0p25_tau_m0p20` | 0.009029 | 0.024951 | 0.109879 | +0.084928 | 0.025153 | 0.110140 | +0.084987 | 15/1 | WITHIN_AT_2_EXCEEDS_AT_5 |
| `mu_p0p25_tau_p0p00` | 0.010631 | 0.027554 | 0.108001 | +0.080447 | 0.029521 | 0.108639 | +0.079119 | 10/1 | WITHIN_AT_2_EXCEEDS_AT_5 |
| `mu_p0p25_tau_p0p20` | 0.009253 | 0.026109 | 0.105952 | +0.079843 | 0.029974 | 0.116186 | +0.086212 | 0/5 | WITHIN_AT_2_EXCEEDS_AT_5 |
| `mu_p0p50_tau_m0p20` | 0.010477 | 0.028643 | 0.116045 | +0.087402 | 0.031590 | 0.117252 | +0.085661 | 20/2 | WITHIN_AT_2_EXCEEDS_AT_5 |
| `mu_p0p50_tau_p0p00` | 0.010631 | 0.026183 | 0.104257 | +0.078073 | 0.030360 | 0.109375 | +0.079015 | 90/3 | WITHIN_AT_2_EXCEEDS_AT_5 |
| `mu_p0p50_tau_p0p20` | 0.010703 | 0.028513 | 0.105896 | +0.077383 | 0.033125 | 0.111126 | +0.078001 | 80/3 | WITHIN_AT_2_EXCEEDS_AT_5 |

На AP-0 grid все `9/9` geometry были within 10%. При `theta=5 deg` within
осталось `0/9`, exceeding стало `9/9`; зарегистрировано девять переходов
`WITHIN_AT_2_EXCEEDS_AT_5`. Это не определяет threshold-crossing angle:
промежуточные theta не вычислялись и интерполяция не применялась.

### Same-pair gap analysis

AP-1 сохраняет 1026 строк
`screening_pairwise_gap_metrics.csv`: девять geometries, 19 beta и шесть
конкретных соседних пар. Для одной и той же пары `k-(k+1)` записаны `g_k` при
theta 0, 2 и 5 deg, изменения относительно T0 и дополнительное изменение
2→5 deg.

Описательные результаты:

- pooled Spearman `g_pair_T0` versus `abs_delta_g_theta5`: `-0.467406`;
- per-pair: 1–2 `+0.345312`, 2–3 `-0.205777`, 3–4 `-0.555495`,
  4–5 `+0.113905`, 5–6 `-0.364918`, 6–7 `+0.532093`;
- median absolute change в lowest-gap quartile: `0.075242`, в остальных
  observations: `0.032977`;
- opened-pair fraction: `0.945525` в lowest-gap quartile и `0.394018` вне
  него.

Таким образом, pooled finite-grid data показывают descriptive same-pair
tendency к большему изменению при меньшем исходном gap, но результаты сильно
различаются по sorted pair. Наблюдения не считаются независимыми; p-values,
regression law и причинный вывод отсутствуют.

### Quality, performance и candidates

Приняты все 3591 roots, rejected roots отсутствуют. Максимумы determinant и
relative singular residual равны `1.6594e-14` и `1.7552e-12`. Восемь новых
T5 families потребовали `56` global anchors, `76` inventories и `23`
automatic fallbacks; recorded scientific family runtime `94.042 s`.

Candidates следующего gate только записаны:

- stable и nearest-to-threshold: `mu=0,tau=0`,
  `Delta_simpl_5_max=0.104277`;
- sensitive: `mu=0,tau=0.2`, `0.119683`;
- transition candidates: все девять geometries;
- smallest pairwise gap: `mu=0,tau=0,beta=75`, pair 4–5,
  `g_pair_T0=0.00022477`;
- largest pairwise change: `mu=0.5,tau=-0.2,beta=0`, pair 4–5,
  `abs_delta_g_theta5=0.180549`.

### AP-1 outputs и ограничения

```text
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 5 --smoke
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 5 --resume
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 5 --reuse-data
```

AP-1 outputs находятся только в `theta5_small_grid/`; smoke — в отдельном
`theta5_smoke/`. Дополнительно к AP-0 schema записаны pairwise table,
point/geometry theta2-vs-theta5 comparisons и Figures AP1-S1–S3. AP-0
outputs и supervisor Figures 1–12/manifest/report остаются byte-identical.

AP-1 не искал точную theta boundary, не запускал theta=3/4 full grids,
theta/beta refinement, MAC, descendants, forms, energy analysis, FEM/3D FEM,
damping или complex roots. Circular-rod workflows и article workspaces не
использовались; следующий gate не запускался.

## AP-2 — intermediate sampled angles 3 and 4 deg

AP-2 сохраняет геометрию, `beta=0,5,...,90 deg`, fixed
`Lambda_ref`, seven-root inventory, `book_slope_clamp`, `J_book(beta)`,
HMS/DX-209 и quality thresholds AP-0/AP-1. Theta обозначает ориентацию
материальных осей относительно локальной продольной оси; beta обозначает
геометрический угол сопряжения. T3/T4 — более полные одномерные Chapter-2
reference-модели, но не точные 3D-модели.

Для каждой ориентации рассчитаны только восемь новых target families. T0 и
EB0 для всех девяти геометрий переиспользованы из AP-0; baseline
`mu=0,tau=0` точно переиспользован из supervisor Figure 11 при theta=3 и
Figure 12 при theta=4. Reuse audits имеют по 19 PASS-строк. Geometry manifests
theta3/theta4 байт-в-байт равны AP-0 manifest.

### Pointwise и family-level classification

| sampled theta | pointwise within / 171 | pointwise exceeds / 171 | uniformly within families / 9 | exceeding families / 9 |
|---:|---:|---:|---:|---:|
| 2 | 171 | 0 | 9 | 0 |
| 3 | 171 | 0 | 9 | 0 |
| 4 | 171 | 0 | 9 | 0 |
| 5 | 66 | 105 | 0 | 9 |

На prescribed finite grid все configurations удовлетворяют 10%-му критерию
при sampled orientations 2, 3 и 4 deg. При 5 deg превышают 105 из 171
configurations, и каждая из девяти families содержит хотя бы одно превышение.
Это не устанавливает `Delta_simpl(theta)<=0.10` при каждом непрерывном theta
между 2 и 4 deg. Первый превышающий sampled angle не называется критическим
углом.

| mu | tau | Delta simpl max: 2 / 3 / 4 / 5 deg | beta/k: 2 / 3 / 4 / 5 deg | first exceeding sampled theta |
|---:|---:|---:|---|---:|
| 0.00 | -0.20 | 0.032177 / 0.055380 / 0.082501 / 0.114324 | 25/1, 25/1, 20/1, 0/5 | 5 |
| 0.00 | 0.00 | 0.025442 / 0.047743 / 0.071503 / 0.104277 | 0/4, 0/4, 80/4, 0/2 | 5 |
| 0.00 | +0.20 | 0.035693 / 0.058969 / 0.087444 / 0.119683 | 15/5, 10/5, 10/5, 5/5 | 5 |
| 0.25 | -0.20 | 0.025153 / 0.048257 / 0.076960 / 0.110140 | 20/1, 20/1, 15/1, 15/1 | 5 |
| 0.25 | 0.00 | 0.029521 / 0.049671 / 0.076966 / 0.108639 | 35/4, 15/1, 10/1, 10/1 | 5 |
| 0.25 | +0.20 | 0.029974 / 0.053299 / 0.082615 / 0.116186 | 0/5, 0/5, 0/5, 0/5 | 5 |
| 0.50 | -0.20 | 0.031590 / 0.054110 / 0.083786 / 0.117252 | 35/3, 25/2, 25/2, 20/2 | 5 |
| 0.50 | 0.00 | 0.030360 / 0.051821 / 0.078451 / 0.109375 | 35/3, 50/3, 70/3, 90/3 | 5 |
| 0.50 | +0.20 | 0.033125 / 0.054563 / 0.080814 / 0.111126 | 40/3, 50/3, 65/3, 80/3 | 5 |

Все 171 pointwise `Delta_simpl` sequences и все девять family maxima sampled-
nondecreasing при заранее фиксированном absolute tolerance `1e-12`; re-entry
не обнаружен. Median pointwise increments растут от `0.020974` на 2→3 deg к
`0.025741` на 3→4 deg и `0.029689` на 4→5 deg. Для family maxima medians равны
`0.022520`, `0.027295` и `0.032240`. Это описывает постепенный sampled рост с
увеличивающимися приращениями, но не доказывает непрерывную монотонность.
Maximizing sorted position изменяется 72 раза на pointwise interval
comparisons и четыре раза на family-level comparisons; beta family maximum
изменяется 17 раз. Эти изменения относятся к sorted positions, а не к
отслеженным physical modal descendants.

Classification decomposition A/B/C составляет `0/0/171` при theta=3,
`0/0/171` при theta=4 и `76/29/66` при theta=5. Это classification
decomposition, а не аддитивный error budget.

### Sampled same-pair gap diagnostics

Pooled Spearman `g_pair_T0` versus `abs_delta_g` равен `-0.364097`,
`-0.427609`, `-0.443532`, `-0.467406` при theta 2, 3, 4, 5 deg. Pooled median
absolute changes равны `0.012229`, `0.022272`, `0.032795`, `0.040415`.
Доля observations с увеличением absolute gap change составляет `0.870370`,
`0.879142`, `0.873294` на интервалах 2→3, 3→4, 4→5. Sign reversals: 82, 71,
53. Per-pair Spearman сохраняется отдельно в combined summary; знаки
неоднородны по шести sorted pairs. Это descriptive pairwise tendency on the
finite screening grid без p-values, independence assumption, regression law
или причинного вывода.

### AP-2 outputs и команды

```text
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 3 --smoke
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 4 --smoke
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 3 --resume
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 4 --resume
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 3 --reuse-data
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 4 --reuse-data
python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --combine-sampled-theta
```

Target outputs находятся в `theta3_small_grid/` и `theta4_small_grid/`;
combined tables/report и Figures AP2-S3/S4 — в
`theta2_theta3_theta4_theta5_comparison/`. AP2-S1/S2 — target heatmaps.
PDF являются vector, PNG записаны с 300 dpi. AP-0, AP-1 и supervisor outputs
не перезаписываются. Exact theta crossing, additional theta sweep, beta
refinement, MAC, modes/shapes/energy, FEM/3D FEM, damping и circular-rod
workflows не запускались.
