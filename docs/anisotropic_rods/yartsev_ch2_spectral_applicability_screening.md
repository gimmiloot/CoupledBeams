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
