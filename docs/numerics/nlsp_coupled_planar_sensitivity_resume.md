# Ограниченное продолжение пилота чувствительности двух стержней
**Итог: `NUMERICAL_PARTIAL`.** Mixed/sign checks и оба полного finite-δ
контроля проходят. Получено наблюдаемое усиление `Gperp^(4)=1.00989868` для A
и `1.03622225` для B. Методные gates E1–E4 пройдены при p16; отдельная
пространственная qualification сравнения A/B не выполнена в пределах бюджета.
Это первый ограниченный A/B результат, без stability certificate или карты.

Продолжение имеет отдельное разрешение пользователя и собственный bundle
`results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/`.
Исторический [пилот](nlsp_coupled_planar_sensitivity_pilot.md) остаётся
`BLOCKED_BY_TANGENT_CHECK`. Независимая
[крутильная диагностика](nlsp_coupled_torsion_tangent_diagnostic.md) установила
`NUMERICAL_FLOOR_IDENTIFIED`; её результаты здесь переиспользуются.

## Постановка и сохранённые источники

Проверяется развитие малых начальных внеплоскостных возмущений внутри
принятой семиполевой модели `q=(u,w,v,Phi,psi,theta,c)`. Физические kernels,
quartic action, переменная масса, нелинейная инерция, RHS/Jacobian, коэффициенты,
осевые соглашения и узловые условия не менялись. Касательный оператор получен
ограничением полного аналитического Jacobian на инвариантный OUT-блок.
Условие `c1=c2` и прежние внешние value clamps сохранены. Дополнительных
ограничений производных, демпфирования или жёсткостей нет.

Обе конструкции имеют `beta=30deg`, `mu=.20`, `L1=.8`, `L2=1.2`,
`E=rho=1`, `nu=.3`, `kappa=5/6`. В A оба сечения `b=.10,h=.05`.
В B только второе сечение повёрнуто: `b=.05,h=.10`. Площадь и масса
сохраняются; прежние моменты инерции, изгибные коэффициенты и H соответствуют
этим размерам. Общий узел имеет семь общих координат, без дополнительной массы.

Нагрузки и собственные масштабы времени берутся из исторического
`frozen_excitation.json`. Линейный критерий выбора нагрузки
`max|w|/min(h)=.05` не переопределяется после результатов чувствительности.
Собственные нелинейные статические равновесия освобождаются мгновенно,
с нулевыми начальными скоростями основного движения.

| Конструкция | g | Tref | 2Tref |
|---|---:|---:|---:|
| A | 7.820867770536735e-5 | 35.40765544870329 | 70.81531089740658 |
| B | 2.1428767135661873e-4 | 21.110243236568515 | 42.22048647313703 |

Проверены все артефакты обоих parent manifests. SHA256 исходного пилота:
`7a149b76f64165a300e27d33f7cb1f2876a5f3094d5193d72f0c591c57564d80`;
диагностики: `d889f8ef5e090f9b1de11a8ba9dbdaec60963e42e31f5a4750e6857a7249a2f9`.
Физические kernels совпадают с сохранёнными execution-code hashes A.
Сохранённая base trajectory A p16 и её accepted Radau polynomials используются
непосредственно. Исторические статические, изгибные, крутильные, узловые,
геометрические и Jacobian проверки не выполнялись заново.

## Mixed и зеркальная симметрия

Начальное направление — равновесовая по E0 комбинация первого сохранённого
изгибного и крутильного координатных направлений. Это энергетическая
нормировка, не требование чистоты форм. Начальные скорости этого контроля нулевые.
Горизонт `.1TrefA=3.540765544870329`; все 186 OUT-компонентов сохраняются.

| delta | Max relative tangent/full difference | Coordinate K error max | Velocity M error max |
|---|---:|---:|---:|
| 7.958509826438746e-7 | 1.4709922342590958e-7 | 1.262208754726987e-7 | 9.3529077725229e-8 |
| 3.979254913219373e-7 | 3.6775065686660104e-8 | 3.155558278039838e-8 | 2.3382783208538468e-8 |
| −7.958509826438746e-7 | 1.4709922342590958e-7 | 1.262208754726987e-7 | 9.3529077725229e-8 |

Наблюдаемый порядок 1.9999898068; порядок изменения плоского состояния
1.9999980518. Ошибка при меньшем delta уменьшается примерно в четыре раза.
В sampled maximum около31.69% квадрата ошибки приходится на скорости.
Подробные вклады обоих плеч и raw internal/joint groups сохранены в JSON.
Signed group shares учитывают cross terms и не являются отдельными
положительными энергиями безмассового узла.

Во всех сохранённых состояниях `OUT(+delta)+OUT(−delta)=0`, а плоские
состояния совпадают побитово. Это проверка зеркальной симметрии принятой
дискретной модели, без нового физического вывода об устойчивости.
Gate A — PASS; старые threshold `1e-5`, floor `1e-8` и status не ослаблены.

Полные trajectories решены прежним Radau с `rtol=1e-12`, componentwise
`atol=old*.01`, `max_step=old*.5`, как в независимой диагностике.
Для линейной касательной ODE используется DOP853 с тем же полным Jacobian,
`rtol=atol_relative=1e-11`. Короткий контроль против сохранённого Radau
даёт разность `4.170929188639539e-10` в E0 norm при gate `1e-8`.
Математический оператор при смене интегратора не меняется.

## Подпространство и замороженный ресурсный план

Первоначальные шесть входных направлений сокращены до четырёх до длинных
production results: одно низшее изгибное и первое направление с крутильным
участием, каждое в координатной и скоростной форме. Для A используются
сохранённые столбцы `[0,2,3,5]`; для B применяется то же правило отбора.
Каждый входной вектор нормируется в собственной фиксированной E0 norm.
Выход сохраняет все93 координаты и93 скорости OUT при p16.

`||z||E0^2=eta.T Kperp,0 eta+eta_dot.T Mperp,0 eta_dot`.
Матрицы положительно определены и не меняются вдоль base trajectory.
Максимальное усиление определяется через Gram matrix полного выходного
состояния четырёх начальных направлений. Обозначение `Gperp^(4)` относится
именно к этому входному подпространству. Максимум по времени берётся
на сохранённой сетке2001 точка; бесконечномерный supremum не установлен.

`frozen_long_plan.json` сохранён до B и длинных gain results. Прогноз
6602.726s для оставшейся программы при остатке6861.961s — planning estimate,
не гарантия. Часть numerical budget338.039s уже затрачена на mixed preflight.
Сохраняются потолки8 full/16 tangent ODE,7200s,4GiB,один поток; no retry.
Для fit в бюджет заранее исключены A p12 dynamic и длинный tight temporal
control. Предусмотрен один B p12 control с физической проекцией p16 initial
state/subspace и сохранением общих узловых координат. Нового A static solve нет.

Для длинного DOP853 dense output вычисляется только при requested samples.
По установленному SciPy этот вывод требует ещё три RHS evaluations,
но не меняет accepted state и шаги интегратора. Начальная и новая output
версии кода сохранены отдельно; mock regression подтверждает сохранение samples.
План использует консервативное снижение стоимости на10%, а не обещает
пропорциональную экономию. Общая физическая сетка comparison не сдвигается;
velocity берётся из state, без дифференцирования interpolated displacement.

## Фактическое движение B и усиление A/B

B p16 nonlinear Newton выполнен один раз: residual `2.65942e-13`, без
перевыбора нагрузки. Linear max|w|=`.00249999759905576`, nonlinear
max|w|=`.00244614289514882`; энергия освобождения `1.4503105347712061e-9`.
После снятия нагрузки B достигла `42.22048647313703` за25446 accepted steps.
Дрейф энергии `3.19665e-12`. Инвариантность, совместность узла и safety
проверены по actual saved states; на каждом RHS действует прежний safety gate.

| Конструкция | Sampled max Gperp^(4) | Рост нормы относительно 1 | tmax/Tref | tmax | Zero-base Gram error |
|---|---:|---:|---:|---:|---:|
| A | 1.00989868251 | 0.989868% | 1.381 | 48.89797217466 | 1.58917e-11 |
| B | 1.03622224713 | 3.622225% | 1.060 | 22.37685783076 | 8.87568e-12 |

Разность наблюдаемых gain ratios равна `.0263235646261`. B имеет больший
результат при данном протоколе возбуждения и горизонте двух собственных
Tref. Этот вывод ещё не имеет независимой spatial qualification. Сравнение
не является изолированным опытом при одинаковом g или одинаковом physical T:
обе нагрузки заранее выбраны одним amplitude rule, а временные масштабы
собственные. Zero-base norm остаётся1 для обеих линейных конструкций.

В выбранном входном пространстве наиболее чувствительны преимущественно
скоростные изгибные направления, с малой примесью крутильного направления.
Коэффициенты ниже относятся к E0-orthonormal seeds, не к необработанным
Galerkin coefficients. Общий знак вектора не меняет gain.

| Геометрия | Bend coordinate | Torsional coordinate | Bend velocity | Torsional velocity |
|---|---:|---:|---:|---:|
| A | -0.322229118 | -0.017855693 | -0.944149883 | -0.066562517 |
| B | 0.036484479 | -0.034697666 | -0.994853630 | -0.087927302 |

Сохранены physical initial coordinates **и velocities**. Начальные скорости
основного движения нулевые; найденное конечное возмущение включает
ненулевые OUT velocities. Initial mode frequencies A: `.143173920832`,
`1.548126632132`; B: `.101621143957`, `.763090329265`. Крутильные mass fractions
соответственно `.00497198/.78768635` и `.00051936/.52980354`. Это описание
состава двух выбранных линейных направлений, не branch tracking или
требование чистых изгибных/крутильных форм.

| Геометрия | Arm 1 / arm 2 E0 fractions в maximum | Torsional fraction, nonlinear base | Torsional fraction, zero base при том же входе/времени |
|---|---|---:|---:|
| A | 44.6756% / 55.3244% | 10.8317% | 10.9053% |
| B | 35.9346% / 64.0654% | 0.7226% | 0.7287% |

Различие bending/torsion состава отклика в основном присутствует уже в
линейной конструкции. Поэтому долю кручения нельзя автоматически считать
доказательством nonlinear transfer или конкретного parametric resonance.
Формы и энергия двух плеч сохранены отдельно от scalar gain.

## Подтверждение полной семиполевой динамикой

Для каждой геометрии интегрирована одна full434-state trajectory из
`X0(0)+delta*z0`, где z0 — найденное направление. Delta rule заранее
фиксирована: `.03*sqrt(2*Ebase)`; амплитуда не менялась для улучшения
совпадения. При floating-point пересчёте сохранённой энергии A последний
разряд delta отличается от historical value примерно на1e-15 относительно;
физическая величина и правило не изменяются.

| Геометрия | delta | Full/tangent relative mismatch max | Coordinate K error max | Velocity M error max | Full energy drift | Gate |
|---|---:|---:|---:|---:|---:|---|
| A | 7.95850982643876e-07 | 4.4869182e-07 | 3.6047839e-07 | 4.4098406e-07 | 2.83306e-14 | PASS |
| B | 1.6157224274572e-06 | 1.3853658e-06 | 1.2080446e-06 | 1.3729574e-06 | 2.49303e-14 | PASS |

В обоих случаях max relative mismatch меньше прежнего `1e-5`. Согласие
подтверждено по всем OUT координатам и скоростям, не только по одному
перемещению. Математическая асимптотика на полном2Tref заново не измерялась
серией delta: выполнен один finite-amplitude контроль каждой геометрии.

| Геометрия | Max physical v | Max abs Phi, rad | Max abs psi, rad | v / fixed linear w scale .0025 | Plane change E0 norm max |
|---|---:|---:|---:|---:|---:|
| A | 9.1964996e-05 | 0.000228993125 | 0.000147761595 | 3.67860% | 9.818762e-10 |
| B | 0.000262051946 | 0.000305688519 | 0.000439514991 | 10.48208% | 6.269957e-09 |

Здесь максимум берётся по сохранённым2001 times и201 points на каждом плече.
Общий displacement scale `.0025` задан до результатов. Более крупный physical
v у B не тождественен разности gain ratios: unit energy states и собственные
base energies двух конструкций дают различные physical amplitudes.

Saved-state safety (не строгий continuum bound): min mass bound A full
`.999238606`, B base `.992301690`, B full `.992300456`. Max rotation norms
`.00674857/.00968370/.00968364`; max Li*curvature `.0531446/.0719901/.0719894`.
Все saved states finite, mass positive; подробные strain/gradient и joint
checks находятся в `saved_state_safety_joint.json`.

## Выполнение, численная qualification и остановка

Все11 IVPs завершены без solver failure и без retry. Суммарно6 full nonlinear
(3mixed+B base+2optimal finite),5tangent (short+2main+2zero),1B Newton.
Нового A base/static, p24/p48/p64, FEM, geometry B variants или карты нет.

| Actual case | Seconds | Accepted steps |
|---|---:|---:|
| A_short_tangent | 103.040 | 983 |
| A_mixed_0 | 77.353 | 5679 |
| A_mixed_1 | 80.252 | 5969 |
| A_mixed_2 | 77.394 | 5679 |
| B_p16_base | 331.996 | 25446 |
| A_p16_tangent | 1639.627 | 19646 |
| A_zero_tangent | 8.480 | 19646 |
| B_p16_tangent | 1368.341 | 13806 |
| B_zero_tangent | 5.173 | 11713 |
| A_optimal_finite | 1761.680 | 122659 |
| B_optimal_finite | 1066.691 | 74942 |

Total numerical seconds `6520.0275943/7200`; actual scientific process peak
working set `3.632984 GiB/4GiB`, one BLAS thread. Venv launcher working set
в отдельном wrapper log не является scientific memory measurement;
правильный lifetime peak сохранён в `actual_process_memory.json`.
Newton timing учитывается отдельно в `B_static_execution.json`.

После фактических B counters read-only p12 pointwise benchmark повысил
прогноз control до`1316.3856s` (он предполагает не больше accepted steps,
что не сертифицировано). После двух finite controls осталось`679.9724s`.
Поэтому B p12 **NOT_RUN**. `updated_spatial_launch_decision.json` и отдельный
preflight reservation фиксируют отказ до solver call. Директория B_p12_base
содержит только NOT_RUN administrative marker, не trajectory или equilibrium.

| Gate / qualification | Status |
|---|---|
| E1 mixed / sign consistency | PASS |
| E2 B base equilibrium / free motion / safety | PASS |
| E3 одинаково выбранное S4, fixed norm, zero controls, A/B gains | PASS |
| E4 optimal finite perturbations, full seven-field confirmation | PASS при p16 |
| Separate dynamic spatial control | NOT_RUN — resource preflight |
| Independent short tangent temporal control | PASS |
| Independent long tangent temporal control | NOT_RUN; combined temporal qualification PARTIAL |
| Overall | NUMERICAL_PARTIAL |

Read-only A p16→p12 truncation по всем saved base states показывает max
changes c=`2.65679e-6/5.43785e-6` для плеч1/2 и theta=`2.07390e-6/8.79641e-6`.
Это projection sensitivity, **не p12 dynamic solution** и не certification.
Она не используется как numerical bound для G.

Раздельно установлены: работоспособность метода на checked p16 model;
наблюдаемый A/B contrast; отсутствие самостоятельной spatial/long temporal
qualification. Нет основания утверждать глобальную устойчивость или
сертифицированную границу применимости plane model. Перед первой научной
картой нужны ограниченные numerical controls и явно заданный input space;
саму карту этот этап не запускает.

## Figures и underlying data

Не более четырёх основных figures, PDF+PNG, без smoothing или phase fitting:

- [Plane motion A/B PNG](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/planar_A_B.png), [PDF](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/planar_A_B.pdf).
- [Mixed error PNG](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/mixed_tangent_full.png), [PDF](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/mixed_tangent_full.pdf).
- [Restricted gains PNG](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/restricted_gains.png), [PDF](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/restricted_gains.pdf).
- [Optimal full response PNG](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/optimal_full_response.png), [PDF](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/optimal_full_response.pdf).

[Manifest](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/manifest.json)
хеширует raw full trajectories/dense polynomials, norm/basis matrices,
propagator states, actual accepted times, initial/long/postprocessing source
snapshots и execution metadata. Source meshes/native FEM data не создаются.
`gain_histories.csv`, `mixed_error_histories.csv`, `finite_control_histories.csv`,
`optimal_initial_combinations.csv`, `actual_execution.csv` сохраняют numerical
data вместе с physical/normalized times. NPZ содержит полные поля и скорости.
Physical reconstruction из saved coefficients — postprocessing, без нового
ODE, equilibrium или eigenanalysis. Historical source reports также не изменены.

[NLSP-D21](../memory/decisions.md#nlsp-d21) и
[NLSP-K22](../memory/knowledge.md#nlsp-k22) сохраняют решение и bounded outcome.

## Границы интерпретации

Конечновременное усиление внутри выбранного подпространства не означает
асимптотическую неустойчивость, Floquet exponent или критическую амплитуду.
Изменение линейных свойств A/B отделяется от нелинейного плоского движения
контролем при нулевой base trajectory. Физическая точность идеального узла,
включая `c1=c2`, этим сравнением не устанавливается.

`DEFERRED_TWO_PARAMETER_ACTION_EXPANSION` сохранено: независимое разложение
исходного неусечённого действия по A и delta с членами delta²,A delta²,A² delta²
и оценкой более высоких членов ещё не выполнялось. Также отложены постоянные
геометрические несовершенства, полный amplification operator, настоящая
устойчивость периодических орбит, граница применимости и 3D angular joint.
FEM, карты, Floquet, periodic orbit и critical amplitude здесь не запускаются.

## Воспроизведение

Existing CLI использует отдельный scoped режим; нового runnable script нет:

```powershell
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
& 'D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe' scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --resume-sensitivity --compute
& 'D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe' scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --resume-sensitivity --report-only results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536
& 'D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe' scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --resume-sensitivity --plot-only results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536
```

Completed compute/report/plot читают cache и hashes, не выполняя новых IVPs.
Отдельная authorization предотвращает retry через смену code hash/HEAD.
Исторические CLI defaults и failed-attempt guards сохранены.

## Проверки и изменения checkout

75 targeted tests PASS (`7.10s`): новый scoped test file и прежние pilot,
torsion-diagnostic и coupled-kernel regressions. В тестах нет реальных ODE/FEM.
Проверены geometry/joint projection, фиксированный SPD norm, full-output gain,
sign algebra, dense-output sampling, физические оси, atomic preflight reservation
и cached dispatch. Старые scientific thresholds и historical tests не изменены.

```powershell
python -m pytest -q tests/test_nlsp_coupled_sensitivity_resume.py tests/test_nlsp_coupled_sensitivity_pilot.py tests/test_nlsp_coupled_torsion_diagnostic.py tests/test_weakly_nonlinear_coupled_spatial_dynamics.py
git diff --check
```

Actual cached compute/report/plot проверены с запрещёнными scientific functions:
ноль IVPs, equilibria и model reconstruction; manifest и counters не меняются.
Проверены13 local links этого отчёта, оба source bundles, old report hashes
и побитовая сохранность memory prefixes до D21/K22. `git diff --check` PASS.
Полное evidence находится в
[verification.json](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/verification.json).

Checkout остаётся `main`, HEAD `642758754250608f03b18e16ff1874bcce0d071b`,
staging empty. Было24 dirty paths, стало28:12modified+16untracked;
24 исходных пользовательских paths сохранены. Никаких commit/push/pull,
checkout/reset/clean/stash/merge/rebase нет. Новых файлов4; существующие
helpers изменены только для scoped dispatch/численной tangent propagation
и accepted-time archival. Theory/physical kernel files не изменены.

README и CHANGELOG обновлены вместе с navigation, journal и D21/K22.
[Список18 изменённых в этом continuation файлов](../../results/nlsp_coupled_planar_sensitivity_resume/7588614fdb79e536/changed_files.json)
отделён от24 исходных dirty paths. Основные additions — config, scoped helper,
этот report и targeted tests. Historical reports не редактировались.

После данного bounded результата программа остановлена. Параметрические
карты и следующий научный этап автоматически не запускаются.
