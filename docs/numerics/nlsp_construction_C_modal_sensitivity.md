# Модальное возбуждение и внеплоскостная чувствительность конструкции C

**Итог: `C_MODAL_SENSITIVITY_RESOLVED`.** На двух линейных reference
периодах получено `Gperp,disp^(3)=1.05302870`; physical v остаётся близким
к начальному масштабу. Ограниченные spatial/time controls разрешают этот
gain, и полный seven-field control подтверждает tangent prediction.
Spatial certification всех семи полей, в особенности c и velocities,
остаётся `PARTIAL`. Устойчивость на неограниченном времени не установлена.

Отдельный эксперимент 2026-10-11 использует принятую семиполевую модель
двух жёстко сопряжённых стержней. Проверяется конечновременный отклик малых
начальных внеплоскостных деформаций при одном плоском модальном возбуждении.
Это продолжение [A/B методики](nlsp_coupled_planar_sensitivity_resume.md),
без повторения её физических расчётов и без пересмотра исторических статусов.

## Постановка и сохранение источников

Угол сопряжения `beta=30deg`, длины `L1=.8,L2=1.2`. Оба сечения C имеют
`b=.05,h=.10`. Материал `E=rho=1,nu=.3,kappa=5/6`. Внешние концы закреплены
по значениям всех семи полей; производные дополнительно не закрепляются.
Идеальный точечный узел имеет прежние семь общих координат, включая `c1=c2`.
Дополнительных масс, пружин, вязкости и внешней нагрузки нет.

Порядок полей `(u,w,v,Phi,psi,theta,c)` и физические kernels неизменны.
Координаты обоих плеч направлены от внешней заделки к узлу. Используется
`q_global=T7_i q_local`, локальный вектор поворота `(Phi,-psi,theta)` и
`k=-EZ`. Толщина h направлена вдоль n, ширина b — вдоль k. Поэтому
`Bp=4.1666666667e-6`, `Bb=1.0416666667e-6`, `Bp/Bb=4` на обоих плечах.
Плоский изгиб жёстче внеплоскостного. Из этого отношения не следует
неустойчивость плоского движения.

Принятые V0, quartic action, cubic equations, переменная масса, нелинейная
инерция, RHS/Jacobian, базис и узловые условия не изменены. Их hashes
сверяются с execution snapshots исторического пилота. Parent A/B bundle
`7588614fdb79e536` имеет manifest SHA256
`338f25ad4ed37bc9369f567350e02aa7a23e70f945356a86baa1832629e3ed3e`.
Проверяются его артефакты и исходные bundles. Исторические
`BLOCKED_BY_TANGENT_CHECK`, `NUMERICAL_FLOOR_IDENTIFIED` и
`NUMERICAL_PARTIAL` сохраняются.

## Выбор формы до нелинейного решения

Исследованы только восемь низших плоских и двенадцать внеплоскостных
собственных форм constrained K0/M0. Это sorted positions при одной
геометрии, без branch tracking или спектральной карты. Выбрана первая
плоская форма, для которой доля сдвиговой и изгибной потенциальной энергии
больше половины. Доля вычисляется непосредственно из принятого линейного
предела: `sum_i integral[S(w_s-theta)^2+Bp theta_s^2] / (q.T K0 q)`.
Остальная часть включает продольное движение и поле c.

Для выбранной sorted form 1 получены `omega=.31338933639565253`,
`Tref=20.049135619749023` и изгибно-сдвиговая доля `.7238523367966677`.
Вся собственная форма, включая u,theta,c и узловые координаты, масштабируется
до `max_i,s|w|=.005`. Максимум определяется по стационарным точкам
представленных полиномов и концам. Максимальный положительный w находится
на втором плече при `s=.7130415338184252`. Внеплоскостные поля и все начальные
скорости равны нулю. Статической нагрузки и нового равновесия для C нет.

Начальное состояние проходит прежний safety gate: max rotation `.0101823 rad`,
max `Li*curvature=.0787653`, минимальное относительное собственное значение
assembled массы после resting-mass whitening `.99116347`. Начальные joint translation/orientation differences
не превышают `1.11e-16`; contraction difference равна нулю. Это проверка
допустимости представленного дискретного состояния, без доказательства
физической точности всего амплитудного усечения.

## Начальное пространство и показатели

Выбраны две низшие формы с torsional mass fraction менее .5 и первая
форма с долей не менее .5. Это sorted OUT forms 1,2,5, с частотами
`.07773966991,.21565348575,.84179835264`. Их крутильные доли равны
`.00069644,.00313728,.69735729`. Смешанный состав не устраняется.
Все три направления имеют нулевые начальные скорости и удовлетворяют
узловым условиям. Начальное пространство обозначается `S_disp^(3)`.

Фиксированная reference norm равна
`||z||E0^2=eta.T Kperp,0 eta+eta_dot.T Mperp,0 eta_dot`.
Матрицы положительно определены. Basis Gram error `1.26635e-10` меньше
прежнего `2e-10` contract. При p16 выход сохраняет все93 OUT coordinates
и93 velocities. Restriction применяется только к начальным состояниям.

`Gperp,disp^(3)(t)` — корень наибольшего собственного значения Gram matrix
полных выходных состояний трёх seeds. Максимум по времени оценивается на
2001 отсчёте `0...2Tref`, шаг вывода `.001Tref`. Это sampled maximum для
выбранного пространства, не полный amplification operator или stability
certificate. Fixed optimal direction выбирается в момент глобального
sampled maximum; его отклик отличается от переоптимизированной в каждый
момент gain curve.

Для физического отклика один общий множитель задаёт
`max{|v(0)|,Lref|Phi(0)|,Lref|psi(0)|}=.01*Apl=5e-5`, `Lref=1`.
Перемещения выражены в единицах reference length, повороты — в rad.
Все компоненты масштабируются совместно; начальные скорости остаются
нулевыми. `Rv=max_t,i,s|v|/.005` не является критерием устойчивости.
Zero-base trajectory с тем же входом отделяет линейное перераспределение
координат и скоростей от дополнительного влияния плоского движения.

Отношение первой OUT frequency к omega_pl равно `.248061`, четвёртой —
`2.01436`. Эти близости к 1:4 и 2:1 служат только частотной диагностикой.
Механизм параметрического резонанса здесь не установлен.

## Численная программа и воспроизводимость

Математический касательный оператор — OUT-блок полного analytic Jacobian
на принятой плоской траектории. Новый оператор не выводится. Base integration
использует существующий Radau, `rtol=atol_relative=1e-10` и прежнюю
componentwise policy. Касательная ODE использует DOP853 с тем же Jacobian,
`rtol=atol_relative=1e-11`. Max step связан с прежним section cutoff/24;
для C cap около `.00720939`. Dense output требуется только у output samples.
Скорости являются сохранёнными state variables, а не производными графиков.

Один temporal control использует найденный fixed optimal input при
вдесятеро меньших rtol/atol и вдвое меньшем max_step. Это контроль одного
отклика, без нового tight restricted optimization и без отдельного temporal
уточнения base trajectory. Spatial control p12 использует физическую
проекцию того же p16 начального состояния и входного пространства.
Семь общих узловых координат сохраняются. Собственный Tref и возбуждение
для p12 заново не подбираются. Поля сравниваются на общей физической
координате; coefficient vectors разных пространств напрямую не вычитаются.

До ODE сохранён `frozen_resource_plan.json`: прогноз полной исходной
программы `7178.504s`, обязательной части с controls `5962.815s`, budget7200s.
Это extrapolation по actual A/B counters и коротким pointwise RHS/Jacobian
timings, без дополнительной ODE. Используется один поток, memory ceiling4GiB.
Максимумы: две planar trajectories, восемь tangent integrations и один
дополнительный полный nonlinear finite control. После любого реального
отказа нет автоматического retry.

Для полного finite-control расчёта все native accepted coordinates/velocities
записываются последовательно в `accepted_states.bin`. Описание dtype, shape
и state ordering находится в `accepted_state_layout.json`. Сохраняются
все native timestamps, actual final state и2001 full comparison states.
Прогноз архивирования всех cubic dense coefficients с дублированием массивов
даёт5.0835GiB и не соответствует memory ceiling4GiB. Новый output option
shared Radau helper оставляет прежнее dense архивирование по умолчанию.
Метод, допуски, шаги и comparison samples сохраняются; mock regression
подтверждает совпадение sampled states двух вариантов хранения.

Продление до5T не выбрано. В последние .2Tref не появилось нового глобального
maximum gain; первый maximum находится при .745Tref. После numerical
controls оставалось3730.109s. Прогноз дополнительного base/tangent prefix
равен2512.003s, резерв единственного full control на5T —6192.960s.
Эта сумма не укладывается в остаток. Ограничение вывода двумя reference
периодами сохраняется; отсутствие большого роста не переносится на
неограниченное время.

До новых gain results также задан десятикратный planning guide для
отношения `G-1` к наблюдаемым spatial/time changes и прежнему numerical floor.
Эти изменения не являются строгими верхними оценками ошибки. Сам guide
не доказывает устойчивость или физическую применимость плоской модели.

По Script Proliferation Control новый runnable script не создаётся.
Существующий CLI получает preset `--construction-c`; scoped helper
переиспользует прежние full/tangent IVPs, assembly, projection и gain recovery.
Авторизация независима от A/B. Изменение code hash не разрешает повторить
неудачный или незавершённый attempt.

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
& 'D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe' scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --construction-c --compute
```

Конфигурация:
[nlsp_construction_C_modal_sensitivity.json](../../data/input/nlsp_construction_C_modal_sensitivity.json).
Execution bundle:
`results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/`.
Cache compute/report/plot читает неизменяемые результаты без новых ODE,
static Newton, eigenanalysis, Gmsh или CalculiX calls.

## Полученные restricted gains и пространственная чувствительность

| Quantity | p16 | p12, same projected input space |
|---|---:|---:|
| Sampled max Gperp,disp^(3) | 1.05302869908747 | 1.05303167743216 |
| Time at max / Tref | .745 | .746 |
| Physical time at max | 14.9366060367130 | 14.9566551723328 |
| Zero-base Gram error | 1.41081e-10 | 7.84373e-12 |
| Max physical v | 3.1732912151e-5 | 3.1732215326e-5 |
| Rv = max abs v / .005 | .00634658243 | .00634644307 |
| Max abs Phi, rad | 3.7088046738e-5 | 3.7100578558e-5 |
| Max abs psi, rad | 5.1121939290e-5 | 5.1046276218e-5 |

В каждой колонке physical response относится к её fixed optimal direction,
масштабированному одним заранее объявленным правилом. При p16 коэффициенты
в E0-orthonormal basis равны
`(.9969234109,-.01354694065,.07720228792)`: направление преимущественно
содержит первую внеплоскостную изгибную форму. Это коэффициенты в выбранном
initial space, не Евклидова норма произвольных Galerkin coefficients.

Для этого p16 направления `delta_E0=1.5311445173942503e-7`. Начальные
максимумы `(v,Lref*Phi,Lref*psi)` равны
`(3.1326767764e-5,2.7719865435e-5,5e-5)`; все velocities равны нулю.
Для zero-base с точно тем же входом max v=`3.1326570844e-5`,
max Phi=`3.5450811752e-5`, max psi=`4.9999864527e-5`.
Поэтому весь рост отдельных rotational amplitudes нельзя приписывать
нелинейной связи: часть присутствует уже в линейном движении.

Энергетическая норма увеличивается на5.30287%, но physical v остаётся
вблизи начального масштаба. Его максимум составляет0.634658% Apl.
Начальное условие1% относится к совместному максимуму перемещений и
Lref-scaled поворотов; в данном направлении его задаёт psi, не v.

| Observed numerical change | Value | Qualification |
|---|---:|---|
| Difference of max gains p12/p16 | 2.9783446871e-6 | sampled maxima at adjacent output times |
| Max difference of complete gain curves | 9.4576879540e-4 | 1.7835% of measured G−1 |
| (G−1) / max gain-curve change | 56.0694 | exceeds frozen planning factor10 |
| Same-input physical v difference | 5.3391725113e-9 | no phase/amplitude fitting |
| Same-input Phi difference | 1.5044629518e-7 rad | all output modes retained |
| Same-input psi difference | 2.0443536670e-7 rad | all output modes retained |
| Main planar w difference | 3.6829782475e-6 | physical common coordinates |
| Main planar u difference | 3.1664430118e-6 | physical common coordinates |
| Main planar theta difference | 3.9349931482e-5 rad | observed change, not strict error bound |
| Main planar c difference | 2.7112562945e-5 | spatial qualification PARTIAL |
| Main c-velocity difference | 1.1399260527e-3 | spatial qualification PARTIAL |

Близость max gains не означает побитового совпадения всей gain curve:
полная разность достигает9.46e-4. Она остаётся в56 раз меньше измеренного
энергетического усиления. Поле c и некоторые velocities чувствительнее к
пространственному уровню. Поэтому resolved gain относится к данной
наблюдаемой величине и выбранному подпространству; полная пространственная
сертификация семи полей не установлена.

## Полная семиполевая проверка и фактическое исполнение

Один full Radau control использует ровно тот же p16 modal initial state
и найденное fixed optimal OUT direction с `delta_E0=1.5311445173942503e-7`.
Нагрузка отсутствует; начальные скорости всех полей равны нулю. Применена
заранее объявленная policy `rtol=1e-12,atol=base*.01,max_step=base*.5`.
Ни amplitude,ни phase,ни модель не менялись для улучшения согласия.

| Full/tangent quantity | Actual value |
|---|---:|
| Max normalized E0 error | 1.6269796545e-7 |
| Relative error, fixed denominator max tangent norm | 1.5450477806e-7 |
| Existing relative gate | 1e-5: PASS |
| Coordinate K error max | 1.3619398531e-7 |
| Velocity M error max | 1.3516675450e-7 |
| Velocity fraction at total error max | 40.2211% |
| Max physical v error | 7.5133773324e-13 |
| Max physical Phi error | 3.9144461983e-12 rad |
| Max physical psi error | 4.5744170636e-12 rad |
| Full max physical v | 3.1732911416e-5 |
| Full max Phi / psi | 3.7088046494e-5 / 5.1121938056e-5 rad |
| Full mechanical energy relative drift | 8.3187240757e-14 |
| Observed planar-change E0 norm | 1.2695071097e-10 |

Последняя величина включает finite perturbation feedback и различие
временной точности base/full integrations. Эти два вклада отдельно не
разделены. Один finite control подтверждает прогноз на данном delta;
квадратичный asymptotic order заново серией delta для C не измерялся.

В reference E0 norm на velocities приходится83.4693% квадрата norm
найденного отклика в момент global gain maximum. Поэтому gain5.3% нельзя
называть увеличением максимального прогиба на5.3%. Здесь max v увеличился
лишь примерно на1.30% относительно самого большого zero-base v с тем же
входом. Фиксированная E0 perturbation norm также не тождественна отдельной
сохраняющейся nonlinear mechanical energy.

| Actual IVP | Seconds | Accepted steps | Response/input dimensions |
|---|---:|---:|---|
| C p16 planar base | 712.140 | 54275 | 248-state invariant plane subsystem |
| C p16 tangent bank | 962.528 | 9286 | 186 OUT states,3 inputs |
| C p16 zero-base bank | 2.866 | 5562 | 186 OUT states,3 inputs |
| C p16 fixed optimal tight | 935.605 | 11124 | 186 OUT states,1 input |
| C p12 planar base | 438.984 | 38030 | 184-state invariant plane subsystem |
| C p12 tangent bank | 415.753 | 6028 | 138 OUT states,3 inputs |
| C p12 zero-base bank | 2.014 | 5562 | 138 OUT states,3 inputs |
| C p16 full finite control | 2081.951 | 159868 | full434 states |

Actual total `5551.8418467s/7200s`:2 planar trajectories,1 full nonlinear
control и5 matrix tangent IVPs. Пять linear matrix solves сохраняют13
column histories:3+3+1+3+3; это не13 отдельно запущенных интеграторов.
Все jobs последовательны, one thread; static Newton/FEM/Gmsh/CCX calls0.
Все8 IVPs достигли `40.09827123949805` без failure/retry. Full accepted
step sizes лежат между `.0001797049544` и `.0004593957834`.

Scientific engine lifetime peak working sets: base1.103024GiB,
gain.788689GiB, controls.788631GiB, full.337257GiB. Controls phase содержит
три разных IVPs; её peak нельзя приписывать отдельно каждому case.
Venv launcher не является измеряемым scientific engine. Native full archive
содержит159869 actual states, каждая из434 coordinates/velocities;
native times и final state сохранены отдельно. После завершения массивы
проверены на finite values, coverage и соответствие exact initial/final state.

Observed p12/p16 changes нормируются отдельно на sampled p16 physical
maximum: w0.07363%,u0.59151%,theta0.38472%,c16.77874%. Для velocities
соответствующие изменения w3.01612%,u19.12749%,theta17.50842%,c153.36253%.
В частности, `max|c_p16|=.00016158874724` на declared physical grid.
Большая разность c velocities сохраняется; она не скрывается выбором общего
масштаба .005. Поэтому все-field spatial qualification остаётся PARTIAL,
даже при устойчивом restricted gain и малом energy drift.

## Рисунки, данные и cache

Не более четырёх figures, каждая PNG+PDF, без smoothing и fitting:

1. [Complete planar initial mode](../../results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/C_modal_initial_fields.png)
   ([PDF](../../results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/C_modal_initial_fields.pdf)).
2. [Base and tangent OUT motion](../../results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/C_planar_and_OUT_motion.png)
   ([PDF](../../results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/C_planar_and_OUT_motion.pdf)).
3. [Restricted gain and same-input zero control](../../results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/C_restricted_gain_and_zero_control.png)
   ([PDF](../../results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/C_restricted_gain_and_zero_control.pdf)).
4. [Initial OUT shape, energy and full comparison](../../results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/C_optimal_form_and_controls.png)
   ([PDF](../../results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/C_optimal_form_and_controls.pdf)).

Основные CSV: `actual_execution.csv`, `C_p16_gain_observations.csv`,
`C_p12_gain_observations.csv`. NPZ хранит original coefficients, все
comparison states, physical profiles/velocities, complete tangent vectors,
fixed matrices и temporal/spatial comparison arrays. Config, initial/code
snapshots, individual attempts/settings, actual environment, resource decisions,
memory evidence и manifest позволяют проследить каждый результат.

Completed compute/report-only/plot-only проверены с instrumented forbidden
scientific functions:0 новых ODE/static/eigen/FEM calls, manifest unchanged.
Targeted C и предшествующие coupled regressions:99 PASS; native coverage,
source/user-file preservation, navigation и `git diff --check` проверены.
Git branch/main и HEAD6427587 сохранены; staging пуст. README/CHANGELOG,
navigation/journal и [D22](../memory/decisions.md#nlsp-d22)/
[K23](../memory/knowledge.md#nlsp-k23) обновлены без переписывания старых records.

## Временной контроль и границы вывода

В завершённом p16 temporal control максимальное изменение полного
нормированного OUT state равно `1.9335025138e-9` в E0 norm. Отдельные
максимумы координатной и скоростной частей равны `1.8455113332e-9` и
`1.8698829205e-9`. В момент максимума общей разности на скорости приходится
58.19% её квадрата. Максимум нормы fixed optimal response меняется
`1.0530286990872986 → 1.0530286990878033`. Этот контроль подтверждает
устойчивость важнейшего tangent response по времени на общей base trajectory;
отдельный temporal refinement основной плоской динамики не выполнен.

C использует другой протокол возбуждения и другое начальное пространство,
чем исторические A/B. Поэтому значения их restricted gains сами по себе
не изолируют влияние ориентации сечения. Основной контроль C — его собственное
zero-base движение в той же фиксированной норме и с тем же initial input.

Конечновременное усиление не называется математической неустойчивостью.
Контроль одного displacement-only пространства не охватывает все малые
возмущения, скоростные IC, постоянные несовершенства или бесконечномерный
предел. Собственные периоды Tref не объявляются периодами нелинейной орбиты.
Новые геометрии, амплитуды, FEM, maps, Floquet и periodic-orbit search отсутствуют.

`DEFERRED_TWO_PARAMETER_ACTION_EXPANSION`, физическое обоснование
идеального узла и `c1=c2`, независимая 3D проверка углового соединения
и поиск границы применимости сохраняются отложенными. Исторические
FEM/strict-float64 qualifications, LONG CLOSED и EB/RLB-KV PAUSED не меняются.

