# Текущий контекст

## Full-family linear3D FEM-1 — bounded complete / mesh PARTIAL, 2026-10-08

[NLSP-D09](decisions.md#nlsp-d09), [NLSP-K10](knowledge.md#nlsp-k10),
[report](../numerics/nlsp_linear_rectangular_3d_fem_validation.md).
L=1,b=.20,h=.10 выбран до FEM по1D спектру; first axial на позиции8. Все4 families
сопоставлены,3 audited C3D10 meshes/24 complete vectors each; quality/execution/
identification PASS. Mesh convergence и quantitative all-family comparison
PARTIAL:5/8 medium/fine изменений выше заранее принятого.1% критерия.
Observed discrepancies не являются exact continuum/model errors. No extra3D
mode inside window; contraction branch outside. No fitting/четвёртая сетка/
nonlinear calculation. После отчёта остановка, FEM-2/FEM-3 не разрешены.
LONG CLOSED; EB/RLB-KV PAUSED; angular same-clamp UNAVAILABLE; prepared strict
PARTIAL и physical sanity DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS сохраняются.

## 3D FEM infrastructure audit — завершён, 2026-10-08

[NLSP-K09](knowledge.md#nlsp-k09),
[readiness report](../numerics/nlsp_3d_fem_environment_readiness.md).
Gmsh4.15.2 и cached CalculiX2.22 x64 запускаются; linear solid workflow ранее
использовался. Локальная документация подтверждает nonlinear STATIC/DYNAMIC
с NLGEOM, но готовых project workflows для них нет. Для rectangular G20 нужна
небольшая geometry/config адаптация и отдельный execution test. Опциональные
Python FEM packages не установлены и не являются blocker существующего CLI path.

Этот этап — только readiness audit: ни mesh, ни FEM job, ни installation,
ни physical/model/result change. Следующий bounded linear rectangular test
предложен, но не разрешён к запуску. Физическая проверка V0 ещё не выполнена.
NLSP physical sanity остаётся DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS;
LONG CLOSED, EB/RLB-KV PAUSED и angular same-clamp UNAVAILABLE сохраняются.

## Planar physical sanity — diagnostic complete with qualifications, 2026-10-08

[NLSP-D08](decisions.md#nlsp-d08), [NLSP-K08](knowledge.md#nlsp-k08),
[canonical note](../theory/weakly_nonlinear_planar_physical_sanity_checks.md).
Сохранённый p64 large run и единственный новый p64 tight half-amplitude run
достигли общего T1. Наблюдаемые maxima близки 2/4 по ожидаемым field powers;
normalized leading deviations уменьшаются примерно в 4 раза. Symmetry,
formal classical stretching limit, strain/mass/energy/safety checks проходят.

Overall `DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS`: second-order/amplitude
checks PARTIAL из-за отсутствия независимого half-amplitude p/time control;
reaction/momentum PARTIAL сохраняют finite-p residual и explicit rotational
truncation contribution quartic action. Это physical sanity принятого reduced
V0, не experimental/3D validation. Prepared strict certification PARTIAL,
`EXPLORATORY_NOT_CERTIFIED` и `state.admitted=False` остаются; старые zero-u/c
и full nonlinear spatial PARTIAL не повышаются.

Model/action/RHS/Jacobian, four fields, basis, BC и coefficients сохранены;
новый half initial state использует те же profiles/Theta3 с разными amplitude
powers и прежним projection rule. После bounded report остановка, новый
scientific direction не выбран. LONG CLOSED; EB/RLB-KV
PAUSED_FOR_SUPERVISOR_DIRECTION; angular same-clamp out-of-plane reference
UNAVAILABLE. Исторические решения и остановки ниже сохраняются отдельно.


## Исторический prepared one-T1 feasibility — exploratory, spatial PARTIAL, 2026-10-08

[NLSP-D07](decisions.md#nlsp-d07), [NLSP-K07](knowledge.md#nlsp-k07),
[canonical continuation](../theory/planar_prepared_initial_state.md#prepared-one-period-feasibility).
Тот же frozen prepared IVP рассчитан всеми тремя cases p48 tight, p64 tight,
p64 allowed_extra до одного линейного reference T1. Source q0/v0/settings
загружены без projection/MP/BVP/eigen/Theta3 пересчёта. Прежний short prefix,
energy/mass/safety и все 8 temporal component gates проходят.

Spatial p48→p64 остаётся PARTIAL: 7 из 8 проходят, theta_t max relative
около 1.59015e-4 выше 1e-4 при прошедшем L2 gate. Практическая finite-dimensional
feasibility установлена на 0…T1; strict initial PARTIAL и
`EXPLORATORY_NOT_CERTIFIED`, `state.admitted=False` сохранены. Это не full 5T1
или continuous-PDE validation, periodic orbit либо повышение старой zero-u/c
задачи до PASS. V0/model/RHS/Jacobian/basis/BC/coefficients не менялись.

После bounded one-T1 report остановка; дальнейший горизонт, refinement или
новый scientific direction отдельно не выбраны. LONG CLOSED, EB/RLB-KV
PAUSED_FOR_SUPERVISOR_DIRECTION, angular same-clamp out-of-plane reference
UNAVAILABLE. Исторические разделы ниже сохраняют прежние результаты и остановки.


## Исторический prepared-state precision/feasibility — exploratory, spatial PARTIAL, 2026-10-08

[NLSP-D06](decisions.md#nlsp-d06), [NLSP-K06](knowledge.md#nlsp-k06),
[canonical continuation](../theory/planar_prepared_initial_state.md#prepared-precision-feasibility).
Сохранён один общий frozen initial evaluator; новая initial-only
`common_endpoint_constrained_L2` representation проходит критерии 1e-6 для всех полей,
их derivatives 0…2 и endpoint jets при p48/p64. Basis, будущие BC, physics,
V0/mass/RHS, coefficients и сохранённая Theta3 не менялись.

Strict float-Gauss relative action2e-12 остаётся unresolved; independent
MP quadrature evidence локализует numerical storage/moment sensitivity.
`state.admitted=False`; выполненный маршрут явно
`COMPLETED_EXPLORATORY_NOT_CERTIFIED`, а не strict verified PASS.
Ровно 3 short controls достигли 0.1T1: temporal 8/8 PASS, spatial 7/8 PASS,
но theta_t max≈1.4053e-4 превышает1e-4; spatial PARTIAL. Energy/mass/safety PASS.
Короткая finite-dimensional feasibility установлена; full 5T1 и continuous-PDE
validation не заявлены. Старые zero-u/c PARTIAL и D05/K05 history сохранены.

После bounded task остановка; новое refinement/scientific direction не выбрано.
LONG CLOSED; EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION; angular same-clamp
out-of-plane reference UNAVAILABLE. Исторические разделы ниже описывают прежние
задания и не снимаются задним числом новой явно разрешённой работой.


## Исторический prepared initial-state gate — PARTIAL, 2026-10-08

[NLSP-D05](decisions.md#nlsp-d05), [NLSP-K05](knowledge.md#nlsp-k05),
[canonical report](../theory/planar_prepared_initial_state.md).
Stat/harm profiles и derivative jets сходятся на64→96. Общий отдельный
initial state сO2 u,c и одной quinticTheta3 приw3=0 имеет проверенную
acceleration compatibility through cubic order; конечные более высокие
остатки сохраняются, periodic orbit не заявлена.

Полный initial projection PARTIAL: ни32/48, ни48/64 не проходят обеими
сторонамиpreset1e-6. O2 u/cp64PASS не равен four-field projection PASS;
p64 сохраняет theta_ss endpoint error. Дополнительный relative strong/weak
numerical check приp48/p64 также unresolved, fixedgate2e-12 не ослаблен.
New short temporal/spatial NOT_RUN,0 ODE; improvement nonlinear convergence
не установлен. Старые nonlinear PARTIAL и zero-u/c task сохранены.
Следующая numerical/scientific strategy не выбрана, автоматического
p64/96, full5T1, заменыbasis/BC/V0, angular/Floquet исследования нет.
LONG CLOSED; EB/RLB-KV PAUSED; angular same-clamp reference UNAVAILABLE.

## Leading second-order axial diagnostic — COMPLETE / spatial PARTIAL, 2026-10-08

[NLSP-D04](decisions.md#nlsp-d04), [NLSP-K04](knowledge.md#nlsp-k04),
[canonical report](../theory/planar_second_order_axial_response.md).
Ведущий u2,c2 response принятой cubic модели вычислен аналитически по времени
с common continuous bending background и всеми Shen coordinates p16…96.
Derivation/forcing/exact-time checks PASS; диагностическая программа COMPLETE,
но spatial convergence PARTIAL: в64→96 толькоu2 проходит обе1e−3 нормы,
c2 и обе speeds — нет. Sampling refinement и independent expm прошли.
Трудность старогоp24→32 уже присутствует во втором порядке без ODE error;
initialA² boundary mismatch сохраняется, единственная причинность не доказана.
Старые nonlinear comparisons используют только actual saved timestamps;
second-order response не exact full cubic trajectory и не новая production theory.
Полный nonlinear p48 не считался; initial fields/V0/RHS/BC/basis сохранены.
Новых time integrations0; после разрешённогоp96 остановка, p128+ не выбран.



## Planar numerical recovery — PARTIAL, 2026-10-07

[NLSP-D03](decisions.md#nlsp-d03), [NLSP-K03](knowledge.md#nlsp-k03),
[continuation report](../theory/weakly_nonlinear_planar_time_pilot.md#14-адресная-диагностика-и-восстановление-вычислительного-пути-2026-10-07).
Исходные trajectories/criteria и cubic/quartic model сохранены. Physical L2
projection подтверждает и tail, и эволюционную разность общих компонент;
ненулевой initial axial acceleration trace порядка A² подтверждён как
ограничение гладкой совместности, не code/model failure или единственная причина.
Lazy energy/gradient/Hessian optimization эквивалентна исходному RHS/Jacobian;
короткий old/new p32 совпадает, ускорение integration≈1.26×. Ровно3 short
controls выполнены, полный p48 отложен по неизменному900s cost rule:
REFINEMENT_DEFERRED_BY_BUDGET. Full spatial four-field convergence и малый
neighboring-p control unresolved; NLSP_PLANAR_SOLVER_RECOVERY=PARTIAL.
Short p48 до0.1T1 не заменяет full5T1 evidence. Следующее numerical/physical
решение не выбрано; basis/IC/V0/BC не заменялись. История NLSP-D02/K02 ниже
сохранена, её прежняя остановка не является запретом уже выполненного D03.


## Нелинейная пространственная ветвь — первый planar time pilot PARTIAL, 2026-10-07

[NLSP-D02](decisions.md#nlsp-d02), [NLSP-K02](knowledge.md#nlsp-k02),
[численный отчёт](../theory/weakly_nonlinear_planar_time_pilot.md).
Рассчитаны реальные four-field trajectories принятой cubic model:
один прямой G20 fixed-fixed rod, u,w,theta,c независимы, две amplitudes
A/h=.05/.025. Обе finalp32 trajectories доходят до5T1. Linear/action и
time controls PASS, energy/mass safe; spatial convergence u,c и части
velocities не достигнута. Последний neighboring-p small-amplitude case
сохранён до2.648T1 по заранее заданному budget, без ослабления tolerances.
Общий NLSP_PLANAR_TIME_PILOT=PARTIAL; полный nonlinear effect acceptance
не объявлен. Вопросы periodic/Floquet/out-of-plane stability/threshold
не рассчитывались. Новый stage/refinement автоматически не разрешён.

## Предыдущий семиполевой model/action audit — сохранён

Новое явное задание пользователя после LONG-D02 выбрало ограниченную
фиксацию/верификацию семиполевой кинематики и принятого V0:
[NLSP-D01](decisions.md#nlsp-d01), [NLSP-K01](knowledge.md#nlsp-k01),
[canonical note](../theory/weakly_nonlinear_spatial_rod.md).
Это нелинейное пространственное продолжение, не новый Bishop/RL/M-H выбор.
Историческое «следующее направление не выбрано» ниже относится к прежней
memory-sync задаче и не переписывается задним числом в D/K.

Завершены независимая quartic-action/cubic-balance сверка, supplied MD,
boundary/energy/reflection/linear/planar/axial/split controls и проверка
амплитудного порядка на manufactured jets. V0 — наше принятое nonlinear
reduced constitutive assumption, не source-attributed nonlinear Jang/Yartsev
и не полная 3D редукция. c1=c2 сохраняет reduced common-DOF qualification.
На этом audit остановка: nonlinear time evolution, периодические орбиты,
Galerkin/Floquet и критическая амплитуда ещё не рассчитаны/не разрешены
автоматически. LONG закрыт; EB/RLB-KV остаётся PAUSED_FOR_SUPERVISOR_DIRECTION;
упругость/вязкость там уже существуют и сюда не переносились.

## Продольная модель сопряжённых стержней — закрыта, 2026-10-07

[LONG-D01](decisions.md#long-d01) — восстановленная история production выбора;
[LONG-K01](knowledge.md#long-k01) — завершённая цепочка проверок;
[LONG-D02](decisions.md#long-d02) — явное решение пользователя о закрытии.
Основание синхронизации: main, HEAD
`05ccf101fdc7adc6f6a1e87143c5987b0efe351a`, Version 0.5.14.

Для rectangular isotropic in-plane линии выбрана Jang-type reduced
Mindlin–Herrmann axial/contraction + Timoshenko bending, preset
`project_jang_reduced_rectangular`. κ=5/6 — отдельное project значение;
численный κ source Jang не восстановлен и не приписывается публикации.

Finite production rod, beta0 transparency и general-beta rigid-joint gates
прошли; hierarchy/thickness screenings завершены; большие Lambda(beta)
implementation checks прошли в объявленных geometry/grid ranges.
`HIERARCHY_SINGLE_ROD=PARTIAL_PASS` сохраняет qualification: elementary/RL
не имеют independent c и того же resolved c-clamp. Это не отменяет
`MHTIM_SINGLE_ROD_FINITE_SPECTRUM=PASS`.

Closure d1=d2, theta1=theta2, c1=c2 с dual force/moment/R balances —
вариационное common-DOF замыкание reduced 1D frame, поддержанное принятой
published assembly structure, **не прямой вывод из 3D elasticity конечной
области сварного узла**. M-H — comparison reference, не exact truth;
1D проверки не являются 2D/3D validation.

**LONGITUDINAL MODEL QUESTION CLOSED IN THE ADOPTED 1D SCOPE.**
Дополнительные Bishop/RL/M-H comparisons, открытое продолжение толщины,
beta refinement/maps, close-pair analysis и новые implementation checks
не являются текущей задачей. Пользователь решил не исследовать close pair;
исторический кандидат из thickness report не становится следующим этапом.
Возобновление требует отдельного явного решения по новому физическому
вопросу или validation target. Следующее научное направление этим заданием
**не выбрано**.

## EB/RLB rotational spring / Kelvin–Voigt — отдельная ветвь, на паузе

**PAUSED_FOR_SUPERVISOR_DIRECTION**, 2026-09-27:
[RLB-D22](decisions.md#rlb-d22), [RLB-K23](knowledge.md#rlb-k23),
[сводка для руководителя](../laminated_beams/inplane_kelvin_voigt_research_status_for_supervisor.md).
Поворотная упругость и Kelvin–Voigt вязкость **уже введены и исследованы**
в EB/RLB in-plane models. Здесь существуют elastic spring theory,
EB implementation/maps/forms, RLB → EB limit, KV theory, literature
benchmarks, elastic damping-participation screening, complex weak-damping
confirmation и production EB/RLB solver routing ([K05–K23](README.md#темы)).

Это другая model branch. «Впервые добавить упругость и вязкость узла»
не является новым следующим шагом. Перенос/generalization spring/KV
на production M-H/Timoshenko не объявлен выполненным или выбранным:
если потребуется, это отдельное будущее решение. Исторические статусы,
source/Ritz qualifications и пауза D22/K23 сохраняются.

Заключительное сравнение пяти EB/RLB-пар подтверждает прогноз отношения
затуханий в своих исходных границах; новых расчётов по нему не требуется.
Подготовлены [Figure 3](../laminated_beams/figure03_crossing_veering.md),
[Figure 4 revised с сохранённой v1](../laminated_beams/figure04_eb_rlb_shapes.md),
[Figure 5 — огибающие](../laminated_beams/figure05_damping_envelopes.md);
[прежняя столбчатая версия](../laminated_beams/figure05_frequency_damping_comparison.md)
остаётся историей. Новая физическая программа KV не выбрана.

## Circular EB reviewer diagnostic — завершён, 2026-10-05

[D25](decisions.md#rlb-d25) / [K26](knowledge.md#rlb-k26):
[единая научная заметка](../laminated_beams/circular_eb_rotational_spring_rigid_limit.md).
Одна геометрия круглых изотропных EB-стержней: mu=.30, beta=15°,
r=.005 м, l=1 м. Это отдельное продолжение темы поворотного узла,
без переноса результата на RLB, анизотропию, Bishop или вязкость.

Exact RIGID согласован с baseline для первых шести корней и guard 7
(max|Delta Lambda|=3.016e-10). Generic solver разрешил событие около
204.287 Гц как один простой корень позиции 5. Sorted spectrum растёт
1→10→100→RIGID, отличие от RIGID уменьшается для всех шести позиций.

28 основных форм восстановлены с M=1. Ограниченное продолжение по
kappa подтвердило все шесть локальных seeds до RIGID:
1→3.162278→5.623413→10→17.782794→31.622777→100→RIGID.
Min accepted MAC=.968142, margin=.936732; пороги неизменны,
endpoint conflicts и смены sorted position отсутствуют.

### Сохранённые qualifications circular EB diagnostic

- Полная 12-root эквивалентность не установлена: отдельная численная
  SVD-nullity classification на позициях 11–12 около Lambda=18.139529.
- Direct 1→10 для 02/03/05/06 остаётся UNRESOLVED. Успешное
  continuation — отдельный результат, а не переименование direct.
- |Delta psi| и s уменьшаются вдоль проверенного пути у 5/6.
  Seed 06 от 1 до 10 сначала увеличивает их на 2.95634% и 4.53756%.
  В exact RIGID относительный поворот проходит physical gate; s не задана.
- Старые [EB-JOINT-K01](knowledge.md#eb-joint-k01) /
  [K02](knowledge.md#eb-joint-k02) сохраняют остановки прежнего workflow,
  без доказательства физической кратности.
- Локальные raw results игнорируются Git; смысл и ограничения сохраняет
  научная заметка. Ни реальный k_theta, ни C, ни finite-r/l ошибка не найдены.

### Историческая остановка circular EB diagnostic

1D spring/RIGID reviewer diagnostic завершён. D25 сохраняет исторический
выбор темы transmission/equilibrium conditions и границ thin-joint
asymptotics. Это не назначение следующего общего этапа настоящей
синхронизацией. В заметке отдельно записаны литературные
асимптотические условия, баланс на области узла, условный compact-joint
scaling и конечный 1D численный результат.

Локальная FEM-калибровка реального узла отложена пользователем.
Новые корни, kappa, геометрии, r/l sweep и исправление high-spectrum
classification автоматически не разрешены. Текущий этап — документационная
синхронизация без расчётов или изменения моделей.

## Историческая остановка после memory-sync

Продольная линия закрыта в принятом 1D scope; EB/RLB-KV остаётся
`PAUSED_FOR_SUPERVISOR_DIRECTION`; circular EB diagnostic завершён
со своими qualifications. **Следующее научное направление не выбрано
этим заданием.** Автоматическое продолжение любой из этих линий не разрешено.

## Историческая остановка после NLSP audit

Принятая семиполевая модель зафиксирована и проверена в объявленном
математическом/вычислительном scope. Это не PHYSICAL_NONLINEAR_VALIDATION_PASS.
Следующий этап после аудита не выбран автоматически: ни движение, ни порог
выхода из плоскости не рассчитывались. LONG остаётся CLOSED; EB/RLB-KV paused.

## Текущая остановка после planar time pilot

Bounded программа завершена с PARTIAL и сохранёнными trajectories/figures.
Spatial four-field convergence и полный малый-amplitude neighboring-p
контроль unresolved. Дальнейшие вычисления не запускаются автоматически.
LONGITUDINAL MODEL QUESTION CLOSED IN THE ADOPTED 1D SCOPE сохраняется;
EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION сохраняется. Упругость/вязкость
уже реализованы в прежней ветви, не перенесены сюда. Angular same-clamp
out-of-plane reference остаётся UNAVAILABLE; прямой planar pilot его не
закрывает. Следующее научное направление после этого pilot не выбрано.


## Текущая остановка после адресной диагностики

NLSP planar solver recovery PARTIAL; initial low-order smoothness mismatch
подтверждён, умеренное ускорение получено, full p48 DEFERRED_BY_BUDGET.
Небольшое изменение реализации не установило полный four-field PASS.
Новые интегрирования или смена задачи автоматически не разрешены.
LONG — CLOSED; EB/RLB-KV — PAUSED_FOR_SUPERVISOR_DIRECTION; angular out-of-plane
same-clamp reference — UNAVAILABLE. Следующее направление отдельно не выбрано.


## Историческая остановка после exact-time leading контроля

NLSP second-order diagnostic COMPLETE / spatial PARTIAL; full nonlinear spatial
convergence unresolved. Дальнейшая замена basis, numerical scheme или initial
data не выбрана; новые nonlinear integrations/angles/Floquet/thresholds не разрешены
автоматически. LONG CLOSED; EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION;
angular same-clamp out-of-plane reference UNAVAILABLE. Предыдущие PARTIAL
и исторические остановки выше сохранены, следующий этап не выбран.


## Историческая остановка после prepared initial-state gate

Подготовка общего numerical candidate завершена; разрешённые spatial pairs
не допущены к nonlinear trajectories. Overall pilot PARTIAL, new temporal/
spatial NOT_RUN; old nonlinear PARTIAL остаются. Ни новые IC/basis варианты,
ни более высокие nonlinear p автоматически не выбираются. LONG CLOSED,
EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION, angular same-clamp out-of-plane
reference UNAVAILABLE сохраняются. Следующее направление требует отдельного
явного решения пользователя.

## Историческая остановка после bounded precision/feasibility

Projection recovery PASS, strict initial verification PARTIAL; exploratory
short temporal PASS / spatial PARTIAL. Новая 0.1T1 prepared trajectory не заменяет
историческую zero-u/c задачу и не завершает первоначальный 5T1 pilot.
LONG — CLOSED; EB/RLB-KV — PAUSED_FOR_SUPERVISOR_DIRECTION; angular same-clamp
reference — UNAVAILABLE. Ни full-period/p-refinement, ни новая модель/IC/basis,
ни новый scientific stage автоматически не выбраны или не разрешены.
## Историческая остановка после bounded one-T1 feasibility

Все 3 exploratory trajectories достигли T1; prefix, temporal, energy/mass
PASS, spatial PARTIAL сохраняется. Strict initial verification не повышена,
историческая zero-u/c задача и прежний full 5T1 pilot остаются PARTIAL.
LONG — CLOSED; EB/RLB-KV — PAUSED_FOR_SUPERVISOR_DIRECTION; angular same-clamp
reference — UNAVAILABLE. Автоматическое увеличение horizon/p, смена basis/IC/
physics или новый scientific stage не выбраны и не разрешены.
## Текущая остановка после bounded physical sanity checks

Диагностика завершена с qualifications; weakly nonlinear patterns согласованы
в проверенном reduced scope, physical accuracy V0 и strict numerical
certification не установлены. Новых half-amplitude controls или дальнейших
horizon/p/IC/model/angle/stability studies автоматически не выбирается.
LONG — CLOSED; EB/RLB-KV — PAUSED_FOR_SUPERVISOR_DIRECTION; angular same-clamp
reference — UNAVAILABLE. Historical zero-u/c и prepared strict PARTIAL
сохраняются; следующий scientific stage требует отдельного явного решения.

## Текущая остановка после FEM-readiness audit

Readiness установлена в проверенном локальном scope; nonlinear solver support
документирован, но не execution-tested для нашей задачи. Никаких FEM-расчётов
не выполнено. FEM-1 предложен только для отдельного следующего разрешения;
FEM-2/FEM-3 не запускаются автоматически. LONG CLOSED; EB/RLB-KV PAUSED;
angular same-clamp UNAVAILABLE; NLSP physical sanity qualified и strict PARTIAL
сохраняются. Это инфраструктурный аудит, не новая научная модель.


## Текущая остановка после FEM-1

Линейный четырёхсемейный benchmark выполнен, mesh convergence PARTIAL сохранена.
Это не проверка nonlinear terms V0 и не повод автоматически расширять сетку,
геометрию или запускать FEM-2/FEM-3. Существующие scoped statuses неизменны;
следующий расчёт требует отдельного решения пользователя.
