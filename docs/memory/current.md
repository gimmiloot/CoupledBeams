# Текущий контекст

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
