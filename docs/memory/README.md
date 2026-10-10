# Память исследования сопряжённых стержней

## Назначение

Память — вторичный слой навигации и восстановления контекста: что известно,
при каких условиях, откуда это следует и где остановились. Она не заменяет
теорию, model-specific contracts и отчёты и не разрешает новые действия.
Отсутствие записи означает неполноту покрытия памяти, а не отсутствие
исследования в проекте.

Покрыты отдельные scoped линии: продольные модели и production M-H/Timoshenko
(`LONG`), EB/RLB с поворотной пружиной и Kelvin–Voigt (`RLB`), а также
завершённый circular EB reviewer diagnostic. Их решения и остановки не
переносятся между моделями автоматически.

Отдельное явное задание 2026-10-07 добавляет пространственную нелинейную
семиполевую линию (`NLSP`); она не открывает заново выбор продольной модели.

## Как читать

Этот README → [текущий контекст](current.md) → относящиеся к задаче D/K
→ отчёт, теория и данные по их ссылкам. Читать всю историю перед каждой
задачей не требуется. Решение D задаёт разрешённый объём; результат K
фиксирует фактически установленное и его ограничения.

## Файлы

- [current.md](current.md) — изменяемый текущий контекст, основания следующего решения и место остановки.
- [decisions.md](decisions.md) — дополняемая история решений и их происхождения.
- [knowledge.md](knowledge.md) — дополняемая история научных результатов, оснований и ограничений.

## Темы

| Тема | Ключевые записи |
|---|---|
| FEM-3C: limited straight-rod verification complete with qualifications; C1 PASS, full-period illustration, all8/energy PARTIAL | [NLSP-D16](decisions.md#nlsp-d16), [NLSP-K17](knowledge.md#nlsp-k17), [report](../numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3c-numerical-robustness-and-dissertation-verification), [scientific summary](../numerics/nlsp_straight_rod_3d_fem_verification_summary.md) |
| FEM-3B: .25T1 выбран по 1D до FEM; 2 native jobs/502 frames; evolving signal измерен, 3D accuracy и energy qualified | [NLSP-D15](decisions.md#nlsp-d15), [NLSP-K16](knowledge.md#nlsp-k16), [report](../numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3b-nonlinear-correction-evolution-and-longer-horizon) |
| FEM-3AR: two corrected native jobs and p64 references complete .05T1; release/motion pass, energy PARTIAL; dynamic correction uncertified | [NLSP-D14](decisions.md#nlsp-d14), [NLSP-K15](knowledge.md#nlsp-k15), [continuation](../numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3ar-controlled-continuation-after-native-elke-output-path-failure) |
| FEM-3A: preload/release protocol documented, p64 preflight PASS; first native job access violation, no trajectories | [NLSP-D13](decisions.md#nlsp-d13), [NLSP-K14](knowledge.md#nlsp-k14), [pilot](../numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md) |
| FEM-2R: отдельное разрешённое продолжение static comparison; исторический input failure сохранён | [NLSP-D12](decisions.md#nlsp-d12), [NLSP-K13](knowledge.md#nlsp-k13), [continuation](../numerics/nlsp_nonlinear_static_3d_fem_validation.md#fem-2r-controlled-continuation-after-input-serialization-failure) |
| FEM-2: 1D quartic static preflight PASS; first medium input-parser failure, 3D comparison not reached | [NLSP-D11](decisions.md#nlsp-d11), [NLSP-K12](knowledge.md#nlsp-k12), [static report](../numerics/nlsp_nonlinear_static_3d_fem_validation.md) |
| FEM-1R: one .020 mesh; all8 matched and preset mesh criterion accepted; linear only | [NLSP-D10](decisions.md#nlsp-d10), [NLSP-K11](knowledge.md#nlsp-k11), [continuation](../numerics/nlsp_linear_rectangular_3d_fem_validation.md#fem-1r-one-additional-mesh-refinement) |
| Full-family rectangular linear3D FEM-1:3 grids/4 families identified, mesh convergence PARTIAL | [NLSP-D09](decisions.md#nlsp-d09), [NLSP-K10](knowledge.md#nlsp-k10), [report](../numerics/nlsp_linear_rectangular_3d_fem_validation.md) |
| Локальная3D FEM readiness: Gmsh/CCX confirmed, nonlinear support documented; без новых jobs | [NLSP-K09](knowledge.md#nlsp-k09), [readiness audit](../numerics/nlsp_3d_fem_environment_readiness.md) |
| Семиполевая модель: V0 → exact/cubic audit → первый four-field planar time pilot PARTIAL; без stability/порога | [NLSP-D01](decisions.md#nlsp-d01), [NLSP-K01](knowledge.md#nlsp-k01), [NLSP-D02](decisions.md#nlsp-d02), [NLSP-K02](knowledge.md#nlsp-k02), [numerical note](../theory/weakly_nonlinear_planar_time_pilot.md) |
| Адресная planar диагностика: physical L2 split / initial compatibility / lazy RHS recovery; full p48 отложен, PARTIAL | [NLSP-D03](decisions.md#nlsp-d03), [NLSP-K03](knowledge.md#nlsp-k03), [continuation](../theory/weakly_nonlinear_planar_time_pilot.md#14-адресная-диагностика-и-восстановление-вычислительного-пути-2026-10-07) |
| Аналитический по времени leading u2,c2 контроль: все Shen coordinates p16…96, 0 ODE; diagnostic COMPLETE / spatial PARTIAL | [NLSP-D04](decisions.md#nlsp-d04), [NLSP-K04](knowledge.md#nlsp-k04), [note](../theory/planar_second_order_axial_response.md) |
| Prepared initial state: stat/harm и O3 endpoint correction PASS; full projection PARTIAL,0 новых trajectories | [NLSP-D05](decisions.md#nlsp-d05), [NLSP-K05](knowledge.md#nlsp-k05), [note](../theory/planar_prepared_initial_state.md) |
| Prepared numerical feasibility: initial-only endpoint L2 PASS;3 short exploratory trajectories; strict/spatial PARTIAL, temporal PASS | [NLSP-D06](decisions.md#nlsp-d06), [NLSP-K06](knowledge.md#nlsp-k06), [continuation](../theory/planar_prepared_initial_state.md#prepared-precision-feasibility) |
| То же prepared движение до одного linear T1:3 exploratory runs; prefix/temporal/energy PASS, spatial7/8 PARTIAL | [NLSP-D07](decisions.md#nlsp-d07), [NLSP-K07](knowledge.md#nlsp-k07), [continuation](../theory/planar_prepared_initial_state.md#prepared-one-period-feasibility) |
| Продольные модели: Bishop diagnostic → M-H/Timoshenko production selection → joint gates → hierarchy/thickness/Lambda checks → закрыта в принятом 1D scope | [LONG-D01](decisions.md#long-d01), [LONG-K01](knowledge.md#long-k01), [LONG-D02](decisions.md#long-d02) |
| Базовая RLB, эквивалентность ламинатов и переносы | [K01](knowledge.md#rlb-k01), [K02](knowledge.md#rlb-k02), [K03](knowledge.md#rlb-k03), [K04](knowledge.md#rlb-k04) |
| Поворотная пружина: теория, EB-пилот, RLB → EB, карты, формы, шарнир, устойчивость механизма | [K05](knowledge.md#rlb-k05), [K06](knowledge.md#rlb-k06), [K07](knowledge.md#rlb-k07), [K08](knowledge.md#rlb-k08), [K09](knowledge.md#rlb-k09), [K10](knowledge.md#rlb-k10), [K11](knowledge.md#rlb-k11) |
| Отдельный reviewer diagnostic: круглые EB-стержни, SPRING → RIGID | [D25](decisions.md#rlb-d25), [K26](knowledge.md#rlb-k26), [научная заметка](../laminated_beams/circular_eb_rotational_spring_rigid_limit.md); прежние остановки [EB-JOINT-K01](knowledge.md#eb-joint-k01), [K02](knowledge.md#eb-joint-k02) |
| Kelvin–Voigt: теория и внутренний пилот | [K12](knowledge.md#rlb-k12) |
| Внешняя литературная проверка и точность печати | [K13](knowledge.md#rlb-k13), [K14](knowledge.md#rlb-k14) |
| Упругий screening участия демпфера | [K15](knowledge.md#rlb-k15) |
| Адресное слабое демпфирование, диагностика, production routing и завершение | [K16](knowledge.md#rlb-k16), [K17](knowledge.md#rlb-k17), [K18](knowledge.md#rlb-k18), [K19](knowledge.md#rlb-k19) |
| Общая чётность по вязкости | [K20](knowledge.md#rlb-k20), [research note](../laminated_beams/inplane_kelvin_voigt_weak_damping_parity.md) |
| Production RLB+KV и точный EB-предел | [K21](knowledge.md#rlb-k21) |
| Sparse elastic RLB screening и сопоставление с EB | [K22](knowledge.md#rlb-k22) |
| Заключительное complex-сравнение и остановка перед обсуждением | [D22](decisions.md#rlb-d22), [K23](knowledge.md#rlb-k23), [сводка для руководителя](../laminated_beams/inplane_kelvin_voigt_research_status_for_supervisor.md) |
| Рисунки для будущей заметки: crossing/veering и формы | [D23](decisions.md#rlb-d23), [K24](knowledge.md#rlb-k24), [Figure 3](../laminated_beams/figure03_crossing_veering.md) |
| Рисунки: упругие EB/RLB-формы и поворот у узла | [D24](decisions.md#rlb-d24), [K25](knowledge.md#rlb-k25), [Figure 4 revised, v1 сохранена](../laminated_beams/figure04_eb_rlb_shapes.md) |
| Рисунки: огибающие затухания EB/RLB | [K23](knowledge.md#rlb-k23), [Figure 5, R2/R3 при d=.001](../laminated_beams/figure05_damping_envelopes.md); [прежние столбцы](../laminated_beams/figure05_frequency_damping_comparison.md) |

Записи сохраняют ссылки на решения и подробные источники; таблица тем
не заменяет их условия и статусы. Общая навигация:
[исследования](../research_index.md), [результаты](../results_index.md),
[тематический README](../laminated_beams/README.md).

## Иерархия источников

1. Математическая/source theory или публикация, с учётом проверенных формул, допущений и предупреждений источника.
2. Tracked научный отчёт с результатом в объявленных границах.
3. Код, тесты и сохранённые расчётные данные.
4. Краткая запись памяти.

Математический приоритет определяется [AGENTS.md](../../AGENTS.md)
и [правилами проекта](../project_rules.md); эта навигация его не меняет.
Расхождение источников отмечается явно, а не исправляется пересказом.
Для локальных ignored CSV/JSON tracked Markdown-отчёт должен сохранять
достаточный научный смысл и квалификации без наличия локальных файлов.

## Соглашения о новых записях

Компактные ориентиры, не обязательная бюрократическая схема:

- **Decision:** Status / decision; Why; Scope; Stop/revisit condition; Result link when completed.
- **Knowledge:** Status; Claim/result; Scope; Basis; Limitations; Implication; Extends/does-not-invalidate, когда это существенно.

Указывать дату и происхождение, различать разрешение, гипотезу и результат.
Память хранит утверждение, условия и ссылку; большие таблицы корней,
перечни невязок, счётчиков и тестов остаются в отчётах. Новый формат
применяется к новым записям, без переоформления старой истории.

## Правило роста

current.md можно переписывать: это временный рабочий контекст.
D/K остаются append-only с устойчивыми номерами и якорями; исторические
записи меняются только для явной опечатки или битой ссылки. Новое решение
или результат не снимает старые квалификации молча. Косметическая правка
не требует новой научной записи.

Подробные источники — отчёты; raw CSV/JSON служат дополнительным основанием,
а не единственным долговечным изложением. Обновляются только затронутые
записи и текущая остановка в разрешённом объёме. Сейчас файлы не разделяются:
если поиск по decisions/knowledge станет неудобен, тематическое разделение
потребует отдельной миграции с сохранением якорей и истории.
