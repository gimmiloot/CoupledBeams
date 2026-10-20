# Проверка упрощения памяти и навигации — 2026-10-20

Технический отчёт документационной задачи, не новый научный этап.
Новых D/K, вычислений, моделей, статусов исследований и форматов памяти нет.
Текущее состояние направлений читать в [current.md](current.md);
этот отчёт фиксирует проверку на указанную дату и не требует постоянного обновления.

## Исходное состояние и сохранность

Исходный HEAD и локальный `origin/main`:
`fca94541e9ecfd6f7eaddad3b6872cd4e72fa527` (Version 0.6.16).
`git status --porcelain=v1` пуст: незакоммиченных и untracked-файлов не было.
Read-only `git ls-remote origin refs/heads/main` вернул тот же SHA.
Это отличается от состояния предыдущего аудита: новые NSL-документы уже
включены в main. Git add/commit/push/reset/clean/fetch не выполнялись.

До изменений сохранены точные байтовые копии 11 целевых файлов и baseline
SHA256 во временной папке ОС вне репозитория. Местоположение передано
пользователю отдельно. Это страховка конкретной правки, не постоянный архив.
Исходное содержимое `current.md` прочитано полностью и может быть восстановлено
из этой копии; исходный Git-коммит также сохраняет текст.

Исходный SHA256 `current.md`:
`ca7792f3694031bae94ff7120535fc861dc1e202331a9e3c9f4125db8c337ce3`.
Исходные хеши всех целевых файлов приведены ниже.
Защищённый набор включает все исходно отслеживаемые нецелевые файлы
и всё локальное `results/`, всего 56 877 файлов.
Проверка сохранности: **PASS**. Повторный полный SHA256-проход по всем
56 877 файлам дал ноль изменений и ноль отсутствующих файлов. В том числе
проверены 56 009 файлов `results/`, среди них 678 файлов `manifest.json`.
D/K, исторические научные отчёты и исходно отслеживаемый код неизменны.

## Изменённая навигация

- `current.md`: 1297 → 162 строки; семь заголовков первого уровня
  «Текущий контекст» → один. Убраны 11 разделителей «Предыдущий контекст»
  и 14 повторных блоков остановок в конце файла. Их источники перечислены ниже.
- `memory/README.md`: маршрут по независимым направлениям вместо очередной
  хронологической сводки; старые вводные указатели внутри D/K явно
  квалифицированы, сами защищённые файлы не редактировались.
- `AGENTS.md`, `project_rules.md`, `theory/AGENTS.md`: уточнены модельная
  область иерархии источников и ответственность документов. Порядок
  источников базовой модели и проверенные формулы не изменены.
- Главный README и индекс исследований получили короткий начальный маршрут.
  NSL/RLB/NLSP получили паспорта в существующих README; LONG — в индексе.
  Исторические каталоги этапов сохранены и явно отделены от текущего состояния.
- Новые результаты не требуют механически повторять числа во всех индексах.
  Журнал, индекс результатов, assumptions и source index не изменялись:
  новых научных фактов, расчётных материалов или допущений эта работа не добавляет.

## Проверка существенного содержания

Сопоставлены 60 содержательных блоков уровня `##` исходного `current.md`
и отдельно вложенная заметка о переработке отчёта NLSP. Для каждого
научного этапа проверены ссылки на D/K и соответствующий канонический отчёт.
Дополнительно сопоставлены статусные маркеры и длинные численные значения;
различия записи проверены как округления значений из отчётов, не новые результаты.
Это проверка сохранности изложения, не повторная математическая верификация.

Особого внимания потребовали:

- **NSL-1 / NSL-2 / NSL-2R:** печатная квалификация Bachoo отделена
  от концевых невязок NSL-2, устранённых только последующим NSL-2R.
- **NLSP:** численно подтверждённый OUT-рост не снимает PL neutral PARTIAL,
  не доказывает полную p-сходимость и не заменяет прежние отрицательные
  screening-результаты или отдельный Construction C.
- **FEM:** ошибки ввода/вывода FEM-2/FEM-3A отделены от успешных продолжений;
  энергия, c_eff, малые нелинейные поправки и строгие prepared-state
  критерии сохраняют свои ограничения.
- **Circular EB:** результат шести seeds не отменяет high-spectrum,
  direct-matching и немонотонность seed 06; реальный узел не откалиброван.
- **Вложенная редакционная заметка, старые строки 186–195:** описание
  26-страничного отчёта, восьми рисунков, 36 коэффициентных проверок и
  ограничения поля c уже сохранено в
  [README отчёта](../reports/nlsp_beta15_parametric_instability/README.md)
  и [аудите текста](../reports/nlsp_beta15_parametric_instability/model_audit.md).
  Старое упоминание «PDF, 5 страниц» относится к предыдущей версии;
  защищённый технический отчёт ради этой подписи не редактировался.
- Действующие остановки LONG и RLB сохраняются независимо от новых
  отдельных заданий NSL/NLSP. Исторический same-clamp UNAVAILABLE не
  превращён в общую валидацию NLSP через другой линейный узел LONG-K04.

Единственных существенных научных утверждений, требующих переноса
в защищённый отчёт или новой научной записи, не обнаружено.

## Cold-start: ответы только по документам

Это документная симуляция нового читателя, без предыдущей переписки
и без доступа к массивам `results/`. Независимый человеческий тест не проводился.
Маршрут: current → паспорт → прямо указанный итоговый отчёт. Полная история D/K
не требовалась; каждая колонка ниже проверена по этому маршруту.

| Вопрос | NSL | RLB / EB–KV | NLSP |
| --- | --- | --- | --- |
| Физическая модель | Несимметричные плоские FSDT-стержни, совместная ABD-редукция, упругий/KV точечный узел | Симметричная Reddy/FSDT и EB-контроль, поворотный KV-узел | Семь полей, действие V0, кубические балансы, переменная инерция, жёсткий common-DOF узел |
| Канонические уравнения | Single-beam / coupled-joint / KV notes из паспорта NSL | Source contract, rigid/spring notes и KV TeX из паспорта RLB | Weakly nonlinear spatial rod и научный отчёт руководителю |
| Последний принятый итог | NSL-VIS-1 PASS для Ritz-карт, не общей континуальной сертификации | K23: ограниченное сравнение пяти EB/RLB-пар; K26: отдельный circular EB diagnostic | K33: OUT-рост выбранной Type-II ветви, воспроизведение p14, PL PARTIAL |
| Открыто | Печать Bachoo, часть идентичностей форм, общие карты комплексного затухания/FRF | Source-print и исторические спектральные квалификации, физическая калибровка узла | PL neutral gate, полная p-сходимость, 3D-валидация углового узла, deferred action expansion |
| Что устранено позже | NSL-2R: конкретные концевые невязки NSL-2 | Архитектура/восстановление A/C и отдельное завершение weak-damping | FEM-2R/FEM-3AR: ошибки прежнего исполнения; причина PL residual объяснена, строгий gate не закрыт |
| Где и почему остановка | Конечный объём VIS-1 выполнен | PAUSED_FOR_SUPERVISOR_DIRECTION; reviewer diagnostic завершён отдельно | STOP до обсуждения с руководителем, ограничения не сняты |
| Доказательства и данные | Заметки Git; локальные bundles и команды по ссылкам в заметках | Сводка руководителю и complex-confirmation report; локальные результаты laminated_beams | Научный LaTeX/PDF, pre-supervisor technical report; локальный bundle d3f6375a27dba883 |
| Что не разрешено автоматически | NSL-3B, вязкость/FRF, оптимизация и дальнейшие карты | Новая физическая программа, вязкости и калибровка | Новые ODE/FEM/Floquet, уточнения, тип I, насыщение и параметры |

Пробелы остаются локальными и явными: неполное покрытие других направлений
в D/K компенсируется их тематическими README; фактические массивы и кеши
не опубликованы вместе с Git-контекстом. Это не требует чтения прежних
сводок и не является новым научным отказом.

## Git visibility и перечень для следующего коммита

Проверены 112 выбранных навигационных, канонических и непосредственно
связанных файлов: все отслеживаются и присутствуют в дереве `origin/main`
при подтверждённом удалённом SHA. Включены текущие отчёты NSL, RLB, NLSP,
научный LaTeX/PDF NLSP и источники исторических блоков. Не найдено
канонического научного утверждения из удаляемого current, доступного
только в untracked/ignored-файле.

Локально изменены 11 перечисленных ниже отслеживаемых документов;
этот технический отчёт — один новый untracked-файл. Новая редакция навигации станет
доступна другим чатам только после отдельного пользовательского коммита
и публикации. Большие CSV/NPZ добавлять ради этой задачи не требуется.

Рекомендуемый состав одного документационного коммита:

1. `AGENTS.md`
2. `docs/project_rules.md`
3. `docs/theory/AGENTS.md`
4. `docs/memory/README.md`
5. `docs/memory/current.md`
6. `README.md`
7. `docs/research_index.md`
8. `docs/nonsymmetric_laminated_beams/README.md`
9. `docs/laminated_beams/README.md`
10. `docs/numerics/README.md`
11. `CHANGELOG.md`
12. `docs/memory/navigation_review_2026_10_20.md`

## Проверки ссылок и границ задачи

Сопоставимый аудит тех же девяти файлов: **1223 → 1044** локальные ссылки,
ноль отсутствующих файлов и ноль неразрешённых D/K-якорей.
Отдельная проверка 295 новых/изменённых ссылок, включая весь новый current
и этот отчёт, дала ноль ошибок путей и Markdown-якорей.
Сокращение числа ссылок связано с удалением повторных сводок; прежние
файлы и D/K-якоря сохранены. Входящих ссылок на удалённые разделы
`current.md#...` в Markdown проекта не найдено.
Основной маршрут заканчивается в научном отчёте; обратные навигационные
ссылки не требуют обхода истории. Скрипт общего назначения не добавлен.

Научные расчёты, построение графиков, solver/import-regression и FEM не
запускались. Проверки выполняются чтением файлов, SHA256, Markdown-ссылок
и read-only Git-командами. Исторические научные статусы не пересмотрены.
`git diff --check` проходит. Индекс Git и HEAD не изменены; изменения
ограничены 11 существующими документами и этим новым отчётом.

**Итог документационной задачи: PASS.** Основной маршрут укорочен,
доказательства и отрицательные результаты сохранены. Новая бюрократическая
система или обязательный новый шаблон научной записи не введены.
Дальнейшая миграция D/K, архивирование данных и научные продолжения не начаты.

## Карта исходных блоков current.md

Номера строк относятся к исходным 1297 строкам при указанном HEAD и SHA256.
«Исторический этап» означает сохранение результата в источниках, а не его
научную отмену. Повторные остановки в конце сопоставлены с теми же источниками.
Таблица является журналом этой проверки, не новым реестром научных результатов.

| Original lines / block | D/K | Canonical report | Qualification retained |
| --- | --- | --- | --- |
| 3-21 / NSL-VIS-1 — две сравнительные Ritz-карты, 2026-10-20 | [NSL-D07](decisions.md#nsl-d07), [NSL-K07](knowledge.md#nsl-k07) | [source 1](../nonsymmetric_laminated_beams/comparative_spectral_curves.md) | Ritz-map PASS only; no continuum-prefix certificate |
| 22-44 / NSL-3A-ANGLE-C — два комплексных корня и три рисунка, 2026-10-20 | [NSL-D06](decisions.md#nsl-d06), [NSL-K06](knowledge.md#nsl-k06) | [source 1](../nonsymmetric_laminated_beams/nsl3a_angle_complex_confirmation.md) | Two targets; small damping, no forced response |
| 45-69 / NSL-3A-ANGLE — упругий обзор завершён, 2026-10-20 | [NSL-D05](decisions.md#nsl-d05), [NSL-K05](knowledge.md#nsl-k05) | [source 1](../nonsymmetric_laminated_beams/nsl3a_elastic_angle_screening.md) | Finite screening; MAC unresolved at 60/75 degrees |
| 70-96 / NSL-3A — ограниченный комплексный пилот завершён, 2026-10-20 | [NSL-D04](decisions.md#nsl-d04), [NSL-K04](knowledge.md#nsl-k04) | [source 1](../nonsymmetric_laminated_beams/nsl3a_kelvin_voigt_pilot.md) | Three active roots; structural protected zeros |
| 97-120 / NSL-2R — численная диагностика завершена, 2026-10-20 | [NSL-D03](decisions.md#nsl-d03), [NSL-K03](knowledge.md#nsl-k03) | [source 1](../nonsymmetric_laminated_beams/nsl2r_numerical_diagnostics.md) | Recovery repaired; original frequencies and NSL-2 status retained |
| 121-151 / NSL-2 — упругой угловой узел, завершён 2026-10-20 | [NSL-D02](decisions.md#nsl-d02), [NSL-K02](knowledge.md#nsl-k02) | [source 1](../nonsymmetric_laminated_beams/coupled_joint_theory.md) | Historical BC failures and inactive derivative control retained |
| 152-196 / NSL-1 — самостоятельное направление, этап завершён 2026-10-20 | [NSL-D01](decisions.md#nsl-d01), [NSL-K01](knowledge.md#nsl-k01) | [source 1](../nonsymmetric_laminated_beams/single_beam_theory.md), [source 2](../nonsymmetric_laminated_beams/source_comparison.md) | Source-print PARTIAL; model/numerical PASS |
| 197-214 / β15: завершающие проверки и отчёт руководителю — STOP, 2026-10-12 | [NLSP-D32](decisions.md#nlsp-d32), [NLSP-K33](knowledge.md#nlsp-k33) | [source 1](../numerics/nlsp_beta15_pre_supervisor_checks.md), [source 2](../reports/nlsp_beta15_parametric_instability/README.md) | OUT confirmed; PL neutral gate still PARTIAL |
| 217-232 / β15: p14 OUT growth воспроизведён — STOP, 2026-10-11 | [NLSP-D31](decisions.md#nlsp-d31), [NLSP-K32](knowledge.md#nlsp-k32) | [source 1](../numerics/nlsp_beta15_out_growth_p14_convergence.md) | p14 reproduced; two p levels are not full convergence |
| 236-249 / β15: динамически подтверждённый OUT growth — STOP, 2026-10-11 | [NLSP-D30](decisions.md#nlsp-d30), [NLSP-K31](knowledge.md#nlsp-k31) | [source 1](../numerics/nlsp_beta15_dynamic_out_growth_validation.md) | Finite-time validation; no independent physical validation |
| 253-266 / β15: уточнённая периодическая ветвь и OUT spectral candidate — STOP,2026-10-11 | [NLSP-D29](decisions.md#nlsp-d29), [NLSP-K30](knowledge.md#nlsp-k30) | [source 1](../numerics/nlsp_beta15_periodic_orbit_spatial_stability.md) | One HB-gated candidate; three amplitudes NOT_EVALUATED |
| 267-285 / β15: двухмодовый плоский 1:3 — остановка,2026-10-11 | [NLSP-D28](decisions.md#nlsp-d28), [NLSP-K29](knowledge.md#nlsp-k29) | [source 1](../numerics/nlsp_beta15_planar_1to3_internal_resonance.md) | Controlled 1:3 approximation; Type I differs from Type II |
| 290-308 / β15: нелинейная плоская поправка — остановка, 2026-10-11 | [NLSP-D27](decisions.md#nlsp-d27), [NLSP-K28](knowledge.md#nlsp-k28) | [source 1](../numerics/nlsp_beta15_nonlinear_planar_resonance.md) | Single-mode finite-amplitude gate FAIL; M1/M2 unavailable |
| 313-331 / Принципиальная возможность параметрического роста — остановка,2026-10-11 | [NLSP-D26](decisions.md#nlsp-d26), [NLSP-K27](knowledge.md#nlsp-k27) | [source 1](../numerics/nlsp_out_of_plane_instability_feasibility.md) | Formal channel only; direct actual pair remains oscillatory |
| 336-353 / Пять углов β — остановка, 2026-10-11 | [NLSP-D25](decisions.md#nlsp-d25), [NLSP-K26](knowledge.md#nlsp-k26) | [source 1](../numerics/nlsp_beta_resonant_coupling_screening.md) | NO_DIRECT_GROWTH_CANDIDATE; planar identity unresolved |
| 358-376 / Нелинейные модальные связи C — остановка, 2026-10-11 | [NLSP-D24](decisions.md#nlsp-d24), [NLSP-K25](knowledge.md#nlsp-k25) | [source 1](../numerics/nlsp_resonant_spatial_coupling_screening.md) | Coupling present; C growth mechanism unidentified |
| 381-407 / Construction C до низшего OUT периода — остановка, 2026-10-11 | [NLSP-D23](decisions.md#nlsp-d23), [NLSP-K24](knowledge.md#nlsp-k24) | [source 1](../numerics/nlsp_construction_C_longer_sensitivity.md) | Exploratory p12; p16 NOT_RUN by cost; output repair preserved |
| 412-441 / Construction C modal experiment — текущая остановка, 2026-10-11 | [NLSP-D22](decisions.md#nlsp-d22), [NLSP-K23](knowledge.md#nlsp-k23) | [source 1](../numerics/nlsp_construction_C_modal_sensitivity.md) | Restricted gain; c/velocities spatial PARTIAL |
| 444-459 / Coupled sensitivity continuation — текущий итог | [NLSP-D21](decisions.md#nlsp-d21), [NLSP-K22](knowledge.md#nlsp-k22) | [source 1](../numerics/nlsp_coupled_planar_sensitivity_resume.md) | Restricted A/B gain; long time/spatial qualification PARTIAL |
| 463-485 / Torsional diagnostic: NUMERICAL_FLOOR_IDENTIFIED, 2026-10-10 | [NLSP-D20](decisions.md#nlsp-d20), [NLSP-K21](knowledge.md#nlsp-k21) | [source 1](../numerics/nlsp_coupled_torsion_tangent_diagnostic.md) | Numerical floor localized; parent failure retained |
| 486-510 / Двухплечевой planar-sensitivity pilot: BLOCKED_BY_TANGENT_CHECK, 2026-10-10 | [NLSP-D19](decisions.md#nlsp-d19), [NLSP-K20](knowledge.md#nlsp-k20) | [source 1](../numerics/nlsp_coupled_planar_sensitivity_pilot.md) | Historical BLOCKED_BY_TANGENT_CHECK; later tests separate |
| 511-555 / Семиполевое пространственное сопоставление: NUMERICAL_PARTIAL, 2026-10-10 | [NLSP-D18](decisions.md#nlsp-d18), [NLSP-K19](knowledge.md#nlsp-k19) | [source 1](../numerics/nlsp_spatial_nonlinear_3d_fem_verification.md) | NUMERICAL_PARTIAL; corrections/energy/c_eff qualified |
| 556-582 / Spatial-profile diagnostic complete with qualifications, 2026-10-10 | [NLSP-D17](decisions.md#nlsp-d17), [NLSP-K18](knowledge.md#nlsp-k18) | [source 1](../numerics/nlsp_spatial_profile_audit.md) | c_eff is a recovery proxy; no new physical solves |
| 583-626 / FEM-3C — limited straight-rod verification complete, 2026-10-10 | [NLSP-D16](decisions.md#nlsp-d16), [NLSP-K17](knowledge.md#nlsp-k17) | [source 1](../numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3c-numerical-robustness-and-dissertation-verification), [source 2](../numerics/nlsp_straight_rod_3d_fem_verification_summary.md) | Limited straight-rod verification; full-horizon/energy qualified |
| 627-658 / FEM-3B — эволюция поправки измерена, с qualifications, 2026-10-09 | [NLSP-D15](decisions.md#nlsp-d15), [NLSP-K16](knowledge.md#nlsp-k16) | [source 1](../numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3b-nonlinear-correction-evolution-and-longer-horizon) | Evolution measured; original time/mesh/energy limits |
| 659-695 / FEM-3AR — short pilot complete with qualifications, 2026-10-09 | [NLSP-D14](decisions.md#nlsp-d14), [NLSP-K15](knowledge.md#nlsp-k15) | [source 1](../numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3ar-controlled-continuation-after-native-elke-output-path-failure) | Successful continuation; historical native failure retained |
| 696-724 / FEM-3A — BLOCKED_BY_SOLVER, 2026-10-09 | [NLSP-D13](decisions.md#nlsp-d13), [NLSP-K14](knowledge.md#nlsp-k14) | [source 1](../numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md) | BLOCKED_BY_SOLVER; ELKE output-path issue localized |
| 725-750 / FEM-2R — static diagnostic complete with qualifications, 2026-10-09 | [NLSP-D12](decisions.md#nlsp-d12), [NLSP-K13](knowledge.md#nlsp-k13) | [source 1](../numerics/nlsp_nonlinear_static_3d_fem_validation.md#fem-2r-controlled-continuation-after-input-serialization-failure) | Qualified static continuation, not universal V0 validation |
| 751-777 / FEM-2 — 1D preflight complete, first 3D input failure, 2026-10-09 | [NLSP-D11](decisions.md#nlsp-d11), [NLSP-K12](knowledge.md#nlsp-k12) | [source 1](../numerics/nlsp_nonlinear_static_3d_fem_validation.md) | Historical input serialization failure, not physical failure |
| 778-797 / FEM-1R — one refined mesh, all eight preset-accepted, 2026-10-09 | [NLSP-D10](decisions.md#nlsp-d10), [NLSP-K11](knowledge.md#nlsp-k11) | [source 1](../numerics/nlsp_linear_rectangular_3d_fem_validation.md#fem-1r-one-additional-mesh-refinement) | Linear preset accepted; no exact continuum/nonlinear validation |
| 798-811 / Full-family linear3D FEM-1 — bounded complete / mesh PARTIAL, 2026-10-08 | [NLSP-D09](decisions.md#nlsp-d09), [NLSP-K10](knowledge.md#nlsp-k10) | [source 1](../numerics/nlsp_linear_rectangular_3d_fem_validation.md) | Original mesh PARTIAL retained |
| 812-827 / 3D FEM infrastructure audit — завершён, 2026-10-08 | [NLSP-K09](knowledge.md#nlsp-k09) | [source 1](../numerics/nlsp_3d_fem_environment_readiness.md) | Readiness only; not execution/physical validation |
| 828-852 / Planar physical sanity — diagnostic complete with qualifications, 2026-10-08 | [NLSP-D08](decisions.md#nlsp-d08), [NLSP-K08](knowledge.md#nlsp-k08) | [source 1](../theory/weakly_nonlinear_planar_physical_sanity_checks.md) | Physical sanity qualified; strict certification PARTIAL |
| 853-874 / Исторический prepared one-T1 feasibility — exploratory, spatial PARTIAL, 2026-10-08 | [NLSP-D07](decisions.md#nlsp-d07), [NLSP-K07](knowledge.md#nlsp-k07) | [source 1](../theory/planar_prepared_initial_state.md#prepared-one-period-feasibility) | One-T1 exploratory; spatial/strict PARTIAL |
| 875-898 / Исторический prepared-state precision/feasibility — exploratory, spatial PARTIAL, 2026-10-08 | [NLSP-D06](decisions.md#nlsp-d06), [NLSP-K06](knowledge.md#nlsp-k06) | [source 1](../theory/planar_prepared_initial_state.md#prepared-precision-feasibility) | Initial projection recovery; strict float64 PARTIAL |
| 899-917 / Исторический prepared initial-state gate — PARTIAL, 2026-10-08 | [NLSP-D05](decisions.md#nlsp-d05), [NLSP-K05](knowledge.md#nlsp-k05) | [source 1](../theory/planar_prepared_initial_state.md) | Prepared projection PARTIAL; no trajectories in this stage |
| 918-935 / Leading second-order axial diagnostic — COMPLETE / spatial PARTIAL, 2026-10-08 | [NLSP-D04](decisions.md#nlsp-d04), [NLSP-K04](knowledge.md#nlsp-k04) | [source 1](../theory/planar_second_order_axial_response.md) | Exact-time second order; spatial PARTIAL |
| 936-953 / Planar numerical recovery — PARTIAL, 2026-10-07 | [NLSP-D03](decisions.md#nlsp-d03), [NLSP-K03](knowledge.md#nlsp-k03) | [source 1](../theory/weakly_nonlinear_planar_time_pilot.md#14-адресная-диагностика-и-восстановление-вычислительного-пути-2026-10-07) | Physical L2 recovery; refinement deferred by budget |
| 954-967 / Нелинейная пространственная ветвь — первый planar time pilot PARTIAL, 2026-10-07 | [NLSP-D02](decisions.md#nlsp-d02), [NLSP-K02](knowledge.md#nlsp-k02) | [source 1](../theory/weakly_nonlinear_planar_time_pilot.md) | Initial nonlinear time pilot spatial PARTIAL |
| 968-987 / Предыдущий семиполевой model/action audit — сохранён | [NLSP-D01](decisions.md#nlsp-d01), [NLSP-K01](knowledge.md#nlsp-k01) | [source 1](../theory/weakly_nonlinear_spatial_rod.md) | Reduced V0 assumption; not full 3D reduction |
| 988-1015 / Одноосное пространственное KV соединение — линейная проверка завершена, 2026-10-19 | [LONG-D05](decisions.md#long-d05), [LONG-K04](knowledge.md#long-k04) | [source 1](../theory/uniaxial_kelvin_voigt_split.md) | Same-model linear OUT identity; no real-3D-joint claim |
| 1016-1036 / Продольные теории при Kelvin–Voigt — завершено, 2026-10-19 | [LONG-D04](decisions.md#long-d04), [LONG-K03](knowledge.md#long-k03) | [source 1](../theory/longitudinal_kelvin_voigt_comparison.md) | Bounded E/RL/MH-FREE comparison; contraction-BC limits |
| 1037-1056 / Отдельная чувствительность к условиям на c — завершена, 2026-10-19 | [LONG-D03](decisions.md#long-d03), [LONG-K02](knowledge.md#long-k02) | [source 1](../theory/contraction_joint_sensitivity.md) | Bounded c-closure diagnostic; not reopening model choice |
| 1057-1092 / Продольная модель сопряжённых стержней — закрыта, 2026-10-07 | [LONG-D01](decisions.md#long-d01), [LONG-K01](knowledge.md#long-k01), [LONG-D02](decisions.md#long-d02) | [source 1](../theory/mindlin_herrmann_timoshenko_single_rod.md), [source 2](../theory/mindlin_herrmann_timoshenko_rigid_joint.md), [source 3](../theory/coupled_longitudinal_theory_hierarchy_screening.md) | LONG closed; close-pair study declined; common-DOF qualification |
| 1093-1118 / EB/RLB rotational spring / Kelvin–Voigt — отдельная ветвь, на паузе | [RLB-D22](decisions.md#rlb-d22), [RLB-K23](knowledge.md#rlb-k23) | [source 1](../laminated_beams/inplane_kelvin_voigt_research_status_for_supervisor.md), [source 2](../laminated_beams/figure03_crossing_veering.md), [source 3](../laminated_beams/figure04_eb_rlb_shapes.md), [source 4](../laminated_beams/figure05_damping_envelopes.md), [source 5](../laminated_beams/figure05_frequency_damping_comparison.md) | KV already implemented; paused for supervisor direction |
| 1119-1166 / Circular EB reviewer diagnostic — завершён, 2026-10-05 | [RLB-D25](decisions.md#rlb-d25), [RLB-K26](knowledge.md#rlb-k26), [EB-JOINT-K01](knowledge.md#eb-joint-k01), [EB-JOINT-K02](knowledge.md#eb-joint-k02) | [source 1](../laminated_beams/circular_eb_rotational_spring_rigid_limit.md) | Six seeds accepted; high-spectrum/direct-matching/seed06 limits |
| 1167-1173 / Историческая остановка после memory-sync | [LONG-D01](decisions.md#long-d01), [LONG-K01](knowledge.md#long-k01), [LONG-D02](decisions.md#long-d02) | [source 1](../theory/mindlin_herrmann_timoshenko_single_rod.md), [source 2](../theory/mindlin_herrmann_timoshenko_rigid_joint.md), [source 3](../theory/coupled_longitudinal_theory_hierarchy_screening.md) | Repeated stop of stage at old line 1057; no new scientific result |
| 1174-1180 / Историческая остановка после NLSP audit | [NLSP-D01](decisions.md#nlsp-d01), [NLSP-K01](knowledge.md#nlsp-k01) | [source 1](../theory/weakly_nonlinear_spatial_rod.md) | Repeated stop of stage at old line 968; no new scientific result |
| 1181-1192 / Текущая остановка после planar time pilot | [NLSP-D02](decisions.md#nlsp-d02), [NLSP-K02](knowledge.md#nlsp-k02) | [source 1](../theory/weakly_nonlinear_planar_time_pilot.md) | Repeated stop of stage at old line 954; no new scientific result |
| 1193-1202 / Текущая остановка после адресной диагностики | [NLSP-D03](decisions.md#nlsp-d03), [NLSP-K03](knowledge.md#nlsp-k03) | [source 1](../theory/weakly_nonlinear_planar_time_pilot.md#14-адресная-диагностика-и-восстановление-вычислительного-пути-2026-10-07) | Repeated stop of stage at old line 936; no new scientific result |
| 1203-1212 / Историческая остановка после exact-time leading контроля | [NLSP-D04](decisions.md#nlsp-d04), [NLSP-K04](knowledge.md#nlsp-k04) | [source 1](../theory/planar_second_order_axial_response.md) | Repeated stop of stage at old line 918; no new scientific result |
| 1213-1222 / Историческая остановка после prepared initial-state gate | [NLSP-D05](decisions.md#nlsp-d05), [NLSP-K05](knowledge.md#nlsp-k05) | [source 1](../theory/planar_prepared_initial_state.md) | Repeated stop of stage at old line 899; no new scientific result |
| 1223-1230 / Историческая остановка после bounded precision/feasibility | [NLSP-D06](decisions.md#nlsp-d06), [NLSP-K06](knowledge.md#nlsp-k06) | [source 1](../theory/planar_prepared_initial_state.md#prepared-precision-feasibility) | Repeated stop of stage at old line 875; no new scientific result |
| 1231-1238 / Историческая остановка после bounded one-T1 feasibility | [NLSP-D07](decisions.md#nlsp-d07), [NLSP-K07](knowledge.md#nlsp-k07) | [source 1](../theory/planar_prepared_initial_state.md#prepared-one-period-feasibility) | Repeated stop of stage at old line 853; no new scientific result |
| 1239-1248 / Текущая остановка после bounded physical sanity checks | [NLSP-D08](decisions.md#nlsp-d08), [NLSP-K08](knowledge.md#nlsp-k08) | [source 1](../theory/weakly_nonlinear_planar_physical_sanity_checks.md) | Repeated stop of stage at old line 828; no new scientific result |
| 1249-1258 / Текущая остановка после FEM-readiness audit | [NLSP-K09](knowledge.md#nlsp-k09) | [source 1](../numerics/nlsp_3d_fem_environment_readiness.md) | Repeated stop of stage at old line 812; no new scientific result |
| 1259-1266 / Текущая остановка после FEM-1 | [NLSP-D09](decisions.md#nlsp-d09), [NLSP-K10](knowledge.md#nlsp-k10) | [source 1](../numerics/nlsp_linear_rectangular_3d_fem_validation.md) | Repeated stop of stage at old line 798; no new scientific result |
| 1267-1277 / Текущая остановка после FEM-1R | [NLSP-D10](decisions.md#nlsp-d10), [NLSP-K11](knowledge.md#nlsp-k11) | [source 1](../numerics/nlsp_linear_rectangular_3d_fem_validation.md#fem-1r-one-additional-mesh-refinement) | Repeated stop of stage at old line 778; no new scientific result |
| 1278-1289 / Текущая остановка после FEM-2R | [NLSP-D12](decisions.md#nlsp-d12), [NLSP-K13](knowledge.md#nlsp-k13) | [source 1](../numerics/nlsp_nonlinear_static_3d_fem_validation.md#fem-2r-controlled-continuation-after-input-serialization-failure) | Repeated stop of stage at old line 725; no new scientific result |
| 1290-1297 / Текущая остановка после FEM-3A | [NLSP-D13](decisions.md#nlsp-d13), [NLSP-K14](knowledge.md#nlsp-k14) | [source 1](../numerics/nlsp_nonlinear_dynamic_3d_fem_pilot.md) | Repeated stop of stage at old line 696; no new scientific result |

## Original SHA256 of editable files

| File | SHA256 before edits |
| --- | --- |
| `AGENTS.md` | `b210971eee789f78b8038026e58d8143f5fd8a6eaf2c6f30fdbb18aba6a98430` |
| `docs/project_rules.md` | `cbcb8b35cdf1a0ffd56a0ec509748fcf689d90c6e179c556ec8a4e7565d242c4` |
| `docs/theory/AGENTS.md` | `6dd4baa97568c16107e2c30c4a5f5f137f38d92a687538fac3d87e2c1a9c1cc9` |
| `docs/memory/README.md` | `94af20ac5a9a4037507f3fc568651fc90956da6af1456cfdea13c7d7ff21446f` |
| `docs/memory/current.md` | `ca7792f3694031bae94ff7120535fc861dc1e202331a9e3c9f4125db8c337ce3` |
| `README.md` | `18e4f02831a92d95521747a6ee3008659675dceb6cda2de4e42dfdab478d34ae` |
| `docs/research_index.md` | `a752d5d1f08957838e944efa96b884b91cd307a53e3f1fedf868c29a27b75351` |
| `docs/nonsymmetric_laminated_beams/README.md` | `bcfd166c91c209f4afed4a74f76729bd55eb5e73a1aeb996472933a691374433` |
| `docs/laminated_beams/README.md` | `9df4200792694ea468f2e04f95428cc34ed37cae5e98b6a8b6ca93ce61adf130` |
| `docs/numerics/README.md` | `6ee5b782af1acb3704bb091b5a5fbf4cbf018b7ff4ebc3c5fc141ca6cbc395d1` |
| `CHANGELOG.md` | `a81198d2f71738b73ff7dd391e05afa0c873d42765ce93899c201a15a11e3a56` |

Protected D/K SHA256 (unchanged criterion):

- `docs/memory/decisions.md`: `80764b5bc693f6ffd65c48422289969ee83fb1e6df45d4c832198486b0f9971e`.
- `docs/memory/knowledge.md`: `535e87ecd96a0aa8c94862f40483885c243ff90169e468b18b654797fc33cb47`.
