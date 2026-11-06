# Текущий контекст

Навигация обновлена **2026-11-02**; научный срез — **2026-10-28**,
с отдельным дополнением RR-префикса **2026-11-02** ниже.
История остаётся в [D](decisions.md) и [K](knowledge.md), правила — в
[README памяти](README.md), направления — в [индексе](../research_index.md).
Сокращение сводки не меняет результаты, прежние статусы и разрешения.

## NSL — Nonsymmetric Laminated Beams

[Паспорт](../nonsymmetric_laminated_beams/README.md).
Плоская FSDT-модель PROJECT_REDDY_COUPLED с изгибно-продольной связью
и полной модальной массой; отдельна от симметричной RLB.
Теория: [одиночный стержень](../nonsymmetric_laminated_beams/single_beam_theory.md),
[отражение, оси и узел](../nonsymmetric_laminated_beams/coupled_joint_theory.md).

- **Итог:** [NSL-VIS-1](../nonsymmetric_laminated_beams/comparative_spectral_curves.md),
  [D07](decisions.md#nsl-d07)/[K07](knowledge.md#nsl-k07) — PASS для
  исследовательских Ritz-карт BASE/MIRROR/SAME_LOCAL. Это sorted-позиции,
  без доказательства полноты непрерывного спектра на всей сетке.
- **Подтверждено:** перенос оси, отражение и защищённый класс MIRROR;
  SAME_LOCAL не получает такую защиту автоматически. [Два комплексных
  подтверждения при 15°](../nonsymmetric_laminated_beams/nsl3a_angle_complex_confirmation.md)
  не превращают остальные упругие предикторы в фактическое затухание
  и не доказывают сильную диссипацию.
- **Устранено адресно:** [NSL-2R](../nonsymmetric_laminated_beams/nsl2r_numerical_diagnostics.md),
  [K03](knowledge.md#nsl-k03), закрыл концевые невязки без изменения частот.
  Исторический NSL-2 PARTIAL_PASS сохранён. NSL-1 PARTIAL_PASS относится
  к [печати Bachoo](../nonsymmetric_laminated_beams/source_comparison.md),
  а не к открытому численному отказу проектной модели.
- **Открыто / STOP после VIS-1:** часть соответствий форм, продолжение ветвей,
  crossing/veering, NSL-3B и FRF; новая вязкость и оптимизация не разрешены.

## RLB / EB–KV — симметричные ламинаты и поворотное соединение

[Паспорт](../laminated_beams/README.md): симметричная Reddy/FSDT с раздельными
продольной и изгибной подсистемами; EB — контроль.
[Пружина](../laminated_beams/inplane_rotational_spring_joint.md) и
[KV](../laminated_beams/inplane_kelvin_voigt_joint_theory.tex) уже исследованы;
вязкость сосредоточена только в узле.

- **Итог:** [сводка руководителю](../laminated_beams/inplane_kelvin_voigt_research_status_for_supervisor.md),
  [D22](decisions.md#rlb-d22)/[K23](knowledge.md#rlb-k23):
  структурная защита, упругие предикторы и пять комплексных EB/RLB-пар.
  Состояние — **PAUSED_FOR_SUPERVISOR_DIRECTION**; рисунки не снимают паузу.
- [Архитектура восстановления](../laminated_beams/inplane_kelvin_voigt_solver_architecture.md)
  и [завершение A/C](../laminated_beams/inplane_kelvin_voigt_targeted_weak_damping_completion.md)
  не отменяют прежние отказы K16, отдельный singular flag и source/Ritz-квалификации.
- Отдельный [circular EB diagnostic](../laminated_beams/circular_eb_rotational_spring_rigid_limit.md),
  [D25](decisions.md#rlb-d25)/[K26](knowledge.md#rlb-k26): шесть seeds доведены
  до RIGID; позиции 11–12, прямое назначение 1→10 и немонотонность seed 06
  остаются квалифицированными. Реальная жёсткость и FEM-калибровка не установлены.
- **Открыто:** универсальность пилотов и перенос на LONG/NLSP. Новые спектры,
  вязкости и калибровка требуют отдельного решения.

## NLSP — семиполевая модель и шестиполевое сравнение

[Паспорт](../numerics/README.md#nlsp), [семиполевая теория](../theory/weakly_nonlinear_spatial_rod.md):
шесть движений сечения и c, принятое действие V0, кубические балансы,
переменная инерция; идеальный жёсткий узел с c₁=c₂. Это не NSL/KV-узел.

- **Последний итог:** [шестиполевой OUT-этап](../numerics/nlsp_six_vs_seven_out_stability.md),
  [D35](decisions.md#nlsp-d35)/[K36](knowledge.md#nlsp-k36). При p14/A=.005
  четыре Hill-контроля подтверждают положительный рост соответствующего режима
  и отрицательный multiplier; σ ниже семиполевой на .2961%. Разность
  разрешена по наблюдаемому индикатору, не строгой границе ошибки. **STOP перед временным подтверждением и p-уточнением.**
- Основания: [линейное сравнение](../numerics/nlsp_six_vs_seven_linear_comparison.md)
  сохраняет EA₀ и u первого порядка, сопоставляет PL2/PL4 и совпадающие OUT-блоки;
  [шестиполевая орбита](../numerics/nlsp_six_vs_seven_planar_orbit.md),
  [D34](decisions.md#nlsp-d34)/[K35](knowledge.md#nlsp-k35), проверяет действие
  степени ≤4, полный residual/tangents и HB11/HB13 типа II, Ψ=0.
- **Семиполевой итог:** [проверки](../numerics/nlsp_beta15_pre_supervisor_checks.md),
  [D32](decisions.md#nlsp-d32)/[K33](knowledge.md#nlsp-k33),
  [научный отчёт](../reports/nlsp_beta15_parametric_instability/README.md):
  PRE_SUPERVISOR_CHECKS_COMPLETE_WITH_PL_PARTIAL. OUT-рост при β=15°, μ=.1,
  A=.005 подтверждён p12 и воспроизведён p14. OUT2-отклонение без скоростей
  возбуждает рост, но прогиб за 5T составляет лишь .2864% плоской амплитуды.
- **Не закрыто:** объяснение PL neutral residual не отменяет строгий PARTIAL.
  Отсутствие найденного ненейтрального PL-роста не доказывает общей устойчивости;
  два p не доказывают полной пространственной сходимости.
- [Прямой стержень / FEM](../numerics/nlsp_straight_rod_3d_fem_verification_summary.md):
  FEM-2R/FEM-3AR устранили отдельные ошибки ввода/вывода, сохранив прежние отказы.
  [Пространственное сравнение](../numerics/nlsp_spatial_nonlinear_3d_fem_verification.md)
  остаётся NUMERICAL_PARTIAL для малых поправок, c/скоростей и нативной энергии;
  c_eff — proxy. [Prepared-state](../theory/planar_prepared_initial_state.md)
  не закрыл strict float64/все пространственные критерии.
- Отрицательный screening и [Construction C](../numerics/nlsp_construction_C_longer_sensitivity.md)
  не заменяются результатом β15. Одномодовое приближение около 1:3 не прошло
  критерий; ветвь типа II не доказывает результат для движения типа I.
- **Открыто / обсуждение с руководителем:** DEFERRED_TWO_PARAMETER_ACTION_EXPANSION,
  физика идеального узла/c₁=c₂, независимая 3D-проверка угловой конструкции,
  p-сходимость, тип I, насыщение и наблюдаемость. Новые ODE/FEM/Floquet не разрешены.

## LONG — разные продольные модели и податливость связей

### RL_FULL / Timoshenko: угловая чувствительность и горизонтальная опора

[Каноническая двухпружинная модель](../theory/isotropic_two_springs_rl_timoshenko.md):
μ=.01, b/h=4, J_RL=ρν²I_P, исходная ε=.005. Это отдельная линия от M–H ниже.

Отдельная [RR-проверка 8+guard9](../numerics/rl_full_rr_prefix_verification.md):
PASS для четырёх регрессий 2+guard3 и шести углов при μ=.01,τ=0,ε=.005;
полные формы, count, Ritz и прямой предел. Старые исходники/архивы сохранены.
[C0/C1: индивидуальные сечения и ограниченный μ-пилот](../numerics/rl_full_rr_geometry_verification.md):
PASS, 20 новых K8-точек с отражением/прямым пределом/Ritz; C1 использует также
пять принятых A–B-точек. Только μ=0,.01,.1, пять β и четыре τ≠0 контроля.
Индивидуальные ветви между редкими углами не установлены; события записаны без уточнений.
Этап D по ε, полный τ-sweep и 3D-проверка не выполнены; старые отчёты и D/K сохранены.

[Настраиваемые карты Λ(β), Λ(μ), Λ(τ)](../numerics/rl_full_rr_frequency_maps.md):
fast_plot 8+guard9, пять предварительных сеток PASS (44 точки, 13 архивных включений).
Один конфигурируемый CLI; CSV и перерисовка без решений. ε=.005, sorted_positions;
tracking и научная интерпретация новых карт не выполнялись.

[Квадратные сечения и расширенные карты](../numerics/rl_full_rr_square_maps.md):
a=1, ε=.005; контроль a=4/прямого предела/отражения/Ritz PASS, затем все
191 точки пяти сеток PASS. Старый физический код сохранён; тот же CLI,
настраиваемое a и геометрическая отметка L/h≥10. Это не 3D-валидация.

[Научная сводка](../theory/isotropic_two_springs_scientific_synthesis.md) и
[критический аудит](../theory/isotropic_two_springs_open_questions.md),
[D23](decisions.md#long-d23)/[K22](knowledge.md#long-k22), — вход в механизм
роста, работу реакций, критерии перехода и ограничения.

**Последний полный расчёт горизонтальной опоры — [HS-B](../theory/isotropic_two_springs_horizontal_geometry_map.md)**,
[D47](decisions.md#long-d47)/[K46](knowledge.md#long-k46), PASS: 30 новых
геометрических вариантов и четыре эталона [HS-A](../theory/isotropic_two_springs_horizontal_angle_map.md).
Все сеточные максимумы при 90°, обе нормировки совпадают. Изменения ε и
знака μ заметно сдвигают пороги; эффект aspect мал, его разность для A2
численно не разрешена. Комбинированный прогноз точен до .111% для A1,
но ошибки A2 достигают 74%; у всех A4 нет ведущего корня 95%.
Полные пороги существуют, в том числе поздние. Непрерывная монотонность
и универсальная область применимости предельных моделей не доказаны.
HS-B завершён; новые карты, геометрии и асимптотики автоматически не начинать.

**Предыдущая завершённая линия — HS-LA6**, [§14 канонического отчёта](../theory/isotropic_two_springs_horizontal_leading_spectral_structure.md),
[D45](decisions.md#long-d45)/[K44](knowledge.md#long-k44):
PASS для четырёх порогов H5/H6 при μ=.01, ε=.005, κ_j=.1, q=.90/.95.
Ошибки прямого прогноза +.021–.052%, двухчленного +.069–.123%;
уровень HS-LA5 и коэффициенты зафиксированы до чтения полных порогов.
Частично компенсируются отрицательная ошибка кривой и положительная ошибка уровня.
**HS-L1–HS-LA6 завершена; STOP.** Универсальная точность, непрерывный глобальный
максимум, строгий остаток усечения и перенос на другие параметры не доказаны.

Маршруты к предшествующим результатам (подробные числа и расходы — в источниках):

| Блок | Канонические источники и сохраняемая квалификация |
| --- | --- |
| Горизонтальная опора и пилот | [Постановка](../theory/isotropic_two_springs_horizontal_spring.md), [HS-P1](../theory/isotropic_two_springs_horizontal_pilot.md), [HS-L1](../theory/isotropic_two_springs_horizontal_leading.md): нет внешней направляющей/поворотной пружины; точные нули сохранены. Исходный PARTIAL_PASS связан с прежним хеш-контрактом. Максимумы HS-P1 сеточные; малая частотная ошибка HS-L1 не означала малой ошибки порогов |
| Спектральная структура, насыщение, внутренний поворот | [HS-LA2–HS-LA4, §§1–12](../theory/isotropic_two_springs_horizontal_leading_spectral_structure.md): коэффициенты без подгонки; отказ binary64 сохранён и отдельно проверен высокой точностью тех же матриц. Континуальный хвост и глобальная монотонность порогов по κ не доказаны |
| Прямой угол | [HS-LA5, §13](../theory/isotropic_two_springs_horizontal_leading_spectral_structure.md): различие массовых форм объясняет основную часть разницы пределов; конечный остаток сохраняется, отдельные поправки полной модели не выделены |
| Геометрическая асимметрия | [Длины](../theory/isotropic_two_springs_length_asymmetry.md), [толщины](../theory/isotropic_two_springs_thickness_asymmetry.md), [проверка ведущей функции](../theory/isotropic_two_springs_thickness_asymmetry_leading_comparison.md), D25–D27/K24–K26: конечные сравнения завершены; знак асимметрии существенен при заданной внешней опоре. Сходимость при ε→0 и идентичность близких форм не доказаны |
| Продольная внешняя связь | [Опора](../theory/isotropic_two_springs_external_axial_spring.md), [слабая жёсткость](../theory/isotropic_two_springs_external_axial_small_stiffness_transition.md), [переход](../theory/isotropic_two_springs_external_axial_transition_mechanism.md): свободный предел без механизма; общий ведущий оператор не делает внешнюю и узловую пружины эквивалентными при конечном угле. Два знака малой кривизны UNRESOLVED; непрерывная траектория максимума не установлена |
| Поперечная опора и внешний поворот | [Полная PASS / ведущая STOP](../theory/isotropic_two_springs_external_transverse_spring.md) → [диагностика квадратуры](../theory/isotropic_two_springs_transverse_leading_diagnostic.md) → [отдельное завершение](../theory/isotropic_two_springs_external_transverse_leading_validation.md), D29–D31/K28–K30. Допуск не ослаблен; исторический STOP сохранён. [Свободный внешний поворот](../theory/isotropic_two_springs_external_rotational_replacement.md) — отдельная конструкция, D32/K31 |
| Высокочастотные уровни | [Основные 90/95%](../theory/isotropic_two_springs_high_frequency_angle_ranges_90_95.md), [95/98%](../theory/isotropic_two_springs_high_frequency_angle_ranges.md), [продольная опора до χ_u=.01](../theory/isotropic_two_springs_external_axial_high_frequency_ranges.md): только проверенные узлы; почти равные максимумы при χ_u≥.01 не переносятся на более слабые опоры. Нулевой вход в диапазон не означает большой абсолютной частоты |

**Общие незакрытые границы RL-линии:**
- [Малоугловой переход](../theory/isotropic_two_springs_angle_asymptotics.md)
  формален; свойства ведущей задачи не доказывают общей сходимости полной.
  [Проверка частот](../theory/isotropic_two_springs_leading_validation.md)
  не закрывает сходимость истинной кривизны: архивный оцениватель зависит от окна.
- Закон β_pl≈εB требует контроля производных, порога и хвоста до 90°;
  точное подобие и монотонность B по κ не следуют из min–max.
  [Карта толщин](../theory/isotropic_two_springs_thickness_grid.md) сохраняет
  чувствительность к схеме/опорному окну. Геометрическое колено, β_tol,
  β_pl и диапазоны 90/95% различаются; универсального допуска нет.
- [Опыт с узловой податливостью](../theory/isotropic_two_springs_longitudinal_slip.md)
  поддерживает роль передачи продольного усилия в двух парах;
  изменение толщины само по себе не изолирует EA. Энергии и формы могут
  меняться после рабочего перехода; малый зазор не доказывает veering.
- [Запаздывающий рост](../theory/archived_delayed_frequency_growth.md),
  [D24](decisions.md#long-d24)/[K23](knowledge.md#long-k23), сохраняет
  контрпримеры универсальности начального окна; физические модели не смешиваются.
  Полная история D06–D45/K05–K44 доступна через [тематический индекс](../research_index.md#long)
  и [таблицу памяти](README.md#темы), без повторения этапов здесь.
- Реальный узел/3D, универсальные критерии, новые карты, оптимизация,
  MAC/FEM и расширение асимптотики не разрешены автоматически.

### M–H / Timoshenko: выбранная продольная теория и отдельные KV-контроли

[Паспорт и источники](../research_index.md#long): reduced Jang-type
project_jang_reduced_rectangular, axial/contraction + Timoshenko.
Коэффициент 5/6 — проектный, не восстановленный из публикации Jang.
[LONG-D02](decisions.md#long-d02) закрывает выбор **в принятой 1D области**;
[K01](knowledge.md#long-k01) сохраняет finite-spectrum PASS и hierarchy PARTIAL
из-за различающихся внешних условий сокращения; M–H не является точным 3D-эталоном.

Отдельно завершены [условия на c](../theory/contraction_joint_sensitivity.md),
[сравнение продольных теорий при KV](../theory/longitudinal_kelvin_voigt_comparison.md)
и [одноосный пространственный split](../theory/uniaxial_kelvin_voigt_split.md),
[D05](decisions.md#long-d05)/[K04](knowledge.md#long-k04).
Неизменность OUT доказана для конкретных изотропии, заделок, FREE и отсутствия
эксцентриситета. Common-DOF c₁=c₂ и FREE — разные идеализации;
это не валидация реального узла или нелинейной NLSP.
**STOP:** новые области применимости, close-pair, пространственная податливость
и расширенные карты требуют отдельного задания.

## Другие сохранённые направления

Отсутствие D/K означает неполное покрытие памяти, а не отсутствие исследования.

| Направление | Состояние и вход |
| --- | --- |
| Базовая EB и ветви | Замороженные [уравнения](../theory/equations.tex), [код](../../src/my_project/analytic/formulas.py), [правила](../project_rules.md); sorted position не равна branch identity |
| Разные толщины, EB/Timoshenko | [Паспорт](../thickness_mismatch/README.md): safe-prefix закрыт; lower-envelope контрпримеры и отсутствие выгоды Rule S сохранены. Прерванная epsilon-карта не готова для статьи |
| Анизотропные стержни Ярцева | [Паспорт](../anisotropic_rods/README.md#current-boundary): 1D UT COMPLETE с PARTIAL_PASS; sampled AP-2 не даёт критического угла. 3D-A0 BLOCKED_MATERIAL_AND_SOLVER относится только к этой модели |
| Veering / локальная теория | [Критерии](../veering/README.md), [теория](../theory/veering_bridge/README.md): условные результаты не доказывают прикладной случай или veering на новом рисунке |
| Исторический OUT EB/кручение | [Теория](../theory/out_of_plane_eb_torsion.md); иные теория/заделки не дают same-BC эталона LONG или NLSP |

## Продолжение

Сначала выбрать модель и её канонический источник; выводы и разрешения
между направлениями не переносятся. PASS/PARTIAL относится к конкретному этапу.
[Индекс результатов](../results_index.md) ведёт к локальным доказательствам;
отсутствие массива в Git не отменяет отчёт и не разрешает пересчёт.
Ни исторический STOP, ни план не являются разрешением новой работы.
