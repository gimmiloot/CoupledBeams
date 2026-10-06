# Память исследования EB/RLB и поворотного узла

## Назначение

Память — вторичный слой навигации и восстановления контекста: что известно,
при каких условиях, откуда это следует и где остановились. Она не заменяет
теорию, model-specific contracts и отчёты и не разрешает новые действия.
Отсутствие записи означает неполноту покрытия памяти, а не отсутствие
исследования в проекте.

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
