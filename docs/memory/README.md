# Память исследования сопряжённых стержней

## Назначение

Память — вторичный слой навигации: что установлено, при каких условиях,
почему принято решение и где остановились. Она не заменяет каноническую
теорию, модельные соглашения, научные отчёты и данные; не разрешает новые действия.
Отсутствие записи означает неполное покрытие памяти, а не отсутствие исследования.

Покрыты независимые направления NSL (несимметричные ламинаты), RLB/EB–KV
(симметричная Reddy/FSDT и контроль Эйлера–Бернулли с поворотным узлом),
NLSP (нелинейная пространственная семиполевая модель), LONG (продольные
модели M–H/Timoshenko) и отдельный circular EB reviewer diagnostic.
Их выводы, ограничения и разрешения не переносятся между моделями автоматически.
Остальные направления доступны в [индексе исследований](../research_index.md).

## Как читать

Новый чат: [краткий текущий контекст](current.md) → паспорт выбранного
направления из таблицы ниже → прямо указанный итоговый научный отчёт.
Это основной маршрут; последовательное чтение всей истории не требуется.
Для причин решения и происхождения результата открыть соответствующие D/K.
Для воспроизведения — входной контракт и команды в отчёте, затем локальные данные.

## Файлы

- [current.md](current.md) — изменяемое действующее состояние, ограничения и остановки.
- [decisions.md](decisions.md) — дополняемая история научных решений и происхождения.
- [knowledge.md](knowledge.md) — дополняемая история результатов, оснований и ограничений.

Вводные заголовки и указатели «последние записи» внутри D/K отражают прежний
охват. Для актуальной навигации использовать этот README; ради косметики
исторические файлы не переписывать. Номера D и K не обязаны совпадать.

## Темы

| Направление | Паспорт и актуальное основание | Исторические основания |
| --- | --- | --- |
| NSL | [README](../nonsymmetric_laminated_beams/README.md), [NSL-D07](decisions.md#nsl-d07), [NSL-K07](knowledge.md#nsl-k07) | [NSL-D01](decisions.md#nsl-d01), [NSL-K01](knowledge.md#nsl-k01); устранение конкретных невязок: [NSL-K03](knowledge.md#nsl-k03) |
| RLB / EB–KV | [README](../laminated_beams/README.md), пауза [RLB-D22](decisions.md#rlb-d22), итог [RLB-K23](knowledge.md#rlb-k23) | [RLB-K01](knowledge.md#rlb-k01), [теория узла](../laminated_beams/inplane_rotational_spring_joint.md); рисунки [RLB-K24](knowledge.md#rlb-k24), [RLB-K25](knowledge.md#rlb-k25) |
| Circular EB reviewer diagnostic | [Отчёт](../laminated_beams/circular_eb_rotational_spring_rigid_limit.md), [RLB-D25](decisions.md#rlb-d25), [RLB-K26](knowledge.md#rlb-k26) | [EB-JOINT-K01](knowledge.md#eb-joint-k01), [EB-JOINT-K02](knowledge.md#eb-joint-k02) |
| NLSP | [Паспорт](../numerics/README.md#nlsp), шестиполевой OUT-рост [NLSP-D35](decisions.md#nlsp-d35) / [NLSP-K36](knowledge.md#nlsp-k36) | [Семиполевая модель](../theory/weakly_nonlinear_spatial_rod.md), [NLSP-D01](decisions.md#nlsp-d01) / [NLSP-K01](knowledge.md#nlsp-k01); семиполевой итог [NLSP-D32](decisions.md#nlsp-d32) / [NLSP-K33](knowledge.md#nlsp-k33); шестиполевая база [NLSP-D34](decisions.md#nlsp-d34) / [NLSP-K35](knowledge.md#nlsp-k35) |
| LONG | [HS-LA4: PASS, section 12](../theory/isotropic_two_springs_horizontal_leading_spectral_structure.md), [LONG-D43](decisions.md#long-d43) / [LONG-K42](knowledge.md#long-k42); [HS-LA3: PASS, section 11](../theory/isotropic_two_springs_horizontal_leading_spectral_structure.md), [LONG-D42](decisions.md#long-d42) / [LONG-K41](knowledge.md#long-k41); [HS-LA2: PASS](../theory/isotropic_two_springs_horizontal_leading_spectral_structure.md), [LONG-D41](decisions.md#long-d41) / [LONG-K40](knowledge.md#long-k40); [HS-L1: PASS](../theory/isotropic_two_springs_horizontal_leading.md), [LONG-D40](decisions.md#long-d40) / [LONG-K39](knowledge.md#long-k39); [HS-P1: PASS](../theory/isotropic_two_springs_horizontal_pilot.md), [LONG-D39](decisions.md#long-d39) / [LONG-K38](knowledge.md#long-k38); [one horizontal spring](../theory/isotropic_two_springs_horizontal_spring.md), [LONG-D38](decisions.md#long-d38) / [LONG-K37](knowledge.md#long-k37); [angular transition](../theory/isotropic_two_springs_external_axial_transition_mechanism.md), [LONG-D37](decisions.md#long-d37) / [LONG-K36](knowledge.md#long-k36); [weak axial support: PASS](../theory/isotropic_two_springs_external_axial_small_stiffness_transition.md), [LONG-D36](decisions.md#long-d36) / [LONG-K35](knowledge.md#long-k35); [axial-support 90%/95% ranges: PASS](../theory/isotropic_two_springs_external_axial_high_frequency_ranges.md), [LONG-D35](decisions.md#long-d35) / [LONG-K34](knowledge.md#long-k34); [90%/95% primary ranges: PASS](../theory/isotropic_two_springs_high_frequency_angle_ranges_90_95.md), [LONG-D34](decisions.md#long-d34) / [LONG-K33](knowledge.md#long-k33); [95%/98% angle ranges: PASS](../theory/isotropic_two_springs_high_frequency_angle_ranges.md), [LONG-D33](decisions.md#long-d33) / [LONG-K32](knowledge.md#long-k32); [free outer rotation: PASS](../theory/isotropic_two_springs_external_rotational_replacement.md), [LONG-D32](decisions.md#long-d32) / [LONG-K31](knowledge.md#long-k31); [transverse leading validation: PASS](../theory/isotropic_two_springs_external_transverse_leading_validation.md), [LONG-D31](decisions.md#long-d31) / [LONG-K30](knowledge.md#long-k30); [strong residual: DIAGNOSTIC_PASS](../theory/isotropic_two_springs_transverse_leading_diagnostic.md), [LONG-D30](decisions.md#long-d30) / [LONG-K29](knowledge.md#long-k29); [kappa_w: full PASS / leading STOP](../theory/isotropic_two_springs_external_transverse_spring.md), [LONG-D29](decisions.md#long-d29) / [LONG-K28](knowledge.md#long-k28); [chi_u](../theory/isotropic_two_springs_external_axial_spring.md), [LONG-D28](decisions.md#long-d28) / [LONG-K27](knowledge.md#long-k27); [Паспорт](../research_index.md#long), [F_full/F_lim](../theory/isotropic_two_springs_thickness_asymmetry_leading_comparison.md), [LONG-D27](decisions.md#long-d27) / [LONG-K26](knowledge.md#long-k26); [A1/A2](../theory/isotropic_two_springs_thickness_asymmetry.md), [LONG-D26](decisions.md#long-d26) / [LONG-K25](knowledge.md#long-k25); [L1/L2](../theory/isotropic_two_springs_length_asymmetry.md), [LONG-D25](decisions.md#long-d25) / [LONG-K24](knowledge.md#long-k24); [запаздывающий рост: архив](../theory/archived_delayed_frequency_growth.md), [LONG-D24](decisions.md#long-d24) / [LONG-K23](knowledge.md#long-k23); [научная сводка](../theory/isotropic_two_springs_scientific_synthesis.md) и [аудит](../theory/isotropic_two_springs_open_questions.md), [LONG-D23](decisions.md#long-d23) / [LONG-K22](knowledge.md#long-k22); [проверка ведущей модели](../theory/isotropic_two_springs_leading_validation.md), [LONG-D22](decisions.md#long-d22) / [LONG-K21](knowledge.md#long-k21); [ведущая асимптотическая задача](../theory/isotropic_two_springs_angle_asymptotics.md), [LONG-D21](decisions.md#long-d21) / [LONG-K20](knowledge.md#long-k20); [толщина и карта угла](../theory/isotropic_two_springs_thickness_grid.md), [LONG-D20](decisions.md#long-d20) / [LONG-K19](knowledge.md#long-k19); продольная податливость [LONG-D19](decisions.md#long-d19) / [LONG-K18](knowledge.md#long-k18) | Масштабирование ε [LONG-D18](decisions.md#long-d18) / [LONG-K17](knowledge.md#long-k17); угловая вариация [LONG-D17](decisions.md#long-d17) / [LONG-K16](knowledge.md#long-k16); формы и энергии [LONG-D16](decisions.md#long-d16) / [LONG-K15](knowledge.md#long-k15); сетка 5×5 [LONG-D15](decisions.md#long-d15) / [LONG-K14](knowledge.md#long-k14); локальное уточнение [LONG-D14](decisions.md#long-d14) / [LONG-K13](knowledge.md#long-k13); относительный наклон [LONG-D13](decisions.md#long-d13) / [LONG-K12](knowledge.md#long-k12); ранний рост [LONG-D12](decisions.md#long-d12) / [LONG-K11](knowledge.md#long-k11); двумерная сетка [LONG-D11](decisions.md#long-d11) / [LONG-K10](knowledge.md#long-k10); критерий плато [LONG-D10](decisions.md#long-d10) / [LONG-K09](knowledge.md#long-k09); геометрическая диагностика [LONG-D09](decisions.md#long-d09) / [LONG-K08](knowledge.md#long-k08); конечные карты [LONG-D08](decisions.md#long-d08) / [LONG-K07](knowledge.md#long-k07); четыре предела [LONG-D07](decisions.md#long-d07) / [LONG-K06](knowledge.md#long-k06); RL-верификация [LONG-D06](decisions.md#long-d06) / [LONG-K05](knowledge.md#long-k05); выбор модели [LONG-K01](knowledge.md#long-k01), закрытие [LONG-D02](decisions.md#long-d02); отдельные диагностики [LONG-K02](knowledge.md#long-k02), [LONG-K03](knowledge.md#long-k03), пространственное разделение [LONG-D05](decisions.md#long-d05) / [LONG-K04](knowledge.md#long-k04) |

Таблица даёт вход в направление, а не заменяет область действия отдельной записи.
[Индекс результатов](../results_index.md) указывает вычислительные материалы.

## Иерархия источников

Сначала установить физическую модель по паспорту направления.
Для модели базового `equations.tex` его прежний приоритет сохраняется.
Отдельные NSL/RLB/NLSP/LONG используют собственную указанную каноническую
теорию, допущения, координаты и реализацию. Иерархия внутри модели
не меняется; см. [AGENTS.md](../../AGENTS.md) и [правила проекта](../project_rules.md).

Научные отчёты сохраняют выводы и квалификации, код/тесты/данные — их
проверяемые основания, память — краткие утверждения со ссылками.
Публикация учитывается вместе с зафиксированными предупреждениями о печати.
Расхождение источников отмечается явно, не исправляется пересказом.
Новизна даты сама по себе не делает результат каноническим для другой модели.

## Git, локальные доказательства и история

- **Git-visible scientific context:** модель, параметры области исследования,
  принятые результаты, ограничения, причина остановки и ссылки в Markdown/TeX.
  Новый читатель должен понимать научный итог без игнорируемых файлов.
- **Local scientific evidence:** CSV, NPZ, графики, кеши и manifests в `results/`.
  Ссылки на них обозначают локальные материалы; Git не обязан содержать массивы.
  Наличие файла локально не доказывает доступность через GitHub. Некоторые
  команды требуют уже сохранённых bundles; отсутствие данных не разрешает пересчёт.
- **Historical records:** D/K, научные отчёты, исходные квалификации и результаты.
  Их сохранность регулируется [archive policy](../archive_policy.md).

При передаче контекста проверять отслеживаемость документов, незакоммиченные
изменения и доступную версию `origin/main`; указывать, проверена ли удалённая
версия или только локальная remote-tracking ссылка. Никакой автоматической
публикации, миграции результатов или настройки резервного копирования это не требует.

## Соглашения о новых записях

Компактные ориентиры, не обязательная бюрократическая схема:

- **Decision:** решение, причина, область, условие остановки/возврата, ссылка на результат.
- **Knowledge:** утверждение, область, основание, ограничения, следствие; связь с прежним результатом, когда существенна.

Указывать дату и происхождение; различать разрешение, гипотезу и результат.
Память хранит утверждение, условия и ссылку. Большие таблицы частот, полные
невязки и перечни тестов остаются в отчёте. Навигационная или косметическая
правка не создаёт научного D/K record. Старые записи не переоформляются.

## Правило роста

Каждый документ обновляется только при изменении информации, за которую он
отвечает: [единая таблица ответственности](../project_rules.md#documentation-responsibilities).
Новый результат не требует повторения чисел во всех README и индексах.

`current.md` можно переписывать: ориентир 100–200 строк, точность важнее лимита.
Оставлять текущую модель, итог, открытые вопросы, остановку и прямые ссылки.
Не добавлять в конец или начало полные прежние сводки. Перед сокращением
прочитать удаляемое содержание и установить для каждого существенного
утверждения источник, включая отказы, контрпримеры, границы и причины остановки.
Единственное описание результата сохраняется до появления допустимого источника;
при конфликте с неизменностью требуется явная квалификация, а не молчаливое удаление.

D/K остаются append-only с устойчивыми номерами и якорями. Исторические записи
меняются только для явной опечатки или битой ссылки, если конкретная задача
и проверки сохранности это допускают. В переработке 2026-10-20 они полностью
неизменны. Данные некоторых bundle проверяют байты прежних D/K; механическая
миграция файлов и исправление заголовков без проверки зависимостей недопустимы.

Текущая сводка различает подтверждённое, нерешённое, устранённое последующим
этапом и неисследованное обычным текстом, без нового набора статусов.
Например, NSL-2 сохраняет PARTIAL_PASS, а NSL-2R отдельно закрывает его
конкретные концевые невязки. Исследовательский PASS не становится общей валидацией.
Разделение D/K сейчас не требуется; возможная миграция — отдельная задача.
