# Source Index

Индекс ниже собран под задачу `CoupledBeams`: связь стержней под углом, собственные частоты, условия сопряжения, граничные условия и сравнение аналитики с FEM.

Citation keys синхронизированы с `docs/literature/bibliography.bib`.

## `tao_2023_wave_coupled_beams`
- PDF: `docs/literature/pdf/Wave-basedin-planevibrationanalysisofmultiplecoupledbeamstructureswitharbitraryconnectionangleandelastic__boundaryrestraints.pdf`
- Тип: статья.
- Роль: основной.
- Что важно для CoupledBeams: wave-based постановка для нескольких сопряженных балок, произвольный угол соединения, упругие граничные закрепления, отражение и прохождение волн на стыке.
- Обозначения: волновые амплитуды, reflection/transmission matrices, connection angle, boundary stiffnesses; точные символы нужно поднимать уже по полному тексту статьи.
- Критично смотреть: abstract; разделы с выводом reflection/transmission matrices на угловом стыке и на упругой границе; параметрические графики по углу и жесткости закрепления; весь диапазон pp. 5250--5269.
- Метаданные: проверены по странице SAGE.

## `albarracin_2005_restrained_frames`
- PDF: `docs/literature/pdf/Vibrations_of_elastically_restrained_fra.pdf`
- Тип: статья.
- Роль: основной.
- Что важно для CoupledBeams: точная постановка для frame с упругими ограничениями на концах и в промежуточной точке; useful as reference for coupling axial and transverse motions through joint and boundary conditions.
- Обозначения: длины `l_1`, `l_2`; rotational springs `r_{1(1)}`, `r_{2(1)}`, `r_{1(2)}`; translational springs `t_{1(1)}`, `t_{2(1)}`, `t_{3(1)}`, `t_{4(1)}`, `t_{1(2)}`; frequency coefficients `\lambda_i`.
- Критично смотреть: раздел `Variational derivation of the boundary and eigenvalue problem`; раздел `Determination of the exact solutions`; таблицы и раздел `Results and discussion`; весь короткий текст pp. 467--476.
- Метаданные: проверены по ScienceDirect.

## `umar_2020_jib_crane_frame_vibration`
- PDF: `docs/literature/pdf/Vibration_Analysis_of_a_Jib_Crane_using.pdf`
- Тип: статья.
- Роль: вспомогательный.
- Что важно для CoupledBeams: пример сведения конструкции из двух стержневых элементов к frame-модели с Euler--Bernoulli beams и последующим расчётом собственных частот и форм.
- Обозначения: mast/jib splitting, assumed-mode amplitudes, natural frequencies and mode shapes; точные символы быстро не извлечены.
- Критично смотреть: sections with governing equations, assumed-mode reduction, and numerical evaluation; для этой работы practically важен весь текст pp. 71--80.
- Метаданные: `NEEDS_CHECK` только по DOI; остальные поля читаются из первой страницы.

## `ouisse_2003_connecting_angle`
- PDF: `docs/literature/pdf/2003JSVb.pdf`
- Тип: статья.
- Роль: основной.
- Что важно для CoupledBeams: напрямую про connecting angle и coupled beams/plates; полезно для понимания чувствительности собственных частот к углу сочленения и к типу сопряжения.
- Обозначения: угол соединения, coupled modes, чувствительность частот; локальный PDF — авторский manuscript, поэтому для точной нотации лучше сверять опубликованную версию.
- Критично смотреть: sections with formulation of the connecting-angle model, frequency sensitivity, and comparison cases; для проекта важна практически вся статья pp. 809--850.
- Замечание по корректности: в printed finite-beam matrix `T` формулы `(39)` знаки у beam-2 bending block в первых двух кинематических строках не согласуются с конечномерным аналогом условий `(10)--(11)`. В manuscript напечатано `(+,+,-,-)` для членов `\sin k_2L_2 \cos\alpha`, `\sinh k_2L_2 \cos\alpha`, `\sin k_2L_2 \sin\alpha`, `\sinh k_2L_2 \sin\alpha`, тогда как после переноса кинематических условий в левую часть получается противоположный sign pattern `(-,-,+,+)`. При переписывании determinant эти четыре элемента нужно проверять вручную.
- Метаданные: страница DOI и авторский архив подтверждают библиографию; локальный PDF не является издательским финальным layout.

## `quintana_2010_restrained_timoshenko_beams`
- PDF: `docs/literature/pdf/1.pdf`
- Тип: статья.
- Роль: основной.
- Что важно для CoupledBeams: Timoshenko beam с общими упругими ограничениями и промежуточными elastic constraints; полезно как reference для general boundary restraints beyond the simplest Euler--Bernoulli setting.
- Обозначения: beam length `l`, area `A`, inertia `I`, elastic moduli `E`, `G`, density `\rho`, translational/rotational restraint parameters, dimensionless frequency parameter; в аннотации отдельно отмечены Ritz and Lagrange multiplier methods.
- Критично смотреть: постановку с intermediate elastic constraints, sections с Ritz/Lagrange multipliers и новые результаты для end conditions; статья целиком pp. 117--125.
- Метаданные: проверены по репозиторной записи, указывающей на публикацию SAGE.

## `nikolai_1926_bent_rod_oscillations`
- PDF: `docs/literature/pdf/lfmo8.pdf`
- Тип: статья.
- Роль: вспомогательный.
- Что важно для CoupledBeams: исторический и очень близкий по геометрии источник про согнутый стержень; полезен как предшественник задачи о связи двух прямых участков через угол.
- Обозначения: в кратком описании Math-Net фигурируют два прямолинейных отрезка длины `2l`, соединённые под углом `2\delta`.
- Критично смотреть: вся статья pp. 77--88; особенно постановку геометрии согнутого стержня и вывод частотного условия.
- Метаданные: проверены по карточке Math-Net.

## `starshin_2015_rod_bend_vibrations`
- PDF: `docs/literature/pdf/vgsa2015-3-8.pdf`
- Тип: статья.
- Роль: вспомогательный.
- Что важно для CoupledBeams: модель стержня с изломом и контактной пружины; полезно как nearby engineering case for broken-geometry beam modeling.
- Обозначения: в OCR видны длины двух участков и постановка через контактную пружину; точная система символов требует ручной сверки с исходным PDF.
- Критично смотреть: весь короткий текст; особенно постановку контактной пружины, вывод уравнений свободных колебаний и описание геометрии излома.
- Метаданные: `NEEDS_CHECK` — автор, номер выпуска и страницы пока восстановлены только по OCR.

## `obradovic_2020_planar_serial_frames`
- PDF: `docs/literature/pdf/работа1.pdf`
- Тип: статья.
- Роль: основной.
- Что важно для CoupledBeams: planar serial frame structures, rigid joints, coupled axial/bending vibrations through boundary and joint conditions; близкий reference для структуры уравнений и переноса граничных условий.
- Обозначения: Euler--Bernoulli beams с variable cross-section and axially functionally graded material; transfer of boundary conditions to a Cauchy initial-value problem.
- Критично смотреть: постановку PDE/ODE после separation of variables, перенос boundary conditions, and numerical example; весь диапазон pp. 221--239.
- Метаданные: проверены по DOI Serbia.

## `ratazzi_2013_internal_elastic_hinge`
- PDF: `docs/literature/pdf/работа2.pdf`
- Тип: статья.
- Роль: основной.
- Что важно для CoupledBeams: exact free in-plane vibrations for two orthogonal beam members with an internal elastic hinge and elastic boundary conditions; очень близко к задаче о joint flexibility and connection angle.
- Обозначения: internal elastic hinge flexibility, boundary stiffnesses, Euler--Bernoulli assumption, Hamilton principle, separation of variables.
- Критично смотреть: постановку через Hamilton's principle, derivation of the exact frequency equation, and comparison with FEM/experiment; для короткой статьи полезен весь текст pp. 1--10.
- Метаданные: DOI и article ID извлечены надёжно.

## `shayna_2022_rod_systems_localized_features`
- PDF: `docs/literature/pdf/Диссертация_Шайна_Е.А..pdf`
- Тип: диссертация.
- Роль: вспомогательный.
- Что важно для CoupledBeams: математическое моделирование стержневых систем с локализованными особенностями; полезно для cases with singular supports, discontinuities, and nonclassical localized effects.
- Обозначения: по введению и оглавлению заметны generalized functions / spectral viewpoint / mixed boundary-value formulations; точная нотация зависит от выбранной модели и требует целевого чтения.
- Критично смотреть: главу 1 про математическую модель малых колебаний стержневой системы с особенностями (по оглавлению примерно pp. 13--66) и главу 2 про адаптацию метода конечных элементов для моделей 4-го порядка (примерно pp. 69--101).
- Метаданные: `NEEDS_CHECK` — библиография восстановлена по локальной титульной странице и OCR; перед внешней ссылкой лучше проверить официальный репозиторий ВГУ.

## `mai_2025_nonstationary_thinwalled_bodies`
- PDF: `docs/literature/pdf/02.-Dissertatsiya_M.-K.-CHien-_21.03.25_.pdf`
- Тип: диссертация.
- Роль: вспомогательный.
- Что важно для CoupledBeams: broader generalized-continuum and moment-elasticity context; может быть полезно, если проект уйдёт от классической балки к более богатым моделям тонкостенных и моментных упругих тел.
- Обозначения: моментные упругие среды, оболочки, пластины и стержни; из оглавления явно видны разделы по продольным колебаниям и изгибу моментного упругого стержня.
- Критично смотреть: главу 3 `Начально-краевые задачи моментных упругих пластин и стержней`, особенно 3.3--3.9 (примерно pp. 53--106), где есть уравнения движения стержней, продольные колебания и изгиб.
- Метаданные: `NEEDS_CHECK` — локальный PDF выглядит как диссертационная рукопись; имя автора и финальный статус публикации нужно сверить по официальной карточке МАИ/ВАК.

## `bauer_2025_coupled_rods`
- PDF: `docs/literature/pdf/Статья-Дорофеев-2025.pdf`
- Тип: статья.
- Роль: основной.
- Что важно для CoupledBeams: это, по сути, наиболее прямое попадание в тему репозитория: thin elastic rods coupled at an angle, dimensionless frequency equation, conditions of coupling, comparison with COMSOL, asymptotics for small angle and thickness parameter.
- Обозначения: dimensionless variables, coupling angle, thickness parameter, frequency equation, axial and transverse vibrations of rods; exact symbol names уже стоит поднимать из опубликованной версии.
- Критично смотреть: весь текст pp. 73--81; особенно введение, безразмерную постановку, условия сопряжения, частотное уравнение и асимптотики для малого угла.
- Метаданные: local PDF похож на draft build, но библиографические поля проверены по официальному выпуску журнала.
