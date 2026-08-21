# Source Index

Индекс ниже собран под задачу `CoupledBeams`: связь стержней под углом, собственные частоты, условия сопряжения, граничные условия и сравнение аналитики с FEM.

Citation keys синхронизированы с `docs/literature/bibliography.bib`.

См. также отдельную заметку по новым источникам для Timoshenko theory, shear
coefficient и circular rods:
`docs/literature/timoshenko_shear_sources.md`.

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

## `li_2012_two_beams_arbitrary_angle`
- PDF: `docs/literature/pdf/s0894-91662960007-x.pdf`
- Тип: статья.
- Роль: основной, прямой источник по двум сопряженным балкам.
- Что важно для CoupledBeams: free vibrations of two beams elastically coupled at an arbitrary angle; directly relevant to coupling-angle geometry, elastic coupling at the joint, and comparison language for the two-beam determinant line.
- Обозначения: coupled beams, arbitrary angle, elastic coupling, free vibrations; exact symbols should be checked against the PDF before importing notation.
- Критично смотреть: abstract/introduction, model formulation for the two elastically coupled beams, boundary/joint conditions, numerical examples over connection angle; pp. 61--72.
- Метаданные: recovered from local PDF XMP; DOI `10.1016/S0894-9166(12)60007-X`.

## `berkolaiko_2022_3d_elastic_beam_frames`
- PDF: `docs/literature/pdf/2104.01275v2.pdf`
- Тип: arXiv/preprint and journal article.
- Роль: вспомогательный для 3D frame/joint formulation.
- Что важно для CoupledBeams: rigorous variational/differential formulation of 3D elastic beam frames with rigid joint conditions; useful if the project later audits joint conditions or extends from planar rods to spatial frames.
- Обозначения: 3D elastic beam frames, rigid joint conditions, variational formulation, differential formulation; do not import notation into the current verified planar determinant without an explicit theory audit.
- Критично смотреть: abstract/introduction and sections deriving rigid joint conditions in variational and differential form.
- Замечание по применимости: spatial/frame formulation, not a replacement for the current Euler--Bernoulli two-rod determinant.
- Метаданные: local PDF is arXiv:2104.01275v2; DOI `10.1111/sapm.12485`.

## `perkins1986`
- PDF: `docs/literature/pdf/perkins1986.pdf`
- Тип: статья.
- Роль: основной общий источник по различению crossing и curve veering в задачах на собственные значения.
- Что важно для CoupledBeams: показывает, что veering возникает и в точной непрерывной задаче, а не только как артефакт дискретизации; формулирует критерии различения crossing и veering и связывает сближение собственных значений с быстрым изменением собственных векторов.
- Обозначения: eigenvalue loci, curve veering, crossing, continuous and discretized eigenvalue problems; обозначения общей спектральной задачи не следует переносить в проектный determinant.
- Критично смотреть: abstract/introduction, общую постановку задачи на собственные значения, критерии crossing/veering и непрерывный пример.
- Замечание по применимости: источник близок по спектральному механизму, но не по геометрии жёстко соединённых стержней; сам по себе не доказывает наличие veering в текущем `mu`-sweep.
- Метаданные: подтверждены по локальному журнальному PDF и издательской записи; DOI `10.1016/0022-460X(86)90191-4`. Сводная оценка: [literature assessment](../veering/literature_assessment.md#perkins1986).

## `pierre1988`
- PDF: `docs/literature/pdf/pierre1988.pdf`
- Тип: статья.
- Роль: основной общий источник по связи mode localization и eigenvalue-loci veering.
- Что важно для CoupledBeams: связывает малые нерегулярности в почти периодических слабосвязанных системах с сильной локализацией форм и veering близких собственных значений; рассматривает системы связанных осцилляторов и многопролётную балку.
- Обозначения: disordered structures, mode localization, eigenvalue loci veering, nearly periodic and weakly coupled systems; параметр disorder не является проектным `mu`.
- Критично смотреть: abstract/introduction, perturbation analysis, примеры связанных осцилляторов и многопролётной балки, выводы о совместном появлении localization и veering.
- Замечание по применимости: полезен по механизму и терминологии, но почти периодическая геометрия и слабая связь не являются прямым аналогом текущего жёсткого углового стыка.
- Метаданные: подтверждены по локальному журнальному PDF и издательской записи; DOI `10.1016/0022-460X(88)90226-X`. Сводная оценка: [literature assessment](../veering/literature_assessment.md#pierre1988).

## `liu2002`
- PDF: `docs/literature/pdf/liu2002.pdf`
- Тип: статья.
- Роль: основной методический источник по спектральным производным около veering/localization.
- Что важно для CoupledBeams: предлагает характеризовать veering и localization через вторую производную собственного значения и первую производную собственного вектора, а также операционально связывает эти признаки с близкими собственными значениями.
- Обозначения: derivatives of eigenvalues/eigenvectors, close eigenvalues, curve veering, mode localization; производные берутся по параметру конкретной общей задачи и не тождественны автоматически производным по проектному `mu`.
- Критично смотреть: определения диагностических производных, критерий close eigenvalues, пример слабосвязанных пружин и выводы о совместной интерпретации eigenvalue/eigenvector sensitivity.
- Замечание по применимости: полезен как общий диагностический источник; применение его критериев к CoupledBeams требует отдельного расчёта производных и не заменяет continuation-based `branch_id` и анализ форм.
- Метаданные: подтверждены по локальному журнальному PDF и издательской записи; DOI `10.1006/jsvi.2002.5010`. Сводная оценка: [literature assessment](../veering/literature_assessment.md#liu2002).

## `nair_1973_quasi_degeneracies`
- PDF: `docs/literature/pdf/nair1973.pdf`
- Тип: статья.
- Роль: вспомогательный.
- Что важно для CoupledBeams: источник по vocabulary and interpretation of quasi-degeneracy, true frequency crossing, and rapid modal/nodal-pattern changes in vibration spectra.
- Обозначения: `frequency crossing`, `transition`, `quasi-degeneracy`, `symmetry group`, close eigenvalues/eigenfunctions, nodal patterns; exact operators and symmetry notation are plate-specific.
- Критично смотреть: pp. 975--976 summary, terminology, and notation; analytical discussion of symmetry-group crossing rules; conclusion around pp. 985--986; rectangular/skew-plate examples and figures for rapid nodal-pattern changes.
- Замечание по применимости: useful by mechanism and terminology, not by geometry. Do not transfer the plate symmetry-group criteria directly to the two-rod CoupledBeams system.
- Метаданные: recovered from the local PDF title page.

## `manconi_2017_veering_strong_coupling`
- PDF: `docs/literature/pdf/vib_139_02_021009.pdf`
- Тип: статья.
- Роль: основной для veering-линии по механизму.
- Что важно для CoupledBeams: ключевой общий theoretical source для различения rapid veering under weak coupling и slow evolution under strong coupling; вводит `uncoupled-blocked system`, `skeleton` и `critical points`.
- Обозначения: mode veering, weak/strong coupling, uncoupled-blocked system, skeleton, critical point, eigenvector rotation; малый параметр coupling order не является проектным `mu`.
- Критично смотреть: p. 021009-1 abstract/introduction; Sec. 2 and Eq. (3) for weak coupling; Sec. 2.2 and Eqs. (17)--(19) near a critical point; Figs. 1, 3, 5 for skeleton/eigenvector rotation; Fig. 6 for strong coupling and gradual evolution; Sec. 5.3/Fig. 15 for continuous examples.
- Замечание по применимости: очень близко по spectral mechanism, но не по геометрии. Использовать как главный источник для осторожной формулировки `not strict veering, possibly slow evolution / modal-character reorganization`.
- Метаданные: title/authors/pages recovered from local PDF; DOI `10.1115/1.4035109` cross-checked because it is not exposed clearly in the local PDF metadata.

## `ehrhardt_2018_clamped_beam_veering`
- PDF: `docs/literature/pdf/ehrhardt2018.pdf`
- Тип: статья.
- Роль: основной/сильный вспомогательный для beam-like veering analogy.
- Что важно для CoupledBeams: лучший близкий beam-аналог по механизму: symmetry-preserving crossing vs symmetry-breaking veering, eigenvector correlation/self-MAC, mode-shape mixing in a beam assembly.
- Обозначения: clamped-clamped cross-beam, bending/torsion LNMs, movable tip masses, `self-MAC`, linear normal mode veering, nonlinear normal modes; tuning variables are mass positions/asymmetry, not `mu`.
- Критично смотреть: pp. 1--2 introduction; Sec. 2/Fig. 1 system and model; Sec. 3.1/Fig. 3 linear crossing vs veering; Sec. 3.2/Figs. 4--5 nonlinear crossing/veering analogue; Secs. 4--5/Figs. 6--8 for forced/experimental comparison; Sec. 6 conclusion.
- Замечание по применимости: близко по beam mechanism, но не по геометрии CoupledBeams; nonlinear parts are secondary for the current linear `mu` question.
- Метаданные: verified from local PDF XMP and title page.

## `lacarbonara_2005_imperfect_beams_veering`
- PDF: `docs/literature/pdf/lacarbonara2005.pdf`
- Тип: статья.
- Роль: вспомогательный, второй эшелон для текущей линейной veering-задачи.
- Что важно для CoupledBeams: useful beam/nonlinear source linking veering, one-to-one internal resonance, nonlinear stretching, bifurcations, frequency islands, and mode localization.
- Обозначения: imperfect/shallow beam, torsional spring constant `k`, rise `b`, natural-frequency veering, internal resonance detuning, mode localization.
- Критично смотреть: pp. 987--988 abstract/introduction; Sec. 2 and Eqs. (1)--(7) formulation and boundary conditions; Sec. 2.1/Figs. 2--4 natural frequencies and veering/crossing; Sec. 3 perturbation analysis; Sec. 4/Figs. 6--15 bifurcation/localization; Sec. 5 conclusion.
- Замечание по применимости: полезно для nonlinear beam/localization background, but not a primary source for strict linear veering under `mu`.
- Метаданные: recovered from local PDF title page, outline, and DOI metadata.

## `fontanela_2021_nonlinear_localisation_coupled_beams`
- PDF: `docs/literature/pdf/s11071-020-05760-x.pdf`
- Тип: статья.
- Роль: вспомогательный для localization in two-beam systems.
- Что важно для CoupledBeams: источник по nonlinear vibration localisation in a symmetric system of two weakly coupled beams; useful for localization vocabulary and arm-wise response intuition.
- Обозначения: vibration localisation, symmetry breaking bifurcation, clearance nonlinearity, piecewise linear stiffness, in-phase/out-of-phase modes, localized state, coupling stiffness `k_c`.
- Критично смотреть: pp. 3417--3418 abstract/introduction; Sec. 2.1 and Eqs. (1)--(6) two-DOF model; Figs. 3--4 backbone/bifurcating localized branches; Sec. 3.1/Fig. 7 test setup; Figs. 9--11 measured localized states; Sec. 4 summary.
- Замечание по применимости: близко по two-beam geometry, but the mechanism is nonlinear contact/clearance localization, not strict linear veering.
- Метаданные: verified from local PDF XMP/title page.

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

## `zheltkov_chan_2008_spatial_rod_spectrum`
- PDF: `docs/literature/pdf/opredelenie-spektra-svobodnyh-kolebaniy-prostranstvennoy-sistemy-pryamyh-odnorodnyh-sterzhney.pdf`
- Тип: статья.
- Роль: вспомогательный/близкий по пространственным стержневым системам.
- Что важно для CoupledBeams: determination of free-vibration spectrum for a spatial system of straight homogeneous rods; useful background for multi-rod spectral equations and spatial rod-system boundary/joint setup.
- Обозначения: spatial system of straight homogeneous rods, free-vibration spectrum; exact symbols require manual PDF reading because local metadata is only a CyberLeninka wrapper.
- Критично смотреть: постановку пространственной стержневой системы и вывод спектрального условия; pp. 58--65.
- Метаданные: `NEEDS_CHECK` -- bibliographic fields inferred from filename/local CyberLeninka metadata and secondary search; verify author initials and journal record before external citation.

## `pavlov_2019_rod_package_vibrations`
- PDF: `docs/literature/pdf/Dissertatsiya-Pavlov-A.M.pdf`
- Тип: диссертация.
- Роль: вспомогательный для пакетов стержней, симметрии и forced/free vibration background.
- Что важно для CoupledBeams: собственные и вынужденные колебания пакета стержней; useful for interpreting rod bundles/packages, modal classification, and symmetry-based decomposition when moving beyond two rods.
- Обозначения: rod package, free and forced vibrations, symmetry/classification language; exact notation should not be imported into the current model without targeted reading.
- Критично смотреть: title/abstract/introduction and chapters on free-vibration classification of rod packages; local PDF metadata gives specialty 01.02.04.
- Метаданные: `NEEDS_CHECK` -- author/title/specialty decoded from local PDF metadata; verify official defense organization and page count before external citation.

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
- Additional local copy from the Timoshenko/shear-source import: `docs/literature/pdf/A2-4.pdf`
- Тип: статья.
- Роль: основной.
- Что важно для CoupledBeams: exact free in-plane vibrations for two orthogonal beam members with an internal elastic hinge and elastic boundary conditions; очень близко к задаче о joint flexibility and connection angle.
- Обозначения: internal elastic hinge flexibility, boundary stiffnesses, Euler--Bernoulli assumption, Hamilton principle, separation of variables.
- Критично смотреть: постановку через Hamilton's principle, derivation of the exact frequency equation, and comparison with FEM/experiment; для короткой статьи полезен весь текст Article ID 624658, 9 pages.
- Метаданные: DOI и article ID извлечены надёжно; локальный `A2-4.pdf` пишет "9 pages", поэтому страницу/объём перед внешним цитированием лучше держать как Article ID 624658, 9 pages.

## `shayna_2022_rod_systems_localized_features`
- PDF: локальная копия отсутствует; ранее ожидавшийся файл не найден в текущем worktree.
- Тип: диссертация.
- Роль: вспомогательный.
- Что важно для CoupledBeams: математическое моделирование стержневых систем с локализованными особенностями; полезно для cases with singular supports, discontinuities, and nonclassical localized effects.
- Обозначения: по введению и оглавлению заметны generalized functions / spectral viewpoint / mixed boundary-value formulations; точная нотация зависит от выбранной модели и требует целевого чтения.
- Критично смотреть: главу 1 про математическую модель малых колебаний стержневой системы с особенностями (по оглавлению примерно pp. 13--66) и главу 2 про адаптацию метода конечных элементов для моделей 4-го порядка (примерно pp. 69--101).
- Метаданные: `NEEDS_CHECK` — библиография восстановлена по локальной титульной странице и OCR; перед внешней ссылкой лучше проверить официальный репозиторий ВГУ.

## `mai_2025_nonstationary_thinwalled_bodies`
- PDF: локальная копия отсутствует; ранее ожидавшийся файл не найден в текущем worktree.
- Тип: диссертация.
- Роль: вспомогательный.
- Что важно для CoupledBeams: broader generalized-continuum and moment-elasticity context; может быть полезно, если проект уйдёт от классической балки к более богатым моделям тонкостенных и моментных упругих тел.
- Обозначения: моментные упругие среды, оболочки, пластины и стержни; из оглавления явно видны разделы по продольным колебаниям и изгибу моментного упругого стержня.
- Критично смотреть: главу 3 `Начально-краевые задачи моментных упругих пластин и стержней`, особенно 3.3--3.9 (примерно pp. 53--106), где есть уравнения движения стержней, продольные колебания и изгиб.
- Метаданные: `NEEDS_CHECK` — локальный PDF выглядит как диссертационная рукопись; имя автора и финальный статус публикации нужно сверить по официальной карточке МАИ/ВАК.

## `bauer_2025_coupled_rods`
- PDF: `docs/literature/pdf/Статья-Дорофеев-2025.pdf`
- Тип: статья.
- Роль: основной источник по геометрии двух жёстко сопряжённых стержней; не источник материальной модели композита.
- Что важно для CoupledBeams: изотропные круглые упругие стержни, сопряжённые под углом; безразмерное частотное уравнение, условия сопряжения, сравнение с COMSOL и асимптотики для малого угла и параметра толщины.
- Обозначения: dimensionless variables, coupling angle, thickness parameter, frequency equation, axial and transverse vibrations of rods; exact symbol names уже стоит поднимать из опубликованной версии.
- Критично смотреть: весь текст pp. 73--81; особенно введение, безразмерную постановку, условия сопряжения, частотное уравнение и асимптотики для малого угла.
- Замечание по применимости: использовать только для общей геометрии жёсткого углового сопряжения и частотного уравнения. Не смешивать круглое изотропное сечение этой статьи с направлением `anisotropic_rods` и прямоугольными композитными стержнями.
- Метаданные: local PDF похож на draft build, но библиографические поля проверены по официальному выпуску журнала.

## `kramer_2024_plane_frames_timoshenko`
- PDF: `docs/literature/pdf/A2-1.pdf`
- Тип: статья.
- Роль: основной для Timoshenko/frame-направления.
- Что важно для CoupledBeams: modern plane-frame formulation based on Timoshenko-Ehrenfest beam theory, useful for future coupled-beam/frame assembly and for keeping axial/transverse frame DOFs explicit.
- Обозначения: density, `E`, `G`, area `A`, inertia `I`, Timoshenko-Ehrenfest shear coefficient; the paper uses `k = 5/6` for a rectangular cross-section.
- Критично смотреть: Sec. 2 governing equations and frame structure, boundary/interface conditions, numerical assembly technique, and shifted/deflated Newton solution.
- Метаданные: recovered from local PDF XMP; circular-rod coefficient source: no.

## `howson_1973_axially_loaded_timoshenko_frames`
- PDF: `docs/literature/pdf/A2-2.pdf`
- Тип: статья.
- Роль: основной historical frame/Timoshenko reference.
- Что важно для CoupledBeams: dynamic stiffness method for plane frames whose members include axial load, rotary inertia, and shear deflection; useful as a benchmark for Timoshenko members in frames.
- Обозначения: dynamic stiffness matrix, axial load, rotary inertia, shear deflection, parameters `p`, `r`, and `s` in the title-page abstract.
- Критично смотреть: title-page abstract, dynamic member stiffness derivation, and H-frame theory/experiment comparison.
- Метаданные: recovered from scanned first page; DOI not found in the local PDF.

## `gladwell_1964_vibration_frames`
- PDF: `docs/literature/pdf/A2-3.pdf`
- Тип: статья.
- Роль: вспомогательный frame-background source.
- Что важно для CoupledBeams: early method for free vibration of plane frames using assumed modes and Rayleigh--Ritz style matrix formulation; useful for frame-method history, not for Timoshenko shear correction.
- Обозначения: assumed modes, inertia and stability matrices, rectangular plane frames.
- Критично смотреть: discussion of frame-vibration methods and the comparison against exact solutions.
- Метаданные: recovered from scanned first page and PDF metadata PII; DOI not found in the local PDF.

## `diaz_de_anda_2012_timoshenko_predictions`
- PDF: `docs/literature/pdf/Т1.pdf`
- Тип: статья.
- Роль: основной для experimental validation, critical frequency, and second Timoshenko spectrum.
- Что важно для CoupledBeams: experimental study of Timoshenko beam theory predictions for cylindrical rods and rectangular beams, with 3-D FEM comparison; useful for deciding the safe diagnostic frequency range of future Timoshenko corrections.
- Обозначения: Timoshenko shear coefficient, critical frequency `f_c`, first/second TBT spectra, cylindrical rods, rectangular beams, free-free boundary conditions.
- Критично смотреть: abstract, Sec. 2 Timoshenko beam theory, experimental/FEM comparisons, and conclusions on the second spectrum and valid range.
- Метаданные: recovered from local PDF XMP/title page; DOI `10.1016/j.jsv.2012.07.041`.

## `diaz_de_anda_2005_locally_periodic_timoshenko_rod`
- PDF: `docs/literature/pdf/Т2.pdf`
- Тип: статья.
- Роль: основной для circular/cylindrical rods and baseline circular shear coefficient.
- Что важно для CoupledBeams: locally periodic aluminum rods of circular cross-section are modeled with Timoshenko beam theory and a transfer matrix method, then compared with EMAT measurements.
- Обозначения: Timoshenko shear coefficient `k`, transfer matrix, unit cells, circular rods, EMAT measurements; the local PDF gives `k = (6 + 12*nu + 6*nu^2)/(7 + 12*nu + 4*nu^2)` and uses `k = 0.925` for `nu = 0.3`.
- Критично смотреть: Sec. II transfer matrix method, the coefficient choice near Fig. 6, and the experiment/theory comparison.
- Метаданные: recovered from local PDF title page; DOI `10.1121/1.1880732`.

## `franco_villafane_2014_best_shear_coefficient`
- PDF: `docs/literature/pdf/Т3.pdf`
- Тип: arXiv/preprint.
- Роль: основной для non-uniqueness, best-fit shear coefficient, and critical-frequency caution.
- Что важно для CoupledBeams: explicitly treats the shear coefficient as an adjustment/modeling parameter and compares one-coefficient, two-coefficient, below-critical, and above-critical choices against experimental data.
- Обозначения: `kappa`, `kappa_1`, `kappa_3`, critical frequency `f_c`, best-fit coefficients, first/second TBT spectra.
- Критично смотреть: introduction on coefficient non-consensus, Table 1 of coefficient choices, and conclusions about different coefficients below/above `f_c`.
- Метаданные: local PDF is arXiv:1405.4885v2, submitted to Elsevier May 26, 2014; DOI not found in the local PDF.

## `stephen_2002_check_timoshenko_accuracy`
- PDF: `docs/literature/pdf/Т4.pdf`
- Тип: статья.
- Роль: основной caution source for shear-coefficient accuracy and second-spectrum interpretation.
- Что важно для CoupledBeams: short note on why the Timoshenko shear coefficient is not a unique universal constant, with discussion of Cowper, Hutchinson, two-coefficient theory, and second-spectrum behavior.
- Обозначения: shear coefficient, wavelength/beam-depth ratio, Rayleigh surface wave limit, second spectrum.
- Критично смотреть: introduction, coefficient comparison, comments on the second spectrum, and references to Cowper/Hutchinson.
- Метаданные: recovered from local PDF metadata/title text; DOI not found in the local PDF.

---

**Направление: `anisotropic_rods / rectangular composite rods`.** Ниже
собраны журнальные источники литературной основы отправленной статьи. Они не
изменяют verified theory, формулы или результаты проекта.

## `miller_1975_orthotropic_beam_resonances`
- PDF: `docs/literature/pdf/miller_1975_orthotropic_beam_resonances.pdf`
- Тип: статья.
- Роль: основной источник литературного фона.
- Что важно для CoupledBeams: общий ортотропный стержень, изгибные и крутильные собственные частоты и одно из ранних исследований влияния материальной анизотропии на связанный спектр.
- Обозначения: generally orthotropic beam, flexural and torsional resonances, fibre-orientation angle, characteristic matrix; обозначения статьи не переносятся в текущую модель автоматически.
- Критично смотреть: pp. 433--435 для постановки и происхождения изгибно-крутильной связанности; pp. 435--443 для аналитического решения; pp. 443--447 для граничных условий, характеристической матрицы и численного примера; pp. 447--448 для выводов.
- Замечание по применимости: один прямой стержень, а не два стержня с жёстким угловым сопряжением; источник фона, не прямой источник текущего частотного определителя.
- Метаданные: подтверждены по первой странице и PII локального журнального скана; DOI `10.1016/S0022-460X(75)80107-6`. Подробности: [source note](notes/miller_1975_orthotropic_beam_resonances.md).

## `teh_huang_1980_fibre_orientation_composite_beams`
- PDF: `docs/literature/pdf/teh_huang_1980_fibre_orientation_composite_beams.pdf`
- Тип: статья.
- Роль: основной и особенно близкий источник по постановке вопроса об угле армирования.
- Что важно для CoupledBeams: непосредственное исследование влияния ориентации волокон на свободные колебания композитной балки, собственные частоты, изгибно-крутильную связанность и изменение форм.
- Обозначения: fibre-orientation angle, predominantly flexural/torsional modes, bending and twisting moments; локальная нотация относится к одному консольному стержню.
- Критично смотреть: pp. 327--328 для постановки; раздел 2 для уравнений движения; Figs. 2--4 для частот и мер изгибно-крутильного влияния; Figs. 5--11 и conclusion для изменения форм при варьировании ориентации волокон.
- Замечание по применимости: геометрия одного стержня, а не двух сопряжённых стержней.
- Метаданные: подтверждены по первой странице локального журнального PDF; DOI `10.1016/0022-460X(80)90616-1`. Подробности: [source note](notes/teh_huang_1980_fibre_orientation_composite_beams.md).

## `chandrashekhara_1990_composite_beam_free_vibration`
- PDF: `docs/literature/pdf/chandrashekhara_1990_composite_beam_free_vibration.pdf`
- Тип: статья.
- Роль: основной источник по уточнённой теории композитной балки.
- Что важно для CoupledBeams: свободные колебания симметрично слоистых композитных балок с учётом поперечного сдвига и инерции вращения сечений; обосновывает обращение к уточнённой теории Тимошенко для композитов.
- Обозначения: first-order shear deformation, rotary inertia, laminated beam resultants and arbitrary end conditions; не подменяют обозначения моноклинной модели проекта.
- Критично смотреть: pp. 269--271 для постановки и первого порядка сдвиговой кинематики; математическую формулировку и точное решение; таблицы влияния сдвига, анизотропии и граничных условий на собственные частоты.
- Замечание по применимости: статья не является источником теории обобщённого кручения Фойгта--Лехницкого и не задаёт жёсткий угловой стык.
- Метаданные: подтверждены по первой странице локального журнального PDF; DOI `10.1016/0263-8223(90)90010-C`. Подробности: [source note](notes/chandrashekhara_1990_composite_beam_free_vibration.md).

## `han_1999_four_beam_theories`
- PDF: `docs/literature/pdf/han_1999_four_beam_theories.pdf`
- Тип: статья.
- Роль: обзорно-методический источник по линейным теориям балки.
- Что важно для CoupledBeams: сопоставляет Эйлера--Бернулли, Рэлея, сдвиговую и Тимошенко-модели, включая уравнения движения, граничные условия, частотные уравнения и формы.
- Обозначения: Euler--Bernoulli, Rayleigh, shear and Timoshenko beam theories; частотные параметры относятся к однородной поперечно колеблющейся балке.
- Критично смотреть: обзор четырёх моделей, вывод через принцип Гамильтона, частотные уравнения для основных граничных условий и численный пример для непротяжённой балки; pp. 935--988.
- Замечание по применимости: материал статьи не является специальной теорией моноклинного композитного стержня.
- Метаданные: подтверждены по первой странице локального журнального PDF и издательской записи; DOI `10.1006/jsvi.1999.2257`.

## `labuschagne_2009_linear_beam_theories`
- PDF: `docs/literature/pdf/labuschagne_2009_linear_beam_theories.pdf`
- Тип: статья.
- Роль: обзорно-методический источник по применимости линейных теорий балки.
- Что важно для CoupledBeams: систематическое сравнение Эйлера--Бернулли, Тимошенко и двумерной теории упругости по собственным частотам и формам; подчёркивает роль поперечного сдвига и вращательной инерции.
- Обозначения: cantilever beam, natural frequencies and modes, Euler--Bernoulli, Timoshenko and two-dimensional elasticity.
- Критично смотреть: постановку трёх моделей, их спектральное сравнение и выводы об области практической применимости; pp. 20--30.
- Замечание по применимости: общий источник фона; не является прямым доказательством десятипроцентного критерия текущей работы.
- Метаданные: подтверждены по первой странице локального журнального PDF и издательской записи; DOI `10.1016/j.mcm.2008.06.006`.

## `banerjee_williams_1996_composite_timoshenko_dynamic_stiffness`
- PDF: `docs/literature/pdf/banerjee_williams_1996_composite_timoshenko_dynamic_stiffness.pdf`
- Тип: статья.
- Роль: основной методический источник.
- Что важно для CoupledBeams: точная динамическая матрица жёсткости композитной балки Тимошенко, материальная изгибно-крутильная связанность, поперечный сдвиг, инерция вращения и вычисление собственных частот алгоритмом Wittrick--Williams.
- Обозначения: dynamic stiffness matrix, composite Timoshenko beam, bending--torsion coupling, Wittrick--Williams algorithm.
- Критично смотреть: pp. 573--574 и introduction; Sec. 2 для теории; Sec. 3 для применения динамической матрицы; Sec. 4 для сопоставлений; Sec. 5 для границ метода.
- Замечание по применимости: методический источник не заменяет постановку Ярцева и условия жёсткого углового сопряжения проекта.
- Метаданные: подтверждены по первой странице, встроенным метаданным и издательской записи; DOI `10.1006/jsvi.1996.0378`. Подробности: [source note](notes/banerjee_williams_1996_composite_timoshenko_dynamic_stiffness.md).

## `song_librescu_1993_anisotropic_thinwalled_beams`
- PDF: `docs/literature/pdf/song_librescu_1993_anisotropic_thinwalled_beams.pdf`
- Тип: статья.
- Роль: основной источник по связанным колебаниям анизотропных тонкостенных стержней.
- Что важно для CoupledBeams: анизотропные композитные толстостенные и тонкостенные однозамкнутые стержни, связанные собственные колебания, поперечная сдвиговая податливость и неравномерное кручение.
- Обозначения: anisotropic composite thin-walled beam, closed cross-section contour, non-uniform torsion, transverse shear flexibility.
- Критично смотреть: p. 129 для состава модели и заявленных неклассических эффектов; разделы с динамической постановкой и анализом cantilever beam; численные результаты о роли анизотропии и других неклассических эффектов; pp. 129--147.
- Замечание по применимости: геометрия и теория значительно шире сплошного прямоугольного стержня текущей статьи; не является прямым источником её определителя.
- Метаданные: заголовок, авторы, том и страницы подтверждены визуально по первой странице защищённого журнального PDF; DOI `10.1006/jsvi.1993.1325` подтверждён издательской записью. Подробности: [source note](notes/song_librescu_1993_anisotropic_thinwalled_beams.md).

## `piovan_2008_tapered_shear_flexible_composite_beams`
- PDF: `docs/literature/pdf/piovan_2008_tapered_shear_flexible_composite_beams.pdf`
- Тип: статья.
- Роль: основной методический источник по связанным композитным балкам.
- Что важно для CoupledBeams: связанные свободные колебания, сдвигово-деформируемые тонкостенные композитные балки переменного сечения, открытые и замкнутые профили и точные решения методом степенных рядов.
- Обозначения: tapered thin-walled composite beam, CUS/CAS laminations, shear flexibility, power-series solution, coupled modes.
- Критично смотреть: Sec. 2 для модели и граничных условий; Sec. 3 для метода степенных рядов; Sec. 4 для численных сравнений, coupled mode labels и влияния сужения; Sec. 5 для выводов.
- Замечание по применимости: полезен как источник по связанным композитным балкам, но не как прямой источник текущего частотного определителя.
- Метаданные: подтверждены по локальной издательской Article-in-Press версии с финальными томом и страницами; DOI `10.1016/j.jsv.2008.02.044`. Подробности: [source note](notes/piovan_2008_tapered_shear_flexible_composite_beams.md).

## `ryabov_yartsev_2016_box_beams_part1`
- PDF: `docs/literature/pdf/ryabov_yartsev_2016_box_beams_part1.pdf`
- Тип: статья; локальный PDF — русский оригинал официальной английской переводной публикации.
- Роль: основной источник линии Рябова--Ярцева по математической модели коробчатого стержня.
- Что важно для CoupledBeams: вывод модели затухающих связанных колебаний анизотропного тонкостенного коробчатого стержня, вариационный подход и сведение трёхмерных соотношений к одномерной системе.
- Обозначения: шесть осевых перемещений и поворотов, депланация, комплексные модули, функционал Гамильтона; нотация относится к замкнутому тонкостенному профилю.
- Критично смотреть: русские pp. 221--229: геометрию и кинематику, вариационный вывод, уравнения движения и построение частотного условия.
- Замечание по применимости: это не та же геометрия, что сплошной прямоугольный моноклинный стержень главы 2 монографии.
- Метаданные: каноническая запись объединяет оригинал и перевод: English version 49(2), 130--137, DOI `10.3103/S1063454116020126`; русский DOI `10.21638/11701/spbu01.2016.206` сохранён только как дополнительная информация.

## `ryabov_yartsev_2016_box_beams_part2`
- PDF: `docs/literature/pdf/ryabov_yartsev_2016_box_beams_part2.pdf`
- Тип: статья; локальный PDF — русский оригинал официальной английской переводной публикации.
- Роль: основной источник физической интерпретации связанных форм.
- Что важно для CoupledBeams: численный эксперимент по влиянию ориентации армирующих волокон на частоты, потери и взаимную трансформацию связанных форм для HMS/DX-209.
- Обозначения: symmetric/asymmetric box beams, partial and coupled frequencies, mechanical loss factors, mode transformation regions.
- Критично смотреть: русские pp. 429--439, особенно описание HMS/DX-209, зависимости от угла армирования и обсуждение взаимной трансформации мод.
- Замечание по применимости: источник терминологии и физической интерпретации, а не прямое подтверждение всех результатов текущей статьи.
- Метаданные: каноническая запись объединяет оригинал и перевод: English version 49(3), 260--268, DOI `10.3103/S1063454116030110`; русский DOI `10.21638/11701/spbu01.2016.311` сохранён только как дополнительная информация.

## `ryabov_yartsev_2021_monoclinic_strip`
- PDF: `docs/literature/pdf/ryabov_yartsev_2021_monoclinic_strip.pdf`
- Тип: статья; локальный PDF — русский оригинал официальной английской переводной публикации.
- Роль: наиболее близкий журнальный источник по материальной модели прямоугольного моноклинного стержня и терминологии статьи.
- Что важно для CoupledBeams: прямоугольная моноклинная полоса, уточнённая теория изгиба Тимошенко, обобщённое кручение Фойгта--Лехницкого, определения `E_x`, `G_{xy}`, `G_{xz}`, изгибно-крутильная связанность, ортотропные пределы при `theta=0°` и `90°` и экспериментальная проверка.
- Обозначения: `w`, `Phi`, `E_x(theta)`, `G_{xy}(theta)`, `G_{xz}(theta)`, mutual-influence coefficients, shear coefficient `k`, `I_y`, `I_p` and generalized torsional rigidity.
- Критично смотреть: русские pp. 696--697 для уравнений (1)--(5), определений модулей, граничных условий и ортотропного предела; pp. 697--698 для частотного решения; pp. 698--700 для экспериментальной проверки; последующие разделы для влияния угла и длины.
- Замечание по применимости: модель одного стержня наиболее близка материально, но статья не выводит условия жёсткого сопряжения двух стержней.
- Метаданные: каноническая запись объединяет оригинал и перевод: English version 54(4), 437--446, DOI `10.1134/S1063454121040166`; русский DOI `10.21638/spbu01.2021.415` сохранён только как дополнительная информация. Подробности: [source note](notes/ryabov_yartsev_2021_monoclinic_strip.md).

## `ryabov_yartsev_2023_composite_wing_coupling`
- PDF: `docs/literature/pdf/ryabov_yartsev_2023_composite_wing_coupling.pdf`
- Тип: статья; локальный PDF — русский оригинал официальной английской переводной публикации.
- Роль: вспомогательный источник общего физического контекста.
- Что важно для CoupledBeams: управление связанностью колебаний композитного крыла, разделение упругой и инерционной связанности, материал HMS/DX-209 и влияние ориентации армирующих слоёв.
- Обозначения: elastic and inertial coupling coefficients, bending and torsional components, fibre-orientation angle; коэффициенты определены для энергии формы крыла.
- Критично смотреть: русские pp. 344--356, особенно определение упругой и инерционной связанности, Sec. 3 с HMS/DX-209 и угловыми зависимостями и итоговое обсуждение.
- Замечание по применимости: не переносить определения коэффициентов связанности крыла в текущую статью без отдельного вывода.
- Метаданные: каноническая запись объединяет оригинал и перевод: English version 56(2), 252--260, DOI `10.1134/S1063454123020152`; русский DOI `10.21638/spbu01.2023.214` сохранён только как дополнительная информация.

## `carrera_2012_advanced_beam_joined_wings`
- PDF: `docs/literature/pdf/carrera_2012_advanced_beam_joined_wings.pdf`
- Тип: статья.
- Роль: дополнительный фон; не входит в 14 источников отправленной статьи.
- Что важно для CoupledBeams: higher-order beam formulations for conventional and joined wings, coupled bending/torsion modes, cross-section warping and comparison with shell/solid FEM.
- Обозначения: Carrera Unified Formulation, higher-order cross-section expansion, conventional/joined wings and mixed vibration modes.
- Критично смотреть: abstract and introduction for model hierarchy; numerical joined-wing cases and comparisons of classical and higher-order beam theories; pp. 282--293.
- Замечание по применимости: joined-wing geometry and refined finite-element kinematics are relevant background, but this paper is not a source for the submitted article's material model or determinant.
- Метаданные: локальный издательский PDF; metadata verified from its first page, XMP and the ASCE record; DOI `10.1061/(ASCE)AS.1943-5525.0000130`.

## `yartsev_2024_coupled_composite_structures`
- Локальные PDF-фрагменты одной монографии: `docs/literature/pdf/Глава 1_compressed.pdf`, `docs/literature/pdf/Глава 2_compressed.pdf`, `docs/literature/pdf/Глава 3 Часть 1.pdf`, `docs/literature/pdf/Глава 3 Часть 2_compressed.pdf`, `docs/literature/pdf/Глава 4_compressed.pdf`, `docs/literature/pdf/Применения_compressed.pdf`, `docs/literature/pdf/Литература_compressed.pdf`. Фрагменты являются локальными и могут отсутствовать в публичном clone; это не семь самостоятельных публикаций.
- Тип: монография.
- Роль: основной литературный источник направления `anisotropic_rods`; не заменяет verified theory и baseline изотропной модели.
- Что важно для CoupledBeams: глава 1 вводит линейную упругость и вязкоупругость, комплексные модули, преобразование характеристик однонаправленного слоя, коэффициенты взаимного влияния и параметры материалов в табл. 1.2. Глава 2 задаёт моноклинный стержень прямоугольного сечения, связанный изгиб по Тимошенко и обобщённое кручение, изгибно-крутильное материальное взаимодействие, безопорный и консольный случаи и экспериментальное определение характеристик. Глава 3 рассматривает многослойные пластины и относится к общей теории связанных колебаний композитов, но не является исходной моделью первого стержневого этапа. Глава 4 развивает пространственную теорию тонкостенных стержней замкнутого профиля с продольным движением, двумя изгибами, сдвигом, кручением и депланацией; это потенциальный более общий будущий этап, а не прямая замена модели главы 2. Глава 5 содержит инженерные применения и примеры управления связанностью.
- Основные обозначения: комплексные модули `M*`, упругие и сдвиговые характеристики `E_i`, `G_ij`, коэффициенты Пуассона `nu_ij`, плотность `rho`, угол армирования `theta`, геометрия прямоугольного стержня `a`, `b`, `L`, коэффициент сдвига `k`, изгибное перемещение `w`, поворот `psi` и кручение `Phi`. Не переносить эти обозначения в verified isotropic theory без отдельного consistency audit.
- Критично смотреть: глава 1, печатные стр. 24--25 для комплексных модулей и (1.32), (1.34); стр. 30--31 для определений модального коэффициента потерь и (1.41), (1.42); стр. 40--46 для преобразования податливостей, (1.50)--(1.56) и таблиц материалов. В главе 2 критичны стр. 52--55 с (2.1)--(2.18), стр. 56--57 с безопорным образцом и рис. 2.2, а также стр. 64--68 с консольной задачей и рис. 2.8. Для будущего более общего этапа глава 4 начинается на стр. 137.
- Замечания по корректности и статус: в буквально напечатанной (2.1) отсутствует множитель `I_y` при инерционном члене; внутренняя согласованность, размерности и воспроизведение рис. 2.2 требуют `rho * I_y * psi_tt`. Напечатанные после (2.16) знаки `d0` и `f0` не воспроизводят независимую положительную крутильную часть спектра при `theta = 0°` и `90°`. Варианты `state_corrected` с восстановленным `I_y` и `eliminated_corrected` с положительными `d0`, `f0` совпали по первым восьми положительным корням с максимальным относительным расхождением около `1.4e-9` и воспроизвели расчётные сплошные кривые рис. 2.2 со статусом `PASS_WITHIN_GRAPH_RESOLUTION`. Сопоставление сохранённых частот и коэффициентов потерь с рис. 2.8 дало `BOOK_SLOPE_CLAMP_CONFIRMED`: source-faithful консоль использует `w=0`, `w'=0`, `Phi=0`. Буквально напечатанный sign variant остаётся только диагностическим. Source-faithful внешний clamp не задаёт автоматически условия будущего внутреннего жёсткого узла; их необходимо вывести отдельно. Экспериментальные точки рис. 2.2 отдельно не оцифровывались.
- Подтверждённые метаданные: Б. А. Ярцев; *Связанные колебания композитных конструкций*; монография; Санкт-Петербург: ФГУП «Крыловский государственный научный центр», 2024; 216 с.: ил.; ISBN `978-5-6048511-4-2`; УДК `534.12:678.067`; ББК `35.719`; язык русский. Данные прочитаны с первой библиографической страницы локального фрагмента главы 1.
- Подробности: [source note и карта скана](../anisotropic_rods/source_note_yartsev_2024.md); [free-free reproduction note](../anisotropic_rods/yartsev_ch2_single_rod_reproduction.md); [cantilever reproduction note](../anisotropic_rods/yartsev_ch2_cantilever_reproduction.md); локальный generated [free-free report](../../results/anisotropic_rods/yartsev_ch2_free_free/single_rod_reproduction_report.md).
