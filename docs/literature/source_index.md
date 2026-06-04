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
