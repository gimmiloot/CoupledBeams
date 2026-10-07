# Source Index

Индекс ниже собран под задачу `CoupledBeams`: связь стержней под углом, собственные частоты, условия сопряжения, граничные условия и сравнение аналитики с FEM.

Citation keys синхронизированы с `docs/literature/bibliography.bib`.

См. также отдельную заметку по новым источникам для Timoshenko theory, shear
coefficient и circular rods:
`docs/literature/timoshenko_shear_sources.md`.

Регистрация 2026-10-05: 11 новых локальных публикаций, включая одну главу
книги. Тематическая навигация: [продольные модели и контрольные задачи](longitudinal_rod_models_sources.md),
[внутриплоскостные и внеплоскостные движения](nonlinear_inplane_outofplane_sources.md).
Для этих 11 записей полный текст доступен локально, метаданные проверены,
а фактически прочитанные места указаны отдельно. Формулы независимо не
проверялись, расчёты не воспроизводились в рамках регистрации; она не является обзором
или подтверждением результатов. Имена PDF сохранены. Дубликатов этих
публикаций по SHA256, DOI, названию и авторам в существующем каталоге не найдено.

Последующий отдельный этап 2026-10-05: [воспроизведение Marais и
Попова–Садовского](../theory/bishop_literature_reproduction.md).
Статусы транскрипции, численной проверки и совпадения печати разделены.

Следующая регистрация 2026-10-05: три новых PDF (Rucka, Jang–Park–Lee,
Liu et al. 2021) и дополнительное чтение ранее зарегистрированного
Banerjee–Ananthapuvirajah. [M-H + Timoshenko: сравнительная карта,
constitutive reductions и ограничения](mindlin_herrmann_timoshenko_sources.md).
Это документационный этап, без расчётов и реализации. Полные тексты
Mindlin–Herrmann и Martin–Gopalakrishnan–Doyle локально не найдены;
их статус — cited / full text unavailable, без новых BibTeX-записей.

## Subsequent Fernandes registration and Jang project decision, 2026-10-06

После предыдущего finite source audit пользователь добавил настоящий
Fernandes PDF. Прежний `RECTANGULAR_MH_PRESET_MAPPING_UNRESOLVED` ниже
сохраняется как результат отдельного вопроса о переносе Ng-type normal
block к текущей energy convention. Новый источник подтверждает повторное
использование correction formula, не равенство normal blocks. Отдельно
выбран [Jang-type project closure и finite single-rod gate](../theory/mindlin_herrmann_timoshenko_single_rod.md#14-production-formulation-decision).

### `fernandes_2022_spectral_tower_cable`

- [PDF](pdf/energies-15-07725-v2.pdf), publisher article, 26 страниц;
  SHA256 `27cb6de0c2038fd849ed153c8e78ffc773bfaea5ef362ef6bed217943c35a2fc`.
  File suffix v2 сохранён; дата специальной revised version не выведена
  из имени файла.
- Yanne Marcela Soares Fernandes; Marcela Rodrigues Machado; Maciej
  Dutkiewicz, *The Spectral Approach of Love and Mindlin-Herrmann Theory
  in the Dynamical Simulations of the Tower-Cable Interactions under the
  Wind and Rain Loads*. Energies 15(20) (2022), article 7725;
  DOI `10.3390/en15207725`. PDF 1: received 2022-08-18, accepted 2022-10-17,
  published 2022-10-19. Issue 20 адресно подтверждён
  [издателем](https://www.mdpi.com/1996-1073/15/20/7725); остальные metadata
  — PDF. Широкого литературного web-поиска не было.
- Прочитаны abstract p.1, §2.4 pp.6–8 в части полей, (19)–(23),
  описание роли rod в §5.1 p.13 и ref.32 p.25. Визуально сверены PDF
  1,6,7,25. Source psi — transverse contraction, independent field.
  I назван inertia of the cross-section; в этих местах нет `I=bh³/12`,
  polar Ip или aspect ratio. Для tower заданы A=.01 m², без b,h.
  Поэтому статью не объявляем проверкой конкретной rectangular geometry.
- (21), p.7: `K_r1=12/pi²`,
  `K_r2=K_r1*((1+nu)/(.87+1.12nu))²`. (19),(20),(23) помещают первый
  factor в GI-gradient stiffness, второй — в rhoI inertia. Это прямое
  повторное использование формул Ng, не square/root/reciprocal mapping.
  Factors заданы formula, specimen fit этих factors не описан.
- Normal block сохраняет `(2mu+lambda)A,lambda*A`,
  `lambda=nu*E/[(1+nu)(1−2nu)]`, `mu=G`. Следовательно, это отдельный
  Ng/Fernandes-type literature variant, **не** selected Jang normal closure
  `C=EA/(1−nu²),nu*C`. Prescription 12/pi² не стала production default.
  Два harmonic expressions (19) напечатаны без `=0`; (23) явно задаёт
  homogeneous zero system. PDF не исправлен; solver этой статьи не добавлен.
- Ref. **32**, p.25 — Doyle 1997. Книга не прочитана/проверена независимо.
  Source problem: M-H axial spectral elements tower, cable и их response
  на wind/rain. Tower-cable calculations, damping, loading, FEM comparisons
  и source FRF benchmarks не воспроизводились. Не источник наших joint BC.

## Rectangular M-H prescription audit, 2026-10-06

Проверены два фактически новых PDF. Первый — Ng 2014, второй — препринт
Elishakoff–Tharu, **не** ожидаемый Fernandes–Machado–Dutkiewicz 2022,
DOI из задания `10.3390/en15207725`. Поиск по титульным страницам и
содержимому 92 локальных PDF проекта, включая untracked/ignored, не нашёл
эту статью Energies. Её полный текст не прочитан, metadata/key/hash не
выдуманы. [Mapping и hard gate](../theory/mindlin_herrmann_timoshenko_single_rod.md#13-production-rectangular-m-h-correction-prescription)
остаются unresolved; повторное употребление prescription двумя доступными
rectangular sources не подтверждено. Имена и байты PDF сохранены.

### `ng_2014_bayesian_guided_wave_damage`

- [PDF](pdf/hdl_85788.pdf), **accepted version** с листом репозитория,
  41 PDF-страница; SHA256
  `d755a482f4e9b1760ba7385dc423324b0aedaa9694a6db452808cdf7f95ae8ae`.
- Ching-Tai Ng, *Bayesian model updating approach for experimental
  identification of damage in beams using guided waves*. Structural Health
  Monitoring 13(4) (2014), 359–373; DOI `10.1177/1475921714532990`.
  Metadata — PDF 1–2. Дата листа репозитория 2014-10-02 не является
  установленной online-first/accepted date. Строка copyright Taylor &
  Francis на листе соседствует с SAGE permissions; издатель по этой строке
  не переопределяется. Для формул используются **PDF pages**, без
  приписывания accepted manuscript журнальной пагинации.
- Прочитано: PDF 9–11, подразделы Mindlin–Herrmann и начало frequency-domain
  formulation, (1)–(7); PDF 33, ref. 37. Визуально проверены PDF 1–2,9–11,33.
  (1): `u_j≈bar(u_j)(x,t)`, `v_j≈y*bar(phi_j)(x,t)`; phi здесь contraction
  scaling, не Timoshenko rotation, несмотря на словесное rotational angle.
  `A_j=b_j*h_j`, `I_j=b_j*h_j^3/12`: прямоугольник, planar I, не polar Ip.
- PDF 10, (2) и следующий абзац прямо содержат `S1=12/pi^2` и
  `S2,j=S1*((1+nu_j)/(0.87+1.12*nu_j))^2`. S1 умножает `mu_j*I_j*phi_j,xx`,
  S2,j — `rho_j*I_j*phi_j,tt`; `mu_j=G_j`. Следовательно, роли
  `H=S1*GI` и `j=S2*rhoI` совпадают с project convention **без** корня,
  квадрата или обращения самих S. b/h в этих формулах отсутствует; S2
  зависит только от nu. Они заданы как model formula перед inverse
  specimen identification; fit этих MH factors не описан. Constants .87
  и 1.12 — напечатанные decimal constants, не точные физические рационалы.
- **Constitutive warning:** в (2) normal diagonal — `(2mu+lambda)A`,
  internal cross coefficient — `lambda*A`; напечатано
  `lambda=nu*E/((1+nu)*(1-2nu))`. Это не текущие `C=EA/(1-nu^2)` и `nu*C`.
  Совпадение correction-factor roles не доказывает равенства полной
  энергии/низкочастотного предела при тех же физических E,nu.
  Эти строки не исправлены на plane-stress reduction и не реализованы.
- Дополнительное печатное предупреждение: PDF 11, нижняя правая запись
  (6) содержит `-mu_j*I_j*S1` без k_j², хотя (2),(4) имеют вторую
  производную. Для mapping использована (2); characteristic matrix (6)
  в новый solver не переносилась.
- Citation chain: ref. **37**, PDF 33 — Doyle, *Wave propagation in
  structures spectral analysis using fast discrete Fourier transforms*,
  2nd ed., Springer, 1997. Книга отдельно не прочитана и prescription
  по ней независимо не проверена. Эксперименты/Bayesian damage fit Ng
  не воспроизводились, article-wide correctness не утверждается.

### `elishakoff_tharu_ssrn_5985611`

- [PDF](pdf/ssrn-5985611.pdf), 100 страниц, **SSRN preprint / not peer
  reviewed**; SHA256
  `49336ea7b0eeaad7062e9251a9bf41f2fb10a6c81f19374669d7ccacb6f3a48d`.
- Isaac Elishakoff; Janak Kumar Tharu, *Sixteen Refined Theories of
  Longitudinal Vibration of Rods: A Review, From Lord Rayleigh to Modernity*.
  Title/authors/identifier 5985611 — PDF 1. Журнал, volume/issue,
  publication year и DOI не установлены по локальному PDF; BibTeX `misc`
  не выдумывает эти поля. File creation 2025-12-23 не является publication
  date. Это **другая публикация**, не статья Fernandes et al.
- Ограниченное чтение: §2.2, PDF 14–21, и references PDF 95. Визуально
  проверены PDF 1,17,19–21. Осесимметричный круглый rod: radius a, radial r,
  longitudinal z; source u — radial surface displacement, w — axial.
  (37),(38), PDF 20, вводят `kappa^2` и `kappa_1^2` в stress-displacement
  relations; это не прямое rectangular S1/S2 inertia mapping. Извлечённый
  полный текст не содержит требуемую формулу `.87+1.12nu`; прочитанная
  M-H постановка не даёт второго подтверждения этой prescription.
- Ref. **10**, PDF 95, цитирует Graff, *Wave Motion in Elastic Solids*,
  1976 (как напечатано); ref. **15** — Doyle, *Wave Propagation in
  Structures*, 1997. Книги независимо не проверены. Circular derivation,
  все шестнадцать теорий и численные comparisons препринта не валидированы
  в этом этапе; никаких production coefficients из него не принято.

## Junction sources, 2026-10-05

Семь локальных публикаций ниже зарегистрированы при документационном
аудите [circular EB joint review](../laminated_beams/circular_eb_rotational_spring_rigid_limit.md).
Метаданные — `VERIFIED_LOCAL_PDF`, без web-поиска; имена файлов сохранены.
Прочитаны только указанные места. Доказательства независимо не проверялись,
численные результаты источников не воспроизводились. Это чтение обосновывает
разграничение asymptotic junction conditions и конечной 1D-диагностики,
не валидацию наших шести уравнений при конечной толщине.

### `leugering_2019_junction_two_elastic_beams`

- [PDF](pdf/The_asymptotic_analysis_of_a_junction_of_two_elast.pdf), издательская статья;
  SHA256 `52b6440bff9380150a0846d4c18df10f39f0b58d307f8bbd3c76fe6f71aec163`.
- G. Leugering, S. A. Nazarov, A. S. Slutskij, *The asymptotic analysis of a
  junction of two elastic beams*. ZAMM 99 (2019), e201700192;
  DOI `10.1002/zamm.201700192`. Год выпуска 2019 подтверждён первой
  страницей и How to cite; copyright/accepted 2018 не заменяет его.
- Прочитано: PDF pp. 1–4, §1.1–1.2, (1.4),(1.7); текст §2.4
  около (2.28),(2.33). Плоская статическая задача, тонкие балки,
  непрерывность вращения и Kirchhoff-type transmission conditions.
  Подвижность элементов и внешнее закрепление существенны; возможны
  алгебраические неизвестные и нелокальные условия. Это не самостоятельное
  доказательство нашего круглого 3D-спектра или finite-thickness точности.

### `kerdid_2026_multi_rod_modes`

- [PDF](pdf/art_01.pdf), издательская статья;
  SHA256 `cd5e2f6a2661e3edf318e72c4ef7fd1dcfe929a3658d5850f95b49255aed57e2`.
- Nabil Kerdid, Mohammed Messaoudi, *Asymptotic analysis of flexural,
  torsional, and stretching modes in a multi-rod structure*.
  Journal of Applied Mathematics and Computational Mechanics 25(2) (2026),
  5–27; DOI `10.17512/jamcm.2026.2.01`. Метаданные — первая страница.
- Прочитано: pp. 5–6 (§1), p. 12 (Lemma 2, (32)–(33), Remark 1),
  p. 26 (§5). Пространственная линейная упругость, два перпендикулярных
  тонких стержня с закреплением обоих внешних концов. Низшие предельные
  формы изгибные; условия (33) сохраняют перпендикулярность проекций осей.
  Это source-specific thin-domain limit, не наша круглая геометрия beta=15°
  и не универсальное описание всех высокочастотных мод.

### `kerdid_1997_multi_rod_vibrations`

- [PDF](pdf/M2AN_1997__31_7_891_0.pdf), оцифрованная журнальная статья с
  обложкой NUMDAM; SHA256 `df7cb20ffec0d7263f4d19dfb11a443c1c6bbbc8b186fc1f7173fa4009a4ca19`.
- N. Kerdid, *Modeling the vibrations of a multi-rod structure*.
  M2AN — Modélisation mathématique et analyse numérique 31(7) (1997),
  891–925. Метаданные — обложка и p. 891; DOI по памяти не добавлен.
- Прочитано: обложка, начало §0, текст введения о предельных формах и
  текст около Lemma 7 / (3.42). Сходимость eigenpairs 3D elasticity к
  1D–1D модели, с различием изгибных и крутильных перемещений;
  один внешний конец закреплён. Не отождествлять с двухзаделочной
  постановкой Kerdid–Messaoudi 2026 или нашим конечным spring-узлом.

### `nazarov_2016_l_shaped_junction`

- [PDF](pdf/Asymptotic_Analysis_of_an_L-Shaped_Junction_of_Two.pdf),
  издательская статья с добавленным листом условий использования;
  SHA256 `98d95ca443f75b0241bd06d49f4833692fb835c87d142519b1ebe39f3f29a024`.
- S. A. Nazarov, A. S. Slutskij, *Asymptotic Analysis of an L-Shaped
  Junction of Two Elastic Beams*. Journal of Mathematical Sciences
  216(2) (2016), 279–312; DOI `10.1007/s10958-016-2901-3`.
  Метаданные — p. 279 и последняя страница статьи 312.
- Прочитано: pp. 279–280, abstract и §1.1–1.2. Две плоские тонкие
  балки под прямым углом; передаточные условия зависят от внешних
  закреплений, узловой пограничный слой участвует в выводе.
  Не выдавать этот источник за общий finite-thickness спектральный тест.

### `nazarov_2002_plane_anisotropic_beams`

- [PDF](pdf/tm295.pdf), русский журнальный текст;
  SHA256 `6de01a4b55f116260b472927da64187255ee8f6807cad78b843fcb5ddf2eb0d1`.
- С. А. Назаров, А. С. Слуцкий, *Произвольные плоские системы
  анизотропных балок*. Труды Математического института им. В. А. Стеклова
  236 (2002), 234–261; метаданные — первая страница.
- Прочитано: pp. 234–235, abstract и §1.1–1.2. Классификация
  закреплённых, малоподвижных и подвижных элементов; при подвижных
  балках возникают алгебраические уравнения и нелокальные условия.
  Это отдельная публикация, не отсутствующий текст *Asymptotics of Natural
  Oscillations of Elastic Junctions with Readily Movable Elements*.

### `kolpakov_2014_joined_elastic_beams`

- [PDF](pdf/2014ZAMMAndrianov.pdf), издательская статья;
  SHA256 `3d5e753668e326d6b4e7c80f45d982c79f92929079694deb6765abb2a87fda23`.
- Alexander G. Kolpakov, Igor V. Andrianov, *Asymptotic decomposition
  in the problem of joined elastic beams*. ZAMM 94(10) (2014), 818–836;
  DOI `10.1002/zamm.201200278`. Метаданные — p. 818;
  published online 2013-06-17 отличается от года выпуска 2014.
- Прочитано: pp. 818–819, abstract, Introduction и начало §1.
  «Normal type» ограничен узлами размера порядка диаметра балок и
  упругих свойств того же порядка; мягкие слои и глубокие выточки отделены.
  Асимптотическое разделение глобальной 1D и локальной 3D задач не
  означает отсутствия локальных напряжений или универсальной точности
  point joint. Коэффициент k_theta для нашего углового узла не вычислен.

### `dockerty_1981_jointed_beam_stress`

- [PDF](pdf/0020-74032990052-7.pdf), скан журнальной статьи;
  SHA256 `b40caede3f453567875fc8d9f36590a3b55673fc9e571c511a7a11f05dbafe1a`.
- G. J. Dockerty, C. M. Leech, *Stress propagation through jointed beam
  systems using finite element theory*. International Journal of Mechanical
  Sciences 23(8) (1981), 457–471. Метаданные — p. 457;
  DOI не выводился из непрозрачного имени файла.
- Прочитано: pp. 457–458 и Conclusions на p. 471. Timoshenko/FEM
  с тремя моделями узла (rigid block, filament, flexible) и экспериментом;
  отмечена неточность малых transmitted shear signals. Фоновый пример
  зависимости результата от модели узла, не асимптотическое доказательство
  `SPRING -> RIGID` и не разрешение FEM-этапа в проекте.

## `rucka_2010_l_joint_guided_waves`

- PDF: `docs/literature/pdf/j.jsv.2009.12.004.pdf`, 20 PDF-страниц (1760–1779); SHA256
  `d15786c6aa2bf71731c63b16bf594739a626e6c45461e63689d90a36d472ce5f`.
- Тип: журнальная статья, издательский PDF с финальными томом и страницами.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_PUBLISHER`, 2026-10-05. Magdalena Rucka;
  *Experimental and numerical study on damage detection in an L-joint using guided wave propagation*;
  Journal of Sound and Vibration 329(10) (2010), 1760–1779;
  DOI `10.1016/j.jsv.2009.12.004`. Автор, название, том, страницы и
  available online 2009-12-30 — PDF p. 1; выпуск 10 и дата выпуска
  2010-05-10 уточнены у [издателя](https://doi.org/10.1016/j.jsv.2009.12.004).
- Роль: прямой M-H-Tim frame precedent с экспериментом: straight rod,
  затем steel L-joint, квадрат 6×6 mm, intact/notched specimens.
- Прочитано: §2.1–2.3, pp. 1761–1764 (PDF 2–5); §3, p. 1764;
  §4, pp. 1766–1768; §5, pp. 1768–1775; §6, pp. 1776–1778;
  refs. [34],[35], p. 1779. PDF 3–5,7 сверены по изображениям.
- Структура: independent contraction **ψ**, не латинская c; source φ —
  Timoshenko rotation. Четыре поля (u,ψ,v,φ); D/E/μ блочные, (8)–(13).
  Time-domain SEM, GLL quadrature; (26)–(28) преобразуют local в global
  и собирают frame. §5.2 p. 1770 и §5.3 p. 1774: mode conversion
  longitudinal→flexural и flexural→longitudinal соответственно.
- **SOURCE-SPECIFIC / NOT A PROJECT DEFAULT:** §2.1–2.2, p. 1762,
  K₁ᴹᴴ=1.1, K₂ᴹᴴ=2.1 — совместный least-squares fit axial wave velocities
  при 100 и 120 kHz; K₁ᵀⁱᵐ=.95 — fit flexural velocities при тех же
  частотах; K₂ᵀⁱᵐ=12K₁ᵀⁱᵐ/π² выбран по совпадению cutoff с Lamb modes,
  не отдельный экспериментальный fit. §4 p. 1766 также подбирает damping
  по амплитудам отражений. Это не independently calibrated low-frequency
  validation и не значения для будущей модели CoupledBeams.
- Ограничения: локальная блочность совместима с глобальной конверсией;
  source frame assembly не выбирает наши условия углового узла.
  Подробности и разграничение обозначений — в [тематической заметке](mindlin_herrmann_timoshenko_sources.md).

## `jang_2014_timoshenko_composite_patch_guided_waves`

- PDF: `docs/literature/pdf/j.compositesb.2013.12.050.pdf`, 13 PDF-страниц (248–260); SHA256
  `8d945ae828f5a794090df67b968ca84901ff7e3441e6a2d5a4435c025b5c1929`.
- Тип: журнальная статья, издательский PDF.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_PUBLISHER`, 2026-10-05. Injoon Jang,
  Ilwook Park, Usik Lee; *Guided waves in a Timoshenko beam with a bonded composite patch:
  Frequency domain spectral element modeling and analysis*;
  Composites Part B: Engineering 60 (2014), 248–260;
  DOI `10.1016/j.compositesb.2013.12.050`. Авторы, название, том, страницы,
  available online 2014-01-03 — PDF p. 1. Полное имя журнала и том без
  отдельного номера выпуска — [издатель](https://www.sciencedirect.com/science/article/abs/pii/S1359836813007804).
- Прочитано: §2.1–2.2, pp. 249–250 (PDF 2–3); Appendix A, pp. 258–259
  (PDF 11–12), (A1)–(A10); §5.1, pp. 254–255, Table 1, Figs. 6,8;
  описание §5.2 p. 255; Conclusions p. 258. PDF 2,11,12 сверены визуально.
- Роль: четыре поля isotropic base beam (1): u_b0, w_b0, θ_b (rotation),
  ψ_b (independent lateral contraction); `u_b=u_b0−zθ_b`, `w_b=w_b0+zψ_b`.
  По (A5)–(A7) strains и energy делятся на axial/M-H и bending/Timoshenko
  contributions с **разными reduced constitutive relations**.
- Appendix A: M-H обнуляет σ_yy, τ_xy, τ_yz; C₁₁*=E/(1−ν²),
  C₁₃*=Eν/(1−ν²), C₅₅*=G. Timoshenko дополнительно обнуляет σ_zz;
  Q₁₁*=E, Q₅₅*=G. κ_b включён в оба shear terms (5) и resultants
  (A8)–(A9); его численный default в проверенном тексте не установлен.
  Нельзя описывать это как применение единого полного 3D law без редукций.
- **Source warnings:** (A3), p. 258: справа напечатано σ_xx^(b) вместо
  ожидаемой strain ε_xx^(b), что несовместимо с размерностью Q₁₁*=E;
  (A10), p. 259: напечатано `A_b=b∫ z dz=b h_b` при симметричных пределах,
  лишний z несовместим с определением площади. В (A7) нет явного κ_b,
  в (5),(A8),(A9) он есть. Печать подтверждена изображениями; PDF сохранён,
  программные исправления не делались. [Точные записи и основания](mindlin_herrmann_timoshenko_sources.md#предупреждения-к-печати-appendix-a).
- Validation: SEM сравнивается с собственным 1D FEM и ANSYS 2D plane-stress
  на задачах metallic base + composite patch. Coupling шести полей после
  constraints (3) не доказывает обязательную axial-bending coupling bare
  isotropic symmetric beam. Это не validation нашей геометрии и не
  отмена результата предыдущего Bishop audit.

## `liu_2021_multibody_beams_rigid_bodies`

- PDF: `docs/literature/pdf/Journal Paper 122-accepted version.pdf`, 37 PDF-страниц
  (репозиторная обложка и accepted manuscript); SHA256
  `bcf7aa91672df9407bd9d624053c296dcdad6f4e166658a8ed6abdf169efe7a6`.
- Тип: журнальная статья; accepted manuscript с обложкой City Research Online.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_REPOSITORY`, 2026-10-05. Xiang Liu,
  Chengli Sun, J. Ranjan Banerjee, Han-Cheng Dan, Le Chang;
  *An Exact Dynamic Stiffness Method for Multibody Systems Consisting of Beams and Rigid-Bodies*;
  Mechanical Systems and Signal Processing 150 (2021), article 107264;
  DOI `10.1016/j.ymssp.2020.107264`. Имена и название — PDF p. 2,
  финальный том/год/article number — обложка. [Репозиторий](https://openaccess.city.ac.uk/id/eprint/26468/)
  отдельно указывает online 2020-09-16 и выпуск 2021-03-31; отдельный
  номер issue не указан. Не смешивать с `liu_2022_longitudinal_dynamic_stiffness`.
- Прочитано: title/abstract; §2.1 (PDF 8–9), описание §2.2–2.3 (PDF 9–16),
  §3.1 representative validations (PDF 19–21), conclusion (PDF 25–26),
  Appendix (PDF 26–30). Matrix (2) на PDF 9 сверена визуально.
- Роль: general DS assembly beams + rigid bodies; local axial/bending
  могут быть независимы, structural coupling возникает через coordinate
  transformations, относительные положения и assembly.
- Существенное разграничение: обзор §2.1 PDF 8 упоминает classical,
  Rayleigh–Love, **Rayleigh–Bishop** и higher-order bending. Но demonstration
  PDF 9, §3.1 PDF 21 и Appendix PDF 26–30 используют только axial
  **classical/Rayleigh–Love** и bending **Euler–Bernoulli/Timoshenko**.
  Статья не является прямой validation Bishop + Timoshenko или M-H.
- Validation у авторов: conventional DSM, опубликованные результаты и
  ANSYS FEM; здесь не воспроизводилась. Например k=1 в §3.1 — source
  input, не проектный default. Accepted manuscript не принят как готовый
  проверенный набор формул/API или условий нашего углового узла.

## `marais_2015_rayleigh_bishop_cylindrical_rod`

- PDF: `docs/literature/pdf/s13370-014-0286-3.pdf`, 12 PDF-страниц; SHA256
  `19fec24bfb5dda35d27141a09661eaf68226e59fe08153a3f1a1c441551e5739`.
- Тип: журнальная статья; издательский online-first PDF 2014 года без финальной пагинации.
- Роль: источник опубликованной контрольной задачи Рэлея—Бишопа для отдельного будущего воспроизведения.
- Что важно для CoupledBeams: линейные продольные колебания прямого изотропного стержня из цилиндрических участков; в §4 один конец закреплён, другой свободен, между участками заданы условия сопряжения.
- Обозначения / terminology: Rayleigh–Bishop, longitudinal/lateral displacements, Green function; эти термины не задают модель углового узла проекта.
- Критично смотреть: §2–3, PDF pp. 2–8, для постановки; §4, PDF pp. 9–10, (22)–(25) и Fig. 2 для двухсекционного примера. `123` внизу листа — издательская метка, не номер страницы.
- Прочитано: первая страница, abstract, начало introduction; постановка и численный пример §4. Визуально просмотрена PDF p. 9.
- Замечание по применимости: прямой осевой стык не определяет условия нашего углового узла. В (23) визуально подтверждено повторение `cos`; возможные опечатки (23)–(24) и согласованность (25) требуют отдельной проверки, исправления здесь не предлагаются. [Карта вопросов](longitudinal_rod_models_sources.md#вопросы-для-отдельной-проверки).
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. J. Marais, I. Fedotov, M. Shatalov; Afrika Matematika 26(7–8) (2015), 1549–1560; DOI `10.1007/s13370-014-0286-3`. Авторы и название — по PDF; том, выпуск, страницы и online 2014-11-27 — по [Crossref](https://api.crossref.org/works/10.1007/s13370-014-0286-3). Citation year — 2015.
- Последующая проверка 2026-10-05: PDF pp. 3–5, 7–10 сверены по изображениям.
  На p. 9 в (23) повтор cos заменён независимым sin; в (24) восстановлен
  квадрат q² под внутренним радикалом, q=(EA−Jω²)/(2H). Основание —
  характеристический многочлен (11), размерность и подстановка, а не числа
  публикации. В (2),(3), PDF p. 3, отсутствуют пространственные интегралы;
  интегральные энергии согласованы с (18)–(21), pp. 7–8. Γ′y=−N по (12),
  p. 5; C-F и непрерывность U,U′,N,P подтверждены по (25), p. 9.
  Независимые расчёты проходят, но с четырьмя значащими цифрами на p. 10
  совпадает только мода 2 из пяти (при округлении и усечении). Причина
  остальных расхождений не установлена, параметры не подбирались.
  [Исходная печать, реконструкция, результаты и ограничения](../theory/bishop_literature_reproduction.md).

## `popov_sadovsky_2021_longitudinal_models_experiment`

- PDF: `docs/literature/pdf/grigory_ne,+08.pdf`, 12 PDF-страниц (270–281); SHA256
  `5b6f65855cceb13d5aa00f370f940d4ee78f5dc4c5bae206ba546665009d7979`.
- Тип: журнальная статья; издательский русский оригинал с английской аннотацией.
- Роль: источник экспериментальных частот и вопросов идентификации параметров продольной модели.
- Что важно для CoupledBeams: линейные продольные колебания длинного цилиндрического стержня; сопоставление волновой модели, поправки Релея и модели Бишопа с экспериментом, оценка скорости волн и коэффициента Пуассона.
- Обозначения / terminology: `c`, `nu`, частотные отношения; свободные концы, жёсткое защемление и свободный упор в стенки различаются.
- Критично смотреть: таблица на p. 273; §5, pp. 274–275, идентификация; §6, pp. 276–278, (10)–(15) и рис. 5.
- Прочитано: титульная страница, аннотация, начало введения; таблица p. 273 и текст о подборе параметров и условиях на pp. 275–278.
- Замечание по применимости: условия (10), (11) и условия перед формой (14) нельзя подменять друг другом; напечатанная точность таблицы не является подтверждённой точностью исходных измерений. Ссылки на номера формул требуют проверки, см. [тематическую заметку](longitudinal_rod_models_sources.md#вопросы-для-отдельной-проверки).
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. А. Л. Попов, С. А. Садовский; Вестник Санкт-Петербургского университета. Математика. Механика. Астрономия 8 (66), №2 (2021), 270–281; DOI `10.21638/spbu01.2021.207` — по блоку цитирования PDF. [Английский перевод](https://doi.org/10.1134/S1063454121020114): Vestnik St. Petersburg University, Mathematics 54(2) (2021), 162–170; его PDF отсутствует. Это связанная версия той же работы, не второй независимый эксперимент и не второй canonical key.
- Последующая проверка 2026-10-05: с. 272–278 сверены по изображениям.
  Все 30 строк таблицы с. 273 сохранены, включая `10.998` при n=11;
  это число не признано опечаткой и не сглажено. f30=38557.99 Hz на с. 274
  хранится отдельно от произведения округлённого отношения на f1.
  (5),(9)–(15) проверены независимо собранной граничной системой; минимум
  (6) по округлённой таблице ν=.336842781730 округляется до .337.
  В (5) A в fwn трактуется как c по (1),(3) и с. 273, не как площадь F.
  Ссылки «уравнение (6)» перед (14) и «алгоритм (7)» на с. 277 не
  соответствуют содержимому этих формул. На с. 278 ссылка (11) в
  утверждении о близости не подтверждается для C-C: при ν=.34 разности
  с (15) 2.36–69.80 Hz. Для F-F (10) максимум .5175 Hz; замена номера
  ссылки остаётся интерпретацией, авторского erratum нет. Ранжирование
  теорий воспроизведено, отдельные столбцы Fig. 5 визуально расходятся;
  цифровые данные рисунка и неокруглённые измерения недоступны.
  [Полный отчёт](../theory/bishop_literature_reproduction.md).

## `liu_2022_longitudinal_dynamic_stiffness`

- PDF: `docs/literature/pdf/Manuscript DSM of rods.pdf`, 30 PDF-страниц (обложка + рукопись pp. 1–29); SHA256
  `ecd7ac4df12f236b1ea7006bb5b0d7a913af029ef84fa300f6e26062d2cd2851`.
- Тип: журнальная статья; accepted manuscript с обложкой City Research Online.
- Роль: источник контрольных продольных задач и метода динамической жёсткости.
- Что важно для CoupledBeams: линейные свободные колебания стержней и ферм по классической, Rayleigh–Love, Rayleigh–Bishop и Mindlin–Herrmann теориям; несколько типов концевых условий и ступенчатые стержни.
- Обозначения / terminology: DS/DSM, Wittrick–Williams, `J0`; слово exact относится к принятой в статье стержневой теории.
- Критично смотреть: §2.3, рукопись pp. 7–9; §3.2, p. 12; §4.1 и Table 1, pp. 14–15 (PDF pp. 15–16); §4.3–4.5 для условий и сборок. Пагинация рукописи не совпадает с журнальной.
- Прочитано: обложка, abstract, начало §1; описание контрольного сравнения и Table 1 на p. 15 рукописи.
- Замечание по применимости: продольные элементы и шарнирная ферма не задают изгибно-продольное сопряжение нашего углового узла; значения Table 1 ещё не воспроизведены.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. Xiang Liu, Yaxing Zhao, Wei Zhou, J. Ranjan Banerjee; Applied Mathematical Modelling 104 (2022), 401–420; DOI `10.1016/j.apm.2021.11.023`. Год и страницы по репозиторной обложке и Crossref; online 2021-12-24 по [City Research Online](https://openaccess.city.ac.uk/id/eprint/28102/). Заголовок рукописи `01 (2021) 1–29` не используется как финальные выходные данные.

## `shatalov_2011_longitudinal_rod_theories`

- PDF: `docs/literature/pdf/InTech-Longitudinal_vibration_of_isotropic_solid_rods_from_classical_to_modern_theories.pdf`, 30 PDF-страниц (28 страниц главы + библиографическая страница книги + лицензия); SHA256
  `603d87922f24bd145536981412c4c6078cec10f6dd15447abc677c6255c8ce57`.
- Тип: глава 10 в коллективной книге; издательский PDF.
- Роль: навигационный источник по продольным теориям изотропных стержней.
- Что важно для CoupledBeams: иерархия линейных моделей осесимметричного движения от классической до уточнённых, поперечная инерция и сдвиг; постановки распространения волн и колебаний рассматриваются в разных разделах.
- Обозначения / terminology: unimodal/multimodal theories, Rayleigh–Love, Rayleigh–Bishop, Mindlin–Herrmann, Pochhammer–Chree.
- Критично смотреть: §2, p. 190; §3.1, p. 192; §3.2, p. 200; §4, p. 205; §5, p. 207. Первая страница главы не имеет печатного номера; диапазон 187–214 определяется последовательной пагинацией, начиная со следующей страницы 188.
- Прочитано: заголовок, начало §1 (отдельного abstract нет), библиографическая страница книги; заголовки тематических разделов просмотрены для навигации.
- Замечание по применимости: теория прямого изотропного стержня; конкретные условия и область применимости каждой модели нужно проверять отдельно, без переноса на угловой узел по названию теории.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_PUBLISHER`, 2026-10-05. Michael Shatalov, Julian Marais, Igor Fedotov, Michel Djouosseu Tenkam; *Advances in Computer Science and Engineering*, ed. Matthias Schmidt, InTech, 2011, ch. 10, pp. 187–214; ISBN `978-953-307-173-2`. Данные книги и online 2011-03-22 — PDF p. 29; DOI `10.5772/15662` подтверждён [издателем](https://www.intechopen.com/chapters/14403) и Crossref. Полное имя четвёртого автора сохранено по PDF: Crossref опускает `Tenkam`.

## `fedotov_2010_rayleigh_bishop_rod`

- PDF: `docs/literature/pdf/Michel.pdf`, 6 PDF-страниц (609–614); SHA256
  `531c8643fd023f9a041a8d92f3ee13a815adc2039cbc2d3d66ccbfd6118a7597`.
- Тип: журнальная статья; издательский английский перевод.
- Роль: источник продольной модели Рэлея—Бишопа и опубликованных контрольных примеров.
- Что важно для CoupledBeams: линейные свободные и вынужденные продольные колебания изотропного цилиндрического и конического стержней; на концах в (2) заданы `u = u_xx = 0`.
- Обозначения / terminology: Rayleigh–Bishop rod, Green function, two types of orthogonality; `fixed ends` читать вместе с (2).
- Критично смотреть: постановка, pp. 609–610, (1)–(7); Example 1, p. 612, и Example 2, p. 613.
- Прочитано: первая страница и введение без отдельного abstract, постановка на p. 610; просмотрены расположение и описание примеров.
- Замечание по применимости: `fixed ends` в этой работе не означает автоматически `u = u_x = 0`; это осевая задача прямого стержня. Формулы ортогональности и примеры независимо не проверялись.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. I. A. Fedotov, A. D. Polyanin, M. Yu. Shatalov, H. M. Tenkam; Doklady Physics 55(12) (2010), 609–614; DOI `10.1134/S1028335810120062`. Первая страница указывает русский оригинал: Доклады Академии наук 435(5) (2010), 613–618; его PDF отсутствует. Оригинал и перевод связаны одной canonical записью.

## `banerjee_2019_rayleigh_love_timoshenko`

- PDF: `docs/literature/pdf/Rayleigh_Love accepted version.pdf`, 30 PDF-страниц (обложка + рукопись pp. 1–29); SHA256
  `26db797736decf81428d3f990afac5a9b48a5ee60f06cfc31eb6b52b052fcd48`.
- Тип: журнальная статья; accepted manuscript с обложкой City Research Online.
- Роль: методический источник по сочетанию продольной и изгибной динамической жёсткости.
- Что важно для CoupledBeams: линейные свободные колебания однородных и ступенчатых стержней и плоской рамы; продольная часть Rayleigh–Love, изгибная — Timoshenko.
- Обозначения / terminology: dynamic stiffness, Wittrick–Williams, Rayleigh–Love bar; не отождествлять с Rayleigh–Bishop по сходству названий.
- Критично смотреть: §2.1, рукопись p. 6; §2.3, p. 15; §3, p. 16; §4.1, p. 19 (PDF p. 20), закреплённый с двух концов и консольный стержни; §4.2–4.3, pp. 21–25, ступенчатый стержень и рама.
- Прочитано: обложка, abstract, начало introduction, вводное описание трёх примеров и §4.1 на p. 19 рукописи.
- Замечание по применимости: источник по Rayleigh–Love, не готовое обоснование дополнительных условий Bishop; условия нашего узла остаются невыбранными.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. J. R. Banerjee, A. Ananthapuvirajah; International Journal of Mechanical Sciences 150 (2019), 337–347; DOI `10.1016/j.ijmecsci.2018.10.012`. Финальные данные по обложке и Crossref; online 2018-10-10 по [City Research Online](https://openaccess.city.ac.uk/id/eprint/20916/). Журнальные и рукописные страницы различаются.
- Дополнительное чтение для M-H/Timoshenko navigation, 2026-10-05:
  конец §2, рукопись p. 5 (PDF 6); (1),(2), p. 6 (PDF 7);
  (22)–(28), pp. 9–10 (PDF 10–11); §2.3 (65),(66), p. 15 (PDF 16);
  §4.3, pp. 24–25 (PDF 25–26). Axial Rayleigh–Love и Timoshenko bending
  явно приняты uncoupled, treated independently; отдельно полученные
  DSM объединяются simple superposition в 6×6 plane-frame matrix.
  Роль — **published precedent for modular assembly**, не источник M-H
  equations. В frame example k=2/3 — заданный source parameter.
  Citation key, PDF и результаты прежнего Bishop audit сохранены.
  [Сравнение с M-H formulations](mindlin_herrmann_timoshenko_sources.md).

## `georgiades_2017_nonlinear_l_shaped_beams`

- PDF: `docs/literature/pdf/j.euromechsol.2017.03.007.pdf`, 32 PDF-страницы (91–122); SHA256
  `875ea5d8ac0f57563d40fc7cb7e0ac9d06977e43511838ce323b61c523dea895`.
- Тип: журнальная статья; издательский PDF.
- Роль: источник постановки пространственной нелинейной задачи для L-образной балки; направление проекта отложено до обсуждения с руководителем.
- Что важно для CoupledBeams: изотропные нерастяжимые балки Euler–Bernoulli, один внешний конец закреплён, другой свободен; внутриплоскостной изгиб, внеплоскостной изгиб и кручение. Заявлен вывод уравнений и условий до второго порядка с вращательной инерцией; это работа о постановке, не выбранный вынужденный режим нашего проекта.
- Обозначения / terminology: inextensional, global displacements, rotary inertia, second-order nonlinearity.
- Критично смотреть: §2, p. 92 и Fig. 1, для допущений; §3, начиная с p. 96, для уравнений; Table 1 и §6 на p. 108 для структуры результатов.
- Прочитано: abstract, introduction и начало §2, pp. 91–92; заголовки дальнейших разделов просмотрены.
- Замечание по применимости: нерастяжимая консольная L-схема не заменяет растяжимую модель двух стержней проекта. Утверждения автора о полноте нелинейной модели здесь не проверены. [Уточнение об инерции](nonlinear_inplane_outofplane_sources.md#проверенные-описания-и-границы-проверки).
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. Fotios Georgiades; European Journal of Mechanics A/Solids 65 (2017), 91–122; DOI `10.1016/j.euromechsol.2017.03.007`; online 2017-03-18 напечатана на первой странице.

## `warminski_2008_autoparametric_beam_structure`

- PDF: `docs/literature/pdf/j.jsv.2008.01.048.pdf`, 23 PDF-страницы (486–508); SHA256
  `ed8ac4ac78d02fbacbbcbb2c68eda57b51a87393e82d06360433096454bb3c2b`.
- Тип: журнальная статья; издательский Article-in-Press PDF с финальными томом и страницами.
- Роль: источник модели и эксперимента по пространственным автопараметрическим движениям сопряжённых балок.
- Что важно для CoupledBeams: две стеклоэпоксидные прямоугольные балки в L-схеме, заделка B, соединение C и присоединённые массы C/A; нерастяжимость, изгиб в двух плоскостях и кручение. Заявлены уравнения до третьего порядка; эксперимент включает случайное и гармоническое возбуждение.
- Обозначения / terminology: autoparametric, HOT (higher order terms), torsional inertia, Galerkin mode shapes.
- Критично смотреть: §2–3, pp. 488–491, (3)–(16); условия (33)–(44), pp. 493–494; §4 и Table 2/Fig. 5, pp. 496–498; §5, pp. 499–501, (56)–(58) и Fig. 8.
- Прочитано: abstract, начало introduction; описание модели и инерционных допущений на pp. 488–491, пояснение HOT на p. 494, описание FEM/эксперимента на p. 498 и выбора форм на pp. 500–501.
- Замечание по применимости: крутильная инерция сохраняется; нельзя пересказывать упрощение как полное отсутствие вращательной инерции. HOT скрывает члены второго и третьего порядков в напечатанных условиях; это ограничение явной записи, не установленная ошибка. Присоединённые массы и композит отличаются от текущей изотропной схемы. [Подробности](nonlinear_inplane_outofplane_sources.md#проверенные-описания-и-границы-проверки).
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. J. Warminski, M. P. Cartmell, M. Bochenski, Ivelin Ivanov; Journal of Sound and Vibration 315(3) (2008), 486–508; DOI `10.1016/j.jsv.2008.01.048`. Первая страница проверена визуально; в ней фамилии напечатаны без диакритики. Выпуск — Crossref; online 2008-03-05 — PDF.

## `bux_roberts_1986_coupled_beam_interactions`

- PDF: `docs/literature/pdf/0022-460x2990304-4.pdf`, 24 PDF-страницы (497–520); SHA256
  `55b0ff224a6e9d31779ac60bdb5f17fe09d6f0afd8fd20fbaf9ab78e8579236b`.
- Тип: журнальная статья; скан журнальных страниц с OCR-слоем.
- Роль: источник сокращённой модели и эксперимента по нелинейным взаимодействиям сопряжённых балок.
- Что важно для CoupledBeams: пара балок под прямым углом с поперечным гармоническим возбуждением основной балки; внутриплоскостной изгиб основной балки и изгибно-крутильное движение присоединённой, квадратичная нелинейная связь и внутренние резонансы.
- Обозначения / terminology: autoparametric interaction, finite degree-of-freedom model, combination internal resonance.
- Критично смотреть: Fig. 1, p. 498; §2 `System equations of motion`, начиная с p. 500, особенно §2.1 `Kinematics`.
- Прочитано: первая страница визуально, abstract и introduction; начало §2.1, p. 500.
- Замечание по применимости: движение основной балки ограничено плоскостью, вторичная балка упрощена дискретной инерцией; это не полная пространственная распределённая модель двух растяжимых стержней.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. S. L. Bux, J. W. Roberts; Journal of Sound and Vibration 104(3) (1986), 497–520; DOI `10.1016/0022-460X(86)90304-4` подтверждён [Crossref](https://api.crossref.org/works/10.1016/0022-460X(86)90304-4); прочие поля — первая страница PDF.

## `ho_scott_eisley_1976_nonplanar_free_motions`

- PDF: `docs/literature/pdf/0022-460x2990943-3.pdf`, 7 PDF-страниц (333–339); SHA256
  `28a5e241a3f9ba9b4e8e9fb91a75d8cc6b076985df5cb005ee9324b0cf58934f`.
- Тип: журнальная статья; скан с OCR-слоем.
- Роль: источник по свободным нелинейным колебаниям с движением в двух плоскостях.
- Что важно для CoupledBeams: одиночная прямая балка с близкими главными моментами инерции сечения; изгиб в двух плоскостях, продольная инерция и эффекты Пуассона отброшены, продольное перемещение исключено.
- Обозначения / terminology: non-planar free motions, whirling, Galerkin approximation.
- Критично смотреть: §2, pp. 333–334, особенно допущения, simply supported ends перед (3)–(4) и сокращение к (7)–(8).
- Прочитано: первая страница визуально, abstract, introduction и §2 на pp. 333–334.
- Замечание по применимости: слово `fixed` в abstract не следует переводить как изгибную заделку: на p. 334 явно используются шарнирные изгибные условия. Внешнее возбуждение и демпфирование в рассматриваемой свободной задаче обнулены; углового узла нет.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. C.-H. Ho, R. A. Scott, J. G. Eisley; Journal of Sound and Vibration 47(3) (1976), 333–339; DOI `10.1016/0022-460X(76)90943-3` подтверждён [Crossref](https://api.crossref.org/works/10.1016/0022-460X(76)90943-3); остальные поля — первая страница.

## `pai_nayfeh_1990_nonplanar_cantilever`

- PDF: `docs/literature/pdf/0020-74622990012-x.pdf`, 20 PDF-страниц (455–474); SHA256
  `8a34a1e9f21ad11d88fb1cfbe9043027be43efbb9c9ca72e61e8f7dda29493b3`.
- Тип: журнальная статья; скан с OCR-слоем.
- Роль: источник по вынужденным нелинейным пространственным движениям консольной балки.
- Что важно для CoupledBeams: однородная прямоугольная консоль с гармоническим поперечным движением основания; два изгибных перемещения, кубические геометрические и инерционные члены, вязкое демпфирование в принятой модели.
- Обозначения / terminology: lateral base excitation, one-to-one internal resonance, inextensibility constraint.
- Критично смотреть: §2 `Problem formulation`, pp. 456–457, допущения (a)–(d), Fig. 1, уравнения (1) и условия (2).
- Прочитано: первая страница визуально, abstract, introduction и постановка §2 на pp. 456–457.
- Замечание по применимости: распределённой крутильной инерцией пренебрегают при оговорённом разделении частот (§2, p. 456); речь об одиночной консоли, а не сопряжении двух растяжимых стержней.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-05. Perngjin F. Pai, Ali H. Nayfeh; International Journal of Non-Linear Mechanics 25(5) (1990), 455–474; DOI `10.1016/0020-7462(90)90012-X` подтверждён [Crossref](https://api.crossref.org/works/10.1016/0020-7462(90)90012-X). Диапазон страниц проверен визуально: OCR искажает 474.

## `tan_ko_2004_connection_dampers`

- PDF: `docs/literature/pdf/tan_ko_2004_connection_dampers.pdf`, 24 страницы
  (журнальные 707–730), SHA256
  `f2b4ca729f82633f785dc1f688126f865319200c4aafc42dce1d83ae1db681de`.
- Тип: журнальная статья; локальный PDF — издательская версия с некорректно
  извлекаемым текстовым слоем, поэтому титульные данные прочитаны визуально.
- Роль: потенциальный источник по демпферам в соединениях длиннопролётных балок.
  Метка будущего чтения: `service-vibration control`.
- Что важно для CoupledBeams: рассматриваются соединения балки с колоннами,
  содержащие вязкоупругие демпферы, экспериментальная постановка и модель
  с дробными производными для вертикальных колебаний.
- Обозначения / terminology: beam–column connection, viscoelastic damper,
  fractional derivative model.
- Критично смотреть: abstract и §1 `Introduction`, pp. 707–708;
  §2 `Experiments on a Beam with Various Connections`, начиная с p. 708,
  особенно §2.1 `Description of the Beam–Column Connections`.
- Замечание по применимости: речь о вертикальной вибрации балки между колоннами;
  дробная вязкоупругая модель не тождественна локальному Kelvin–Voigt-узлу проекта.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_PUBLISHER`, 2026-10-01. X. M. Tan,
  J. M. Ko; *Vibration Control of Long-Span Beams: Experimental and Analytical
  Study of Beam Structures Incorporated with Connection Dampers*;
  Journal of Vibration and Control 10(5) (2004), 707–730;
  DOI `10.1177/1077546304040132`. Первая и последняя страницы локального PDF;
  номер выпуска подтверждён [SAGE](https://journals.sagepub.com/doi/10.1177/1077546304040132).

## `hsu_fafitis_1992_viscoelastic_connections`

- PDF: `docs/literature/pdf/hsu_fafitis_1992_viscoelastic_connections.pdf`,
  16 страниц (журнальные 2459–2474), SHA256
  `4c030c71d6848342d410332d6927333b009c66e0d571464fc2de6655db7ce418`.
- Тип: журнальная статья; локальная копия — скан журнальных страниц с OCR
  и отметкой загрузки ASCE. Фамилия Fafitis сверена с изображением первой страницы.
- Роль: потенциальный источник по вязкоупругим соединениям рам.
  Метка будущего чтения: `seismic application`.
- Что важно для CoupledBeams: представлены эластомерное устройство в соединении,
  его экспериментальное описание и модель типа Kelvin–Voigt для анализа рам.
- Обозначения / terminology: connection isolator, elastomeric pad,
  Kelvin–Voigt-type model, durometer hardness.
- Критично смотреть: abstract и `Introduction`, p. 2459;
  `Experimental Data`, `Connection Modeling` и Fig. 1, p. 2460.
- Замечание по применимости: рассматривается сейсмическое возбуждение рам
  с эластомерными соединениями; соответствие устройства идеализированному
  вращательному узлу двух стержней требует отдельного чтения.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_PUBLISHER`, 2026-10-01. Sheng-Yung Hsu,
  Apostolos Fafitis; *Seismic Analysis Design of Frames with Viscoelastic
  Connections*; Journal of Structural Engineering 118(9) (1992), 2459–2474.
  Заголовок, авторы, выпуск и страницы — локальный PDF;
  DOI `10.1061/(ASCE)0733-9445(1992)118:9(2459)` подтверждён
  [ASCE](https://ascelibrary.org/doi/10.1061/%28ASCE%290733-9445%281992%29118%3A9%282459%29)
  и Crossref, на титульной странице не напечатан.

## `song_hong_2007_nonconservative_joints`

- PDF: `docs/literature/pdf/song_hong_2007_nonconservative_joints.pdf`,
  15 PDF-страниц: статья на журнальных pp. 15–28 (PDF pages 1–14), затем
  рекламная страница Hindawi; SHA256
  `02f6793419426c94400dbf7a32d197e833d84db1a76a03771ff7f8415a395270`.
- Тип: журнальная статья; локальный PDF — издательская верстка IOS Press
  с добавленной рекламной страницей, не отдельная версия статьи.
- Роль: потенциальный источник по диссипативным соединениям балочных сетей.
  Метка будущего чтения: `mathematical joint model`.
- Что важно для CoupledBeams: рассматриваются три пары пружин и демпферов
  в плоском узле и передача колебательной энергии между балками под углом.
- Обозначения / terminology: non-conservative joint, spring–dashpot model,
  energy flow analysis (EFA), wave transmission.
- Критично смотреть: abstract и §1 `Introduction`, pp. 15–16;
  §2 `Wave transmission analysis of beam networks with compliant and
  dissipative joints` и Fig. 1, начиная с p. 16.
- Замечание по применимости: основной предмет — энергия и интенсивность
  колебаний в среднем и высоком частотных диапазонах; это не готовый
  источник низкочастотной модальной задачи CoupledBeams.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-01. Jee-Hun Song,
  Suk-Yoon Hong; *Development of non-conservative joints in beam networks
  for vibration energy flow analysis*; Shock and Vibration 14(1) (2007), 15–28.
  Первая страница и конечная пагинация — локальный PDF; DOI
  `10.1155/2007/273472` и номер выпуска — [Crossref](https://api.crossref.org/works/10.1155/2007/273472).
  Год 2007 сохранён по журнальному заголовку и copyright. Дата `published`
  в Crossref, 2005-03-30, совпадает с напечатанной датой `Received` и
  не использована как год публикации. Это самостоятельная работа,
  отличная от `hong_kim_1999_damped_timoshenko_joints`.

## `zeng_2026_piecewise_space_truss`

- PDF: `docs/literature/pdf/zeng_2026_piecewise_space_truss.pdf`, 17 страниц,
  SHA256 `906fc9a717f5a6404cafab639dabdb8446348aaac5a5284d00393f6572dc69c0`.
- Тип: журнальная статья, Research Article; локальный PDF — издательская
  open-access версия Wiley.
- Роль: потенциальный источник по сегментированным космическим фермам
  с податливыми диссипативными соединениями. Метка будущего чтения: `beam/frame application`.
- Что важно для CoupledBeams: рассматривается эквивалентный кусочно-однородный
  стержень с сосредоточенными осевыми пружинами и демпферами в соединениях.
- Обозначения / terminology: large space truss structures (LSTS),
  piecewise elastic bar, axial spring–damper joint.
- Критично смотреть: abstract и §1 `Introduction`, pp. 1–2;
  §2 `Model Development`, §2.1 `Physical Model and Assumptions`, p. 3.
- Замечание по применимости: модель ограничена продольным движением
  эквивалентного стержня; изгиб, кручение и их связь с продольным движением
  не входят в описанную постановку.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-01. Xianhong Zeng,
  Yin Zhang, Hongcheng Chen, Han Wu, Zundi Huang; *Axial Vibration
  Characteristics of Piecewise Large Space Truss Structures*;
  International Journal of Aerospace Engineering 2026(1), Article ID 2258617;
  DOI `10.1155/ijae/2258617`. Первая страница, пагинация `1 of 17`–`17 of 17`
  и отметка Wiley в PDF; номер выпуска также подтверждён
  [Crossref](https://api.crossref.org/works/10.1155/ijae/2258617).

## `xu_2023_right_angle_viscoelastic_damper`

- PDF: `docs/literature/pdf/xu_2023_right_angle_viscoelastic_damper.pdf`,
  22 страницы, SHA256
  `8f7159b26e65e684658b0b31da6c11cb6443b04e51d575b1cc3c6561e0f2561f`.
- Тип: журнальная статья, Research Article; локальный PDF — издательская
  open-access версия Hindawi/Wiley.
- Роль: потенциальный источник по конструкции и испытаниям узлового демпфера.
  Метка будущего чтения: `experimental joint/device`.
- Что важно для CoupledBeams: рассматривается угловой вязкоупругий
  демпфер с полиуретаном для соединения балки с колонной; описаны испытания
  материала и устройства, а также рамная постановка.
- Обозначения / terminology: right-angle viscoelastic damper (RVD),
  polyurethane, dynamic mechanical analysis (DMA).
- Критично смотреть: abstract и §1 `Introduction`, pp. 1–2;
  Fig. 1 и §2 `Dynamic Thermodynamic Testing of Polyurethane Rubber`, p. 2.
- Замечание по применимости: работа посвящена конкретному устройству
  и свойствам материала, включая нелинейное поведение; сведение к линейному
  Kelvin–Voigt-узлу проекта отдельно не проверялось.
- Метаданные: `VERIFIED_LOCAL_PDF`, 2026-10-01. Jun-Hong Xu, Zhe-Yu Zhu,
  Guang-Dong Zhou, Hao Wang, Ai-Qun Li; *Dynamic Characteristics of a Novel
  Right-Angle Viscoelastic Damper (RVD) Using Polyurethane Damping Materials*;
  Structural Control and Health Monitoring, Volume 2023, Article ID 2568963,
  22 pages; DOI `10.1155/2023/2568963`. Все основные поля и дата
  `Published 8 February 2023` напечатаны на первой странице.

## `beshara_keane_1997_dissipative_beam_joints`

- PDF: `docs/literature/pdf/beshara_keane_1997_dissipative_beam_joints.pdf`,
  19 страниц (журнальные 321–339), SHA256
  `37313cc9c70a04d58fb15cf4d2d9c2661fce363749fbb7bb289a173dcc2e44c5`.
- Тип: журнальная статья; локальная копия — скан журнальной верстки
  Academic Press без извлекаемого текста.
- Роль: потенциальный источник по податливым диссипативным узлам балочных сетей.
  Метка будущего чтения: `mathematical joint model`.
- Что важно для CoupledBeams: рассматривается плоская сеть упругих балок
  с тремя парами пружин и демпферов в узле и углом между соединяемыми балками.
- Обозначения / terminology: compliant and dissipative joint, receptance,
  vibrational energy flow, spring–dashpot model.
- Критично смотреть: аннотацию на p. 321 и §1 `Introduction`, pp. 321–322,
  для исходных предположений об узле и формулировки задачи передачи энергии.
- Замечание по применимости: анализ потоков энергии через узлы не тождествен
  задаче о комплексных собственных частотах двух закреплённых стержней.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_CROSSREF`, 2026-10-01. M. Beshara,
  A. J. Keane; *Vibrational energy flows in beam networks with compliant
  and dissipative joints*; Journal of Sound and Vibration 203(2) (1997), 321–339.
  Титульные данные и конечная страница проверены визуально;
  DOI `10.1006/jsvi.1996.0889` подтверждён
  [Southampton ePrints](https://eprints.soton.ac.uk/id/eprint/21080) и Crossref.
  Репозиторная карточка помечает файл как `Accepted Manuscript`, однако
  локальный скан содержит журнальные колонтитулы, пагинацию и copyright;
  здесь зафиксирован непосредственно наблюдаемый вид копии.

## `attarnejad_pirmoz_2014_damped_semirigid_frames`

- PDF: `docs/literature/pdf/attarnejad_pirmoz_2014_damped_semirigid_frames.pdf`,
  9 страниц (журнальные 165–173), SHA256
  `58662093bac8a4f8011ee3b8c1c693b54ff4c883ac95ec17172c347e47feb724`.
- Тип: журнальная статья; локальный PDF — издательская версия Elsevier.
- Роль: потенциальный источник по моделированию рам с податливыми соединениями
  и вращательными демпферами. Метка будущего чтения: `beam/frame application`.
- Что важно для CoupledBeams: рассматриваются балки Эйлера–Бернулли,
  нелинейные вращательные пружины и параллельные вращательные демпферы
  с учётом взаимодействия момента и поперечной силы в соединении.
- Обозначения / terminology: partially restrained (PR) connection,
  moment–shear interaction (MVI), rotational damper.
- Критично смотреть: abstract и §1 `Introduction`, p. 165;
  §2 `Basic assumptions`, pp. 165–166;
  §4 `Modeling of the flexible connections` и Fig. 1, p. 166.
- Замечание по применимости: предварительно нагруженные рамы и нелинейная
  характеристика соединения отличаются от линейного узла текущей модели;
  научное сопоставление ещё не выполнено.
- Метаданные: `VERIFIED_LOCAL_PDF`, 2026-10-01. Reza Attarnejad, Akbar Pirmoz;
  *Nonlinear analysis of damped semi-rigid frames considering moment–shear
  interaction of connections*; International Journal of Mechanical Sciences
  81 (2014), 165–173; DOI `10.1016/j.ijmecsci.2014.02.016`.
  Первая страница и XMP; номер выпуска не указан ни в PDF, ни в Crossref.
  Это самостоятельная работа, отличная от Failla 2014.

## `xu_zhang_2001_connection_dampers`

- PDF: `docs/literature/pdf/xu_zhang_2001_connection_dampers.pdf`,
  12 страниц (журнальные 385–396), SHA256
  `81f58cd207d4fc0dfb8cc9c149ead13b2d19588778e0177b5ce85aa94df289b4`.
- Тип: журнальная статья; локальный PDF — издательская версия Elsevier.
- Роль: потенциальный источник по рамным моделям с демпфированием в соединениях.
  Метка будущего чтения: `seismic application`.
- Что важно для CoupledBeams: рассматриваются балочные элементы с
  вращательными пружинами и демпферами на концах, комплексный модальный
  анализ и отклик стальной рамы на сейсмическое возбуждение.
- Обозначения / terminology: connection damper, rotational spring,
  modal-damping ratio, generalized pseudo-excitation method.
- Критично смотреть: abstract и §1 `Introduction`, pp. 385–386;
  §2 `Formulation of element matrices`, Figs. 1–2, начиная с p. 386.
- Замечание по применимости: рамная система и вынужденный сейсмический
  отклик отличаются от собственной модальной задачи для двух стержней;
  сопоставление параметров узла оставлено следующему этапу.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_PUBLISHER`, 2026-10-01. Y. L. Xu,
  W. S. Zhang; *Modal analysis and seismic response of steel frames with
  connection dampers*; Engineering Structures 23(4) (2001), 385–396.
  Первая страница и конечная пагинация — локальный PDF; DOI
  `10.1016/S0141-0296(00)00062-6` и номер выпуска подтверждены
  [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0141029600000626).
  Журнальный год — 2001; copyright 2000 не использован как год публикации.

## `failla_2014_viscoelastic_discontinuous_beams`

- PDF: `docs/literature/pdf/failla2014.pdf`, 12 страниц (журнальные 52–63),
  SHA256 `0d10c20796b1d35368e11aaad57764ae68b239b3f493a01361ea01782e3521e2`.
- Тип: журнальная статья Giuseppe Failla, *On the dynamics of viscoelastic
  discontinuous beams*, Mechanics Research Communications 60 (2014), 52–63.
  DOI `10.1016/j.mechrescom.2014.06.001`; номер выпуска в PDF не указан.
- Роль: внешний численный reference для EB, локальных rotational Kelvin–Voigt
  joints, translational supports и complex modal analysis. Независимый от K12
  опубликованный результат; общий численный корректор не является второй
  независимой реализацией всего CoupledBeams.
- Обозначения: `u,theta,M,S`, `xi=x/L`, `psi=U/L`, `mu=M*L/EI`, `T=S*L²/EI`.
  Source `psi` — прогиб, а не проектный поворот; source `mu` — момент, не разность длин.
  Из (1),(6) скачок S равен **минус V**, хотя V назван shear-force discontinuity;
  из (7)–(9) скачок theta равен `-M/(k_theta+c_theta*p)`.
- Время: `u=U*exp(i*varpi*t)`, (12); `omega_F²=varpi²*m*L⁴/EI`.
  При `t_ref=L²*sqrt(m/EI)` проектный `z=i*omega_F=-q_F+i*p_F`.
  Damping ratio — (34b), `q_F/sqrt(p_F²+q_F²)`.
- Benchmark: только §6.1, Fig. 2, Table 1, p.59; Fig.3,p.60 — качественный
  контроль непрерывности прогиба и скачков поворота/силы. Три совмещённых TS/RJ
  при xi=.25,.5,.75: kappa_u=100,gamma_u=.1,kappa_theta=10,gamma_theta=.1.
  Пять опубликованных eigenvalues и damping ratios; mode 4 не возбуждает TS/RJ.
- Проверяет EB-комплексную задачу, закон RJ/TS, перевод временного соглашения,
  комплексные корни и неактивный локальный демпфер. Не проверяет Timoshenko/RLB,
  продольное движение, слоистую редукцию или двухплечевой узел под углом beta.
- Метаданные: `VERIFIED_LOCAL_PDF`, 2026-09-15, первая страница и пагинация;
  формулы/таблица дополнительно просмотрены в изображениях PDF.
  [Транскрипция и benchmark](../laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md#failla-source).
- Фактическое воспроизведение: Table1 — LITERATURE_MISMATCH по части
  печатных разрядов; неактивная mode4 подтверждена. См.
  [численное сравнение](../laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md#failla-table-1).
- Дополнение 2026-09-26: отдельный анализ печати подтверждает несовместимость
  mode4 с обычным округлением `(4*pi)^2` и отображённого ratio mode3 с
  напечатанными p,q; конкретная причина не установлена.
  [Второй проход](../laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md#second-pass-source-precision-audit-and-hong-table-3).

## `hong_kim_1999_damped_timoshenko_joints`

- PDF: `docs/literature/pdf/hong1999.pdf`, 20 страниц (журнальные 787–806),
  SHA256 `e8622f7d407ffcd353fea8578cd40f840b530c33be272a9d1a81af47ac015088`.
- Тип: журнальная статья S.-W. Hong, J.-W. Kim, *Modal analysis of multi-span
  Timoshenko beams connected or supported by resilient joints with damping*,
  Journal of Sound and Vibration 227(4) (1999), 787–806.
- Роль: внешний reference для Timoshenko shear deformation, rotary inertia,
  Laplace-domain state и exact dynamic matrix с damped resilient joints.
- Обозначения и знаки: (1)–(3), `Psi=[u*,phi*,F*,M*]`,
  `u'=phi-F/(kAG)`, `phi'=M/(EI_d)`, `F'=-rho*A*s²*u`,
  `M'=F+rho*I_d*s²*phi`. Относительно `[w,psi,Q,M]` проекта:
  `Psi=diag(1,-1,-1,-1)*y_transverse`. Простая замена F=Q,phi=psi неверна.
- Время: Laplace variable `s`, нулевые начальные условия; (2),(17),(22).
  Полюсы Table 3 `lambda_k=sigma_k+j*omega_k` сравниваются с проектным p
  напрямую, без умножения на i. Частоты Table 2 и eigenvalues Table 3 размерные.
- Benchmark: Numerical example 1 (§4.1), Fig.2 и Tables1–3,p.796,797,799.
  Две **supporting translational** KV-опоры на концах, `k_t=2e6 N/m`,
  `c_t=20 Ns/m`; повороты свободны, концевые моменты нулевые.
  Table2: по пять положительных частот hinged–hinged/free–free без опорных
  элементов; нулевые rigid-body modes исключены. Table3: пять complex roots,
  target — Proposed method, не FEM. Пример не использует connecting rotational joint.
- Ограничения: не проверяет продольное движение, ламинатную редукцию,
  внутренний rotational KV-узел и angled two-arm геометрию. E,G,nu из Table1
  сохраняются как отдельные benchmark inputs; G не пересчитывается из E,nu.
- Метаданные: `VERIFIED_LOCAL_PDF_AND_PUBLISHER`, 2026-09-15. PDF содержит
  Article No. jsvi.1999.2385; полный DOI `10.1006/jsvi.1999.2385` подтверждён
  [издательской страницей](https://www.sciencedirect.com/science/article/pii/S0022460X99923854).
  [Транскрипция и signed mapping](../laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md#hong-source).
- Фактическое воспроизведение: Table2 — LITERATURE_MISMATCH при согласии
  матрицы с формулами сноски; Table3 — NOT_RUN_B1_GATE. Демпфированная
  численная проверка этим запуском не выполнена; см.
  [результат](../laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md#hong-table-2).
- Дополнение 2026-09-26 по D13: printed FAIL Table2 сохранён, отдельный
  equation-level PASS разрешил Table3. Все пять новых complex roots проходят
  прежний printed gate и solver checks. Исторический NOT_RUN не перезаписан;
  [данные второго прохода](../laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md#second-pass-source-precision-audit-and-hong-table-3).

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
