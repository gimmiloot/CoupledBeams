# β15: пространственная проверка внеплоскостного роста при p14

**P14_OUT_GROWTH_REPRODUCED**, 2026-10-11.
Та же плоская периодическая ветвь типа II, Ψ=0, при A=0.005 даёт
положительный OUT-показатель при p14: σ=0.003411315508.
Его наблюдаемое отличие от p12 составляет 0.00277%.
Сохраняются отрицательный вещественный Floquet multiplier и физическая
форма растущего направления. Это сильное согласование двух разрешений,
но не доказательство полной p-сходимости или физическая 3D validation.

OUT2-отклонение с нулевыми скоростями имеет разрешимую ненулевую проекцию
на растущий режим. Новая временная траектория такого отклонения не вычислялась.
Дополнительный PL spectral control — **PARTIAL**: ненейтрального роста в
рассмотренном спектре нет, но строгий time-shift neutral-residual gate не пройден.

## Источники и неизменность модели

Использованы [принятое действие](../theory/weakly_nonlinear_spatial_rod.md),
[ветвь 1:3](nlsp_beta15_planar_1to3_internal_resonance.md),
[p12 HB/Hill](nlsp_beta15_periodic_orbit_spatial_stability.md) и
[p12 dynamic validation](nlsp_beta15_dynamic_out_growth_validation.md).
Исторические bundles:

- `results/nlsp_beta15_planar_1to3_internal_resonance/b59944c2ef43d0b7/`;
- `results/nlsp_beta15_periodic_orbit_spatial_stability/6b93bc782316296d/`;
- `results/nlsp_beta15_dynamic_out_growth_validation/49dc29161a5dda78/`.

Их manifests, PDF и остальные файлы сохранены по исходным SHA256.
Новый bundle: `results/nlsp_beta15_out_growth_p14_convergence/2f77d2952671122b/`.
Все p12 частоты, коэффициенты, спектры и динамические результаты переиспользованы.
Новые eigenproblems относятся только к p14; повторного p12 HB/Hill/ODE нет.

L₁=0.9,L₂=1.1,β=15°,μ=0.1; оба b=0.05,h=0.10;
E=ρ=1,ν=0.3,κ=5/6. Сохранены семь полей u,w,v,Φ,ψ,θ,c,
S≤4, cubic mass-form residuals, variable mass, nonlinear inertia,
constitutive coefficients, Galerkin family, заделки, local/global mappings,
идеальный жёсткий узел и c₁=c₂. Physical kernels не изменены.
Модельные длина/время нормированы; углы измеряются в радианах, Lref=1.
Forcing, damping и несовершенства ненапряжённой формы отсутствуют.

При p14 — 189 coordinates: 108 PL и 81 OUT. Полные линейные bases
нужны для full-space HB и canonical projection; широкий поиск резонансов
или изменение физической геометрии не выполнялись.

## Межпространственный перенос и линейная идентичность

Перенос построен массовой Galerkin-проекцией физических полей на общей
материальной Gauss quadrature. Используются те же семь resting-mass weights,
endpoint lifts, whitening и узловые transformations. Коэффициенты разных
дискретизаций не отождествляются по индексам.

Для transfer L14,12 получены relative errors: сохранение M₀ — 3.33e−14,
K₀ — 1.01e−13; physical profiles и их производные порядков0/1/2 — 5.56e−16.
Сопоставление mass-weighted MAC выполнено после физического переноса.
Знаки активных p14 eigenvectors согласованы с сохранёнными p12 формами.

| Форма в этой фиксированной геометрии | p14 frequency | p12/p14 mass MAC |
|---|---:|---:|
| PL2 | 0.434820604535 | 0.999999999936 |
| PL4 | 1.306328966049 | 0.999999999272 |
| OUT2 | 0.220127282806 | ≈1 |

Все три соответствия однозначны; MAC margins близки к единице.
Mass orthogonality residuals PL/OUT≤2.89e−15; individual low-mode
eigen-residuals≤1.35e−11. Sorted positions используются только как местные
обозначения уже сопоставленных форм, без утверждения о продолжении по β.
При p14 ωPL4−3ωPL2=0.001867152445.
Для уточнённой орбиты Ω14−2ωOUT2=−0.005139107309.
Обе прежние частотные особенности сохраняются; эти расстройки сами по себе
не устанавливают неустойчивость.

## Периодическая ветвь типа II

Сохранённые physical HB11 profiles p12 перенесены и разложены в полном
mass-orthonormal PL basis p14. Прежний FullPlanarHB решает исходные
mass-form equations для всех Fourier coefficients и неизвестной Ω.
Косинусная постановка сохраняет time reversal и нулевые начальные скорости.
Четвёртая форма не задаётся дополнительно; A — точный spatial max|w(s,0)|.

HB11 прошёл gate после одного Newton update. Единственный контроль HB13
сохранил то же движение: начальная enlarged-system residual уже ниже frozen
projected target1e−10, поэтому дополнительных Newton updates не потребовалось.
Нулевую наблюдаемую разность Ω11/Ω13 нельзя считать нулевой точной ошибкой.
Полная off-grid residual проверена независимо и остаётся контролируемой.

| Проверка p14 | Результат |
|---|---:|
| Ω11=Ω13 в принятой solver accuracy | 0.435115458303789 |
| Actual initial max abs(w) | 0.005000000000016 |
| Independent full relative residual, 61 shifted phases | 9.010882e−9 |
| Absolute full mass-form residual | 5.662472e−13 |
| Полная energy variation | 4.289650e−11 |
| Минимальная sampled mass eigenvalue | 0.985693057 |
| Joint compatibility | 1.110223e−16 |
| Scaled Jacobian condition, HB11 / HB13 | 5732.4 / 6180.4 |
| Изменение активных координат от перенесённого seed | 4.105454e−5 |
| Physical waveform mass-norm difference p12/p14 | 8.579825e−6 |

U₁,₂ и U₃,₄ положительны, Ψ=0; переключения на Ψ=π нет.
Сохранены средняя составляющая и все четыре PL поля. Relative position-E₀
norms средней/второй/третьей гармоник к основной при p14 равны
0.02545690 / 0.03990328 / 0.09528999; при p12 —
0.02545690 / 0.03990383 / 0.09528617. Это нормы коэффициентов, не energy fractions.
Последние retained HB11 harmonics имеют relative norm1.12e−8.

В начальной полной квадратичной reference energy p14 доли PL2/PL4/прочих
форм равны98.73660% / 1.03855% / 0.22485%.
Это другое определение, чем доля PL4 в **ведущих активных** U₁,₂/U₃,₄,
приведённая ниже для прямого сравнения с историческим отчётом.
Ни одно такое разбиение не утверждает отдельного сохранения энергии каждой формы.

## OUT Hill analysis и numerical gates

Оператор получен прежними производными quartic action и analytic Jacobian.
Канонический переход Π=Mη̇+Cη, Z=S(t)z включает ṠS⁻¹.
Учитываются среднее состояние, третья и остальные гармоники орбиты;
stiffness-only surrogate или прежний single-harmonic ansatz не используются.
PL/OUT cross blocks точно нулевые. Canonical Hamiltonian, duality,
action/Jacobian и Fourier-conjugacy checks проходят frozen operator gate2e−9.
На12 независимых фазах Fourier reconstruction residual=1.09e−14.

Подстановка z=exp(λt)Σzₙexp(inΩt) даёт Hₙₘ=Fₙ₋ₘ−inΩδₙₘI.
Использованы прежние alias/edge conventions и fundamental strip.

| OUT modes | Hill N | σ |
|---:|---:|---:|
| 24 | 7 | 0.003366258955 |
| 48 | 7 | 0.003366908679 |
| 72, дополнительный контроль | 7 | 0.003411210156 |
| Все81 | 7 | 0.003411315507701 |
| Все81 | 9 | 0.003411315507706 |

Первоначальная 48/81 разность4.440683e−5 превышала frozen gate5e−6.
Этот непрошедший результат и его расходы сохранены отдельно. До нового solve
объявлена одна дополнительная 72-mode проверка; HB и прежние спектры не повторялись.
Порог не менялся. Разность72/81=1.053521e−7, N7/N9=4.379743e−15.
Observed uncertainty indicator=max(last modal difference, temporal difference,
orbit comparison, floor1e−8)=1.053521e−7. Это не строгая error bound.
Сравнение HB11/HB13 в24 модах не обнаружило изменения σ при данном solver target.

Last full Hill local eigen-residual=1.80e−14; independent convolution
согласуется с dense action до6.01e−16. Quartet defect=1.10e−14,
edge weight3.67e−7 проходит заранее принятый gate1e−6.
Physical state MAC для N7/N9=0.999999999999996;
для72/all81=0.999999485469. Прямое p12/p14 state-E₀ MAC=0.999999999759.
Сопоставлены физические перемещения **и скорости** на той же основной фазе.
Согласован только знак eigenvector; amplitude/phase fitting основной орбиты нет.
Все81 modes входят в оператор, но квалифицирован только спектр interior replicas
при N7/9. Это не полный спектральный аудит всех высокочастотных Floquet aliases.

## Основное сравнение

| Величина | p12, saved | p14 | p14−p12 |
|---|---:|---:|---:|
| Periodic Ω | 0.435115586723 | 0.435115458304 | −1.284195e−7 |
| PL4 leading active reference energy share, % | 0.899484850 | 0.899556351 | +0.000071501 процентного пункта |
| Full PL relative residual | 8.855097e−9 | 9.010882e−9 | +1.557854e−10 |
| OUT σ, 1/model time | 0.003411410067 | 0.003411315508 | −9.455889e−8 |
| σT=log abs(ρ) | 0.049261672671 | 0.049260321754 | −1.350917e−6 |
| Negative-real multiplier ρ | −1.050495200663 | −1.050493781532 | +1.419131e−6 |
| e-folding time / T | 20.299757312 | 20.300314013 | +0.000556701 |
| OUT phase-zero bending/torsion, % | 97.299643 / 2.700357 | 97.299683 / 2.700317 | практически сохранены |

Observed relative differences σ и σT равны0.00277184% и0.00274233%.
Они меньше last modal uncertainty indicator; малую разность не следует
трактовать как точно определённую ошибку p12. Знак роста и физическая идентичность
согласуются значительно сильнее заранее объявленного ориентира5%.
В Fourier energy-coordinate composition доминирует OUT2:93.38713%,
с OUT4≈5.77380% и OUT5≈0.78509%. Это mixed bending/torsion state, не pure bending.

При p14 multiplier имеет imaginary residual1.90e−10, а reciprocal partner
равен−0.951933288497 с imaginary residual1.72e−10. Их знаки и reciprocal/conjugate
relations проходят проверки. Raw λ=0.003411315508+1.958019562355i — Fourier replica;
положительный σ не зависит от ветви complex logarithm. Форма меняет знак через T.

## Начальное отклонение с нулевыми OUT-скоростями

До проекции выбрана сопоставленная линейная OUT2, с η̇(0)=0 и физическим
масштабом5e−5. Начальное PL состояние — прежнее HB13 Type-II turning point.
Maxima (v,Φ,ψ)=(1.328526e−5,1.429647e−5,5e−5); velocities exactly0,
joint compatibility сохранена. Перебора начальных форм нет.

Сохранена reciprocal stable Floquet eigenvector из **того же** eigensolve.
Для canonical real vectors rᵤ,rₛ используется symplectic dual
lᵤᵀ=rₛᵀJ/(rₛᵀJrᵤ). Это следует из reciprocal relation симплектической
эволюции и не является евклидовым скалярным произведением правых eigenvectors.
Проверены lᵤᵀrᵤ=1, lᵤᵀrₛ≈0, phase invariance и physical canonical mapping.

Initial OUT E₀-norm=1.926794e−7; коэффициент единичного E₀-growing direction
αᵤ=−1.877948e−7, то есть αᵤ/norm(initial)=−0.974649105.
Этот коэффициент **не является долей энергии** неортогональных Floquet modes.
Dual E₀-norm=1.143366; N7/N9 relative coefficient change=1.63e−11.
Unstable/stable pointwise canonical residuals≤5.10e−4 проходят frozen0.002 gate;
Wronskian relative variation=4.41e−7. Консервативный observed projection
uncertainty indicator с учётом этих остатков=1.12e−10, существенно меньше abs(αᵤ).
Это диагностический indicator, не строгая adjoint error bound.

Следовательно, механизм доступен для данного displacement-only initial state.
Фактический рост полной суммы stable/neutral/unstable компонентов на коротком
интервале и большой физический выход из плоскости здесь не проверялись.
Short tangent control не запускался: корректная dual projection получена без ODE.

## Дополнительный PL control

All108 PL modes, N5/N7 дают максимальные ненейтральные real parts
2.69e−15 / 6.23e−15. Нового разрешённого ненейтрального PL роста не найдено.
Однако нейтральный кластер содержит real splitting≈±2.62e−8.
Эти значения находятся в объявленной neutral band5e−6, но выше floor1e−8.
Full time-shift convolution с N13 имеет residual7.59e−8 и не проходит
строгий operator gate2e−9. Поэтому контроль отмечен PARTIAL;
neutral splitting не объявляется ни доказанным ростом, ни полностью разрешённым нулём.
Порог не повышался. Новый HB/Fourier refinement ради вторичного PASS не выполнялся.
Это конкретное оставшееся препятствие для полной PL qualification,
а не отказ основного OUT spatial comparison.

## Расходы, tests и воспроизводимость

До production сохранён прогноз≈500s при бюджете1800s, включая основные
HB/OUT controls и вторичный PL control. Матрицы OUT3078 и PL3240 имеют
144.6/160.2MiB соответственно, без eigenvector/workspace allocations.

Два p14 linear eigenproblems, два HB solves, восемь Hill eigenproblems.
HB11/HB13 заняли2.23/1.92s. OUT spectral solves, включая HB11 comparator
и triggered72-mode control, около75s; PL N5/N7 —19.54/48.10s.
Все production/scientific diagnostics вместе156.83s.
Regression runs17.53s+4.47s+4.10s консервативно включены в budget:
всего **182.93s/1800s**.
Новых ODE/FEM=0; p16, новые геометрии/амплитуды/ветви не выполнялись.

61 scoped tests PASS:20 новых p14 и41 прежних HB/Hill/dynamic regressions.
После дополнительной saved-mode MAC diagnostics новые20 tests снова прошли.
Финальный выборочный cache/history run:4 tests PASS.
Проверены transfer на независимых материальных точках, complete HB residual
на дополнительных фазах, coefficients, identity, mass/joint/clamps,
Ṡ/canonical/reflection identities, Hill signs/alias/residuals, multiplier pairing,
resolution controls, symplectic dual и physical zero-velocity initialization.
Cache compute/report/plot проверен на temporary copies с scientific calls запрещёнными.
Исторические source hashes сохранены; `git diff --check` пройден.

Existing CLI:
`scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --beta15-p14-convergence`.
Completed `--compute` читает cache; `--report-only <copy>` и `--plot-only <copy>`
не выполняют научных расчётов. Нового runnable script или physical solver нет.
Расширен прежний scoped helper; reciprocal vector сохраняется опционально,
с прежними default outputs для старых режимов. Runtime source snapshots
и supplemental policy объясняют технические изменения controller; physics frozen.

Bundle содержит p12 references, transfer и MAC matrices, p14 eigenpairs,
HB coefficients/diagnostics, Fourier operators, eigenvalues/vectors,
independent Hill action, reciprocal pair, zero-velocity initial state,
projection controls, PL neutral evidence, costs, provenance и manifest.

## Рисунки и остановка

1. [PL ветвь и модальный состав](../../results/nlsp_beta15_out_growth_p14_convergence/2f77d2952671122b/p12_p14_planar_branch.png):
   physical midpoint w, проекция обеих орбит в один согласованный PL basis,
   leading active PL4 reference share. Основные фазы не выравнивались.
2. [OUT growth и физическая форма](../../results/nlsp_beta15_out_growth_p14_convergence/2f77d2952671122b/p12_p14_OUT_growth_mode.png):
   две пространственные точки, observed uncertainty indicators, profiles v(s/L)
   при одинаковом комбинированном физическом масштабе1%A. Bars не являются CI/bounds.
   В bundle также PDF; никакой карты геометрической неустойчивости нет.

Gates A/B/C/D/E PASS; F PARTIAL. Основной статус P14_OUT_GROWTH_REPRODUCED.
Предыдущий p12 dynamic status не пересмотрен; p14 подтверждение спектральное,
без новой nonlinear dynamic validation. Два p — существенное свидетельство
пространственной устойчивости результата, не теорема p-convergence.
Для строгой количественной p-qualification дополнительное разрешение могло бы
быть отдельным следующим этапом, но немедленный p16 для выяснения знака не требуется.
PL neutral qualification имеет отдельный временной/орбитальный остаток;
его нельзя автоматически приписывать недостатку p или устранять запуском p16.

DEFERRED_TWO_PARAMETER_ACTION_EXPANSION, physical ideal joint/c₁=c₂,
independent 3D validation, further p-refinement, Type-I dynamics,
large OUT amplitude и saturation остаются отложенными. **STOP после p14.**
