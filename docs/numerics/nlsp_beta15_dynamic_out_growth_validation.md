# β15: динамическое подтверждение внеплоскостного роста

**DYNAMIC_OUT_GROWTH_CONFIRMED_P12**, 2026-10-11.
Прямое интегрирование касательной системы и исходной нелинейной семиполевой
системы воспроизводит растущую Floquet-компоненту сохранённого движения типа II.
За пять периодов её усиление составляет соответственно 1.279293578 и 1.279292325.
Результат относится к одной геометрии, одной амплитуде и одной начальной форме
возмущения в p12-модели. Большого физического выхода из плоскости за 5T нет.

## Постановка и источники

β=15°, μ=0.1, L₁=0.9, L₂=1.1; оба сечения b=0.05, h=0.10;
E=ρ=1, ν=0.3, κ=5/6. Внешние заделки, идеальный жёсткий узел,
c₁=c₂, локальные оси, Galerkin space и все семь полей сохранены.
Принятое [quartic action](../theory/weakly_nonlinear_spatial_rod.md),
variable mass, nonlinear inertia, cubic residuals и analytic Jacobian не изменены.
Физические kernels совпадают по SHA256 с исходными данными предыдущего этапа.
Внешние силы после освобождения равны нулю; damping и несовершенства
ненапряжённой конфигурации не вводились.

Основной источник — [HB/Hill-аудит](nlsp_beta15_periodic_orbit_spatial_stability.md):
`results/nlsp_beta15_periodic_orbit_spatial_stability/6b93bc782316296d/`,
manifest SHA256 `0ec30aad088ed9c3d85464bdcacfe024ec88ff9f6b8f901a43d34f7637721d3b`.
Предшествующая [ветвь 1:3](nlsp_beta15_planar_1to3_internal_resonance.md):
`results/nlsp_beta15_planar_1to3_internal_resonance/b59944c2ef43d0b7/`,
manifest SHA256 `dd1f5d5300783fc67f0386081715a0f2a64865d1dabe8e2f7ce0362f06157df6`.
Все файлы обоих исторических bundles, включая PDF и manifests, сохранены побайтово.

Новый bundle:
`results/nlsp_beta15_dynamic_out_growth_validation/49dc29161a5dda78/`.
Полная дискретизация содержит 161 coordinate: 92 PL и 69 OUT.
HB11 coefficients задают основное начальное состояние типа II, Ψ=0,
с фактическим max|w(s,0)|=0.004999999999999998 и нулевыми PL-скоростями.
Это специально подготовленная периодическая ветвь, а не начальное возбуждение
одной линейной формы типа I. Другие амплитуды и ветви не вычислялись.

Сохранённая частота Ω=0.4351155867232466 даёт T=14.440267135674869.
Период использован с полной доступной точностью и не подгонялся к ODE.
Длины и время приведены в единицах нормированной модели; углы — в радианах.
Lref=1. Нормировка E=ρ=1 не задаёт физический материал и единицу секунды.

## Проверка плоской орбиты

Интегрировалась существующая инвариантная PL-система с полной нелинейной
массой и инерцией. Radau использует прежний analytic Jacobian. На первом
периоде выполнены исходный и более строгий расчёты; уточнённое конечное
физическое состояние продолжено сначала до 3T, затем до 5T.
Префиксы не пересчитывались. Restart начинается с точного принятого состояния,
но не воспроизводит внутренний adaptive solver path побайтово.

Фиксированная reference-energy norm имеет квадрат
qᵀK₀q+q̇ᵀM₀q̇, отдельно для PL и OUT. Нормы не менялись во времени.

| Проверка | Результат |
|---|---:|
| Относительная E₀-невязка возвращения за T | 6.417501e−10 |
| Абсолютная максимальная coefficient-невязка возвращения | 1.067426e−14 |
| Возвращение: position K₀-norm / velocity M₀-norm | 8.396397e−14 / 2.932485e−14 |
| Максимальная относительная E₀-разность ODE/HB11 за T | 7.430075e−10 |
| Та же разность до 5T | 7.833104e−10 |
| Относительный drift полной энергии за T | 4.084257e−13 |
| Разность двух temporal tolerance levels, E₀ относительно PL | 8.518757e−13 |
| Разность при половинном max_step, E₀ относительно PL | 8.382123e−13 |
| Минимальное собственное значение PL mass, sampled | 0.985729697 |
| Максимальный остаток joint compatibility | 2.220446e−16 |
| OUT RHS при точно нулевых OUT state/velocities | 0 |

PL/HB comparison использует все поля и скорости. Фазы ODE output не совпадают
с узлами прежнего HB solve. Нормировка PL-разностей — E₀-норма исходного
состояния 1.385860187e−4. Safety и joint diagnostics проверены на сохранённых
состояниях; точная совместность обеспечивается неизменной assembled space.
Новая длинная seven-field zero-OUT trajectory не потребовалась.

## Физическая Floquet-мода и касательная проверка

Взяты сохранённые all-69-OUT/N9 eigenvector и R из предыдущего bundle.
Новых HB/Hill solves нет. Использованы оригинальные Fourier indices −9…9
и raw λ=0.003411410066594682−1.958020140242108i. Переход к физическим
перемещениям и скоростям выполнен через S(t)⁻¹R, где
S=[[I,0],[C(t),M(t)]]. Сопряжённый импульс не отождествляется со скоростью.
Прежний канонический оператор включает Ṡ; direct physical Jacobian используется
без повторного добавления этого преобразования в ODE.

Multiplier за T равен −1.050495200663−1.89645e−10i. Поэтому вещественная
форма меняет знак через каждый период. Это ожидаемая Floquet-структура.
Imaginary alias нельзя заменять без соответствующего сдвига Fourier indices.
Real/imag Gram имеет eigenvalues 1.075041e−9 и 1.000854742; отношение
1.074123e−9 соответствует одному разрешённому вещественному направлению.
Выбрана его главная Gram-комбинация до ODE, с фиксированным знаком формы.
Малая вторая компонента не объявляется независимой неустойчивой модой.
Kinematic reconstruction error=2.226449e−6; saved eigenpair residual=6.562052e−15.

Начальное направление d нормировано в OUT E₀-norm к единице. Масштаб
δE=4.338799640750e−7 даёт
max(|v|,Lref|Φ|,Lref|ψ|)=5e−5=0.01A.
Начальные максимумы (v,Φ,ψ) равны (1.502877e−5,3.298738e−5,5e−5),
а M₀-норма OUT-скоростей — 3.674726e−7. Эти скорости сохранены.
Фаза и направление не подбирались по результату полного нелинейного расчёта.
Физическая форма смешанная; в исходном спектральном аудите преобладает OUT2,
с меньшим участием OUT4/OUT5. Sorted labels относятся только к этой геометрии.

Касательная система интегрировалась около **реальной nonlinear PL ODE trajectory**
по её accepted cubic dense output, а не около повторяемого HB массива.
На стробоскопических временах вычислен cₖ=dᵀE₀z(kT). Сравниваются знак,
полный физический вектор и |cₖ/c₀|; σdyn — наклон log|cₖ/c₀| по kT.
При продолжении не нормировались заново ни d, ни tangent state на 3T.

| Величина | Saved Hill | Tangent ODE | Full nonlinear ODE |
|---|---:|---:|---:|
| Спектральный/линейный growth exponent | 0.003411410067 | 0.003411405495, fit 0…5T | NOT_ESTIMABLE как nonlinear Floquet exponent |
| Finite-window log-projection rate, 0…5T | — | — | 0.003411392055 |
| Усиление компоненты за 3T | 1.159263649 | 1.159263270 | 1.159262825 |
| Усиление компоненты за 5T | 1.279294000 | 1.279293578 | 1.279292325 |
| Начальный физический OUT-масштаб | 5e−5 | 5e−5 при сравнении | 5e−5 |
| Максимальный физический abs(v) до 5T | 3.673504e−5 | 3.673501e−5 | 3.673498e−5 |
| Максимальный abs(Φ) до 5T, rad | 1.118983e−4 | 1.118976e−4 | 1.118974e−4 |
| Максимальный abs(ψ) до 5T, rad | 1.527965e−4 | 1.527960e−4 | 1.527959e−4 |
| Квалификация | Исторический spectral candidate | Dynamic tangent confirmation | Finite-perturbation confirmation |

Tangent fit 0…3T даёт 0.003411403496, разность с Hill 6.570199e−9;
разность fit 0…5T — 4.571716e−9. Максимальная стробоскопическая ошибка
полного физического вектора относительно reconstructed Floquet prediction
до 5T равна 2.389892e−5. Знаки cₖ чередуются, как предсказывает multiplier.
Нормы отдельных перемещений могут осциллировать: их монотонность не требовалась.

Допуски PL rtol=1e−10, control rtol=1e−11; component atol задаются
сохранённым field-scaled правилом. Tangent DOP853 rtol/relative atol=1e−11,
tolerance control — в десять раз строже. Max_step≈0.007209393.
Поскольку оба tolerance levels упёрлись в одинаковый max_step, отдельно выполнен
контроль одного периода с max_step/2 при тех же строгих допусках.
Разность unit-E₀ tangent response=1.219205e−12; разность одно-периодного
lograte=6.735064e−16. Tolerance-control разность lograte=5.121771e−16.
Это наблюдаемые различия, а не строгая оценка global error.

## Полная семиполевая динамика

После прохождения предыдущих gates запущен прежний physical seven-field Radau
solver с тем же PL initial state и добавленным физическим Floquet-возмущением.
Использованы rtol/relative atol=1e−12 и max_step≈0.0036046965.
Проверка до 3T продолжена с сохранённого native endpoint до 5T.

Относительная E₀-ошибка OUT nonlinear/tangent составляет 1.236632e−6
до 3T и 2.415089e−6 до 5T. Абсолютные максимальные ошибки соответственно
6.220025e−13 и 1.340519e−12. Сравнение включает OUT-скорости.
PL backreaction до 5T имеет E₀-norm 7.008741e−10, или 5.057322e−6
относительно исходной PL-нормы. Полная энергия сохраняется с относительным
drift 3.019807e−14. Минимальная sampled mass eigenvalue=0.985729684,
joint compatibility=2.220446e−16; assembled generalized nodal balance≤1.905824e−21.
RHS на PL-компонентах этих состояний при обнулённых OUT равен нулю в OUT block.

Максимальный физический выход из плоскости равен 0.007346996A, то есть
0.7346996% плоской амплитуды. Максимумы вращений указаны в таблице.
Рост Floquet-компоненты примерно на 28% не означает большой прогиб.
Отношение максимального прогиба к его начальному значению также не является
чистым экспоненциальным gain: начальные скорости ненулевые, carrier осциллирует.
Временные максимумы sampled на сетке T/256; пространственные — точные экстремумы
представленных полиномов с учётом корней производных. Сглаживания и fitting нет.

Один короткий nonlinear control при OUT-масштабе 2.5e−5 выполнен до T,
при той же PL-базе, направлении, фазе и solver policy. Нормированные ошибки
OUT response равны 3.987678e−7 и 9.969167e−8; наблюдаемый порядок 2.000004.
PL backreaction уменьшается с 3.201310e−10 до 8.003277e−11, порядок 2.000000.
Однако квалификация OUT order — **NUMERICAL_FLOOR_LIMITED**: меньшая ошибка
чуть ниже frozen 10×floor=1e−7. Порог не снижался ради PASS.
Норма положительного component-atol vector относительно δE=2.446951e−8
служит лишь proxy масштаба tolerances; это не достигнутая ошибка и не её bound.
Главное усиление 0.279 существенно выше этих масштабов, а tangent growth
подтверждён независимо. Ограниченная проверка порядка не отменяет основной результат.

## Gates, вычислительные расходы и воспроизводимость

До production сохранены config, pre_ODE_policy и исходные source snapshots.
Frozen gates: PL orbit/closure relative E₀≤5e−6, temporal state difference≤5e−7,
energy drift≤1e−8, joint≤2e−10, Floquet reconstruction/vector≤0.002,
Hill/dynamic exponent difference≤5e−5, temporal exponent difference≤5e−6,
full/tangent relative error≤1e−5, floor=1e−8. Gates A/B/C/D/E — PASS.
Контроль половинного max_step добавлен из-за одинаковых фактических шагов
двух tolerance levels; прежние пороги не менялись.

| Новые ODE calls | Число | Runtime, s |
|---|---:|---:|
| PL primary/refined/half-step, каждый T | 3 | 98.702490 |
| PL continuation T→3T и 3T→5T | 2 | 97.732789 |
| Tangent primary 3T и continuation 3T→5T | 2 | 583.100634 |
| Tangent tolerance/half-step controls, каждый T | 2 | 348.331643 |
| Full nonlinear 1% OUT до 3T и continuation до 5T | 2 | 333.003092 |
| Full nonlinear 0.5% OUT, T | 1 | 66.766990 |
| **Все ODE** | **12** | **1527.637638** |

Итого восемь nonlinear calls (пять PL и три full spatial segments), четыре
tangent calls; failed/retried integrations=0. HB/Hill solves=0, FEM=0.
Production, preflight, pointwise diagnostics и scientific postprocessing:
1541.182446 s. Два regression runs, включая первоначальный отказ одного
побайтового assertion, консервативно добавляют 34.61 s.
Суммарно около **1575.79 s из 3600 s**. Сохранены расходы каждого случая,
native time grids, nfev/njev/nlu, mass-factorization counters и endpoints.

Перед full spatial solve профиль десяти различных сохранённых состояний дал
RHS≈0.001167 s, Jacobian≈0.006234 s, safety≈0.000403 s на evaluation.
Консервативный прогноз по историческому C/p16 nfev density — 2413.24 s
при остатке 2754.09 s и reserve=60 s. Actual full3T занял 199.02 s:
12019 accepted steps,84135 nfev. Исторический C control имел гораздо больше
шагов; геометрия, p и подготовка начального состояния различаются.
Поэтому причина разности стоимости не приписывается одному физическому фактору.
Solver не оптимизировался. Прогноз продолжений до 5T=529.81 s; resource gate PASS.

Минимальные технические правки writers сохраняют native accepted times и
точные final state/matrix также при t_start=0. Они не меняют RHS, Jacobian,
integrator или adaptive path. Сохранены первоначальные, pre-ODE и финальные
версии scoped controller: whole-wrapper hashes между подготовительными
правками различаются, frozen physical-kernel hashes остаются неизменными.

80 scoped regression tests PASS. Проверены source/initial identity, physical
scaling и velocities, canonical mapping, отрицательный multiplier и real rank,
fixed norm, reflection/joint/clamps, temporal controls, native restart,
full/tangent comparison, бюджет и сохранение всех historical artifact hashes.
Compute/report/plot replay проверен на временных копиях с scientific functions
запрещёнными. Первоначальный exact-array parity assertion обнаружил BLAS
roundoff 2.71e−20; проверка теперь требует точной чётности Fourier basis и
выведенного floating-point summation bound. Scientific acceptance gates не изменены.
`git diff --check` пройден; staging и исходный HEAD сохранены.

CLI: `scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --beta15-dynamic-out-growth`.
После completed compute читает cache.
Для повторного оформления копии bundle разрешены `--report-only <copy>`
и `--plot-only <copy>`. Новый runnable script не создан.
В manifest включены исходные состояния, HB/Floquet source identities, фиксированные
norm matrices, dense PL archive, все actual solver grids, checkpoints, tangent/full
responses, physical observations, scalar projections, cost forecasts и verification.

## Рисунки

1. [PL ODE/HB11 и полная E₀-разность](../../results/nlsp_beta15_dynamic_out_growth_validation/49dc29161a5dda78/planar_orbit_validation.png):
   прогибы середин обоих плеч в model-length units; справа относительная state norm.
2. [Стробоскопическое tangent/Floquet усиление](../../results/nlsp_beta15_dynamic_out_growth_validation/49dc29161a5dda78/tangent_Floquet_validation.png):
   знак формы меняется каждый T, изображён модуль projected component.
3. [Полная пространственная динамика и tangent prediction](../../results/nlsp_beta15_dynamic_out_growth_validation/49dc29161a5dda78/nonlinear_spatial_validation.png):
   слева max|v|/A, справа fixed E₀-norm относительно начальной нормы; до 5T.
   В bundle сохранены также PDF. Кривые не сглажены и не выровнены по фазе.

## Научный вывод и остановка

Достигнуты level 2 — независимое временное подтверждение Hill prediction —
и level 3 — finite-perturbation подтверждение полной семиполевой динамикой
в той же p12-модели. Это даёт первое адресное динамическое подтверждение
предсказанного начального внеплоскостного роста данного движения типа II.
Нет основания объяснять результат только numerical floor или PL orbit drift.
Условное e-folding time Hill/tangent≈293.134 model-time, или 20.300T;
до столь длинного горизонта интегрирование не проводилось.

Finite-window nonlinear log-projection rate не является показателем Флоке
нелинейной траектории или доказательством её асимптотического поведения.
Точная периодическая орбита не доказана теоремой существования.
Spatial p-convergence эффекта, физическая корректность идеального узла,
реальное соединение и независимая 3D/экспериментальная validation не установлены.
Level 4 не достигнут. Тип I, другие направления возмущения, критическая амплитуда,
большой выход из плоскости, насыщение и длительная динамика не исследовались.

Минимальный следующий численный шаг для более надёжного пространственного
вывода — отдельно разрешённая проверка той же Type-II ветви и OUT-эффекта
при p-refinement, с физическим сопоставлением форм и той же нормировкой.
Она не заменит независимую 3D/экспериментальную проверку.
Сохранены отложенные DEFERRED_TWO_PARAMETER_ACTION_EXPANSION,
physical ideal joint/c₁=c₂, independent 3D validation, Type-I dynamics,
spatial p-refinement и large-amplitude/saturation. **STOP:** новых расчётов нет.
