# β15: две активные плоские формы около внутренней близости 1:3

**PLANAR_1TO3_RESONANT_REDUCTION_VALIDATED**, 2026-10-11.
Получены контролируемые в пределах зафиксированных p12 numerical gates
приближённые периодические движения **типа II**. Четвёртая форма входит
в активный блок; её резонансное возбуждение не делится на малый знаменатель.
Устойчивость этих движений и их возникновение из начального состояния типа I
не установлены. OUT-операторы, M1/M2 и Floquet здесь не вычислялись.

## Источники и неизменность постановки

β=15°, μ=0.1, L₁=0.9,L₂=1.1; оба b=0.05,h=0.10;
E=ρ=1,ν=0.3,κ=5/6,p12. Заданы ровно четыре A=.005,.006,.008,.010,
где A=max|w(s,0)|. Зафиксированы прежние заделки, идеальный узел/c₁=c₂,
local/global maps, quartic action, variable mass, nonlinear inertia и kernels.
Нет новых ODE,FEM,геометрий,амплитуд или исходного двухпараметрического вывода.

Использованы [действие](../theory/weakly_nonlinear_spatial_rod.md),
[его полиномиальное разложение](../theory/weakly_nonlinear_spatial_rod_expansion_generated.md),
[канонический аудит](nlsp_out_of_plane_instability_feasibility.md) и
[предыдущий PL audit](nlsp_beta15_nonlinear_planar_resonance.md).
Parent: `results/nlsp_beta15_nonlinear_planar_resonance/4faf1f36cd1017cd/`,
manifest SHA256 `2365f2b0c8d566808c4ebbdcba372f1f51f31171e501fa86f35b7f424b3f2562`.
Новый bundle: `results/nlsp_beta15_planar_1to3_internal_resonance/b59944c2ef43d0b7/`.

Переиспользованы все 92 сохранённые planar eigenvectors; новых спектров нет.
Их общая mass orthogonality error 2.66e−15, eigen-equation residual 2.08e−15.
Для отдельных active форм residual соответственно1.14e−11 и6.99e−13;
агрегированная норма не используется для сокрытия ошибки низших форм.

| Локальная sorted форма | ω, из общего saved p12 eigenbasis | Изгибно-сдвиговая доля |
|---|---:|---:|
| 2 | 0.4348207304902722 | 0.960266266 |
| 4 | 1.3063293949106531 | 0.977523101 |

δ₁₃=ω₄−3ω₂=0.00186720343983660. Разность с прежней отдельно вычисленной
ω₂=0.434820730491217 составляет 9.45e−13; коэффициенты используют один
общий eigenbasis. Обе формы включают u,w,θ,c и удовлетворяют joint mappings.
φ₂ ориентирована по прежнему initial shape; у φ₄ знак наибольшего |w| положителен.
Продолжение по β не заявляется.

## Каноническое взаимодействие из принятого действия

При φᵢᵀM₀φⱼ=δᵢⱼ положено q=ΣφⱼQⱼ,Π=M₀ΣφⱼPⱼ,

\[
a_j=(\sqrt{\omega_j}Q_j+iP_j/\sqrt{\omega_j})/\sqrt2,
\quad \{a_j,\bar a_k\}=-i\delta_{jk},\quad I_j=|a_j|^2.
\]

Legendre expansion сохранённого действия даёт

\[
H_3=V_3-T_3(q,M_0^{-1}\Pi),\qquad
H_4=V_4-T_4+\tfrac12(T_{3,v})^TM_0^{-1}T_{3,v}.
\]

Последний член необходим для переменной массы. Это аналитическое
приближение Hamiltonian, не изменение physical RHS или его обращения массы.
Прямое ограничение на две формы не принимается как окончательный коэффициент.

Устраняемые кубические monomials имеют частотную комбинацию D≠0.
Использовано W₃=iH₃/D, так что {H₂,W₃}=−H₃ и

\[
H_{4,\mathrm{eff}}=H_4+\tfrac12\{H_3,W_3\}.
\]

Учтены все 90 остальные плоские формы. Для active quartic terms достаточно
H₃ monomials с нулём/одной passive переменной: остальные не дают четыре
active фактора после одного Poisson bracket. Минимальный кубический
частотный знаменатель 0.041045354753; это не δ₁₃.
Близкий quartic monomial a₂³bar(a₄) **не устраняется**.

Сохранённая каноническая система:

\[
H_{\rm res}=\omega_2I_2+\omega_4I_4+
\alpha_{22}I_2^2+\alpha_{24}I_2I_4+\alpha_{44}I_4^2
+\gamma a_2^3\bar a_4+\bar\gamma\bar a_2^3a_4.
\]

| Коэффициент | Прямой | Непрямой | Сумма |
|---|---:|---:|---:|
| α₂₂ | 20949.431499 | −14048.339560 | 6901.091939 |
| α₂₄ | 54448.211439 | −46823.625544 | 7624.585895 |
| α₄₄ | 26697.669919 | −16940.071212 | 9757.598707 |
| γ | −5912.352883 | +3053.768713 | **−2858.584170** |

Это коэффициенты при указанной mass/canonical normalization, не величины
в единицах 1/time, автоматически сопоставимые с прежним OUT coupling g.
γ имеет отрицательный знак при выбранном знаке φ₄; изменение фазы формы
меняет представление γ, сохраняя Hamiltonian и физическое взаимодействие.

Плоские T₃=j∥cθ̇²,T₄=j∥c²θ̇²/2 сохранены. В V₃ участвуют,
например, −Cνcθ²/2,−Cu_sθ²/2+Su_sθ² и
(Cνc+Cu_s−Su_s)θw_s. V₄ содержит θ⁴,θ³w_s,
θ²u_s²,θ²w_s² и θ²cu_s. Они дают смешанные modal projections;
узел меняет перекрытия через прежнюю вариационную сборку, не добавляя энергии.
Для γ прямой V₄ вклад −5912.352892, прямые rotational Legendre terms
+8.93e−6; непрямой V₃–V₃ +3053.703718, T₃–T₃ −0.001398,
mixed terms +0.066394. В данном случае геометрическая энергия доминирует;
инерция проверена и сохранена, а не заранее отброшена.

При последовательном включении passive sorted modes γ равен
−6025.51(8),−4959.64(16),−3450.19(32),−2947.14(64),−2858.58(92).
Высокие плоские моды существенны. Итог использует полное p12 пространство;
это не доказательство сходимости по p.

Reality defects H₃,W₃,H₄eff ≤1.96e−15. Коэффициентная projection H₃
и её полный passive gradient совпадают с прямой action assembly до 1.12e−15.
Из α₂₂ получено Ω₂=11.4382901441 в прежней physical-amplitude convention;
разность с PL limit 2.13e−10. Проверка scalar variational limit независимо
подтверждает знаки Lie correction; artificial Duffing model не вводилась.

## Периодическое движение и модуляция — разные начальные состояния

При aⱼ=√Iⱼexp(−iθⱼ),Ψ=3θ₂−θ₄, резонансный член равен
2γI₂^(3/2)√I₄cosΨ. Из канонических уравнений следуют

\[
\dot I_2=6\gamma I_2^{3/2}\sqrt{I_4}\sin\Psi,
\qquad \dot I_4=-2\gamma I_2^{3/2}\sqrt{I_4}\sin\Psi,
\quad J=I_2+3I_4=\mathrm{const}.
\]

Самостоятельные и перекрёстные diagonal frequency shifts определяются
2α₂₂I₂+α₂₄I₄ и α₂₄I₂+2α₄₄I₄. Диагональная расстройка

\[
\delta_{13,\mathrm{diag}}=\omega_4-3\omega_2+
(\alpha_{24}-6\alpha_{22})I_2+(2\alpha_{44}-3\alpha_{24})I_4
\]

дополняется фазовым членом γ(I₂^(3/2)/√I₄−9√(I₂I₄))cosΨ.
Stationarity требует sinΨ=0 и равенства полной поправленной частоты формы4
утроенной частоте формы2. Проверены обе real phase choices Ψ=0,π.
Не утверждается исчерпывающий поиск всех периодических ветвей.

**Тип I:** physical initial Q₄=0, состояние по φ₂ с нулевыми скоростями.
Для чистого canonical a₄=0 имеется ненулевой source
dot(a₄)=−iγa₂³. Это поддерживает медленный обмен/модуляцию, но не
назначает ему конкретную trajectory. Обратная quadratic canonical map
показывает малые initial a₄ и passive amplitudes даже при physical Q₄=0.
Свободные passive transients и ошибки более высокого порядка записаны
как ограничения. Type-I modulation не интегрировалась и не подтверждена
невязкой полного модулированного решения.

**Тип II:** совместное фазово согласованное движение с ненулевой формой4
и согласованными NR corrections. Только оно ниже прошло numerical gates.
Оно не возникает автоматически из специально не подготовленного состояния типа I.

## Алгебраическая проверка периодических кандидатов

Normal-form stationary equations дали начальные guesses. Их вспомогательный
root-sum-square amplitude scale **не принимался** за physical A.
Затем решались три скалярных уравнения по x,y,Ω:
две резонансные проекции full planar mass-form equations и точная
пространственная нормировка max|w(s,0)|=A.

\[
q(\tau)=x\varphi_2\cos\tau+y\varphi_4\cos3\tau+q_{\rm nr}(\tau).
\]

Все остальные p12 Fourier coefficients вычислялись как passive algebraic
corrections через spectral resolvents, с восемью фиксированными итерациями.
Active slots (mode2,harmonic1) и (mode4,harmonic3) никогда не делились
на линейные denominators. Это три nonlinear root unknowns, не новая
глобальная многомодовая normal form. Основное усечение — гармоники0…9,
контроль —0…7; 64 фазы для coefficients, 53 другие фазы со сдвигом0.173
для полной невязки. Скорости при τ=0 точно нулевые; closure за2π следует
из cosine representation. Не запускалась IVP или periodic-orbit time search.

| A | Ω, Ψ=0 | x, mass-normalized Q₂ | y, mass-normalized Q₄ | E₄/(E₂+E₄), начальная reference | Full residual | Статус |
|---:|---:|---:|---:|---:|---:|---|
| .005 | .435115586714 | 3.19171259e−4 | 1.01213918e−5 | 0.8995% | 1.213e−7 | PASS |
| .006 | .435238797040 | 3.85702215e−4 | 2.36343098e−5 | 3.2779% | 1.410e−6 | PASS |
| .008 | .435471177662 | 5.04860330e−4 | 9.85364159e−5 | 25.5855% | 6.469e−5 | PASS |
| .010 | .435686014899 | 5.90632808e−4 | 2.00570935e−4 | 51.0007% | 3.995e−4 | PASS |

Reference energies Eⱼ=ωⱼ²Qⱼ²/2 относятся к ведущим компонентам
при нулевых initial velocities. Это не полная nonlinear energy partition
и не отдельно сохраняемые энергии форм. NR coefficient norm ratio
√Σ_n(q_nr,nᵀK₀q_nr,n) / corresponding active norm =4.73…6.54%.
Норма гармонических коэффициентов не выдаётся за усреднённую энергию.

При .008 также найден Ψ=π: Ω=.435717443713,E₄ reference fraction43.18%;
при .010 — Ω=.436164214738,fraction7.903%. Все шесть найденных кандидатов
проходят numerical checks; ни один не назван устойчивым.
Отсутствие других ветвей при .005/.006 в ограниченном наборе guesses
не является доказательством их отсутствия.

Для Ψ=0 max|w| leading mode4=.0001504,.0003511,.0014638,.0029795.
Её K₀ norm относительно mode2 составляет9.53%,18.41%,58.64%,102.02%.
Оставшаяся **нерезонансная** третья гармоника имеет max|w|
лишь1.06e−6…3.64e−6. Следовательно, большой третий harmonic content
не исчез: он отделён как реальное active участие формы4, а не ошибочно
назван малой обычной поправкой. Среднее смещение max|w|=.0001429… .0006348,
вторая гармоника=.0001030,.0001405,.0001698,.0001556. Данные не сглажены.

Новые Ψ0 частоты отличаются от прежних **формальных**, не принятых PL
значений на −1.046e−5,−2.822e−5,−1.673e−4,−4.477e−4.
Это существенная поправка к исходной OUT-гипотезе; прежний M0 прогноз
нельзя переносить на эту другую пространственную форму основного движения.

## Контроль точности и пригодность к следующей задаче

До результатов зафиксированы gates: полный relative mass-form residual≤1e−3,
relative energy variation≤1e−4, observed frequency control≤2e−5,
NR coefficient norm ratio≤0.2 и action/projection gate2e−10.
В full residual normalization используется линейный инерционный масштаб
Ω²√(x²+81y²) в resting-mass coordinates.

На независимых фазах полный energy variation для основной ветви
8.02e−10…6.87e−6. Direct physical mass-form/independent action difference
≤1.18e−14 относительно суммы норм членов уравнения.
Относительная разность **двух почти нулевых residuals** сохранена отдельно;
она не используется как operator gate. Joint compatibility<1.12e−16,
масса в проверенных phase states имеет минимальное eigenvalue≈0.96748.

7/9 гармоник: максимальная |ΔΩ| основной ветви7.03e−8,
relative change y≤1.77e−4. Линеаризованный индикатор ошибки от одного
дополнительного passive update≤3.93e−10 для этой ветви.
Никакой новый производственный solve с изменённым числом итераций не проводился.
Scaled stationary Jacobian condition основной ветви688…1446;
вторая ветвь достигает4763. Невязки корней≤4.0e−13 в scaled equations.

Это наблюдаемые numerical indicators, не строгая апостериорная оценка.
Но они значительно меньше прежнего OUT margin2.675e−4; частота
приближения разрешена для последующей **постановки**, без вывода о знаке OUT роста.
Все90passive modes включены, поэтому within-p12 modal truncation γ отсутствует.
Сходимость по p, строгая existence настоящей periodic orbit и её stability
не проверены. Leading quartic normal form и уточнённый full-cubic HB waveform
различаются по порядку приближения; NF detuning через quadratic inverse map
маркирована как дополнительная диагностика, не точная finite-amplitude величина.

Пригодна конкретная проверенная **Type-II** ветвь Ψ=0. Для исследования
исходного Type-I движения остаются modulation/passive transients;
его нельзя молча заменить периодической Type-II орбитой.
OUT-усреднение, tangent trajectories и полный nonlinear контроль не запускаются.

## Воспроизводимость, тесты и остановка

Научные вычисления под таймером **135.056 s из600 s**,80initial canonical
evaluations, затем checks/decomposition; всего182canonical evaluations,
1803HB forcing evaluations и232746homogeneous residual evaluations.
Новые ODE,FEM,frequency eigensolves,OUT/Floquet calculations — **0**.
39 targeted tests PASS: 14новых и25 прежних action/PL/canonical checks.
Проверены независимая projection, Lie signs,conjugacy,J invariant,
normalization,symmetry,residual/harmonic controls,history preservation и cache.

Cache test предыдущего PL отчёта первоначально перезаписал historical PDF
CreationDate и manifest. Массивы не менялись. Исходные PDF bytes и manifest
**точно восстановлены по исходному SHA256**; изменённый тест теперь работает
на временной копии. Evidence сохранена в `historical_cache_restore.json` и
`cache_history_test_failure/`. Initial Type-I array writer также исправлен:
первый seed-only output архивирован, теперь сохранён полный canonical inverse-map
state. Ни один periodic solve не повторялся из-за этих writer/cache исправлений.
Numerical gates не ослаблялись.

[Config](../../data/input/nlsp_beta15_planar_1to3_internal_resonance.json),
[helper](../../scripts/lib/nlsp_planar_1to3_internal_resonance.py),
[tests](../../tests/test_nlsp_planar_1to3_internal_resonance.py),
[manifest](../../results/nlsp_beta15_planar_1to3_internal_resonance/b59944c2ef43d0b7/manifest.json).
Существующий CLI получил режим `--planar-1to3-internal-resonance`.
Новый runnable solver не создаётся; scoped library contract отличается от
одномодового PL и переиспользует его action/mass-form/source helpers.

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
python scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --planar-1to3-internal-resonance --compute
python scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --planar-1to3-internal-resonance --report-only results/nlsp_beta15_planar_1to3_internal_resonance/b59944c2ef43d0b7
python scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --planar-1to3-internal-resonance --plot-only results/nlsp_beta15_planar_1to3_internal_resonance/b59944c2ef43d0b7
```

[Частоты и линейная reference detuning](../../results/nlsp_beta15_planar_1to3_internal_resonance/b59944c2ef43d0b7/two_mode_frequencies.pdf),
[участие форм](../../results/nlsp_beta15_planar_1to3_internal_resonance/b59944c2ef43d0b7/two_mode_participation.pdf).
Две фигуры показывают только четыре заданных A; это не stability map.
Cached compute/report/plot не вызывает scientific routines.

[NLSP-D28](../memory/decisions.md#nlsp-d28),[K29](../memory/knowledge.md#nlsp-k29).
`DEFERRED_TWO_PARAMETER_ACTION_EXPANSION`, физическое обоснование идеального
узла/c₁=c₂ и независимая3Dvalidation остаются отложенными. После отчёта — остановка.
