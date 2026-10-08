# Аналитический по времени контроль ведущего продольного отклика

## 1. Область и исходная остановка

Дата: 2026-10-08. Исходный checkout: main,
`b608aa118247819dfd118ddafb2a0a6155512c80`, Version0.6.1, clean/staging unchanged.
Это отдельная диагностика после `NLSP_PLANAR_TIME_PILOT=PARTIAL` и
`NLSP_PLANAR_SOLVER_RECOVERY=PARTIAL`. Она не заменяет полную cubic IVP.

Источники: [принятая spatial модель](weakly_nonlinear_spatial_rod.md),
[generated expansion](weakly_nonlinear_spatial_rod_expansion_generated.md),
[первый pilot и solver recovery](weakly_nonlinear_planar_time_pilot.md),
[finite M-H/Timoshenko reference](mindlin_herrmann_timoshenko_single_rod.md).
Crespo da Silva1988 — зарегистрированный вариационный antecedent, не источник
конкретных forcing coefficients принятого V0. Новых источников нет.

Защищены V0/T4/V4, cubic residuals, nonlinear RHS/Jacobian, Shen basis,
четыре independent fields, material/geometry/κ, initial pair, BC и старые bundles.
Один прямой G20 fixed-fixed rod: E=ρ=1, ν=.3, b=.20, h0=.05, L=1, κ=5/6.
У u2,c2 essential values и initial values/velocities нулевые. Slopes,
нерастяжимость, static correction и c=-νu_s не вводятся.

## 2. Особая амплитудная специализация начальной задачи

Используется один непрерывный analytic Timoshenko eigenpair старого reference:
max|w_hat|=1; theta_hat делится на тот же signed midpoint factor. Сохранённые
ω1≈.317474290788 и T1≈19.7911625902 одинаковы для всех p; новых roots нет.

$$\epsilon_a=A/h_0,\qquad W=h_0\widehat w,\qquad
\Theta=h_0\widehat\theta.$$

W имеет размерность длины, Θ безразмерна. Первый порядок:

$$w_1=W\cos\omega_1t,\qquad \theta_1=\Theta\cos\omega_1t.$$

$$w=\epsilon_aw_1+\cdots,\quad\theta=\epsilon_a\theta_1+\cdots,
\quad u=\epsilon_a^2u_2+\cdots,\quad c=\epsilon_a^2c_2+\cdots.$$

Это разложение конкретного решения с заданным planar bending initial state,
не изменение общего порядка полей в семиполевой теории. Коэффициенты m,jp,C,H,S
прочитаны из audited model, без EA scaling по εa. Физический leading response
u_asym=εa²u2, c_asym=εa²c2 и аналогично velocities. Одна решённая u2,c2 система
даёт обе исходные amplitudes .05/.025; четверть scaling — алгебра, не новый
эмпирический amplitude law.

## 3. Независимое выделение второго порядка

При jet weights u,c→2; w,theta→1 и исключённых out-of-plane jets выделены
εa² terms из обоих точных Fraction residual representations A/B. Независимо
varied weighted quartic action; все coefficient differences нулевые.
Символы u,c в извлечённых Polynomial обозначают их second-order coefficients.

$$m u_{2,tt}-C u_{2,ss}-\nu C c_{2,s}=\partial_s Z_t,$$

$$j_p c_{2,tt}-Hc_{2,ss}+C(c_2+\nu u_{2,s})
=-\nu C(\theta_1w_{1,s}-\theta_1^2/2)+j_p\theta_{1,t}^2,$$

где

$$Z_t=(C-S)\theta_1w_{1,s}+(S-C/2)\theta_1^2.$$

Знак inertial source +jp*theta1_t² следует из variation kinetic energy и
обоих residuals. Forcing не содержит неизвестные u2,c2; слева тот же linear
M-H block. Источники для w2,theta2 точно нулевые. Соответствующий homogeneous
linear bending IBVP с zero initial/essential data даёт w2=theta2=0;
это вывод данного порядка, не удаление полей из полной модели.

Размерности: m — kg/m, jp — kg*m, C,S — N, H — N*m²; Z — N;
axial strong source ∂sZ — N/m, contraction source — N. u2 — длина,
c2 — безразмерная coordinate. All-term dimension tests проходят.

## 4. Постоянная и 2ω1 части источника

Положим Ωd=2ω1,

$$Z=(C-S)\Theta W_s+(S-C/2)\Theta^2,\qquad D=\Theta W_s-\Theta^2/2.$$

$$F_{u0}=F_{u2}=Z_s/2,$$

$$F_{c0}=(-\nu C D+j_p\omega_1^2\Theta^2)/2,\qquad
F_{c2}=(-\nu C D-j_p\omega_1^2\Theta^2)/2.$$

Тогда F=F0+F2*cosΩdt. Проверены прямые harmonic substitutions, включая t=0
и T1/4; отдельно сохранён противоположный знак inertial части в Fc2.
Это внутреннее действие prescribed first-order bending, не внешнее нагружение.

## 5. Слабая M-H система и квадратура

Каждое поле сохраняет независимые p−1 Shen coefficients:
B_n=P_n−P_(n+2), ξ=2s/L−1, n=0,…,p−2. Essential endpoints нулевые;
дополнительных slope constraints нет. Обратимое M0 whitening повторяет
existing helper с matrix nq=2p+1 и не меняется при forcing quadrature refinement.

$$M\ddot x+Kx=f_0+f_2\cos\Omega_dt.$$

Для физических Bu,Bc:

$$M_{uu}=m\int B_u^TB_u\,ds,\quad M_{cc}=j_p\int B_c^TB_c\,ds,$$

$$K_{uu}=C\int B_{u,s}^TB_{u,s}\,ds,\quad
K_{uc}=\nu C\int B_{u,s}^TB_c\,ds,\quad K_{cu}=K_{uc}^T,$$

$$K_{cc}=C\int B_c^TB_c\,ds+H\int B_{c,s}^TB_{c,s}\,ds.$$

$$f_{u0}=f_{u2}=-\frac12\int B_{u,s}^TZ\,ds,\qquad
f_{c0,2}=\int B_c^TF_{c0,2}\,ds.$$

Минус axial weak source получен интегрированием по частям; граничная часть
равна нулю по test functions, а не вследствие искусственного зануления Z_s.
Strong/weak assembly сравнивается на одинаковых quadrature/function sets.
Independent M,K совпадают с frozen planar helper MH submatrix; SPD/symmetry
и mass orthogonality проверены. Ширина/section factors не перенормированы.

Matrices полиномиальны; analytic W,Θ во forcing неполиномиальны. Поэтому
nq=2p+1 не объявляется точной forcing quadrature. Проверены последовательности
2p+1→3p+3→4p+5; coefficient changes и strong/weak differences должны быть≤2e−11,
значительно меньше принятого spatial gate1e−3. Никакой reduced integration нет.

## 6. Exact-time evaluator всех координат

$$KV=MV\operatorname{diag}(\omega_j^2),\qquad V^TMV=I,
\qquad b_{0,2}=V^Tf_{0,2},\qquad x=Vy.$$

Сохранены все 2(p−1) columns. Это вычисление полной finite-dimensional matrix
function, не modal truncation. Acoustic/contraction coordinates, высокие
частоты и малые amplitudes не удаляются. Одна primary generalized eigensolve
на p; для optional extension прежние spectral matrices восстанавливаются
с hash/coefficient checks без нового eigh.

$$y_j=b_{0j}H(\omega_j,0,t)+b_{2j}H(\omega_j,\Omega_d,t),$$

$$H(\omega,\Omega,t)=\frac{\cos\Omega t-\cos\omega t}{\omega^2-\Omega^2}
=\frac{t^2}{2}\operatorname{sinc}_u(at)\operatorname{sinc}_u(bt),$$

где a=(ω+Ω)/2, b=(ω−Ω)/2, sinc_u(x)=sinx/x, sinc_u(0)=1. NumPy pi-normalized
sinc не подменяет эту функцию. Без division по малой detuning:

$$H_t=\frac t2[\cos(at)\operatorname{sinc}_u(bt)+
\operatorname{sinc}_u(at)\cos(bt)],$$

$$H_{tt}=\cos(at)\cos(bt)-(a^2+b^2)H.$$

Прямое differentiation даёт H_tt+ω²H=cosΩt. При ω=Ω>0
H=t*sinωt/(2ω); при t→0 H≈t²/2. Близкая ненулевая detuning не округляется
до нуля. Synthetic tests включают Ω=0, exact/near resonance, малое t и разные
frequency scales; также ω=Ω=0 даёт t²/2,t,1.

Zero ICs проходят; M*x_tt(0)=f0+f2. Независимый augmented constant-matrix
exponential включает state generators z0'=0,zcos'=−Ωd*zsin,zsin'=Ωd*zcos
с z0=zcos=1,zsin=0. Проверены несколько отдельных времён p16; ODE integrator
не используется. Экспонента имеет собственную numerical error, не machine-zero
requirement. Detunings сохранены как diagnostics, не resonance study.

## 7. Forced power и начальная совместность

$$E_2=\tfrac12\dot x^TM\dot x+\tfrac12x^TKx,\qquad
\dot E_2=\dot x^T(f_0+f_2\cos\Omega_dt).$$

Energy E2 не обязана быть постоянной: prescribed bending передаёт мощность
этой выделенной подсистеме. Полная cubic модель и её conservation не менялись.
Power identity используется для verification, не energy modal classification.

На концах Θ=0, поэтому continuous initial trace

$$u_{2,tt}|_{\partial}=(C-S)\Theta_sW_s/m.$$

Для G20: W_s=(+.004135580879,−.004135580879),
Θ_s=(1.369151063,1.369151063); u2_tt=(+.004407417320,−.004407417320).
Multiplication εa² при A=.0025 даёт прежние ±1.101854330e−5.
Initial mismatch сохраняется. Galerkin acceleration нулева на концах по
basis, что не зануляет предельный continuous strong source. Weak solution,
clamp reactions и требование гладкости до endpoints различаются.
Ни forcing, ни initial bending не исправлялись; причинность всех spatial
разностей только этим mismatch не доказана.

## 8. Sampling, нормы и comparison semantics

Основной interval0…5T1, отдельно0…0.1T1. Integration step отсутствует; есть
sampled times аналитического evaluator. Для каждой пары используется общий
grid с16 samples/period самой высокой сохранённой частоты. Итоговая пара
повторена с32 samples/period. Maxima refinement — numerical diagnostic,
не certificate непрерывного supremum; target относительного изменения1%.

Physical L2 вычисляется точно через Legendre Gram weights L/(2n+1), после
преобразования обратно из whitening. Max — прежние100 Gauss spatial nodes.
Даны max-in-time L2/max-in-space-time, absolute values, собственный ненулевой
reference scale, тот же relative floor coefficient1e−10 и unchanged gate1e−3.
Фазы/амплитуды/periods не подбираются; края и начальный interval не исключаются.

Обработка time blocks2048; все modal coordinates участвуют в каждом evaluation.
Сохраняются spectral representations, streaming maxima/integrals, компактные
norm samples и selected full1001-point profiles. Полные2001-point observations
на5T1 являются coarse observations и не certificate разрешения быстрых полей.
Для response figures сохранён отдельно resolved short grid по максимальной ωp.

Physical L2 projector старого recovery переиспользуется для p24→32,
p48→64 и p64→96; проверяются tail/common split, Pythagoras и orthogonality
для u2,c2 и обеих velocities. Это approximation diagnostic, не energy fractions.

Old nonlinear comparison вычисляет analytic response точно на сохранённых
временах, без interpolation. Разделены same-p comparison, difference-of-p
comparison и amplitude normalization. Same-p remainder включает следующие
амплитудные порядки; second-order control не является exact cubic IVP truth.
Старый nonlinear IC — finite-p projection continuous pair; наш forcing — один
continuous pair. Projection/frequency discrepancy отдельно qualified;
ω1 не заменяется на ω1,p.

## 9. Results и границы вывода

Ниже приведены фактические результаты завершённого bounded run. LONG CLOSED, EB/RLB-KV PAUSED и angular same-clamp
out-of-plane reference UNAVAILABLE сохранены. Full nonlinear spatial
convergence остаётся unresolved независимо от статуса leading diagnostic.


### Пространственные разности на0…5T1

В ячейках relative L2 / relative max. Последняя строка использует сгущённую
32-samples grid; прежний gate1e−3 в обеих нормах не менялся.

| Пара p | u2 | c2 | u2_t | c2_t |
|---|---:|---:|---:|---:|
| 16→24 | 0.0228621 / 0.0337617 | 0.403213 / 0.446831 | 0.0923198 / 0.148072 | 0.740068 / 0.674832 |
| 24→32 | 0.00271001 / 0.00384253 | 0.030589 / 0.0455912 | 0.0219773 / 0.0335837 | 0.0575671 / 0.077636 |
| 32→48 | 0.00365501 / 0.00466146 | 0.149938 / 0.18139 | 0.0420199 / 0.0507992 | 0.281066 / 0.284942 |
| 48→64 | 0.000376484 / 0.000601205 | 0.00684992 / 0.0125914 | 0.00503974 / 0.00818417 | 0.0145532 / 0.0226349 |
| 64→96 | 0.00010286 / 0.000183729 | 0.00159964 / 0.00294958 | 0.00169123 / 0.00347779 | 0.0041746 / 0.00658215 |

В p48→64 и p64→96 только u2 проходит обе нормы. Последний p96 не является
точным continuum solution и не разрешает автоматически p128 или выше.
Уточнение32→48 увеличивает c/velocity differences; это немонотонность
конечных пространственных сравнений, не основание переставить/фазово выровнять решения.

### Общий начальный участок0…0.1T1

| Пара p | u2 | c2 | u2_t | c2_t |
|---|---:|---:|---:|---:|
| 16→24 | 0.00758779 / 0.0119397 | 0.0406602 / 0.0764294 | 0.0418003 / 0.0857706 | 0.0720315 / 0.134843 |
| 24→32 | 0.0014265 / 0.00228501 | 0.0131633 / 0.0195431 | 0.0141092 / 0.0237378 | 0.0246685 / 0.0351665 |
| 32→48 | 0.00272302 / 0.00302403 | 0.0273192 / 0.0320993 | 0.031407 / 0.0365396 | 0.051688 / 0.0607951 |
| 48→64 | 0.000155251 / 0.000266117 | 0.00240749 / 0.00479046 | 0.00230775 / 0.00493349 | 0.00546085 / 0.0115202 |
| 64→96 | 4.88118e-05 / 9.85268e-05 | 0.000884649 / 0.0016574 | 0.000998981 / 0.00266531 | 0.00244304 / 0.00442895 |

Здесь reference scales берутся на коротком interval; они не тождественны
полным5T1 denominators. В итоговой короткой паре также только u2 проходит
обе нормы; отдельный L2 PASS для других components не означает общего PASS.

### Проверки и сохранённые координаты

| p | Независимых координат / сохранённых eigenvectors |
|---|---:|
| 16 | 30 / 30 |
| 24 | 46 / 46 |
| 32 | 62 / 62 |
| 48 | 94 / 94 |
| 64 | 126 / 126 |
| 96 | 190 / 190 |

Фактические maxima: forcing refinement8.75e-14, strong/weak6.07e-14,
eigenpair residual1.01e-15, mass orthogonality3.33e-15,
forced equation2.32e-12, forced power6.69e-13.
Independent p16 augmented-expm comparison maxrelative=6.24e-12.
SPD и все zero IC/initial acceleration checks прошли. Матрицы/ω/V/f/b
сохранены для каждого p; ничего не отбрасывалось.

Final sampling:799581→1599159 actual times; max relative change measured
peak/reference=2.71e−5, target1%=PASS. Это sampling evidence, не continuous
supremum proof. Projection scaled identities в итоговой паре≤2.90e−19;
tail fractions squared-L2 time integral: u2=3.953%,c2=10.414%,
u2_t=7.471%,c2_t=12.756%; common-space evolution dominates again.

### Сопоставление с сохранёнными nonlinear histories

| Case | Фактический interval/T1 | L2 u2 | L2 c2 | L2 u2_t | L2 c2_t |
|---|---:|---:|---:|---:|---:|
| p24_large | 5 | 0.0240204 | 0.0118378 | 0.0303991 | 0.00275262 |
| p32_large | 5 | 0.0240306 | 0.0118416 | 0.0304272 | 0.00275252 |
| p24_small_prefix | 2.64813 | 0.00263143 | 0.00135252 | 0.00230396 | 0.000210013 |
| p32_small | 5 | 0.00602393 | 0.00296965 | 0.00762976 | 0.000690598 |
| p48_short | 0.1 | 0.000516263 | 0.000259394 | 0.000123043 | 7.66365e-06 |

Это norm разности NL/εa² и same-p leading response, с leading reference scale.
p24 small заканчивается только на2.648132T1; p48 — только на0.1T1.
Сравнения делаются по точным timestamps, duplicate final snapshots не
превращаются в более поздние данные. Уменьшение норм при второй исходной
amplitude наблюдается на одинаковом full p32 interval, без нового power-law fit.

Для difference-of-p24/32 при большой amplitude L2 space-time correlations
с εa² leading difference:
u2=0.999890462, c2=0.999999609, u2_t=0.999990700, c2_t=0.999999938.
Relative L2 remainder на масштабе самой NL difference:
0.0236751 / 0.00154438 / 0.00655289 / 0.000573315.
Это сильное evidence, что значительная часть spatial difficulty уже присутствует
в leading forced system без ODE error; весь nonlinear remainder не объявлен
объяснённым по correlation. Representative difference profiles сохранены.

Из independent background audit: projectionL2 relative порядка1e−14,
ω discrepancy порядка1e−12, sampled full5T1 linear-background difference≤2.9e−10.
Исторический preflight и independent audit дают согласованные малые bounds;
они не делают continuous forcing exact historical semidiscrete asymptotic coefficient.

### Provenance, стоимость и cache

Final bundle `results/planar_second_order_axial_response/b3ea4eb6ac95d6e1/`.
Primary source execution `4ba7abd6e2e8b064` сохранён; optional p96 впервые
сохранил spectral matrices в checkpoint `78c83d20c48e1d17` до metadata error
от relative parent path. После path fix p64 и p96 восстановлены с checks
и0 повторных eigensolves; все p96 arrays побитово совпадают с checkpoint.
Метаданные и archive исполнявшихся CLI/helper отделены от численных данных.
Исторические c972/054874 manifests/hashes неизменны, current code hash к ним
не предъявлялся. Отдельные exact Fraction и background audits входят в bundle.

Charged primary+conditional work373.599s/1200s;
учтены prior20s и failed metadata attempt2.283s. Main primary MH eigh6 —
один на каждый p16/24/32/48/64/96; cached restores2, augmented expm6,
analytic evaluator block calls8064, ODE/PDE integrations0, continuum roots0.
Polynomial model получен один раз на процесс и переиспользован для всех p;
из-за primary/metadata-failure/resume процессов derivation calls3, что записано
отдельно. Validation calls/tests не выдаются за primary eigensolves; independent
background audit имел4 Timoshenko validation eigh calls (2 утраченны при
сериализации,2 сохранены) и не запускал ODE. Фактический numerical work
остаётся значительно ниже fixed1200s; отмены из-за близкого прогноза нет.

Три PDF/PNG figures: `leading_response`, `spatial_convergence`,
`nonlinear_comparison`. Response/overlay показывают resolved initial0…0.1T1
при s/L=.25; convergence figure использует полный5T1 interval. Без галереи,
energy labels, phase alignment или обработки full cubic feedback.

```powershell
python scripts/analysis/verify_planar_second_order_axial_response.py --compute
python scripts/analysis/verify_planar_second_order_axial_response.py --report-only results/planar_second_order_axial_response/b3ea4eb6ac95d6e1
python scripts/analysis/verify_planar_second_order_axial_response.py --plot-only results/planar_second_order_axial_response/b3ea4eb6ac95d6e1
```

Matching default compute использует conditional coverage только при совпадении
current code/config/source identity и полном наличии primary grid. Report/plot
получают явный путь и делают0 eigendecompositions, analytic evaluations,
derivations и time integrations. Code identity включает section factory.
`--extend-p96 <primary bundle>` — явно условный optional step; checkpoint
option — восстановление уже сохранённых spectral matrices после metadata
failure, не новая resolution series. На matching final cache он также не решает
матрицы заново. После bounded report остановка.

### Раздельные статусы

| Статус | Результат |
|---|---|
| NLSP_SECOND_ORDER_DERIVATION | PASS |
| NLSP_SECOND_ORDER_FORCING_ASSEMBLY | PASS |
| NLSP_SECOND_ORDER_EXACT_TIME_EVALUATOR | PASS |
| NLSP_SECOND_ORDER_SPATIAL_CONVERGENCE | PARTIAL |
| NLSP_SECOND_ORDER_NONLINEAR_COMPARISON | COMPLETE_WITH_ACTUAL_PREFIX_QUALIFICATIONS |
| NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC | COMPLETE |

Диагностика COMPLETE означает выполненный объявленный контроль; spatial
convergence PARTIAL остаётся самостоятельным результатом. На принятом gate
нет evidence достаточности p48/p64 для всех leading fields/velocities; даже
p64→96 сохраняет FAIL для c2 и обеих speeds. Последний p не exact truth.
У full cubic solution дополнительные feedback/amplitude effects; требуемый
p полной задачи этим контролем не установлен. Initial smoothness mismatch
сохранился и согласуется с наличием пространственной трудности, но не
доказан её единственной причиной. V0/RHS/IC/BC/basis не изменены.
Новый nonlinear p48, Floquet, angles, out-of-plane perturbations и thresholds
не считались; старые PARTIAL не повышены. Ни p128+, ни новый solver/initial
correction следующий этап автоматически не выбираются.


### Targeted regressions и границы проверки

220 уникальных checks PASS:215 combined для нового exact-time helper/CLI,
старого recovery/planar helper, spatial action, finite M-H и Bishop, ещё3
rectangular Timoshenko coefficient/convention checks и2 новых source-cache checks.
52 новых targeted cases:50 входят в215,2 отдельно проверяют source identity. Final runs13.94s и.78s;
дополнительный ранний targeted run2.80s, source-cache checks.44s и rendering2.69s учитываются отдельно
в verification budget charge (общий conservative charge≈396.92s<1200s).

Три старых tests намеренно deselected во всех совместных regressions:

- `test_full_semidiscrete_linear_reference_exact_time_against_short_integrator`;
- `test_pure_axial_short_time_control_matches_full_linear_reference`;
- `test_tiny_planar_nonlinear_smoke_keeps_all_four_fields`.

Они вызывают time integrators и нарушили бы нулевой ODE contract. Новые
cache/report/plot tests запрещают derive/eigh/compute paths; loader tests
запрещают eigensolve и проверяют corrupted matrices/forcing rejection.
Validation fixtures и старые linear regressions имеют собственные небольшие
eigen/expm calls; число6 относится к primary diagnostic cases, не ко всем
вызовам linear algebra во время разработки/tests. Production/primary counters
и additional validation evidence в manifest/verification различены.

Links, append-only memory prefixes, protected file/source hashes, staging и
`git diff --check` проверены. [NLSP-D04](../memory/decisions.md#nlsp-d04) и
[NLSP-K04](../memory/knowledge.md#nlsp-k04) дополняют прошлую историю, не
переписывают D02/K02/D03/K03. README/CHANGELOG/navigation обновлены;
assumptions/equations.tex/source index/BibTeX и frozen scientific code/results
не менялись. Specialization εa относится только к этому диагностическому
решению и описана здесь, не promoted в новую production physics.


Final cache identity дополнена прямыми hashes audited coefficient result и
linear reference result/manifests. Actual optional execution408f79b9042e7b64
сохранена отдельно; finalb3ea4eb6ac95d6e1 копирует её numerical artifacts.
`cache_identity_revision.json` и exactdiff подтверждают, что кроме identity
весь CLI AST неизменён;0 новых analytic/eigen/ODE evaluations. Historical
source manifests проверены повторно, source result hashes совпадают. Это
исправление provenance bookkeeping, не ещё один научный расчёт.

Контроль u2/c2 не устанавливает convergence theta_t полной nonlinear IVP:
её higher-order bending feedback здесь не моделируется. Нет ни spatial
acceptance полной cubic задачи, ни доказательства единственной причинности
initial mismatch, ни выбора следующего метода/этапа.
