# Первый плоский nonlinear time pilot прямого стержня

Дата: 2026-10-07. Разрешён один ограниченный расчёт задачи Коши для
плоского подпространства [проверенной кубической пространственной
модели](weakly_nonlinear_spatial_rod.md). Это пространственная численная
дискретизация распределённого стержня с четырьмя независимыми полями,
а не двухмодовая модель или поиск nonlinear normal mode.

Новый этап не меняет принятые $V^{(0)}$, cubic equations, source variants,
линейные production solvers и generated appendix. Исторические остановки
`LONG closed` и `EB/RLB-KV paused` сохраняются. Прямой плоский расчёт не
заменяет отсутствующий angular out-of-plane same-clamp reference.

## 1. Задача, геометрия и essential BC

Нормированный G20 control:

| Параметр | Значение |
| --- | --- |
| $E,\rho$ | 1, 1 |
| $\nu$ | 0.3 |
| $b,h,L$ | 0.20, 0.05, 1 |
| $\kappa$ | 5/6, прежнее project value |
| Independent fields | $(u,w,\theta,c)$ |

Одна материальная координата $s\in[0,L]$. На обоих концах
$u=w=\theta=c=0$: закрепляются перемещение, поворот сечения и разрешённая
contraction coordinate. Нулевые slopes $u_s,w_s,c_s$ не назначаются;
`book_slope_clamp` не используется. В задаче нет внутреннего физического
узла, угла beta, forcing, damping, spring или начального prestress.

Restriction $v=\Phi=\psi=0$ берётся из уже доказанного invariant planar
subspace. Эти три fields не интегрируются: сохранение их нулевых значений
не является проверкой пространственной устойчивости.

## 2. Какая энергия дискретизируется

Плоские quartic densities получаются из existing `T4,V4` helper
[weakly_nonlinear_spatial_rod.py](../../scripts/lib/weakly_nonlinear_spatial_rod.py).
Restricted kinetic density точно равна

\[
T_{pl,\le4}=\frac m2(u_t^2+w_t^2)+\frac{j_\parallel}{2}c_t^2
+\frac{j_\parallel}{2}(1+2c+c^2)\theta_t^2.
\]

Restricted potential соответствует степени не выше четырёх в

\[
\frac C2(\Gamma_1^2+2\nu\Gamma_1c+c^2)
+\frac H2c_s^2+\frac S2\Gamma_2^2+\frac{B_\parallel}{2}\theta_s^2,
\]

с planar kinematic measures предыдущего аудита. В time evaluator
используются именно precompiled полиномы, а не полные `sin(theta)` и
`cos(theta)`. Никакая additional closure $c=-\nu u_s$, $\theta=w_s$ или
нерастяжимость не вводится. Все четыре generalized fields независимы.

Коэффициенты $m,j_\parallel,C,H,S,B_\parallel$ остаются прежними,
вычисленными по исходному сечению. У неактивного out-of-plane block
сохранён исходный $C_T$ из accepted audit bundle; он не входит в planar
dynamics. Новые correction factors и fitted coefficients отсутствуют.

## 3. Shen Legendre space и квадратура

При $\xi=2s/L-1$ используются

\[
B_n(\xi)=P_n(\xi)-P_{n+2}(\xi),\qquad n=0,\ldots,p-2.
\]

Каждый $B_n$ обращается в ноль на $\xi=\pm1$, но не имеет imposed
нулевого endpoint derivative. Для каждого field создаётся свой набор
из $p-1$ coefficients. Временное решение включает их все; физические
моды для nonlinear evolution не отбираются.

| Maximum degree $p$ | Coefficients per field | Position coordinates | First-order time state | Gauss points |
| --- | --- | --- | --- | --- |
| 16 | 15 | 60 | 120 | 33 |
| 24 | 23 | 92 | 184 | 49 |
| 32 | 31 | 124 | 248 | 65 |
| 48, единственный разрешённый дополнительный уровень | 47 | 188 | 376 | 97 |

Quartic action содержит произведения пространственной степени не выше
$4p$, включая $c^2\theta_t^2$ и $\theta^4$. Gauss rule с
$n_q=2p+1$ точен до $2n_q-1=4p+1$, поэтому достаточен для энергии и
её дискретных derivatives. Повышение квадратуры на одном состоянии
служит independent aliasing control. Фильтрация и reduced integration
не применяются.

Физические raw Shen coefficients связаны с вычислительными через
constant Cholesky scaling. Если $M_{0,raw}=L_ML_M^T$, то
$d_{raw}=L_M^{-T}d_h$. Resting mass в новых coordinates равна
identity с учётом floating representation. Это обратимое изменение
базиса, а не modal truncation и не замена переменной $M_h(d_h)$
на $M_h(0)$.

## 4. Дискретное действие и переменная инерция

Здесь $d_h$ — вектор всех пространственных coefficients. Эта локальная
численная запись не обозначает rotation vector, амплитуду или угол beta.
Подстановка spatial expansions
и точная для polynomials квадратура дают

\[
L_h=\tfrac12\dot d_h^TM_h(d_h)\dot d_h-V_h(d_h),
\qquad
M_h\ddot d_h+g_h+\nabla V_h=0.
\]

Пусть $B_\theta,B_c$ — уже scaled basis evaluation matrices, $W$ —
диагональная матрица quadrature weights, а $c,\theta_t,c_t$ ниже —
значения fields в quadrature nodes. Постоянные mass blocks соответствуют
$u,w,c$; переменный block

\[
M_\theta=j_\parallel B_\theta^T W\operatorname{diag}((1+c)^2)B_\theta.
\]

Вариация того же kinetic law даёт

\[
g_\theta=2j_\parallel B_\theta^TW[(1+c)c_t\theta_t],\qquad
g_c=-j_\parallel B_c^TW[(1+c)\theta_t^2],\qquad g_u=g_w=0.
\]

Оба contributions сохранены. Система
$M_h\ddot d_h+\nabla V_h=0$ не совпадает с принятым действием.
На каждом state соответствующий mass solve выполняется через Cholesky;
дополнительный Taylor series для $M_h^{-1}$ не вводится.

Polynomial potential, gradient и Hessian precompiled один раз из
audited artifact. Производные, basis arrays, quadrature и constant
matrices не генерируются внутри RHS.

### Analytic RHS Jacobian

Для $f=-M_h^{-1}(g_h+\nabla V_h)$ при variation state

\[
\delta f=-M_h^{-1}\left[(\delta M_h)f+
(\partial_d g_h+\nabla^2V_h)\delta d_h+
\partial_{\dot d}g_h\,\delta\dot d_h\right].
\]

Term $(\delta M_h)f$ сохраняется. Его nonzero theta-row/c-column
contribution имеет вид

\[
2j_\parallel B_\theta^TW\operatorname{diag}[(1+c)\theta_{tt}]B_c.
\]

Прочие derivatives $g_h$ получаются аналитически из expressions выше;
finite differences не используются в production RHS Jacobian. Отдельный
directional finite-difference comparison проверяет этот Jacobian до
начала длинных trajectories.

## 5. Weak projection и дискретное тождество энергии

Action assembly сравнивается с проекцией audited continuum cubic
residuals на те же test functions. Polynomial residual evaluator
независим от сборки $g_h+\nabla V_h$ нового helper. После integrations
by parts endpoint terms исчезают, поскольку test functions равны нулю
на заделках. Strong PDE на endpoints не проверяется как machine-level
identity: clamp reactions там допустимы.

Дискретная энергия определяется тем же quartic action:

\[
E_h=\tfrac12\dot d_h^TM_h(d_h)\dot d_h+V_h,
\]
\[
\frac{dE_h}{dt}=\dot d_h^TM_h\ddot d_h+
\nabla V_h\cdot\dot d_h+
j_\parallel\int(1+c)c_t\theta_t^2\,ds=0.
\]

Последний contribution равен $\dot d_h^Tg_h$. Identity RHS
отдельна от фактического energy drift временного интегратора.
Exact untruncated energy не считается invariant cubic trajectory.
Energy fractions и modal classification не вычисляются.

## 6. Общая непрерывная начальная форма

Из preserved analytic Timoshenko fixed–fixed reference берётся первая
простая bending eigenpair; новые linear root searches не выполняются.
В данном control

\[
\omega_1=0.3174742907880648,\qquad
T_1=2\pi/\omega_1=19.791162590151373.
\]

Все coefficients analytic mode делятся на один signed peak $w(L/2)$.
Stationary peak и global maximum проверяются: $\max_s|\hat w|=1$;
$\hat\theta$ сохраняет тот же множитель. Поэтому $\hat w$
безразмерна, $\hat\theta$ имеет размерность $1/L$.

\[
w(s,0)=A\hat w(s),\quad\theta(s,0)=A\hat\theta(s),\quad
u(s,0)=c(s,0)=0,\qquad \dot d_h(0)=0.
\]

Две заранее заданные amplitudes: $A/h=0.05,0.025$, то есть
$A=0.0025,0.00125$. Одна и та же continuous pair проецируется на все
spatial levels. Projection error сохраняется отдельно. Начальные
условия не подвергаются static correction; $c=-\nu u_s$ не назначается,
быстрый переходный процесс не подавляется. Нулевые initial $u,c$ не
являются условием их дальнейшего сохранения.

## 7. Интервал времени и два linear controls

Для обеих amplitudes общий интервал
$0\le t\le5T_1=98.95581295075687$, observation time $\tau_t=t/T_1$.
Нелинейный период не подгоняется и не меняет reference $T_1$.

Continuous comparison — та же analytic pair, умноженная на
$A\cos(\omega_1t)$, с $u=c=0$. Отдельно полная finite-dimensional
linearized action решается точно по времени generalized
eigendecomposition $K_h\phi=\omega_h^2M_h(0)\phi$. Используются все
coordinates той же initial projection. Таким образом, пространственная
ошибка, time-integration error и nonlinear departure сравниваются
отдельно. Этот linear control не является nonlinear modal reduction.

До main trajectories проверяются первые три frequencies обоих local
blocks — M-H и Timoshenko. Нижний mixed sorted prefix не заменяет
проверку contraction/axial resolution и boundary layers.

## 8. Radau, tolerance scales и быстрые движения

Используется `scipy.integrate.Radau` — тот же неявный метод, который
вызывает `solve_ivp(method="Radau")`. Явный цикл по принятым шагам
позволяет проверять wall-time budget и сохранять доступную часть
траектории. Jacobian аналитический. Процесс запускается с одним BLAS thread; global user
environment и installed dependencies не изменяются. Характерный
contraction frequency и period:

\[
\omega_c=\sqrt{C/j_\parallel},\qquad T_c=2\pi/\omega_c.
\]

Уровни time accuracy зафиксированы до main runs:

| Level | rtol | Relative atol factor | max_step |
| --- | --- | --- | --- |
| coarse | $10^{-8}$ | $10^{-8}$ | $T_c/12$ |
| medium | $10^{-9}$ | $10^{-9}$ | $T_c/18$ |
| tight | $10^{-10}$ | $10^{-10}$ | $T_c/24$ |
| один разрешённый адресный extra | $2\cdot10^{-11}$ | $2\cdot10^{-11}$ | $T_c/36$ |

Atol задаётся в mass-scaled coordinates. Для $n=p-1$ coefficients
per field применяется fixed characteristic scale

\[
d=(A,A,A/L,A/L),\quad m_f=(m,m,j_\parallel,j_\parallel),\quad
s_f=d_f\sqrt{m_fL/n}.
\]

Coordinate atol равен $s_f$ times table factor; velocity atol имеет
дополнительный множитель $\omega_c$. Это analytical L2/mass scaling
per coordinate, а не relaxation по текущей величине field. Малые
возбуждённые $u,c$ не отбрасываются. `max_step` учитывает contraction
scale, а saved output grid дополнительно разрешает самый быстрый
retained semidiscrete linear period; одна bending frequency не служит
единственным критерием шага.

## 9. Заранее заданные convergence gates

Основная amplitude $A/h=.05$: spatial levels $p=16,24,32$ при tight
time accuracy; затем на $p=32$ сравниваются coarse/medium/tight.
Для $A/h=.025$ предусмотрены tight trajectories на $p=32$ и соседнем
$p=24$. Совпадающие runs переиспользуются.

Сравниваются восстановленные physical fields и velocities по spatial
L2 norms на одной time grid и общем physical interval. Phase alignment,
amplitude adjustment, period fitting и across-resolution modal tracking
не выполняются. Нельзя сравнивать коэффициенты разной размерности
напрямую или делить разность на мгновенно нулевое field.

Для каждого field $f$ используются оба отношения:
\[
e_{2,f}=\frac{\max_t\|f_a-f_b\|_{L^2}}{\max_t\|f_b\|_{L^2}},\qquad
e_{\infty,f}=\frac{\max_{t,s}|f_a-f_b|}{\max_{t,s}|f_b|}.
\]
Тот же контроль выполняется для скоростей; PASS требует обоих норм.
Ненулевая reference scale берётся по всему временному интервалу,
а не отдельно при каждом пересечении нуля. Spatial comparison использует
независимые100 Gauss points и ту же output time grid.

| Quantity | Criterion |
| --- | --- |
| $w,\theta$ и их velocities, последние разрешённые trajectories | Relative difference $\le10^{-4}$ |
| $u,c$ и их velocities | Relative difference $\le10^{-3}$ к собственному ненулевому характерному масштабу |
| Максимальный relative drift $E_h$ | $\le10^{-6}$ |
| Final spatial linear frequency comparison | Relative $\le2\cdot10^{-8}$ |
| Algebraic numerical identity | Scaled $\le2\cdot10^{-12}$ |
| Numerical floor qualification | Relative $10^{-10}$ |

Для малых $u,c$ сохраняются absolute differences и numerical-floor
qualification. Сходимость dominant $w$ не заменяет контроля остальных
fields. Допустимы только один additional $p=48$ и один addressed
extra time level в зафиксированном бюджете; automatic indefinite
refinement отсутствует. Tolerances после просмотра результатов не
ослабляются. Если nonlinear difference меньше numerical uncertainty,
используется `EFFECT_NOT_RESOLVED`.

## 10. Safety domain и вычислительный бюджет

Safety gates пилота не являются найденными physical applicability
thresholds. Вдоль trajectories проверяются

\[
\min(1+c)\ge0.9,\quad\max|c|\le0.1,\quad\max|\theta|\le0.1,
\]
\[
\max|u_s|\le0.1,\quad\max|w_s|\le0.1,\quad
L\max|\theta_s|\le0.2,
\]

и minimum eigenvalue относительной mass matrix через $M_h(0)$ не ниже
0.8. Наблюдаются retained axial/shear strain measures, contraction
gradient, energy, mass positivity и conditioning. Negative eigenvalues
не обнуляются; NaN/Inf и safety failure останавливают case с причиной.

Короткий medium smoke, $p=16$, physical duration 0.1, занял
0.2274216000 s при 575 RHS evaluations. Linear wall extrapolation
до $5T_1$ дала около 225 s; это cost estimate, не гарантия длительности
main integration. После smoke зафиксированы 2400 s total и 480 s per
case, максимум 10 cases. Основной план содержит семь trajectories;
последующие разрешённые refinement controls остаются внутри бюджета.
Каждый завершённый case сохраняется до следующего. Budget exhaustion
ведёт к `PARTIAL` с сохранёнными данными и явной причиной.

## 11. Artifacts, cache и figures

Новый workflow состоит из одного
[helper](../../scripts/lib/weakly_nonlinear_planar_dynamics.py), одного
[CLI](../../scripts/analysis/simulate_weakly_nonlinear_planar_rod.py),
[config](../../data/input/weakly_nonlinear_planar_time_pilot.json) и
[targeted tests](../../tests/test_weakly_nonlinear_planar_dynamics.py).
Старые audit/linear files читаются как protected dependencies.

Сохраняются continuous shape/normalization, initial projections,
trajectory coefficients и velocities, observation time series,
selected full spatial snapshots, convergence tables, integrator
settings/stats, mass/domain checks и manifest. Full field grids не
записываются на каждом внутреннем RHS call.

Observations включают $w(L/2)$, $\theta(L/4)$, $u(L/4)$, $c(L/4)$,
$c(L/2)$ и полные spatial norms; нуль в середине symmetric rod не
интерпретируется как отсутствие field. Не более трёх основных figures:
две amplitudes против linear transverse reference; induced $u,c$
с раздельными units; energy drift/convergence.

Cache identity включает accepted model/hash, geometry, basis/quadrature,
initial pair/amplitude, time horizon, integrator/Jacobian/tolerances и
library versions. Повторные complete compute/report/plot запускают
ноль integrations, ноль root searches и ноль symbolic derivations.
Сохранённое integration time отделяется от plot time.

Конфигурация этого workflow описывает только фиксированный G20 case:
coefficients читаются из protected accepted audit bundle, а shape — из
protected single-rod bundle. Изменение material_geometry не является
реализованным здесь способом построения нового физического case.
Timeout сохраняет accepted prefix; при safety/solver exception сохраняются
ранее законченные cases и failure record. Сохранение текущего prefix при
таком exception отдельно не реализовано; в данной серии exceptions не было.

## 12. Фактические результаты

Работа выполнена на local main HEAD
`7b3d667d5418cde247a6ca09fe9186945177537b` (Version0.5.15).
Pre-existing16 modified/8 untracked files завершённого NLSP audit
сохранены; staging не менялся. Использован existing interpreter
`D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe`, Python3.12.4,
NumPy2.1.3, SciPy1.15.2, Matplotlib3.9.2. SymPy отсутствует; новая
зависимость не устанавливалась. Точная алгебра использует уже проверенный
Fraction artifact, не новое ручное переписывание PDE.

### Линейные и вариационные controls

| p | max relative error MH first3 | max relative error Tim first3 | initial-pair max L2 projection error |
| --- | ---: | ---: | ---: |
|16|3.57936e-5|9.33e-13|1.23e-13|
|24|1.82394e-7|2.19e-12|3.26e-14|
|32|1.34866e-10|1.82e-12|2.72e-14|

Slow MH convergence at low p отражает finite contraction-clamp boundary
layers. Численная spatial tolerance здесь отдельна от прежних exact
algebra/source tolerances. M0 eigenvalues отличаются от1 не более8e-15.
Higher-quadrature energy absolute differences составляют не более7.42e-20
на заданных deterministic states; weak/action differences не более1.08e-14,
scaled energy RHS identity не более8.3e-17.

Actual linear time controls на $0\ldots0.1T_1$ сравниваются с full
semidiscrete exact-in-time solution: first bending relative q/velocity
differences3.25e-12/1.63e-11; pure MH axial1.76e-11/1.33e-11.
Zero RHS точно нулевой. Дополнительные unit checks проверяют variable
mass, оба g contributions и analytic Jacobian against finite differences.
Требуемый контроль первых трёх MH frequencies выполнен отдельно.

### Сохранённые nonlinear trajectories и сходимость

Bundle: `results/weakly_nonlinear_planar_time_pilot/c97287772bc461ef/`.

**NLSP_PLANAR_TIME_PILOT=PARTIAL.** Рассчитаны шесть полных траекторий,
включая обе основные amplitudes при p32 до5T1. Последний auxiliary
small-amplitude p24 case остановлен заранее заданным бюджетом при
`t=52.409610446 =2.648132T1`. Prefix сохранён и не заменяет полный5T1
convergence control. Snapshots unfinished case читаются по actual indices/
time: поздние requested targets повторяют последний доступный sample.
Дополнительные p48/time tightening не запускались в исчерпанном бюджете.

| p | A/h | time level | достигнуто t/T1 | integration seconds | RHS calls | max relative energy drift |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
|16|0.05|tight|5.000000|382.732|985510|9.514e-11|
|24|0.05|tight|5.000000|415.646|985510|9.516e-11|
|32|0.05|tight|5.000000|450.525|985510|9.532e-11|
|32|0.05|coarse|5.000000|141.878|311655|2.975e-8|
|32|0.05|medium|5.000000|252.998|554198|1.683e-9|
|32|0.025|tight|5.000000|380.139|828717|5.658e-11|
|24|0.025|tight|2.648132|184.983|438909|3.030e-11|

Полные cases имеют70985 common output samples, dt≈0.001394,
12 samples на период самой быстрой retained linear frequency p32
(omega≈375.5915). Output sampling не заменяет internal-step diagnostics.

### Пространственная сходимость, A/h=.05, tight time

| Компонента | p16→24 L2 relative | p24→32 L2 relative | p24→32 max relative | Критерий | last-pair status |
| --- | ---: | ---: | ---: | ---: | --- |
|u|2.93020e-2|2.73889e-3|3.88764e-3|1e-3|NOT CONVERGED|
|w|1.13846e-4|7.27865e-7|1.06111e-6|1e-4|PASS|
|theta|6.63444e-4|1.19402e-5|1.84444e-5|1e-4|PASS|
|c|4.07653e-1|3.05821e-2|4.57150e-2|1e-3|NOT CONVERGED|
|u_t|1.01767e-1|2.20027e-2|3.38514e-2|1e-3|NOT CONVERGED|
|w_t|4.62138e-3|4.11684e-5|7.71539e-5|1e-4|PASS|
|theta_t|2.70283e-2|1.43245e-3|2.23010e-3|1e-4|NOT CONVERGED|
|c_t|7.40547e-1|5.75033e-2|7.76910e-2|1e-3|NOT CONVERGED|

Последний spatial pair не разрешает u,c и скорости u,theta,c.
Абсолютные max-time L2 differences q(u,w,theta,c):
`[4.40364e-9,1.14890e-9,6.40449e-8,3.30107e-7]`; velocities:
`[2.28235e-7,2.06396e-8,2.44893e-6,2.43844e-5]`.
Ни одна reference scale не floor-limited. Independent exact-Legendre-Gram
reconstruction всех full histories подтверждает эти L2 comparisons.

### Временная сходимость, p32, A/h=.05

| Компонента | coarse→medium L2 relative | medium→tight L2 relative | medium→tight max relative | last-pair status |
| --- | ---: | ---: | ---: | --- |
|u|2.44728e-6|2.61365e-7|1.46370e-6|PASS|
|w|2.66783e-12|1.80827e-13|9.93876e-13|PASS|
|theta|7.39346e-11|4.10763e-12|9.75451e-12|PASS|
|c|5.24310e-5|2.96187e-6|3.44117e-6|PASS|
|u_t|1.32989e-4|1.42885e-5|7.63395e-5|PASS|
|w_t|1.81982e-9|1.47207e-10|7.67245e-10|PASS|
|theta_t|2.26209e-8|1.50527e-9|7.32953e-9|PASS|
|c_t|9.90741e-5|5.72118e-6|5.51278e-6|PASS|

Оба successive time pairs проходят оба norms. Max last-pair L2 relative
time difference1.43e-5, намного меньше spatial differences. Причина
unresolved быстрого c в этом контроле — пространственное разрешение.
Полный small-amplitude neighboring-p контроль не завершён.

### Отличие от линейного движения и масштабы полей

| A/h | max abs u | max abs w | max abs theta | max abs c | w departure / linear L2 scale | theta departure / linear L2 scale |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
|0.05|2.52881e-6|2.50000e-3|7.48109e-3|2.11817e-5|2.005311%|2.005384%|
|0.025|6.37158e-7|1.25000e-3|3.74055e-3|5.30483e-6|0.501714%|0.501874%|

Нормированные max-time L2 departures от semidiscrete linear reference
для q(u,w,theta,c)/A (L=1):
`[0.000643129,0.012661113,0.043026014,0.004317659]` при большой amplitude,
`[0.000323915,0.003167716,0.010767836,0.002158830]` при малой.
Continuous и semidiscrete references сравниваются отдельно. Уменьшение
наблюдается в обеих реально рассчитанных5T1 trajectories. Это descriptive
two-case result, не power-law fit, frequency/period identification или
доказательство заранее заданного отношения4/16.

Ненулевые u,c и отличие w от linear reference заметно больше соответствующих
last-pair absolute differences. Но все четыре distributed fields/velocities
не прошли convergence. Aggregate acceptance/effect flag остаётся
EFFECT_NOT_RESOLVED для полной four-field trajectory. Это не утверждение
«нелинейности нет»: small energy drift и w curves не дают общего PASS.

### Энергия, масса, domain и независимая проверка

У двух основных p32 trajectories min(1+c)=0.999978818/0.999994695.
Relative mass condition upper bounds1.000075626/1.000018988.
Во всей сохранённой серии weighted-Gram eigenvalue bounds:
`[0.999957637,1.000034930]`; sparse explicit generalized
eigensolves дают positive mass. Loewner bounds не заменяют mass matrix
в dynamics и не являются mass lumping. Max|theta|=.00748109,
max|c|=2.11817e-5, max|u_s|≤3.13e-5, max|w_s|≤.007656,
max L|theta_s|≤.06855. Safety neighborhood не нарушена; extrema относятся
к численным spatial samples. Retained strain/contraction-gradient measures
сохранены в case JSON.

Discrete E0 large/small p32:1.26041687e-9/3.14765e-10.
Max drift9.53238e-11/5.65809e-11; максимум всей серии2.97503e-8
при coarse time, ниже1e-6. Independent source T4/V4 integration и exact
Legendre Gram/unwhitening по7 snapshots всех шести полных cases
согласуются со stored energy не хуже4.7e-15 relative. Essential endpoints
точно нулевые; raw→physical snapshot differences≤5.3e-18.
Reader сохранён в targeted tests, не только в OS temporary script.

SciPy сообщает2 Jacobian updates и4 Newton LU на full case. Это standard
factorization reuse внутри adaptive Radau; nonlinear analytic callable
передан, constant Jacobian/M0 substitution не применялись. Force и variable
mass пересчитывались на RHS, ≈985507 mass factorizations в tight large case.

### Статусы и стоимость

| Категория | Статус |
| --- | --- |
|NLSP_PLANAR_DISCRETIZATION|PASS|
|NLSP_PLANAR_LINEAR_TIME_REFERENCE|PASS|
|NLSP_PLANAR_TIME_INTEGRATION|PARTIAL|
|NLSP_PLANAR_SPATIAL_CONVERGENCE|PARTIAL|
|NLSP_PLANAR_TEMPORAL_CONVERGENCE|PASS|
|NLSP_PLANAR_SMALL_AMPLITUDE_LIMIT|PASS|
|NLSP_PLANAR_ENERGY_AND_MASS|PARTIAL|
|NLSP_PLANAR_TIME_PILOT|PARTIAL|

TIME_INTEGRATION и ENERGY_AND_MASS — PARTIAL по coverage: шесть полных
cases безопасны, последний auxiliary имеет prefix. SMALL_AMPLITUDE_LIMIT
PASS означает observed decrease, не physical nonlinear validation.
Spatial и small-amplitude neighboring-p uncertainty остаются unresolved.

Серия: 2425.512s, включая assembly/postprocessing.
Интегрирование остановлено по2400s total scheduling deadline; обязательные
сохранение prefix и comparisons завершились≈25.5s после него.
480s per-case limit применяется к integration, не к postprocessing.
Бюджет/допуски не повышались. Семь attempts:6full+1partial;
5,090,009 RHS,13 Jacobian calls,26 Newton LU,5,089,987 variable-mass
factorizations. Root searches0; одна model derivation в main.
Matching compute/report/plot-only:0 integrations,0 roots,0 derivations.
247 unique targeted/relevant regression tests PASS; git diff --check PASS.
Audited files/reference bundles/staging сохранены; packages не ставились.

### Три figures

Основные p32 full5T1 histories; каждый figure сохранён PDF+PNG:

- `figures/transverse_motion.pdf` — normalized w и linear reference;
- `figures/generated_axial_contraction.pdf` — u(L/4),c(L/4), раздельные scales;
- `figures/energy_drift.pdf` — drift принятой quartic energy.

Графики — evidence pilot с общим PARTIAL. Partial p24 small history
не показана как полный5T1 result.

## 13. Воспроизведение и остановка

```powershell
python scripts/analysis/simulate_weakly_nonlinear_planar_rod.py --check
python scripts/analysis/simulate_weakly_nonlinear_planar_rod.py --smoke
python scripts/analysis/simulate_weakly_nonlinear_planar_rod.py --compute
python scripts/analysis/simulate_weakly_nonlinear_planar_rod.py --report-only results/weakly_nonlinear_planar_time_pilot/c97287772bc461ef
python scripts/analysis/simulate_weakly_nonlinear_planar_rod.py --plot-only results/weakly_nonlinear_planar_time_pilot/c97287772bc461ef
python -m pytest tests/test_weakly_nonlinear_planar_dynamics.py -q
```

`python` обозначает доступный project interpreter. `--report-only` и
`--plot-only` получают явный completed bundle path; окончательная команда
для текущего bundle приведена в результатах.

Это задача Коши, а не построенная периодическая орбита. Все четыре
fields развиваются независимо. Не вычисляются out-of-plane perturbations,
Floquet multipliers, critical amplitude, frequency/amplitude maps или
угловая конструкция. 3D FEM и experiment не используются как truth;
$V^{(0)}$ и coefficients не подбираются. По двум amplitudes не выводится
универсальный power law trajectory error; прошлый четвёртый порядок
full/cubic residual не является таким законом.

Плоский pilot не устанавливает пространственную устойчивость и не
снимает qualification `c1=c2` будущего reduced angular closure.
Следующий nonlinear этап не запускается автоматически после проверки.
