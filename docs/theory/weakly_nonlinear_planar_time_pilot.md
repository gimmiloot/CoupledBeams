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


## 14. Адресная диагностика и восстановление вычислительного пути, 2026-10-07

Это продолжение после Version 0.6.0, main HEAD
`1510d75c106a28a4da899c7eea1a337f11791ce6`. Первоначальные разделы 1–13,
исходный bundle `c97287772bc461ef`, его manifest и семь histories сохранены
неизменными. Новая диагностика не превращает предыдущий PARTIAL в COMPLETE.
V0, quartic action, cubic equations, variable kinetic energy, четыре поля,
Shen–Legendre basis, G20 geometry, initial pair и essential-only clamps не менялись.

### Исторические данные и нормы

Все artifact hashes старого manifest проверены на его собственных условиях;
новый code hash не предъявляется историческому execution identity.
Large-amplitude p16/p24/p32 tight и p32 medium/coarse содержат 70985 одинаковых
фактических времён: 0…98.9558129508 = 5T1, output dt≈0.001394058.
Small-amplitude p24 действительно заканчивается на 52.4096104458 = 2.648132T1;
повторённые последние snapshot indices означают одно и то же время, не 3/4/5T1.
Ни один старый case не пересчитывался.

Воспроизведены прежние max-over-time physical L2 и max-over-time/space нормы,
знаменатели, reporting floor и 100 Gauss comparison nodes. Фазового сдвига,
амплитудного масштабирования и исключения времени/краевых зон нет.
Профили на 1001 равномерной точке служат дополнительной локализацией и
не подменяют прежний max criterion. Обработка coefficients выполняется блоками.

| Компонента | p24→p32, relative L2 | relative max | Прежний допуск | Результат |
|---|---:|---:|---:|---|
| u | 2.738887e-3 | 3.887637e-3 | 1e-3 | FAIL |
| w | 7.278650e-7 | 1.061108e-6 | 1e-4 | PASS |
| theta | 1.194019e-5 | 1.844444e-5 | 1e-4 | PASS |
| c | 3.058208e-2 | 4.571499e-2 | 1e-3 | FAIL |
| u_t | 2.200272e-2 | 3.385145e-2 | 1e-3 | FAIL |
| w_t | 4.116842e-5 | 7.715395e-5 | 1e-4 | PASS |
| theta_t | 1.432452e-3 | 2.230099e-3 | 1e-4 | FAIL |
| c_t | 5.750326e-2 | 7.769096e-2 | 1e-3 | FAIL |

Все восемь исторических relative L2/max значений воспроизведены с разностью 0.
Контекст p16→p24 сохранён отдельно. Основные ограничения дают u,c, u_t,c_t
и theta_t. Превышение fixed L2 gates начинается уже при t/T1≈0.002465 для c_t,
0.005001 для u_t, 0.006762 для c, 0.026203 для u, 0.051490 для theta_t.
К 0.1T1 достигнуто около 43% полного максимума L2-разности c, 51% для u,
64% для u_t. Это не только накопление поздней ошибки.

Максимумы c и c_t по исходным Gauss nodes находятся при
(t/T1,s/L)≈(4.703525,0.507814) и (4.848628,0.972028), соответственно;
точные координаты/времена всех восьми максимумов сохранены в JSON.
Интеграл по времени квадрата разности распределён по всей длине:
71–88% приходится на [.1L,.9L], 12–29% суммарно на две крайние зоны.
При максимуме L2-разности u крайние зоны дают 52.9%.
Это не универсальная локализация всех разностей только у заделок.

### Физическая L2-проекция и разложение ошибки

Для каждого поля и скорости p32 физически восстановлен из собственной
mass-whitening transformation, затем ортогонально спроецирован на p24:

$$ f_{32}-f_{24}=(f_{32}-\Pi_{24}f_{32})+(\Pi_{24}f_{32}-f_{24}). $$

Gram matrices используют physical L2 weights, не необработанные whitened
coefficients. Точное полиномиальное интегрирование выполняется прежней
100-point Gauss quadrature; отдельно проверены ортогональность и Pythagoras.
Максимальные остатки, делённые на квадрат собственного characteristic L2
масштаба поля: 1.90e-17 для Pythagoras и 7.66e-18 для ортогональности.

| Компонента | Доля интеграла квадрата ошибки вне p24, % | В общем пространстве, % |
|---|---:|---:|
| u | 20.484 | 79.516 |
| w | 6.216 | 93.784 |
| theta | 2.545 | 97.455 |
| c | 9.584 | 90.416 |
| u_t | 34.934 | 65.066 |
| w_t | 17.285 | 82.715 |
| theta_t | 1.007 | 98.993 |
| c_t | 9.358 | 90.642 |

На первых 0.1T1 tail fractions u,c,u_t,c_t≈46%,52%,60%,52%.
Таким образом, есть и недостаток spatial detail, и разность эволюции общих
компонент. Второе слагаемое не доказано исключительно фазовой ошибкой.
Это L2 approximation diagnostic, не энергетическая классификация или modal
reduction. Представительные physical Legendre coefficients сохранены;
коэффициенты решателя не фильтровались и не отбрасывались.

### Непрерывная начальная совместность

Из защищённого `V_le4` независимо получен axial flux второго порядка
при u=c=0 и planar bending:

$$ F_u^{[2]}=(C-S)\theta w_s+(S-C/2)\theta^2. $$

Проверка точными Fraction polynomials совпадает с derivative action.
На fixed endpoint theta=0, поэтому

$$ (F_u^{[2]})_s=(C-S)\theta_s w_s,\qquad
   u_{tt}\big|_{\partial}=(C-S)\theta_s w_s/m. $$

Вычисление использует continuous analytic Timoshenko state basis и его
пространственные производные, а не Galerkin acceleration, которая по
конструкции равна нулю на концах. Residual convention: m*q_tt+r=0.
Четыре endpoint body residuals при точных clamp values:

| Поле | r при начальных данных у конца | Вывод |
|---|---|---|
| u | -(C-S)*theta_s*w_s | ненулевой trace порядка A² |
| w | S*(theta_s-w_ss) | 0 по linear eigenpair identity |
| theta | -S*w_s-Bp*theta_ss | 0 по linear eigenpair identity |
| c | 0 | 0 |

Для A=.0025: w_s=(+2.067790439e-4,-2.067790439e-4),
theta_s=(.06845755316,.06845755316), w_ss с теми же значениями,
theta_ss=(-.3181216061,+.3181216061). Axial acceleration traces:
(+1.101854330e-5,-1.101854330e-5). Половинная amplitude даёт четверть этих
значений; это проверка A² expression без новой малой-amplitude траектории.
Другие traces и linear controls равны нулю с roundoff порядка1e-15;
continuous state consistency max3.55e-15. Значения и скорости совместны.

Статус **CONFIRMED_LOW_ORDER_MISMATCH** относится к требованию гладкости
по времени до второго порядка вплоть до неподвижного конца. Он не означает
недопустимость weak initial-boundary problem, ошибку кода, отсутствие slope
constraints или необходимость исправить initial fields. Реакции заделки,
слабая постановка и предельное внутреннее strong equation различаются.
Совместность не является единственной доказанной причиной всех разностей;
более высокие compatibility orders здесь не исследованы. Static correction
или дополнительное w_s=0 не вводились.

### Profiling и минимальная оптимизация

Baseline helper сохранён до изменения отдельно, SHA256
`fa03369ea8aed8f6475679c319b4026eec169d7cad5ea39237f07354d6f37e90`,
совпадает с historical manifest. Выявлено вычисление potential Hessian при
обычном RHS. Теперь energy-only, gradient и Hessian paths разделены;
cache повышает уровень запроса для того же q и инвалидируется при его смене.
Неподвижные monomial selectors подготовлены один раз. Variable M(q), его
factorization при изменившемся c, inertial terms, analytic Jacobian, safety,
quadrature и Radau/tolerances сохранены.

16 разных реальных p32 states, warmup и медиана трёх повторов,
BLAS threads=1 во всех сравнениях:

| Операция | До, ms | После, ms | Ускорение |
|---|---:|---:|---:|
| RHS | .32557 | .23196 | 1.40× |
| potential gradient | .20018 | .12268 | 1.63× |
| potential energy | .17961 | .08318 | 2.16× |
| Hessian request | .59473 | .54012 | 1.10× |
| mass assembly | .03789 | .03713 | 1.02× |
| mass solve including assembly | .05712 | .05470 | 1.04× |
| inertial terms | .01764 | .01743 | 1.01× |
| safety | .04698 | .04578 | 1.03× |
| Jacobian | 1.12449 | 1.08351 | 1.04× |
| physical reconstruction | .01178 | .01167 | 1.01× |

Компоненты timing не складываются как независимые: многие paths включают
другие операции/cache. Отдельный32-state profiler даёт RHS1.32× и gradient1.81×;
это разброс измерений, не гарантированный factor полного расчёта.

На реальных snapshots и малых deterministic states gradient/Hessian/M/inertia,
acceleration/RHS/Jacobian/energy rate/weak residual до и после совпадают
побитово. В fresh16-state измерении V имеет absolute difference≤2.07e-24,
relative≤3.39e-15; неизменный equivalence gate2e-12 пройден. Weak/action
absolute residual≈1.37e-15 проходит2e-12. Отношение к уже компенсированной
итоговой силе≈1.11e-11 одинаково до/после и отдельно сохранено; tests также
проверяют relative identity на масштабе нескомпенсированных local work terms.
Cache controls включают energy→gradient→Hessian, repeat/new q и velocity-only
change. Не используется устаревшая mass factorization при изменённом c.

Старый full p32 tight:140786 accepted steps, dt_min≈.0005580,
dt_max≈.0007029 при max_step≈.0036047; ceiling шага фактически не активен.
nfev985510≈7.00006/step, njev2, Radau nlu4, variable-mass factorizations985507.
Полная последовательность внутренних dt и rejected steps не сохранялась:
median dt или точное число rejected steps не объявляются. Редкий refresh
Jacobian и Radau LU не равен редкому вычислению variable mass/RHS.

### Три коротких контроля и решение о p48

| Контроль 0…0.1T1 | Время, s | Steps | RHS | Jacobian / Radau LU | Drift |
|---|---:|---:|---:|---:|---:|
| old p32 tight | 9.02096 | 2816 | 19720 | 2 / 4 | 3.38052e-12 |
| new p32 tight | 7.13321 | 2816 | 19720 | 2 / 4 | 3.38063e-12 |
| p48 stricter | 13.81326 | 4211 | 29485 | 2 / 4 | 5.64770e-13 |

Old/new p32 на одной общей сетке имеют нулевые L2/max differences для всех
четырёх полей и скоростей; integration speedup1.2646×. Roundoff difference
energy drift не меняет energy gate1e-6. P48 имеет97 quadrature nodes,
исходную форму/амплитуду и stricter `allowed_extra` prescription, не новую
физическую задачу. Safety/mass controls пройдены на этих коротких интервалах.

| Компонента | p32 tight→p48 stricter, short relative L2 | short relative max |
|---|---:|---:|
| u | 2.723316e-3 | 3.024127e-3 |
| w | 1.719204e-7 | 2.248942e-7 |
| theta | 1.094882e-6 | 1.872709e-6 |
| c | 2.731908e-2 | 3.209896e-2 |
| u_t | 3.140829e-2 | 3.654294e-2 |
| w_t | 3.726305e-5 | 4.950172e-5 |
| theta_t | 2.441589e-4 | 4.181399e-4 |
| c_t | 5.168276e-2 | 6.078865e-2 |

Это **SHORT_ONLY** на0…0.1T1: разные p и time prescriptions, отдельные
short-horizon denominators. Таблица не заменяет full spatial comparison и
не доказывает p48 temporal convergence; её нельзя непосредственно трактовать
как улучшение/ухудшение относительно full5T1 таблицы. Полного p48 нет.

Заранее установлены900s profiling+integration, максимум3 short и1 full.
Учтено40.41217s, включая10s консервативного учёта первоначального profiler.
Full forecast=13.8132551×50×1.25=863.32844s превышает остаток859.58783s.
Поэтому **REFINEMENT_DEFERRED_BY_BUDGET**: margin/budget не изменялись,
единственный кандидат full p48 не запускался. Наличие mismatch само по себе
не использовалось как запрет p48. Дополнительный p48 tight/strict temporal
pair отсутствует; smaller-amplitude соседний p-контроль остаётся прежним PARTIAL.

Для p48 fastest retained linear omega≈818.6625 старые output times дают
около5.50 samples/period. Предусмотрена union точных старых times и fine grid
с12 samples/period (dt≈.000639576), без интерполяции исторических данных.
Эта сетка использована в shorts; полная union из225706 samples только
подготовлена, full trajectory на ней не получена.

### Bundle, воспроизведение и фактическая остановка

Final diagnostic bundle:
`results/weakly_nonlinear_planar_recovery/054874a4a4c9c9ff/`.
Исходная новая numerical execution сохранена как `db47d4efb6941bed`.
Final cache содержит те же данные и explicit execution identity, archived
executed CLI, exact code diff и `cache_revision.json`: исправлены три plot
lookup keys и две safety/partial-timestamp проверки неисполненного full path.
Это post-execution revision с0 новых integrations, не повтор старой серии.
Historical `c97287772bc461ef` не изменён. Cache identity содержит оба helper
hashes, CLI, frozen action, исходный config и dependency versions.

Сохранены manifests/hash checks, обе пары старых spatial diagnostics,
нормы/зоны/profiles, projection/Pythagoras, continuous compatibility,
component timings/equivalence,3short histories, cost decision и SHORT_ONLY
comparison. Две диагностические figures, PDF+PNG:
`figures/contraction_projection` и `figures/contraction_localization`.
Они показывают физические differences, не energy classes.

```powershell
python scripts/analysis/diagnose_weakly_nonlinear_planar_rod.py --report-only results/weakly_nonlinear_planar_recovery/054874a4a4c9c9ff
python scripts/analysis/diagnose_weakly_nonlinear_planar_rod.py --plot-only results/weakly_nonlinear_planar_recovery/054874a4a4c9c9ff
```

Matching `--compute` reuse, report-only и plot-only проверены с0 integrations,
0 profiling,0 roots,0 symbolic derivations. `--diagnose` читает старые данные
и выполняет algebra/compatibility diagnostic, без ODE. Первый uncached
`--compute` требует preserved `--baseline-helper`; он является отдельной
bounded numerical программой, не способом просто перерисовать graphs.
Исходный старый CLI после изменения numerical hash может иметь другой cache
identity: повтор его `--compute` в этом задании не выполнялся.

| Статус | Итог |
|---|---|
| NLSP_PLANAR_ERROR_DIAGNOSTIC | PASS: диагностика выполнена; spatial differences сохраняются |
| NLSP_PLANAR_INITIAL_COMPATIBILITY | CONFIRMED_LOW_ORDER_MISMATCH |
| NLSP_PLANAR_RHS_EQUIVALENCE | PASS |
| NLSP_PLANAR_PERFORMANCE | PASS: умеренный измеренный выигрыш |
| NLSP_PLANAR_P48_SPATIAL_CHECK | REFINEMENT_DEFERRED_BY_BUDGET |
| NLSP_PLANAR_SOLVER_RECOVERY | PARTIAL |

Небольшая оптимизация улучшила стоимость вычислений без изменения trajectory,
но не установила full four-field spatial convergence. Дальнейшее решение
требует отдельного выбора пользователя; замена basis, FEM/local grid, modal
reduction, initial correction, новая amplitude/angle/out-of-plane/Floquet
программа не начинались. LONG CLOSED, EB/RLB-KV PAUSED и angular out-of-plane
same-clamp reference UNAVAILABLE сохранены. История memory дополнена
[NLSP-D03](../memory/decisions.md#nlsp-d03)/[NLSP-K03](../memory/knowledge.md#nlsp-k03),
старые D02/K02 не переписаны.


### Проверки продолжения

31 новых targeted tests входят в165 прошедших combined checks для planar helper,
protected spatial action, finite M-H/Timoshenko rod и Bishop literature. Три
старых tests с solve_ivp намеренно deselected: они добавили бы скрытые
integrations сверх разрешённых3short controls. Ещё3 targeted rectangular
Timoshenko coefficient/convention checks прошли; итог168 уникальных PASS.
Новая CLI cache/report/plot проверяется с запрещёнными ODE/profiling/derivation
paths; full-path fail-gate и partial actual-time metadata проверены без ODE.
README/CHANGELOG/navigation обновлены; assumptions/equations/source indices
не менялись, поскольку новые physical assumptions/formulas не вводились.


Дополнительная обработка трёх сохранённых short histories не интегрирует ODE:
`short_safety.json` подтверждает min(1+c)≥.9999827897, relative mass eigenvalue
lower bound≥.9999655797 и condition upper bound≤1.000052435. Sparse explicit
mass eigensolves и исходные small-neighborhood limits пройдены. Эти bounds
относятся только к рассчитанному0…0.1T1, не к отсутствующей full p48 history.
