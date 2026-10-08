# Подготовленное начальное состояние прямого плоского стержня

Дата: 2026-10-08. Этот этап проверяет постоянную и гармоническую части
ведущего продольного отклика, затем определяет **отдельный новый initial
case**. Принятая физическая модель, её энергия и внешние заделки сохраняются.
Историческая задача с initial `u=c=0` и её статусы PARTIAL не изменяются.

Итоговая остановка: periodic profiles разрешены, но ни одна из двух
разрешённых nonlinear spatial pairs не представляет общий initial state
с заданной точностью производных. Новых ODE-интегрирований — **0**.
`NLSP_PREPARED_INITIAL_STATE_PILOT=PARTIAL`; SHORT PASS не заявлен.

## 1. Источники и границы

Основание постановки — [проверенная пространственная модель](weakly_nonlinear_spatial_rod.md),
[generated cubic expansion](weakly_nonlinear_spatial_rod_expansion_generated.md),
[первый planar pilot и recovery](weakly_nonlinear_planar_time_pilot.md),
[exact-time second-order diagnostic](planar_second_order_axial_response.md).
Новая начальная коррекция выводится из нашего принятого действия, а не
приписывается Crespo da Silva или публикации Jang.

Переиспользованы сохранённые spectral representations и historical bundles:

- `results/planar_second_order_axial_response/b3ea4eb6ac95d6e1/`;
- `results/weakly_nonlinear_planar_time_pilot/c97287772bc461ef/`;
- `results/weakly_nonlinear_planar_recovery/054874a4a4c9c9ff/`.

Проверены их собственные manifests/hashes и actual timestamps. Изменение
текущего code hash не требует пересчёта или перезаписи старого bundle.
Сохранённые M,K,V,omega,b0,b2 используются целиком: новые eigenvalues,
модальная фильтрация и миллионные exact-time histories здесь не нужны.

Это diagnostic continuation NLSP. LONG остаётся CLOSED в принятом 1D scope;
EB/RLB-KV — PAUSED_FOR_SUPERVISOR_DIRECTION; angular out-of-plane same-clamp
reference — UNAVAILABLE. Углового узла, новой продольной теории и переноса
spring/KV в эту модель в задаче нет.

## 2. Неизменная задача и два разных cases

Один прямой homogeneous fixed–fixed G20 rod:

| Параметр | Значение |
| --- | --- |
| E,rho | 1,1 |
| nu | 0.3 |
| b,h0,L | 0.20,0.05,1 |
| kappa | 5/6, прежнее project value |
| Поля | независимые u,w,theta,c |
| Essential BC | u=w=theta=c=0 на обоих концах |

Дополнительные slope constraints, book_slope_clamp и изменение c-clamp не
вводятся. Прежние quartic action, cubic residuals, variable mass и inertial
terms заморожены. Shen–Legendre basis остаётся
$B_n=P_n-P_{n+2}$, $\xi=2s/L-1$, n=0,...,p-2.

| Case | Initial preparation | Статус истории |
| --- | --- | --- |
| `linear_bending_zero_axial_ic` | initial u=c=0, линейная bending eigenpair | historical PARTIAL сохранён |
| `prepared_axial_O2_with_cubic_endpoint_compatibility` | общий O2 axial state и однозначно заданная theta3 | отдельный новый candidate; short runs NOT_RUN |

Используются прежние continuous shape и один общий множитель её нормировки:

\[
\omega_1=0.3174742907880648,\qquad T_1=19.791162590151373,
\]
\[
\epsilon_a=A/h_0,\qquad W=h_0\hat w,\qquad\Theta=h_0\hat\theta.
\]

Для разрешённого нового динамического case epsilon_a=.05, A=.0025.
Значения .025 и .0125 используются только в дешёвом algebraic residual audit,
без новых trajectories. Omega_d=2*omega1 одинаково для всех p.

## 3. Периодическая и свободная части ведущего отклика

Для уже проверенной линейной forced системы

\[
M\ddot x+Kx=f_0+f_2\cos(\Omega_dt)
\]

сохранённые mass-normalized eigenvectors дают

\[
x_{stat}=V(b_0/\omega^2),\qquad
x_{harm}=V\left(b_2/(\omega^2-\Omega_d^2)\right),
\]
\[
x_{per}(t)=x_{stat}+x_{harm}\cos(\Omega_dt).
\]

Деление в этих выражениях покомпонентное. Все 2(p-1) spatial coordinates
сохранены; малые coefficients и высокие frequencies не отбрасываются.

Независимые прямые решения проверяют

\[
Kx_{stat}=f_0,\qquad(K-\Omega_d^2M)x_{harm}=f_2.
\]

Минимальный detuning около 2.51108; near-singular preparation не обнаружена.
Condition number динамического оператора растёт от 1.14e3 при p16 до 1.06e6
при p96. Максимальная spectral/direct coordinate difference для выбранных
p32/p96 — 8.50e-12 относительно direct solution; scaled linear-solve
residual не более 1.58e-12. Вязкость и regularization не вводятся.

Прежняя zero-IC задача содержит необходимую свободную часть:

\[
x_{free}(t)=-V\operatorname{diag}(\cos\omega t)
\left(b_0/\omega^2+b_2/(\omega^2-\Omega_d^2)\right).
\]

Проверено $x_{old}=x_{per}+x_{free}$ в нескольких фиксированных временах,
включая values, velocities и accelerations. Свободная часть обеспечивает
старые IC, является частью физического конечномерного отклика и без damping
не обязана затухать. Из старых nonlinear histories она не удалялась.

## 4. Политика profile preparation

Политика зафиксирована до просмотра результатов в
[config](../../data/input/planar_prepared_initial_state.json): цель 1e-6
для L2/max profiles, производных d=0,1,2 и отдельного endpoint audit.

Для каждого field одновременно сохраняются собственные ненулевые масштабы
high-p profile и фиксированные размерные масштабы:

\[
u^{(d)}:\ h_0/L^d,\qquad c^{(d)}:\ 1/L^d.
\]

L2 scale дополнительно умножается на sqrt(L). Numerical floor применяется
согласно config; крайние зоны и производные не исключаются из gate.
L2 вычисляется по physical Legendre Gram; maxima относятся к заданным
spatial samples, а не к доказанному continuous supremum. Endpoint jets
проверяются отдельно и получают те же derivative scales.

Вторая производная вычисляется дифференцированием полинома. Она не
восстанавливается из проверяемой PDE. Это отделяет качество profile
representation от формального уравнения continuous boundary-value problem.

## 5. Сходимость stat/harm profiles

Из existing p16,24,32,48,64,96 восстановлены два physical profiles для u,c,
их первые/вторые derivatives и endpoint jets. Никаких новых p нет.

Последняя пара p64→p96 проходит fixed 1e-6 policy для обоих profiles:

| Проверка последней пары | Максимальное значение |
| --- | ---: |
| own-profile relative max, все derivatives | 1.98e-8 |
| fixed-scaled max, все derivatives | 3.35e-8 |
| own-profile relative max, d=0 | около 7.14e-12 |
| p96 O2 axial endpoint cancellation, scaled | 2.94e-8 |
| p96 O2 contraction endpoint cancellation, scaled | 3.58e-12 |

P48→p64 ещё не проходит фиксированный derivative gate c_ss, хотя profile
values уже очень близки. Разрешён один общий reference p96 из проверенной
последней пары. Это численное приближение, не continuum truth и не доказательство
сходимости прежнего полного time history.

У periodic component допустимы только верхние оценки

\[
\sup_t\|\Delta x_{per}\|\le\|\Delta x_{stat}\|+\|\Delta x_{harm}\|,
\qquad
\sup_t\|\Delta\dot x_{per}\|\le\Omega_d\|\Delta x_{harm}\|.
\]

Они применяются покомпонентно и не называются точными maxima.

## 6. Один общий initial evaluator

Пусть $U_\star,C_\star$ — physical profiles суммы stat+harm в общем reference.
Здесь C_star — contraction field; C без индекса остаётся constitutive
coefficient EA/(1-nu²).

\[
u_0=\epsilon_a^2U_\star,\quad c_0=\epsilon_a^2C_\star,\quad
w_0=\epsilon_aW,\quad
\theta_0=\epsilon_a\Theta+\epsilon_a^3\Theta_3.
\]

Все initial velocities равны нулю. Для каждого spatial p проецируется этот
же evaluator; profiles и Theta3 для него заново не подбираются.

Stat и harm сохранены вместе. Независимые u,c не удерживаются на x_per во
времени; c=-nu*u_s, theta=w_s и inextensibility не вводятся. External forcing
не добавляется в полную автономную cubic RHS. Initial energy нового case
может отличаться от старого: amplitude w сохраняется без energy matching.

## 7. Независимая endpoint compatibility

Из защищённого cubic action при нулевых essential endpoint values и initial
velocities получены четыре traces внутренних сильных уравнений:

\[
m a_u=C u_{ss}+\nu C c_s+(C-S)\theta_s w_s,
\]
\[
m a_w=S(w_{ss}-\theta_s)+(C-S)u_s\theta_s,
\]
\[
j_\parallel a_\theta=B_\parallel\theta_{ss}+S w_s-(C-S)u_s w_s,
\]
\[
j_\parallel a_c=Hc_{ss}-\nu C u_s.
\]

Это internal strong traces и условия гладкой initial compatibility, а не
замена численной Dirichlet row сильной PDE. Reactions заделок не зануляются.
Нулевой Galerkin acceleration endpoint сам по себе не доказывает compatibility.

До theta3 candidate имеет O2 axial/contraction cancellations

\[
CU_{\star,ss}+\nu CC_{\star,s}+(C-S)\Theta_sW_s=0,
\qquad HC_{\star,ss}-\nu CU_{\star,s}=0.
\]

Формальные continuous identities и численные traces independently
differentiated reference представлены раздельно. Остающиеся cubic bending
terms требуют проверки обоих уравнений w,theta; проверки одного axial trace
недостаточно.

## 8. Заданная quintic theta3

Принято w3=0. Из endpoint residuals следует правило

\[
\Theta_3=0,\quad
\Theta_{3,s}=(C-S)U_{\star,s}\Theta_s/S,\quad
\Theta_{3,ss}=(C-S)U_{\star,s}W_s/B_\parallel
\]

на обоих концах. Шесть условий определяют polynomial degree ≤5. Для
локальной координаты eta=s/L преобразование jets имеет множители L и L²; numerical
coefficient fitting по trajectories не применяется.

Для G20 получены коэффициенты в ascending powers eta:

\[
[0,\ 0.018067372420774844,\ 0.04197945490684334,\ -0.34859154383707897,\ 0.48090786084944404,\ -0.1923631443399832].
\]

Determinant малого Hermite system равен 4. First derivative на обоих концах
около 0.01806737242; second derivatives имеют противоположные знаки и модуль
около 0.0839589098. Endpoint values не более 6.25e-17; jets проверены независимо.
Scaled coefficient reflection-parity defect 4.27e-12.

Это одно явно выбранное правило минимальной подготовки, не единственная
физически возможная correction и не построение nonlinear normal mode.
Степень и свободные coefficients не увеличиваются для улучшения результата.

## 9. Through cubic order и finite amplitude

Проверяются endpoint values, нулевые velocities и coefficients ускорений
при epsilon_a¹, epsilon_a², epsilon_a³ для всех четырёх fields. Formal
cancellation с defining profiles и theta3 не означает exact finite-amplitude
compatibility.

После подстановки theta=epsilon_a*Theta+epsilon_a³*Theta3 в cubic PDE могут
появляться степени выше3. Они сохраняются при numerical evaluation. У разных
уравнений не требуется одинаковый ненулевой ведущий порядок остатка.

При epsilon_a=.05 фактический u acceleration trace порядка ±3.596e-10,
w trace порядка 2.388e-11, theta trace порядка ±4.05e-12, c trace порядка4.7e-11. В эти значения входят
реальные higher-order terms и numerical representation residuals; полная
таблица всех четырёх traces и трёх algebraic amplitudes сохраняется в bundle.
Значения не объявлены exact zero.
В условном exact endpoint audit после defining BVP и Hermite cancellations
остаются numerator terms epsilon_a^4*(C-S)*Theta3_s*W_s у u и
epsilon_a^5*(C-S)*U_star_s*Theta3_s у w. У theta,c в этом специальном
идеальном endpoint ansatz ненулевой higher-order coefficient не требуется;
численные профили сохраняют свои измеренные representation residuals.


Статус `ENDPOINT_ACCELERATION_COMPATIBLE_THROUGH_CUBIC_ORDER` относится
только к указанным initial acceleration coefficients. Infinite smoothness,
точная periodic orbit и physical validation не доказаны.

## 10. Projection admission и остановка до ODE

L2 projection одного общего p96 state, без theta3 fitting отдельно по p.
Нужно различать первоначальный two-field O2 audit и окончательный full u,w,theta,c audit.
У полного физического initial state используются прежние characteristic scales
[A,A,A/L,A/L] для [u,w,theta,c], а d-я derivative добавляет множитель L^(-d).
Проверяются relative L2/max и independently differentiated endpoint error
на этих фиксированных scales. Coefficient cancellations используют отдельные
force scales из config. Эти проверки не взаимозаменяемы:

| p | O2 u/c-only audit | Полный four-field projection audit | Admission |
| --- | --- | --- | --- |
|32|c и derivatives не разрешены; axial trace8.62e-3|values/derivatives и coefficient gate не проходят|PARTIAL|
|48|c_ss endpoint error1.94e-5; axial trace1.77e-6|c_ss derivative и coefficient gate не проходят|PARTIAL|
|64|u/c-only jets проходят, worst error5.69e-8|theta_ss endpoint error3.02565e-6 при1e-6 gate|PARTIAL|

Theta3 quintic принадлежит пространству каждого из этих p. Это не гарантирует
точность second derivative проекции continuous Theta. При p64 compatibility
coefficient gate сам по себе проходит (6.79e-7), но theta_ss projection gate
не проходит. Ни эта ошибка, ни O2 u/c ошибки не доказывают invalid weak IVP.
В основной паре p32/p48 и в разрешённой замене p48/p64 не допущены оба участника. Поэтому nonlinear short comparison
не запускается. Не выбирается новая пара p64/p96 и не ослабляется1e-6 policy.

Малый L2 error values p48 не заменяет точность second derivatives/endpoints.
Новых nonlinear integrations —0; spatial/temporal checks нового case NOT_RUN.
Эта остановка определяется quality gate, а не исчерпанием бюджета900s.


### Initial energy и safety

| p | Новая initial quartic energy |
| --- | ---: |
|32|1.2597383472258475e-9|
|48|1.2597383472257464e-9|
|64|1.2597383472257420e-9|

На проверяемых projection states min(1+c) не ниже0.99999539776,
max|c| не выше4.603e-6, max|theta| не выше0.007484,
L*max|theta_s| не выше0.06844. Weighted-Gram Loewner lower bound
relative mass не ниже0.9999907955, condition bound не выше1.000009205.
Энергия и положительность массы проходят эти initial checks, но не заменяют
failed projection gate. Energy drift и temporal uncertainty нового case
не измерены, поскольку новые trajectories отсутствуют.
## 11. Сохранённый OLD short control

Historical recovery cases new_p32 и p48_strict_short имеют4515 одинаковых
actual timestamps на0...1.9791162590151374=0.1T1. P32 использовал tight,
rtol1e-10; p48 — allowed_extra, rtol2e-11. Это qualification: comparison
не является чистым spatial test при одинаковой time prescription.

| Field | OLD absolute L2 | OLD relative L2 | OLD relative max | NEW |
| --- | ---: | ---: | ---: | --- |
|u|4.308e-9|2.723e-3|3.024e-3|NOT_RUN|
|w|2.714e-10|1.719e-7|2.249e-7|NOT_RUN|
|theta|5.873e-9|1.095e-6|1.873e-6|NOT_RUN|
|c|2.950e-7|2.732e-2|3.210e-2|NOT_RUN|
|u_t|3.240e-7|3.141e-2|3.654e-2|NOT_RUN|
|w_t|1.099e-8|3.726e-5|4.950e-5|NOT_RUN|
|theta_t|2.444e-7|2.442e-4|4.181e-4|NOT_RUN|
|c_t|2.185e-5|5.168e-2|6.079e-2|NOT_RUN|

Old histories не пересчитаны. Повторённые snapshots незаконченного small-p24
case не трактуются как более поздние данные. Own-field denominators,
absolute differences и fixed-dimensional scales сохраняются отдельно.
Нельзя заявить improvement только из-за изменения знаменателя нового case.

## 12. Воспроизведение и ограничение claims

Алгебра подготовки реализована в [planar_prepared_initial_state.py](../../scripts/lib/planar_prepared_initial_state.py);
endpoint identities сравниваются с обоими проверенными вариантами frozen residuals.

Focused entry point — [prepare_planar_initial_state.py](../../scripts/analysis/prepare_planar_initial_state.py),
explicit configuration — [planar_prepared_initial_state.json](../../data/input/planar_prepared_initial_state.json).
Он переиспользует protected action, spectral archives и basis machinery.
Итоговый preparation bundle:
`results/planar_prepared_initial_state/5ea8d41faf8ede54/`.
Он содержит manifest/hashes, профили, endpoint audit, initial projection и
actual остановку. Final code/cache identity и execution snapshots входят в manifest. Два ранних
preparation executions f6c06a90f979185c/b392d104d3391942 сохранены отдельно;
они заняли4.087/3.853s. Финальный numerical audit3.857s:6 spectral restores,
4 direct linear validation solves,0 новых M-H/Timoshenko eigensolves и0 ODE.
Исторические nonlinear и million-point analytical histories не воспроизводились.
Независимые read-only audits и validation runs имеют отдельные counters;
они также не интегрируют ODE.

Matching cache/report/plot выполняют0 BVP/eigensolves,0 новых analytic histories
и0 ODE integrations. Ненулевое количество matrix direct checks
для первого audit различается с eigensolves и time integration.

Не объявлено, что подготовка улучшила full nonlinear convergence: new
trajectories отсутствуют. Старая free component физическая, а её отделение
само по себе не доказывает единственную причину прежних spatial differences.
Новый candidate одновременно меняет leading free part и initial compatibility;
возможный последующий эффект также нельзя будет приписывать лишь одному
механизму без отдельного controlled comparison.

Historical NLSP_PLANAR_TIME_PILOT/PARTIAL и SOLVER_RECOVERY/PARTIAL сохранены.
Нет новых amplitudes trajectories, full5T1, basis changes, joint conditions,
Floquet, out-of-plane perturbations, FEM или thresholds. Автоматический переход
к следующему этапу не разрешён.

## 13. Раздельные итоговые статусы

| Статус | Итог |
| --- | --- |
|NLSP_PERIODIC_AXIAL_PROFILES|PASS|
|NLSP_PERIODIC_PROFILE_CONVERGENCE|PASS|
|NLSP_PREPARED_INITIAL_STATE|PASS — общий candidate построен|
|NLSP_INITIAL_COMPATIBILITY_THROUGH_CUBIC_ORDER|PASS — formal и numeric reference checks|
|NLSP_COMMON_INITIAL_PROJECTION|PARTIAL|
|NLSP_PREPARED_SHORT_TEMPORAL_CHECK|NOT_RUN|
|NLSP_PREPARED_SHORT_SPATIAL_CHECK|NOT_RUN|
|NLSP_PREPARED_INITIAL_STATE_PILOT|PARTIAL|

Команды из корня проекта:

```powershell
python scripts/analysis/prepare_planar_initial_state.py --compute
python scripts/analysis/prepare_planar_initial_state.py --report-only <bundle>
python scripts/analysis/prepare_planar_initial_state.py --plot-only <bundle>
```

`<bundle>` — fingerprint path, напечатанный CLI; для завершённого этапа
используется `results/planar_prepared_initial_state/5ea8d41faf8ede54`. Targeted
checks находятся в [test_planar_prepared_initial_state.py](../../tests/test_planar_prepared_initial_state.py).
При matching cache никакой повторный preparation/integration не выполняется.

## 14. Дополнительная numerical action check

На одном общем p96 initial evaluator отдельно проверены floating-point
strong/weak action residual и identity энергии без временного интегрирования.
Для synthetic test velocity (не actual IC) и projections p32/48/64 относительные
невязки на прежнем uncancelled-work scale составляют6.83e-13/1.14e-11/4.31e-11.
При reference-only p96 —3.37e-10. Абсолютный gate2e-12 проходит для всех;
относительный gate2e-12 проходит толькоp32 из этих уровней. Identity энергии
проходит (scaled power ниже5e-18). Сами точные A/B coefficient identities
остаются нулевыми; frozen action/RHS не меняются.

Это дополнительное численное ограничение сохранено в
`auxiliary_weak_identity.json`; не объявляется ошибкой физической модели
или единственной причиной spatial difficulty. Relative gate не ослаблен.
P48/p64 checks обозначены strict XFAIL как нерешённый numerical identity
control; они не становятся PASS. Nonlinear p96 trajectory запрещена и не
запускалась. Обе разрешённые пары уже блокируются initial projection независимо
от этого дополнительного check. Переход к дальнейшей numerical strategy
требует отдельного решения; её настоящий этап не выбирает.

## 15. Verification и preservation

Accepted targeted run:190 PASS,2 strict XFAIL,8 intentionally deselected
за12.32s. Новый test module:83 PASS и те же2 strict XFAIL. Проверены
source Rucka/Jang convention (включая «Jang kappa is never guessed»),
Bishop и rectangular Timoshenko contracts, old planar action/RHS/cache,
recovery projection и28 spatial-action checks. Frozen tests не редактировались.

Три старых integration tests, три frequency-convergence cases, один mixed
runtime/eigen control и один mass-eigenvalue check не запускались в этой
scope. New mass positivity проверена независимыми weighted-Gram bounds;
model eigensolves и ODE routes запрещены regression guards. Стандартное
внутреннее построение Gauss–Legendre nodes этим запретом не подменяется.
Первый слишком широкий validation guard блокировал именно квадратуру;
его failed attempt сохранён отдельно. Mathematical tolerances не менялись.

Новый bundle сохраняет также auxiliary strong/weak data, независимые read-only
audits, execution code/config/test snapshots, OLD/NEW CSV и результаты
cache/preservation/link checks. Report/plot и повторный matching compute
проверены без новых BVP, eigen, analytic-history или time integrations.
Большие исторические bundles, old source inputs, source index/BibTeX,
model/RHS/Jacobian/generated expansion и baseline equations не изменены.
README/CHANGELOG и соответствующая navigation обновлены; добавлены только
NLSP-D05/K05, прежние D/K prefixes остаются неизменными. На этом остановка.


<a id="prepared-precision-feasibility"></a>
## 16. Precision и bounded feasibility continuation, 2026-10-08

Этот раздел фиксирует отдельное решение пользователя после исторического
preparation gate выше: сначала ограниченно проверить практическую вычислимость
того же подготовленного состояния. Исторические §13–15 и bundle `5ea8d41faf8ede54`
не переоценены задним числом. Строгие thresholds, нормы, floors и default admission
сохранены. `EXPLORATORY_NOT_CERTIFIED` обозначает отдельно разрешённое выполнение,
а не изменение `PreparedInitialState.admitted` или объявление numerical certification.

Получены **три полные короткие траектории** до `0.1T1`; `p48→p64` даёт семь
проходящих компонентов из восьми. Временной контроль проходит. Это практическое
подтверждение вычислимости данной конечномерной задачи в заявленном коротком scope.
Требуемая строгая точность пока не установлена: float64 relative strong/weak FAIL
и max-критерий пространственного сравнения `theta_t` FAIL остаются видимыми.

### 16.1. Frozen inputs и ограниченная arithmetic diagnosis

Повторно проверены **собственные** manifests четырёх immutable historical sources:
prepared `5ea8d41faf8ede54`, second-order `b3ea4eb6ac95d6e1`, recovery
`054874a4a4c9c9ff`, pilot `c97287772bc461ef`. Старые BVP, eigensystems,
Theta3 и trajectories не пересчитывались из-за нового code hash.

Один target сохраняет `epsilon_a=.05`, `A=.0025`, numerical p96 U_star/C_star,
сохранённую непрерывную analytic W/Theta pair и прежний quintic Theta3. Входные
binary64 coefficients трактуются как фиксированные числа; более точные операции
над ними не превращают saved numerical reference в точное континуальное решение.
Параметры G20, V0, quartic action, четыре поля, cubic PDE, заделки, variable mass,
все inertial terms и полный Shen trial/test space неизменны.

`mpmath` уже установлен; зависимости не устанавливались. В данной Windows-среде
`np.longdouble` имеет ту же epsilon, что float64. Сверка выполнялась при40/70
decimal digits для проекции и45/70 для independent strong/weak control.

**Проекция.** Exact physical Gram и analytic Legendre moments отделяют ошибку
интегралов/решения коэффициентов от конечного p. Сохранённые U/C coefficients
используются полностью, Theta3 включается напрямую через degree5 moments.
Последующее raw↔whitened преобразование и восстановление проверены уже в том
float64 пути, который использует временной решатель.

* p48: ошибка `c_ss` у unconstrained L2-проекции сохраняется после MP пересчёта
  (endpoint fixed-scaled около `1.938e-5`): здесь есть настоящая finite-p
  approximation error.
* p64: ошибка `theta_ss` около `3.029e-6` сохраняется даже при MP moments.
  Сохранённая аналитическая пара имеет essential-value roundoff порядка
  `5.5e-16` для w и `2.6e-15` для theta; endpoint differentiation L2-проекции
  усиливает эту малую несовместность с точно нулевым Shen trace. Только как
  diagnostic control вычитание linear endpoint lift уменьшает second-jet error
  до `1e-66` и ниже. Этот lift **не применяется** к target или solver input.
* Exact-input dyadic quintic воспроизводится с MP-погрешностью порядка `1e-73`
  в коэффициентах и `1e-66` в jets. Perturbation одного старшего B94 coefficient
  на `1e-12` обнаруживается; source tails не отбрасываются.

**Strong/weak.** Проверены одно saved p48/p64 initial state и прежний synthetic
velocity state при одинаковых q,v,a и whitening transforms. MP sums над прежними
stored arrays и переоценка basis при старых float64 Gauss nodes/weights сохраняют
relative discrepancy. Только refinement узлов/весов вместе с exact polynomial
basis снижает разность:

| Arithmetic path | p48 relative difference | p64 relative difference |
| --- | ---: | ---: |
| MP45, stored arrays | ≈1.136e-11 | ≈4.305e-11 |
| MP45, analytic basis / old Gauss data | ≈1.135e-11 | ≈4.304e-11 |
| MP45, refined Gauss data | ≈6.68e-44 | ≈8.79e-43 |
| MP70, refined Gauss data | ≈3.79e-68 | ≈8.92e-68 |

Таким образом, в проверенных states локализована quadrature-storage/cancellation
ошибка float64, а не установленный дефект знаков, field ordering или inertial
assembly. Frozen-action independent checks не заменены совпадением symbolic A/B.
Все сравнения используют прежний uncancelled-work scale. Условие2e-12 не ослаблено.
Решение `M(q) delta_a = residual_difference` даёт point-state оценки theta
L2/max около `5.73e-12/7.27e-11` (p48) и `2.17e-11/3.46e-10` (p64).
Это не error bound траектории; никакая такая разность не вычитается из RHS.

### 16.2. Одно initial-only numerical representation rule

До всех ODE выбрано **`common_endpoint_constrained_L2`**, одинаковое для p48/p64
и всех четырёх полей: ближайшее L2-представление того же frozen target в прежнем
Shen-пространстве, сохраняющее его известные первые и вторые derivatives на
обоих концах. Null values уже задаются базисом. Решение использует exact Gram
и небольшой Schur complement, без penalty. Это разрешённый способ initial
approximation после измеренной finite-p ошибки, не новая физическая поправка.

Математически он задаётся `a = a_L2 + G^-1 A^T (A G^-1 A^T)^-1 (g-A a_L2)`;
A здесь — матрица endpoint jets, а не площадь сечения. Во время движения
дополнительные derivative constraints отсутствуют; все4(p−1) coordinates
независимы. Theta3 не подбирается заново и не зависит от p. Float64 coefficients
после MP40/70 подготовки совпадают; никакой filtering/truncation не используется.

Исторические comparison quadrature `max(100,2p+1)`, собственные нормы,
fixed scales `[A,A,A/L,A/L]/L^d`, floor1e-10 и tolerance1e-6 сохранены.
Helper-only diagnostics с отдельными O2 physical scales не используются
для переопределения старых admission gates.

| Initial representation check | p48 | p64 | Criterion/status |
| --- | ---: | ---: | --- |
| max own relative L2/max, all fields and d=0,1,2 | 4.912e-8 | 3.951e-9 | ≤1e-6, PASS |
| max endpoint fixed-scaled error | 1.053e-12 | 1.053e-12 | ≤1e-6, PASS |
| formal through-cubic endpoint coefficient | 2.934e-8 | 2.934e-8 | ≤1e-6, PASS |
| float64 strong/weak absolute max | 3.366e-15 | 1.318e-14 | ≤2e-12, PASS |
| float64 strong/weak relative, both states | 1.136e-11 | 4.305e-11 | ≤2e-12, FAIL |

Basic finite q/RHS, exact essential BC, zero initial velocities, positive mass
and original safety bounds pass. Independent precision evidence объясняет
непройденную проверку. Поэтому каждый run имеет явный execution mode
**`EXPLORATORY_NOT_CERTIFIED`**, а state.admitted остаётся False. Strict
initial verification имеет PARTIAL; разрешение exploratory не становится
default strict admission и не позволяет обходить необъяснённый assembly defect.

### 16.3. Три коротких автономных trajectories

`T1=19.791162590151373`, target end `1.9791162590151374=0.1T1`, nq=2p+1.
Output grid общая:13568 actual timestamps, включая прежние recovery times
и16 samples на период conservative retained-frequency bound от K/M0 Gershgorin
и принятого variable-mass lower bound. Это output sampling, не time-error proof.
Dense output не заменяется грубой интерполяцией прежних snapshots. Max/L2 extrema
по времени и spatial max остаются sampled quantities; вся общая сетка сохранена.

| Run | rtol / relative atol | max_step | accepted steps | nfev/njev/nlu | integration seconds |
| --- | --- | ---: | ---: | --- | ---: |
| p48 tight | 1e-10 / 1e-10 | 0.00360469649 | 550 | 3852/1/314 | 3.382 |
| p64 tight | 1e-10 / 1e-10 | 0.00360469649 | 550 | 3852/1/322 | 5.546 |
| p64 allowed_extra | 2e-11 / 2e-11 | 0.00240313100 | 824 | 5770/1/484 | 8.415 |

Все3 достигли target end; finite values/safety сохранялись. Точный вектор atol,
internal dt, RHS/Jacobian/mass-factorization counters и actual times находятся
в каждом case.json/trajectory.npz/internal_steps.npz. Мало Jacobian refresh и
одинаковое число шагов у двух p согласуются с max_step; это не доказательство
ошибки. Число rejected steps не выдумывается.

Сравнение **p48 tight → p64 tight**, без phase/period/amplitude alignment:

| Component | absolute max-time L2 | absolute max-space-time | relative L2 | relative max | Gate |
| --- | ---: | ---: | ---: | ---: | --- |
| u | 2.3925e-13 | 5.5180e-13 | 2.8260e-07 | 4.0261e-07 | PASS |
| w | 3.8335e-11 | 1.1683e-10 | 2.4286e-08 | 4.6758e-08 | PASS |
| theta | 1.4173e-09 | 3.4297e-09 | 2.6422e-07 | 4.5836e-07 | PASS |
| c | 2.5768e-12 | 5.8657e-12 | 5.7320e-07 | 1.2746e-06 | PASS |
| u_t | 1.2175e-11 | 2.9175e-11 | 4.7424e-05 | 7.0192e-05 | PASS |
| w_t | 2.0684e-09 | 6.4487e-09 | 7.0120e-06 | 1.3815e-05 | PASS |
| theta_t | 7.6360e-08 | 1.9599e-07 | 7.6447e-05 | 1.4053e-04 | FAIL |
| c_t | 1.5899e-10 | 3.7632e-10 | 1.1682e-04 | 2.7029e-04 | PASS |

Семь компонентов проходят. `theta_t` проходит L2, но relative max=1.4053e-4
превышает1e-4: **SHORT_SPATIAL_CHECK=PARTIAL**. u,c и их скорости проверяются
по1e-3; w,theta и их скорости по1e-4. Малые поля не исключены, floor не повышен.
Полные таблицы содержат absolute, own-scale relative и fixed-physical-scale
differences для всех8 components; initial projection uncertainty сохранена отдельно.

Временной контроль **p64 tight → p64 allowed_extra**:

| Component | relative L2 | relative max | Gate |
| --- | ---: | ---: | --- |
| u | 4.5441e-11 | 2.4073e-10 | PASS |
| w | 8.9988e-13 | 5.9688e-12 | PASS |
| theta | 5.9839e-11 | 6.4189e-11 | PASS |
| c | 4.1609e-10 | 1.0028e-09 | PASS |
| u_t | 6.4052e-08 | 3.5498e-07 | PASS |
| w_t | 1.4266e-09 | 1.1260e-08 | PASS |
| theta_t | 1.8013e-08 | 9.4839e-08 | PASS |
| c_t | 1.7072e-07 | 1.8392e-06 | PASS |

Все8 проходят; это temporal evidence только для данного short interval и
выбранной задачи. Spatial theta_t difference заметно больше temporal uncertainty.
Ни один short result не объявляется full5T1 convergence или доказательством
сходимости непрерывной PDE.

Для каждого run используется собственная initial semidiscrete quartic-action
energy, около1.25973835e-9; old energy не выравнивается. Max relative energy
 drift:3.129e-13 /3.065e-13 /4.728e-14, против criterion1e-6. Общий lower bound
relative mass≥0.9999907955, min(1+c)≥0.9999953978, condition upper bound≤1.0000092046.
Прежние bounds для c,theta,u_s,w_s,curvature соблюдены; max|c|≈4.603e-6,
max|theta|≈0.007484, max L|theta_s|≈0.06844. Mass positivity подтверждается
weighted-Gram Loewner bounds без новых eigensolves. Energy не используется
для классификации или correction траектории и не доказывает spatial convergence.

### 16.4. Practical assessment, результаты и остановка

Метод практически рассчитывает это подготовленное short движение:3/3 complete
runs, содержательные two-p/time comparisons и невысокая стоимость. Само по себе
это не обосновывает замену базиса или physical model. Оно даёт практическое
основание рассматривать этот метод далее, **если будет отдельно выбрана новая
задача**, но требуемая strict accuracy сейчас не заявляется.

До/после изменён только numerical initial representation и optional explicit-q0
API прежнего runner; его default projection сохранена. Добавлено сохранение
actual prefix при numerical/safety failure и внутренних dt. Model/RHS/Jacobian,
V0, inertia, BC, quadrature и time prescriptions неизменны. Это не очередная
RHS-оптимизация и не корректировка measured residual. Old/new comparisons здесь
не объявляются quantitative speedup или isolated proof одного механизма:
старый recovery использовал другую initial задачу, spatial pair и time settings.

Charged numerical work≈53.49s из900s; local precision/preparation≈21.90s из180s;
три интегрирования суммарно17.342s. Повторный compute, report-only и plot-only
используют сохранённый cache и выполняют0 ODE/BVP/eigen/history evaluations.
Новых integrations в tests нет; real integration evidence исходит ровно из3 runs.

Новый bundle: `results/planar_prepared_feasibility/284a4039177391d1/`.
В нём: manifests/hashes, immutable-source links, precision ladder и reproducible
read-only local probes, до/после projection rows, decimal/raw/whitened coefficients,
pre-run policy/strict table, общий physical target, actual short histories,
internal dt/counters, all8 spatial/temporal CSV/JSON, energy/safety, execution
source snapshots и3 diagnostic PDF/PNG figures:
`initial_projection_before_after`, `prepared_short_trajectories`, `prepared_short_energy`.

```powershell
python scripts/analysis/prepare_planar_initial_state.py --compute --feasibility
python scripts/analysis/prepare_planar_initial_state.py --report-only results/planar_prepared_feasibility/284a4039177391d1
python scripts/analysis/prepare_planar_initial_state.py --plot-only results/planar_prepared_feasibility/284a4039177391d1
```

Новый [config](../../data/input/planar_prepared_feasibility.json) отделён от
старого strict preparation config. Local precision evidence — generated input;
его отсутствие не запускает повтор старых BVP/spectra/histories. Primary claims
фиксируются в этой tracked note, а не только в ignored results.

| Continuation status | Result |
| --- | --- |
| NLSP_PROJECTION_ARITHMETIC_AUDIT | COMPLETED |
| NLSP_NUMERICAL_REPRESENTATION_FIX | PASS |
| NLSP_STRICT_INITIAL_VERIFICATION | PARTIAL |
| NLSP_PREPARED_FEASIBILITY_RUN | COMPLETED_EXPLORATORY_NOT_CERTIFIED |
| NLSP_PREPARED_SHORT_SPATIAL_CHECK | PARTIAL |
| NLSP_PREPARED_SHORT_TEMPORAL_CHECK | PASS |

Нерешены strict float64 identity gate и требуемый spatial max-критерий theta_t;
продолжение на больших временах/других p этим этапом не выбрано. Не объявлена
exact periodic orbit или out-of-plane stability. Исторический zero-u/c PARTIAL,
LONG CLOSED, EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION и angular same-clamp
reference UNAVAILABLE сохранены. Нет новых amplitudes, full5T1, slope constraints,
physical corrections, filtering, FEM, joint/Floquet study или нового V0. Остановка.


### 16.5. Verification и preservation этого continuation

Адресный combined regression: **202 PASS,2 historical strict XFAIL,7 intentionally
 deselected**,16.24s. Три прежних real ODE tests,3 family-frequency recomputations
и1 mixed runtime/eigen control исключены; realODE entry points заблокированы
во всём test run. Новые tests также отдельно проверяют exact-input polynomial
projection/jets, frozen evaluator/Theta3, explicit-q0/default AST equivalence,
policy только для initial coefficients, separate strict/exploratory admission,
own precision manifests и actual3-run metadata. Historical relative XFAIL
остаётся FAIL на неизменённых float64 данных; искусственный PASS не поставлен.
Тестовый prepared helper suite дополнительно имеет49 PASS в новом lightweight
subset; повторных integration controls в tests нет.

Проверены matching feasibility compute, report-only и plot-only с forbidden
ODE/BVP/eigen/preparation routes: все имеют нулевые numerical counters,
deterministic figure hashes сохранены. Уникальное coverage обеих precision
levels, обоих p и двух states проверено независимо. Пройдены731 relative links
и281 fragments изменённых документов; новые NLSP-D06/K06 anchors уникальны.
`git diff --check` проходит.

HEAD `f60f14370713f84b9ada09ce83dae2d1357ec24f` и staging index сохранены.
Проверены683 исходных tracked hashes;666 файлов неизменны,17 изменённых файлов
находятся в разрешённой numerical/documentation области; добавлен один новый
feasibility config. Frozen model/RHS/Jacobian/old inputs/reference bundles и
старые D/K/text prefixes canonical note сохранены. README/CHANGELOG обновлены
из-за нового explicit user mode; source index/BibTeX и assumptions не менялись:
новой физической assumption нет. Обновлённые memory entries только NLSP-D06/K06;
новый научный этап не выбран.
