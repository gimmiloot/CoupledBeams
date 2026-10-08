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
