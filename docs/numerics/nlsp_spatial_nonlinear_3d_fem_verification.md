# Независимое сопоставление нелинейного пространственного движения 1D и 3D

## Научный результат и область вывода

Все шесть одномерных траекторий и четыре последовательных расчёта CalculiX достигли H=0.25T1.
Семиполевая реализация воспроизводит прежнюю четырёхполевую задачу в плоском пределе. При совместном
возбуждении двух изгибных направлений выделена собственная нелинейная эволюция, включая смешанный
отклик. Физическая модель, коэффициенты и нагрузка не подбирались по результатам FEM.

На fine-сетке полные нелинейные перемещения w и v различаются между 1D и 3D на 3.55% и 3.03% общего
масштаба. Для малых нелинейных поправок различия существеннее: 7.55% для w и 47.29% для v; для эволюции
поправок — 13.15% и 38.80%. Сгущение medium→fine сохраняет большое расхождение v. Удовлетворительное
количественное согласие всех пространственных нелинейных откликов не установлено. Итоговый статус —
**NUMERICAL_PARTIAL**.

Сеточные изменения эволюции w,v малы относительно разностей моделей, однако малые поправки ψ,θ
чувствительны к сетке и секционному восстановлению. Ненулевой поворот Φ и статические крутящие реакции
подтверждают наличие отклика, но количественная точность динамической крутильной деформации не
установлена. Поле c, все скорости, один временной уровень 1D/3D и энергетический учёт CalculiX сохраняют
отдельные ограничения.

Этот результат относится к одному прямому стержню с двухкомпонентной нагрузкой. Прежний
[FEM-3C](nlsp_straight_rod_3d_fem_verification_summary.md) проверял плоское подпространство и сохраняет
свой исторический статус. [Аудит c_eff](nlsp_spatial_profile_audit.md) учитывается при интерпретации
эффективного поперечного сокращения. Экспериментальная проверка не проводилась.

## Постановка и неизменный контракт

\[
L=1,\quad b=0.20,\quad h=0.10,\qquad
E=\rho=1,\quad\nu=0.3,\quad\kappa=5/6.
\]

Это нормированный контроль, а не конкретный экспериментальный материал. Частоты выражены в радианах на
единицу нормированного времени; Hz и Lambda не смешиваются. Сохранены порядок полей, локальные оси и
отображение вращений:

\[
q=(u,w,v,\Phi,\psi,\theta,c),\qquad
B=\operatorname{diag}(1,-1,-1),\qquad a=(\Phi,-\psi,\theta).
\]

Все семь полей равны нулю на обоих концах. Условия на производные не добавлены. В 3D закреплены три
поступательные степени свободы двух торцевых поверхностей; боковые поверхности свободны. Полная 3D
заделка не тождественна всем локальным деформационным ограничениям 1D. Сохранены quartic action, cubic
residuals, переменная масса и аналитический Jacobian. Полная кинематика Rodrigues используется в
диагностике ориентаций, не заменяя принятое усечённое действие. C_T=1.759089824002232e−5 взят из
прежнего FEM-1 reference; его не заменяли на GI_p и не пересчитывали редукцию сечения.

Равномерная нагрузка постоянна в исходной глобальной системе:

\[
f_{\rm vol}=\rho(g_n n+g_k k),\qquad q_w=\rho bhg_n,\quad q_v=\rho bhg_k.
\]

Она направлена по отрицательным глобальным Y,Z и не содержит заданного крутящего момента. Основной
кандидат w_lin,max/h=.04, gk/gn=1.25 дал оценочную угловую изгибную деформацию .00924608819<.01;
fallback .03 не нужен. [Frozen
load](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/frozen_load.json) до
новых FEM-результатов фиксирует gn=.0011379800853485065, gk=.001422475106685633,
qw=2.275960170697013e−5, qv=2.8449502133712663e−5.

Протокол: STATIC preload → мгновенное снятие GRAV → свободный DYNAMIC. Линейная и нелинейная системы
начинаются из собственных статических равновесий с нулевыми скоростями; амплитуды начальных состояний не
выравниваются.

\[
T_1=10.37828159055014,\qquad H=0.25T_1=2.594570397637535.
\]

T1 задан первой линейной 1D изгибной частотой, а не найденной нелинейной орбитой.

## Stage A: соответствие реализации принятому действию

[`SpatialGalerkin`](../../scripts/lib/weakly_nonlinear_spatial_dynamics.py) переиспользует Shen–Legendre
basis и компиляцию сохранённых полиномов. Все 7(p−1) координат независимы, nq>=2p+1; масштабирование по
массе покоя обратимо. Действие загружено из неизменного a9cedd4b6de99295: T/V имеют 36/90 членов. T
квадратична по скоростям и зависит от Φ,ψ,θ,c; связанный вращательный блок массы и все кинетические
координатные силы сохранены.

Точные рациональные проверки дают нулевые разности для семи уравнений Euler–Lagrange, энергетического
тождества, плоского ограничения и отражений. При v=Φ=ψ=0 три неактивных остатка исчезают. Для
равномерной нагрузки относительно середины ожидаются чётные w,v,Φ,c и нечётные u,ψ,θ.

Восемь сохранённых 1D частот FEM-1 проверены при p48/p64 через четыре семидискретных линейных блока на
каждом p: максимальная относительная разность 1.644e−11<прежнего 2e−8. Нового поиска континуальных
корней нет. 13 сохранённых состояний FEM-3AR/FEM-3C воспроизводят энергию, gradient/Hessian, массу,
инерцию, RHS/Jacobian и физические поля. Неактивные ускорения и смешанные блоки Jacobian равны нулю.
Историческая плоская статика не пересчитывалась.

A4 независимо компилирует объёмные и потоковые члены из сохранённого действия и проверяет слабую
проекцию до пространственного дифференцирования:

\[
\int_0^L(\delta q^T r_{\rm body}+\delta q_s^T f_{\rm potential})\,ds
-[\delta q^T f_{\rm potential}]_0^L.
\]

Для функций Shen граничная работа нулевая; отдельная проверка функции с ненулевыми концами контролирует
её знак. Это стабильная численная форма проверки, не изменение RHS, физики или квадратуры. Историческая
проекция сильной формы сохранена отдельно. Знаменатель — прежняя сумма норм несокращённых вкладов
потенциальной работы, floor=1e−30; оба допуска 2e−12.

| Проверка | p48 | p64 | Статус |
|---|---:|---:|---|
| Weak-flux relative difference | 3.852693e−16 | 4.931265e−16 | PASS |
| Weak-flux absolute max | 3.787296e−19 | 4.337867e−19 | PASS |
| Strong-projection relative difference | 2.733824e−12 | 1.043188e−11 | PARTIAL |
| Strong-projection absolute max | 2.229816e−15 | 8.871638e−15 | PASS по абсолютному допуску |

Плоский контроль в том же состоянии воспроизводит прежнюю float64-разность до 2.17e−19. Не прошедший
строгий относительный критерий не заменён на PASS. Энергетическое тождество проходит 2e−12; центральная
directional-проверка Jacobian даёт 8.86e−10<прежнего 2e−7. Проверки массы и сборки проходят. [Stage A
evidence](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/stage_a_checks.json):
PASS_WITH_QUALIFICATIONS.

## Stage B: фактические траектории и численная разрешённость

| Возбуждение | p | Принятые шаги | RHS / Jac / LU | Radau, s | Max собственного energy drift |
|---|---:|---:|---|---:|---:|
| Совместное | 48 | 2999 | 21667 / 2 / 538 | 42.41 | 2.4013e−11 |
| Совместное | 64 | 3912 | 28706 / 2 / 1084 | 101.23 | 3.3254e−11 |
| Только w | 64 | 3751 | 27633 / 2 / 1178 | 99.91 | 2.6380e−11 |
| Только v | 64 | 4335 | 32693 / 1 / 2138 | 144.25 | 3.7588e−11 |
| Только w | 48 | 2947 | 21573 / 2 / 778 | 47.38 | 1.7072e−11 |
| Только v | 48 | 2992 | 21774 / 2 / 664 | 46.04 | 4.0732e−11 |

Все шесть расчётов достигли H; с проверками и постобработкой — 550.04 s<3600 s. При совместном p64
min(1+c)=.99998627, относительные границы массы .996721–1.003312, граница условности 1.006612, max нормы
вращения .0114571. Прежние safety gates и энергетический допуск 1e−6 проходят. Tight componentwise
policy сохранена; временная проверка 1D остаётся PARTIAL_SINGLE_TIGHT_LEVEL.

Смешанный поступательный отклик U_NL(joint)−U_NL(w)−U_NL(v) отделён от наличия двух поперечных компонент
и от геометрической неплоскостности оси. Для конечных ориентаций простое сложение вращательных координат
не применяется.

| Показатель на 0…H | w | v |
|---|---:|---:|
| Max эволюции NL−L | 1.989249e−6 | 4.687739e−7 |
| Её изменение joint p48→p64 | 6.391979e−11 | 1.903819e−11 |
| Max эволюции смешанного отклика | 1.598119e−7 | 3.655408e−7 |
| Её изменение p при собственных isolated controls | 1.651136e−11 | 1.586664e−11 |
| Max смешанной нелинейной эволюции при одинаковых IC | 1.280237e−7 | 1.145188e−7 |

Сохранённые полные линейные операторы отделяют распространение статической разности IC от нелинейной
эволюции при одинаковых IC.
[Разложение](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/initial_state_decomposition.json)
воспроизводит тождество до 5.42e−19 без новых ODE/eigen/static расчётов. Для joint p64 собственная
нелинейная эволюция w/v достигает 2.73109e−6/1.57021e−7. Следовательно, смешанный отклик не сводится к
унаследованной статической поправке. Максимумы при разных x,t не являются складываемыми процентными
долями. Изолированные нагрузки для этого разложения рассчитаны только в 1D, не в 3D.

Первичный config и сравнение четырёх расчётов сохранены отдельно. Исходный §12 допускает три типа
возбуждения, §19 требует два p. Агент интерпретировал их совместно как основание для двух isolated
p48 controls, необходимых для проверки смешанного отклика. Нового сообщения пользователя не было.
[Отдельный config](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/additional_1d_controls_config.json)
фиксирует два дополнительных расчёта и общий максимум шесть. Нагрузка, H, tolerances и бюджет 3600 s
сохранены. Первичные данные не переписаны.

All14 p48/p64 использует прежние масштабы, floors и допуски:

| Компонента | Relative max | Relative L2 | Допуск | Статус |
|---|---:|---:|---:|---|
| u | 1.43212e−4 | 9.05461e−5 | 1e−3 | PASS |
| w | 2.44921e−5 | 1.11700e−5 | 1e−4 | PASS |
| v | 7.36077e−5 | 3.17712e−5 | 1e−4 | PASS |
| Φ | 1.70498e−5 | 7.15189e−6 | 1e−4 | PASS |
| ψ | 8.77039e−5 | 5.13749e−5 | 1e−4 | PASS |
| θ | 9.13275e−5 | 5.03051e−5 | 1e−4 | PASS |
| c | 1.16617e−3 | 4.22220e−4 | 1e−3 | PARTIAL |
| u_t | 7.17257e−3 | 4.51575e−3 | 1e−3 | PARTIAL |
| w_t | 3.50930e−3 | 1.22358e−3 | 1e−4 | PARTIAL |
| v_t | 6.02785e−3 | 2.12713e−3 | 1e−4 | PARTIAL |
| Φ_t | 6.53630e−4 | 2.77404e−4 | 1e−4 | PARTIAL |
| ψ_t | 6.21598e−3 | 2.89385e−3 | 1e−4 | PARTIAL |
| θ_t | 1.02761e−2 | 5.28636e−3 | 1e−4 | PARTIAL |
| c_t | 1.04326e−1 | 4.50701e−2 | 1e−3 | PARTIAL |

Эволюция смешанного w также имеет relative max .00010331747>.0001: PARTIAL. Её абсолютное изменение
p=1.65e−11 мало относительно сигнала 1.598e−7. Обнаружение сигнала и строгая относительная сходимость
отвечают разным вопросам. У p64 max|Φ|=2.15278e−4 и max|χ1|=1.14846e−3. Линейная ось уже имеет max
расстояния от лучшей плоскости 7.67e−5, поэтому неплоскостность сама по себе не доказывает нелинейного
взаимодействия.

[Pre-FEM
decision](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/pre_fem_decision.json)
(SHA256 c47ea46f00e4eb7a2c0dbf3bda90a147f9fba7ec5d1893606823c4848780c1b1) записан до новых 3D
результатов. Ограниченный gate использует проходящие w/v координатные критерии и сигнал>10×наблюдаемых
p/output/recovery индикаторов 5e−9/7.78842e−9. Это правило планирования, не строгая оценка ошибки или
критерий физической точности. All14, mixed-w и временные PARTIAL сохраняются.

## Stage C: четыре расчёта и независимое сопоставление

Переиспользованы C3D10 meshes: medium 5649 узлов/3120 элементов, fine 11553/6670. STATIC и DYNAMIC идут
в одном CCX 2.22 job с передачей полного preload. OP=NEW/zeroGRAV/STEP, ALPHA=0, нулевые скорости и
безопасные наборы energy output сохранены; ELKE в линейном STATIC отсутствует. Новых сеток, модальных
или отдельных статических FEM-расчётов нет. Initial/max increments=T1/8000,T1/4000,
minimum=1e−4×initial, output каждый второй increment. Один поток, 4 GiB; medium timeout 1800 s/fine 4200
s, общий CCX budget 14400 s.

| Расчёт | Actual end | DYNAMIC increments / frames | Native s / peak MiB | Статус |
|---|---:|---:|---|---|
| Medium linear | H | 1002 / 501 | 1129.53 / 99.67 | PASS |
| Medium nonlinear | H | 1002 / 501 | 1201.04 / 99.58 | PASS |
| Fine linear | H | 1002 / 501 | 2991.22 / 210.25 | PASS |
| Fine nonlinear | H | 1002 / 501 | 3149.04 / 210.14 | PASS |

Во всех случаях STATIC содержит 1/10 increments для linear/NL, cutbacks=0. CCX execution вместе с
overhead: 8470.86 s; восстановление полей отдельно 319.43/322.62/646.12/660.08 s. Предварительная оценка
8299.27 s не была гарантией. Первый положительный DYNAMIC time=.00259457; подтверждённое статическое t=0
отличается от native dynamic frame. Preload force/moment relative imbalance <=4.914e−8/1.1765e−8
проходит прежний gate 1e−5. Нагрузка снята, external/damping work=0; DAT/FRD displacement
difference=5e−9. Новых повторных попыток нет.

| Fine midspan | Linear | Nonlinear |
|---|---:|---:|
| Initial w | .0038650872414 | .0038598658086 |
| Initial v | .0016204880159 | .0016199686454 |
| Final w | −.0001519772433 | −.0001549081144 |
| Final v | −.0014944527032 | −.0014942060798 |
| Initial Φ | −8.72743e−6 | 1.74268e−4 |
| Final Φ | 1.20871e−6 | −1.43865e−4 |

[Medium](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/comparison_medium/one_d_three_d_comparison.json)
и [fine
metrics](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/comparison_fine/one_d_three_d_comparison.json)
сопоставлены на 201 общей временной точке и 41 материальной координате. Для процентов ниже общий масштаб
A=max(max|1D|,max|3D|) берётся отдельно для каждой величины на всём 0…H. Разность не делится на
мгновенные значения около нуля. Прежние 3D-reference знаменатели JSON сохранены отдельно. L2 —
пространственная норма, затем максимум по времени; все максимумы относятся к сохранённой сетке.

| Величина | Medium max/A, % | Fine absolute max | Fine max-time L2 | Fine A | Fine max/A, % |
|---|---:|---:|---:|---:|---:|
| Linear w | 3.7564 | 1.413923e−4 | 9.523672e−5 | .004000000 | 3.5348 |
| Linear v | 3.2159 | 5.023609e−5 | 3.395749e−5 | .001666074 | 3.0152 |
| Nonlinear w | 3.7692 | 1.417597e−4 | 9.548158e−5 | .003995172 | 3.5483 |
| Nonlinear v | 3.2309 | 5.047833e−5 | 3.410781e−5 | .001665796 | 3.0303 |
| Δw | 6.8484 | 3.940562e−7 | 2.479817e−7 | 5.221520e−6 | 7.5468 |
| Δv | 47.0350 | 2.456626e−7 | 1.536781e−7 | 5.194846e−7 | 47.2897 |
| δ_evol w | 12.8333 | 3.013122e−7 | 1.946944e−7 | 2.290562e−6 | 13.1545 |
| δ_evol v | 38.5780 | 2.972200e−7 | 1.934830e−7 | 7.659938e−7 | 38.8019 |
| Nonlinear Φ | 19.8997 | 4.276700e−5 | 2.786638e−5 | 2.152783e−4 | 19.8659 |
| ΔΦ | 16.0234 | 3.380859e−5 | 2.789591e−5 | 2.152783e−4 | 15.7046 |
| δ_evol Φ | 12.2199 | 4.454841e−5 | 3.997922e−5 | 3.726176e−4 | 11.9555 |
| Nonlinear ψ | — | 1.213597e−4 | 7.505092e−5 | .003412751 | 3.5561 |
| Nonlinear θ | — | 4.256517e−4 | 2.870605e−4 | .010911271 | 3.9010 |

Полные изгибные перемещения близки по масштабу и ходу движения. Это не скрывает существенную разность
малой Δv. Заранее заданного допуска физической разности нет; произвольный порог 10% не вводился. Ни 1D
p64, ни fine FEM не являются точным континуальным эталоном. Коэффициенты, фазы и масштабы не
подгонялись.

## Сетка, восстановление ориентаций и кручение

[Mesh
comparison](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/mesh_comparison.json)
при одной временной policy отделяет изменение сетки от разности моделей:

| Эволюционная поправка | Medium→fine absolute max | Fine 1D/3D gap | Отношение, % | Fine 41→81 endpoint change / model gap, % |
|---|---:|---:|---:|---:|
| w | 1.010890e−8 | 3.013122e−7 | 3.3550 | .7771 |
| v | 4.780302e−9 | 2.972200e−7 | 1.6083 | .8886 |
| Φ | 1.168365e−6 | 4.454841e−5 | 2.6227 | .1810 |
| ψ | 3.018145e−7 | 5.122499e−7 | 58.9194 | 59.1931 |
| θ | 5.433445e−7 | 7.763599e−7 | 69.9862 | 82.7772 |

Последний столбец использует только static anchor и final dynamic frame, не оценку восстановления на
всём интервале. [Endpoint
evidence](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/recovery_sensitivity_endpoints/fine/metrics.json).
Сеточная разность не строгая граница континуальной ошибки. Для δ_evol w/v изменения меньше расхождения
моделей; для ψ/θ такого количественного вывода нет. Часть больших изменений полных полей при 41→81
сокращается в NL−L; это сокращение не сертифицирует все малые поправки.

Основная временная интерполяция линейная; PCHIP сохранён как диагностика. Для fine δ_evol w/v разность
методов 2.609e−11/1.834e−11, .00866%/.00617% соответствующей разности моделей. Ориентации переносятся
как полные R с проекцией на SO(3), а не суммой углов. C(t)=R_L(t)^T R_NL(t), C(0)^T C(t) определяют
инвариантные поворотные поправки. Малость агрегированной разности ориентаций не сертифицирует малые
ψ/θ-поправки. Один временной уровень нового FEM не даёт самостоятельного temporal certificate.

Для fine нелинейный max|Φ_3D|=1.86795e−4; эволюционный поворот 3.28069e−4 превышает прежний planar proxy
1.10627e−5. Это обнаружение отклика, не доказательство точности крутильной кривизны χ1. Её fine-разность
при двух способах дифференцирования ориентаций 4.58387e−4, medium→fine изменение 4.24782e−4, разность
1D/3D 6.62184e−4: количественное кручение NOT_RESOLVED.

[Статические
реакции](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/static_torsional_support_audit.json)
получены как RF_SUPPORT минус независимо интегрированная consistent bodyload. Исходный
moment_about_face_centroid фактически использовал среднее координат закреплённых узлов. Без изменения RF
выполнен перенос момента на геометрическую ось: M_axis=M_node_mean+(node_mean−axis_center)×F_support.
Нелинейные fine Mx на левом/правом торцах −1.935446e−8/−1.933519e−8 против 1D −2.020239e−8; max разность
4.293% общего масштаба. Линейные Mx порядка 1e−11. Это независимое подтверждение ненулевого статического
крутящего отклика при двухосном изгибе. Оно не проверяет все динамические или warping-связи. Внешний
заданный torque отсутствует; сила на деформированной оси может иметь осевой момент. Dynamic RF как
независимая крутящая реакция NOT_TESTED.

## Деформации, энергия и численные ограничения

Во всех сохранённых кадрах значения конечны. Min sampled detF для medium L/NL .9923076/.9923064, fine
.9914564/.9914476; max компонент Green strain .0101546/.0101587 и .0111385/.0111342. Фактические
локальные деформации превышают оценочную corner bound .00924609: она использовалась для выбора нагрузки,
не как строгий предел всех 3D strains. После наблюдения нагрузка не менялась. Положительность detF между
quadrature samples отдельно не доказана.

Во всех четырёх jobs сохранён raw STATIC/DYNAMIC native energy jump +100%. Native drift внутри DYNAMIC
1.3431–1.3481e−7, external/damping work=0. [Независимая кинетическая
энергия](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/energy_diagnostics.json)
из C3D10 mass/quadrature и фактических FRD velocities отличается от DAT до 6.096e−14 (medium) и
2.642e−14 (fine), 1.96e−6/8.46e−7 общего кинетического масштаба. Восстановление внутренней StVK energy
NOT_RUN; ENERGY PARTIAL. Native reference не перенормирован. Малый drift не заменяет time/mesh controls.

c_eff остаётся effective contraction proxy с ограничениями секционного восстановления из D17/K18; он не
тождественен M–H координате c. All14 c/velocity PARTIAL, mixed-w relative PARTIAL и single-level
temporal PARTIAL сохраняются. Совпадение w/v не подтверждает автоматически все вращательные инерционные
члены, зависимую от c массу, nonlinear constitutive closure и угловой узел.

## Фигуры, данные и воспроизводимость

Четыре основных рисунка сохранены в PDF+PNG, без независимой нормировки кривых:

1. [Пространственная ось](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/spatial_free_motion.pdf).
2. [Две изгибные компоненты](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/two_bending_components.pdf).
3. [Нелинейная поправка и её эволюция](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/nonlinear_spatial_response.pdf).
4. [Смешанный отклик и численная диагностика](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/spatial_coupling_diagnostics.pdf).

[Config](../../data/input/nlsp_spatial_nonlinear_3d_fem_verification.json),
[manifest](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/manifest.json),
[summary](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/summary.json),
[all14
CSV](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/all14_spatial.csv), [fine
CSV](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/comparison_fine/comparison_metrics.csv),
[NLSP-D18](../memory/decisions.md#nlsp-d18), [NLSP-K19](../memory/knowledge.md#nlsp-k19), [scripts
guide](../../scripts/README.md). Bundle сохраняет native INP/DAT/FRD/STA/logs, q/v histories, linear
operators, raw R/gradients, NPZ/CSV и code snapshots.

```powershell
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --spatial-verification --compute
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --spatial-verification --postprocess-only results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --spatial-verification --report-only results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0
```

Postprocess-only читает сохранённые поля/градиенты/моменты без CCX/Gmsh/Radau/ static/eigen/BVP calls.
Четыре фактических CLI replay — compute, report-only, plot-only, postprocess-only — дали ноль
scientific/render calls и неизменный manifest. Focused regressions: 218 PASS за 12.34 s в девяти файлах,
включая nested-manifest guard; скрытых FEM/ODE jobs в тестах нет. Checkout main, HEAD
896220778c4f0a9821a20f981928b6173a4a59ec, initial 23 dirty paths, index пуст. Python 3.12.4:
D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe; NumPy 2.1.3, SciPy 1.15.2, matplotlib 3.9.2,
pytest 8.3.4.

Файлы этого этапа: восемь новых helpers в scripts/lib/ —
weakly_nonlinear_spatial_dynamics, nlsp_spatial_verification_checks,
nlsp_spatial_1d_program, nlsp_spatial_comparison, nlsp_spatial_fem_protocol,
nlsp_spatial_native_program, nlsp_spatial_nonlinear_verification,
nlsp_spatial_recovery_endpoints (.py); восемь одноимённых test_*.py в tests/;
новые config и этот отчёт. В существующем resume CLI добавлены четыре строки
dispatch, прежние режимы сохранены. Дополнены 14 tracked paths, уже dirty до
задания: README, CHANGELOG, memory README/current/decisions/knowledge,
numerics README, journal, research/results indexes, scripts README/STATUS,
lib README и CLI. Исходные 23 dirty paths не объявляются целиком изменениями
этого этапа. [Полная проверка Git/source/tests/cache](../../results/nlsp_spatial_nonlinear_3d_fem_verification/e7e9dee6dbf616f0/final_verification.json).

## Итоговые статусы и научная формулировка

| Проверка | Статус |
|---|---|
| Source preservation / frozen model / planar restriction | PASS |
| Stage A | PASS_WITH_QUALIFICATIONS |
| Strict float64 strong projection | PARTIAL, прежний критерий |
| Шесть 1D траекторий / safety | PASS |
| All14 spatial / mixed-w relative / 1D temporal | PARTIAL |
| Stage B | NUMERICAL_PARTIAL |
| Четыре actual CCX jobs / preload / release / recovery | PASS |
| Stage C | COMPLETE_WITH_QUALIFICATIONS |
| Динамическая крутильная деформация χ1 | NOT_RESOLVED |
| Малые эволюционные ψ/θ-поправки / 3D temporal / energy | PARTIAL |
| Общий итог | NUMERICAL_PARTIAL |

Разработанная семиполевая модель сопоставлена с независимым 3D решением пространственного движения
прямого стержня при неизменных коэффициентах и нагрузке. Обе модели воспроизводят две изгибные
компоненты близкого масштаба; нелинейный пространственный отклик численно обнаружен. Вместе с тем малые
нелинейные поправки согласуются неодинаково: для второй изгибной компоненты сохраняется существенное
расхождение. Динамическая крутильная деформация и часть вращательных поправок количественно
не разрешены. Результат даёт ограниченное основание обсуждать характер пространственного
отклика, но не завершает независимую проверку всех нелинейных связей модели.

V0, variable mass, RHS/Jacobian, материал, геометрия, BC и коэффициенты сохранены. Исторические
FEM-1–FEM-3C bundles/manifests/failed ledgers и D17/K18 не переписаны. LONG CLOSED, EB/RLB-KV
PAUSED_FOR_SUPERVISOR_DIRECTION, angular same-clamp UNAVAILABLE, prepared strict/zero-u-c PARTIAL и
прежние qualifications остаются. Нет новых Gmsh meshes, modal jobs, full-period FEM, amplitudes,
p/time/mesh sweeps, angular-joint, Floquet, periodic-orbit или critical-amplitude расчётов. После отчёта
работа остановлена; PHYSICAL_DYNAMIC_VALIDATION_PASS не присвоен.