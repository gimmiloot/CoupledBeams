# Пространственные профили u, θ, c: диагностический аудит после FEM-3C

Зубцы 3D effective contraction присутствуют уже в raw samples. Они зависят
от локального section fit, ширины окна и распределения quadrature samples.
Кубическая интерполяция добавляет экстремумы. Реальная неоднородность
трёхмерной поперечной деформации подтверждена, но физическая устойчивость
отдельных продольных зубцов остаётся **UNRESOLVED**.

Результат — **DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS**, без новых физических
расчётов. Прежний
[FEM-3C](nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3c-numerical-robustness-and-dissertation-verification)
сохраняет статус
STRAIGHT_ROD_NONLINEAR_3D_FEM_VERIFICATION_COMPLETE_WITH_QUALIFICATIONS.

## Источники и неизменная постановка

Основной bundle: results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef/.
Использованы сохранённые nonlinear nodal displacements, medium/fine C3D10
meshes и p48/p64 coordinate histories. Предыдущие FEM-3B/3AR служат
источниками provenance; их результаты не пересчитывались.

\[
L=1,\quad b=0.20,\quad h=0.10,\quad E=\rho=1,\quad\nu=0.3,\quad\kappa=5/6.
\]

Торцы закреплены, боковые поверхности свободны. Прежние g, q и протокол
STATIC preload → полное снятие GRAV → free DYNAMIC неизменны.
\(T_1=10.37828159055014\) — период первой линейной изгибной формы.
Medium имеет 5649 узлов/3120 C3D10; fine — 11553/6670. Fine dynamic data
доступны только до \(0.25T_1\), новые meshes не строились.

## Что именно означает c_eff

Локальный базис \(B=\operatorname{diag}(1,-1,-1)\): s направлена по global X,
толщина η — по −Y, ширина ζ — по −Z. Положительное 1D w соответствует −Y;
положительное θ — повороту вокруг −Z.

Фактическая функция fem2_recover_reference_samples независимо аппроксимирует
перемещение в каждом материальном участке по восьми функциям
\[
1,d,d^2,d^3,\eta/h,\zeta/b,d\eta/h,d\zeta/b,\qquad
d=(x-x_c)/(L/n).
\]
\(x_c\) — weighted центр samples, n — число участков. Веса получены прежней
14-point positive C3D10 reference-volume quadrature; при \(\rho=1\)
mass и volume weighting совпадают. Из coefficients вычисляются \(G=\nabla U\)
и \(F=I+G\). Добавленные нулевые face values не вводят derivative constraints.

fem2_polar_section_orientation использует \(A=[Fe_2,Fe_3]\), Python F[:,1:3],
и правое полярное разложение \(A=Q\,U^{polar}\),
\(U^{polar}=\sqrt{A^TA}\). Это матрица \(3\times2\), не polar decomposition
всей \(F\). Точные определения:

| Показатель | Фактическое определение | Физический смысл |
|---|---|---|
| c_eff | stretch[0,0]−1 | Normal thickness component правой transverse stretch matrix |
| c_small | gradient[1,1] | Линейная локальная thickness deformation |
| Width proxy | stretch[1,1]−1 | Отдельное изменение масштаба ширины |
| Direct small strain | C3D10 gradient[1,1] | Дифференцирование исходного FE displacement interpolant |
| Direct Green strain | \(E_{22}^{GL}\) | Finite strain из того же FE gradient |

При transverse shear \(|Fe_2|^2=(U_{22}^{polar})^2+(U_{23}^{polar})^2\);
c_eff не является главным stretch или просто длиной одного директора.
Для одной и той же матрицы
\[
E_{22}^{GL}
=c_{\mathrm{small}}+\frac12\sum_iG_{i2}^2
=c_{\mathrm{eff}}+\frac12[c_{\mathrm{eff}}^2+(U_{23}^{polar})^2].
\]
Тождество проверено до порядка \(10^{-15}\). При конечном rigid rotation
Green/polar strain равны нулю, но c_small может быть \(O(\theta^2)\).
Разность мер не является сама по себе recovery error.

В принятой 1D кинематике \(x=r+\eta(1+c)d_2+\zeta d_3\).
Ни c_eff, ни width proxy не объявляются M–H coordinate c.
Один scalar c не описывает всё неоднородное 3D transverse field.

## Фактические времена

Обработаны 11 nonlinear состояний. Nodal temporal interpolation не выполнялась.
Для промежуточных фаз полного периода выбраны ближайшие actual native frames.

| Source | Requested t/T1 | Actual t/T1 | Actual t |
|---|---:|---:|---:|
| Full-period medium | 0 | 0, static preload | 0 |
| Full-period medium | 0.25 | 0.249375000034 | 2.588083972 |
| Full-period medium | 0.5 | 0.499374999973 | 5.182654369 |
| Full-period medium | 0.75 | 0.749375000008 | 7.777224767 |
| Full-period medium | 1 | 1 | 10.37828159055014 |
| Medium/fine refined-time | 0 | 0, static preload | 0 |
| Medium/fine refined-time | 0.125 | 0.125187500037 | 1.299231127 |
| Medium/fine refined-time | 0.25 | 0.25 | 2.594570397637535 |

Actual offsets около −0.00648643 для intermediate full-period states и
+0.00194593 для quarter-control midpoint сохранены в metadata.
1D для raw FEM figures восстановлена из accepted saved Radau polynomials
в тех же actual times, без новых временных шагов. Основная 1D диагностика
использует пять точных сохранённых фаз 0/.25/.5/.75/1.

Исторические nominal-time curves сохранены отдельно в historical_figure_profiles.
Повторный raw41 recovery всех семи fields/11 состояний совпал с сохранёнными
native section profiles до max \(3.44319\cdot10^{-15}\).
Static t=0 не переименован в native dynamic frame.
Старые nominal-time figure curves воспроизведены линейной временной
интерполяцией уже восстановленных source fields до \(1.735\cdot10^{-18}\).
Их max c_eff-разности с ближайшими actual native profiles при .25/.5/.75
равны \(1.09039/1.70352/2.29377\cdot10^{-7}\): это разные физические времена,
а не error bars или оценка временной ошибки; при 0/1 разность нулевая.

## A–B: raw recovery, interpolation и ширина окна

В medium static raw41 c_eff имеется 27 extrema до интерполяции.
CubicSpline добавляет 4; max выход за диапазон соседних samples —
\(1.40891\cdot10^{-6}\). Он воспроизводит raw knots до арифметической точности.

| Medium actual phase | Raw41 extrema | Дополнительные cubic extrema | TV raw c_eff |
|---:|---:|---:|---:|
| 0 | 27 | 4 | 2.20828550e−4 |
| 0.249375 | 27 | 0 | 1.43284648e−5 |
| 0.499375 | 27 | 2 | 2.53246456e−4 |
| 0.749375 | 26 | 3 | 3.27499462e−5 |
| 1 | 29 | 2 | 2.28146295e−4 |

Для исходного medium static preload:

| Sections | Volume-weighted mean c_eff | TV | RMS второй производной | Raw extrema | Condition number |
|---:|---:|---:|---:|---:|---:|
| 21 | −2.04201505e−5 | 5.38885621e−5 | 0.00251226 | 11 | 51.98–55.20 |
| 41 | −2.02072165e−5 | 2.20828550e−4 | 0.0216965 | 27 | 51.48–60.41 |
| 81 | −1.95546396e−5 | 7.17013200e−4 | 0.146569 | 61 | 49.67–68.89 |

Сужение окна усиливает неровность; средний уровень гораздо устойчивее
локальных samples. Все fits rank8, minimum sample count462; по всему
аудиту condition numbers 49.67–69.98. Rank loss или большое amplification
float64 не выявлены.

\(TV=\sum_i|y_{i+1}-y_i|\). Для roughness на неравных физических координатах
используются slopes \(a_i=(y_{i+1}-y_i)/(x_{i+1}-x_i)\) и
\(b_i=2(a_{i+1}-a_i)/(x_{i+2}-x_i)\); RMS \(b_i\) взвешена
\((x_{i+2}-x_i)/2\). Forced endpoint zeros не входят в raw roughness.
Сравнения выполняются на общей физической координате без alignment.

41→81 cubic profile difference max/L2:
\(2.83228\cdot10^{-5}/8.15119\cdot10^{-6}\) при t=0 и
\(2.90878\cdot10^{-5}/8.51378\cdot10^{-6}\) при T1.
Линейная interpolation даёт max \(2.70630\cdot10^{-5}/2.77780\cdot10^{-5}\)
для этих двух состояний. Следовательно, основной эффект не создан cubic.

## C: влияние имеющейся FEM mesh

| Одинаковое actual состояние | Mean c_eff medium/fine,41 | TV medium/fine,41 | Max c_eff gap: linear / cubic |
|---|---:|---:|---:|
| Static | −2.02072e−5 / −2.01574e−5 | 2.20829e−4 / 1.83003e−4 | 2.02179e−5 / 2.07163e−5 |
| 0.1251875T1 | −1.06904e−5 / −1.06531e−5 | 1.63507e−4 / 1.33165e−4 | 1.68816e−5 / 1.72861e−5 |
| 0.25T1 | −2.46186e−8 / −3.62397e−8 | 1.39786e−5 / 1.09193e−5 | 2.22445e−6 / 2.25871e−6 |

При 41 fine уменьшает TV примерно на 17% в static и 22% в конце четверти
периода. Локальные значения/положения peaks остаются неустойчивыми.
При 81 static TV слегка возрастает: 7.17013e−4→7.28008e−4.
Монотонное исчезновение зубцов при refinement не подтверждено.

Контрольные max medium/fine gaps u/w/theta в static:
\(3.22867\cdot10^{-7}/1.13914\cdot10^{-5}/3.26811\cdot10^{-5}\);
при 0.25T1: \(1.92549\cdot10^{-8}/7.44118\cdot10^{-6}/2.67286\cdot10^{-5}\).
Все L2/sign/width/control-field metrics сохранены в JSON/CSV.
Fine states после 0.25T1 **NOT_AVAILABLE**.

## D: прямые FEM strains и механизм загрязнения

C3D10 gradients получены прежними shape functions и 14-point quadrature
из actual nodal U, без section polar fit. В medium static direct Green mean
и pointwise polar-stretch mean тоже имеют зубцы:TV 2.23194e−4/2.23237e−4.
Их неровность не доказывает наличие физически устойчивых longitudinal peaks.

Quadrature samples hard-bin по материальной координате. Это не exact
clipped-tet integral: границы участков пересекают элементы, а их части
не интегрируются геометрически отдельно. Medium41 slab volumes равны
0.967145–1.037684 номинального; поперечный sample centroid отклоняется
до 7.25120e−4L. Усреднение сильной неоднородности может дать sampling bias.

Известное гладкое поле \(U_\eta=\eta^2/L\) точно принадлежит C3D10.
Его gradient воспроизведён до 3.10e−14/4.15e−14 на medium/fine.
Среднее transverse gradient симметричного сечения равно 0, но hard-bin
averages дают \(2\langle\eta\rangle/L\); identity error ≤ 2.14e−16.
Исторический affine fit отличается от этих direct means до 0.000935874
на medium41 и 0.000700495 на fine41. Это kinematic postprocessing control,
не новая scientific equilibrium solution.

Один отдельный 11-column WLS probe добавляет η²,ηζ,ζ² к старым 8 функциям.
Он не заменяет recovery и не является continuum truth. Вложенное LS
тождество связывает изменение affine coefficients с проекцией добавленных
квадратичных функций на старый basis; max остаток 1.26316e−17.

| Raw41 state | TV historical/probe | RMS2 historical/probe | Max raw gap |
|---|---:|---:|---:|
| Medium static | 2.20829e−4 / 5.73567e−5 | .0216965 / .00465893 | 1.29946e−5 |
| Medium actual 0.499375T1 | 2.53246e−4 / 7.81598e−5 | .0230829 / .00792807 | 1.30601e−5 |
| Medium T1 | 2.28146e−4 / 5.66955e−5 | .0226097 / .00406841 | 1.44225e−5 |
| Fine static | 1.83003e−4 / 5.34910e−5 | .0166679 / .00463480 | 9.69958e−6 |

В medium static TV снижается примерно на 74%, mass-weighted displacement residual
norm 4.26019e−6→6.67821e−7. Это локализует существенный вклад поперечного fit,
без сглаживания исходного displacement field.

Реальная неоднородность FE interpolant подтверждена в обеих meshes:
max GL transverse std по samples одного участка .00113076/.00110172,
range .00648039/.00715608, при signed mean около 2e−5.
Statistics включают variation по продольной ширине участка; это не
точные central-section integrals. Scalar c не описывает всё 3D поле,
но устойчивость каждого longitudinal зубца остаётся неустановленной.

В medium static mean c_eff/small/Green =
−2.02072e−5/−6.61076e−5/−1.93811e−5. Max fitted c_eff−c_small —
8.73038e−5, около половины периода 9.11578e−5: finite measures/rotation
существенны. Mean width polar/Green:
−1.51903e−5/−1.53311e−5 initially,
−1.56294e−5/−1.57651e−5 при T1. Width strain отдельна от 1D c.

## E: собственная 1D contraction physics

Восстановлены сохранённые поля/derivatives. Диагностические
\[
\Gamma_1=(1+u_s)\cos\theta+w_s\sin\theta-1,\quad
\Gamma_2=-(1+u_s)\sin\theta+w_s\cos\theta,\quad
N=C(\Gamma_1+\nu c)
\]
не заменяют quartic action/cubic RHS сохранённой динамики.
Retained polynomial measures сохранены отдельно; выбранный max
full-versus-retained Gamma1 gap<=5.79e−9.

\(C=.021978021978021983,\ H=5.341880341880343e-6\),
\(\ell_c=\sqrt{H/C}=.01559023911155809L\).
Interior заранее задана как [3ell,L−3ell], примерно [.04677L,.95323L];
edges не исключены из прежних convergence gates.

| t/T1 | Mean c | Mean Gamma1 | \(\|c+\nu\Gamma_1\|_{L2}/\|c\|_{L2}\), interior | Boundary |
|---:|---:|---:|---:|---:|
| 0 | −1.70407349e−5 | 5.84870477e−5 | .007439 | .504783 |
| .25 | −2.46301588e−9 | 7.17413883e−8 | .120664 | .503657 |
| .5 | −1.74624252e−5 | 6.00523609e−5 | .007597 | .510623 |
| .75 | −5.85262577e−8 | 1.89901829e−7 | .074748 | .543595 |
| 1 | −1.82425991e−5 | 6.25017269e−5 | .007178 | .515580 |

При большом изгибе mean Gamma1>0, mean c<0; interior c≈−nuGamma1 до
0.72–0.76% regional L2. Это согласуется с пуассоновским сокращением в reduced law.
У clamps c=0 не требует Gamma1=0; связь нарушается примерно на 50–52%.
При малом изгибе динамическая часть остаётся: это не exact closure,
поскольку c имеет собственную инерцию, Hc_ss и другие взаимодействия.

L/p48=.0208333L и L/p64=.015625L — только bulk polynomial scales.
First Gauss sample distances1.52080e−4L/8.62091e−5L, counts внутри левой
3ell области 13/18; это не uniform-element mesh.
По пяти фазам max/L2 c sensitivity 2.45233e−8/1.34456e−8;
max на фиксированном five-time масштабе 1.9524946753e−5 равен 0.125600%.
c_s max/L2 sensitivity 7.23257e−6/L/1.56042e−6/L, или 0.634706%/1.053717%
собственных fixed five-time scales. Новый порог не вводится.

p64 unfiltered Legendre tail degrees>=48: max 2.00e−9…6.78e−9 в dynamics.
При 0.25T1 он содержит 0.286% c-L2, но 6.704% c_s-L2. Коротковолновая
составляющая существует; гладкий рисунок не означает numerical certification.

| Historical full-T1 gate | Absolute max / L2 | Relative max / L2 | Threshold | Status |
|---|---:|---:|---:|---|
| c | 4.11640e−8 / 1.37172e−8 | .00193586 / .000743918 | .001 | PARTIAL |
| c_t | 3.25824e−6 / 1.00384e−6 | .1136466 / .0627130 | .001 | PARTIAL |

## F: механический смысл u и θ

В пяти p64 состояниях w,c symmetric, u,theta antisymmetric до arithmetic:
max residuals 2.17e−18/1.09e−19/1.52e−20/1.56e−17; endpoints 0.
Saved initial linear theta совпадает с analytic static Timoshenko до 6.07e−16.
Linear \(w_s=\theta+q(L/2-s)/S\); \(\theta=w_s\) не требуется.
Initially nonlinear max|theta−w_s|=.00221906117,
nonlinear-minus-linear theta max 2.56998619e−5. Axis/sign сохранены.

Классический explanatory benchmark
\(u_s\approx\bar\varepsilon-w_s^2/2\),
\(\bar\varepsilon=(2L)^{-1}\int_0^Lw_s^2ds\) не меняет динамические equations.

| t/T1 | Actual u extrema | Max abs(u), actual/classical | Max benchmark difference |
|---:|---:|---:|---:|
| 0 | 4 | 5.54441e−6 / 6.04095e−6 | 5.33594e−7 |
| .25 | 2 | 2.59465e−7 / 1.51659e−8 | 2.74343e−7 |
| .5 | 4 | 5.14392e−6 / 6.05454e−6 | 1.03290e−6 |
| .75 | 2 | 6.33948e−7 / 2.46156e−8 | 6.56920e−7 |
| 1 | 4 | 5.55056e−6 / 5.63697e−6 | 3.48404e−7 |

Геометрическое приближение воспроизводит качественную форму четырёх лепестков
и масштаб при большом изгибе;
extremum positions differ<=~.004L initially и~.0013L при T1.
При 0.25/0.75 actual u значительно больше benchmark; четыре лепестка
не сохраняются на всех фазах движения.
Сдвиг, независимый поворот, c и продольная инерция ограничивают
применимость приближения. Их причинные вклады отдельно не разделены.

Знак u не определяет знак продольной деформации или усилия: вначале u меняет
знак, но Gamma1,N>0 всюду.
При 0.25/0.75 N меняет знак: диапазоны [−2.66595,4.58992]e−8 /
[−8.83545,7.78361]e−8. N — local axial resultant, а global axial flux
\(F_1=N\cos\theta-Q\sin\theta,\ Q=S\Gamma_2\).
Статические диапазоны N/F1=1.25974e−7/2.46771e−10; nonconstant static N
не означает equilibrium failure. При 0.5T1 диапазон F1 около 1.49302e−7:
динамический flux не постоянен.

## Вывод и безопасное представление в диссертации

Подтверждены recovery/sampling effects и mesh sensitivity. Direct strain
averages также загрязняются hard-binning; устойчивые physical longitudinal
teeth не установлены. Это не диагноз solver defect. Real cross-sectional
heterogeneity остаётся важным отличием 3D поля от scalar coordinate c.

w,u,theta profiles пригодны как рассчитанные иллюстрации с прежними
qualifications; знак u интерпретируется через деформацию и усилия.
1D c показывается вместе с p48/p64 sensitivity и spatial PARTIAL.
3D c_eff лучше показывать raw точками + 41/81 recovery variants, отдельно
от c. Эти разности не error bounds или confidence intervals.
11-column probe подписывается как отдельный diagnostic, не corrected truth.
Никакие исходные curves не сглаживались и не заменялись.

Историческая c/c_eff panel пригодна для демонстрации различия model/proxy,
но не количественного доказательства совпадения M–H transverse fields.
FEM-3C, strict float64/energy/all8 PARTIAL и прежние научные статусы неизменны.
V0, mass, RHS/Jacobian, coefficients, BC и Shen space не менялись.
LONG CLOSED, EB/RLB-KV PAUSED_FOR_SUPERVISOR_DIRECTION,
angular same-clamp UNAVAILABLE сохраняются; после отчёта остановка.

## Четыре фигуры и underlying data

- [Raw/interpolated c_eff PNG](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/raw_vs_interpolated_contraction.png), [PDF](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/raw_vs_interpolated_contraction.pdf): actual native times, 41/81 raw samples, spatial interpolation, отдельный quadratic-fit probe и 1D same-time c.
- [Direct transverse strains PNG](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/native_transverse_strain_reconstruction.png), [PDF](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/native_transverse_strain_reconstruction.pdf): different finite/small measures, medium/fine и width strain.
- [1D contraction physics PNG](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/one_d_contraction_physics.png), [PDF](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/one_d_contraction_physics.pdf): c, −nuGamma1, N, spatial sensitivity and unfiltered tail.
- [u/theta mechanics PNG](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/mechanical_consistency_u_theta.png), [PDF](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/mechanical_consistency_u_theta.pdf): u/derivatives/strain/resultants, benchmark, symmetry/shear.

Bundle: results/nlsp_spatial_profile_audit/c19f5a82203c0260/.
[Manifest](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/manifest.json),
[provenance](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/provenance.json),
[scalar JSON](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/numbers.json),
[FEM CSV](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/FEM_profile_metrics.csv).
states/ содержит raw recovery JSON/NPZ/CSV и quadrature arrays;
one_d_profile_audit/actual_time_one_d сохраняют обе 1D time policies.
historical_figure_profiles хранит старые interpolated nominal-time curves.
[Historical time-policy evidence](../../results/nlsp_spatial_profile_audit/c19f5a82203c0260/historical_time_policy.json)
отдельно сохраняет их точное воспроизведение и actual-time differences.

## Provenance, команды и tests

Checkout initially clean main, HEAD 896220778c4f0a9821a20f981928b6173a4a59ec,
index empty; cwd/root D:\PHD\CoupledBeams\CoupledBeams.
Python D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe, 3.12.4,
NumPy 2.1.3, SciPy 1.15.2. [Config](../../data/input/nlsp_spatial_profile_audit.json)
фиксирует отдельную authorization и 3 parent manifest SHA256:

| Source | Manifest SHA256 |
|---|---|
| FEM-3C | 0fa3488d30c1b36de2061894e2f8811443ab80b802f47fddb99fcb30e5677450 |
| FEM-3B | e4bcc291fab5f04a2fb73103c35a2fa6481aecf5f058b8ce6ae4a5661efa0c4e |
| FEM-3AR | 187870ece572dfe1b999b83d99016d01c1fdbe10e8835ed45543d1119d6f36ea |

Проверены 6 собственных historical manifests и 33 используемых выбранных
source artifacts; это адресный source audit, не повторное хеширование
всего multi-GB archive. Старые D/K byte prefixes сохранены точно.
Основной postprocessing 19.59 s, 1D five-time 4.34 s, actual-time 2.62 s,
figures 7.04 s. New scientific calls=0. Matching compute cache 10.70 s
не запускал solvers. Existing CLI reuse не создал нового physics runner.

~~~powershell
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --profile-audit --preflight
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --profile-audit --compute
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --profile-audit --report-only results/nlsp_spatial_profile_audit/c19f5a82203c0260
python scripts/analysis/resume_nlsp_nonlinear_dynamic_3d_fem.py --profile-audit --plot-only results/nlsp_spatial_profile_audit/c19f5a82203c0260
~~~

157 targeted tests прошли за 4.09 s: три новых test files и legacy resume suite.
Использованы synthetic affine/stretch/rigid-rotation fields, saved meshes/arrays
и запрет scientific calls. Новых FEM/ODE jobs внутри tests нет. Итоговые
link/diff/cache checks сохраняются в verification evidence bundle.
[D17](../memory/decisions.md#nlsp-d17), [K18](../memory/knowledge.md#nlsp-k18)
фиксируют scope/result. Теоретические definitions взяты из
[spatial note](../theory/weakly_nonlinear_spatial_rod.md),
[planar note](../theory/weakly_nonlinear_planar_time_pilot.md);
fem2 recovery/polar/sample/strain functions — из
[existing FEM-2 code](../../scripts/analysis/verify_nlsp_nonlinear_static_3d_fem.py).
