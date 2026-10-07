# Семиполевая пространственная модель: действие и кубический аудит

Дата: 2026-10-07. Это отдельное нелинейное продолжение принятой линейной
модели, с воспроизводимой проверкой алгебры и порядка усечения. Оно не
переоткрывает [закрытую продольную линию](../memory/decisions.md#long-d02)
и не заменяет замороженные линейные решатели или `equations.tex`.

**Область результата:** математический и вычислительный аудит принятого
reduced 1D закона. Общий `NLSP_CUBIC_MODEL_AUDIT=PASS`, если выполнен,
не означает `PHYSICAL_NONLINEAR_VALIDATION_PASS`. Критическая амплитуда,
устойчивость плоского движения и точность реального трёхмерного узла
здесь не определяются.

## 1. Область и независимые поля

Рассматриваются два первоначально прямых, незакрученных изотропных
прямоугольных стержня, однородных внутри каждого плеча. Коэффициенты
между плечами могут различаться. Исходная конструкция плоская; угол
`beta` конечный и следует существующему project geometry contract.
Узел идеальный точечный, жёсткий, без массы и эксцентриситета. Внешние
концы неподвижно закреплены. Нет forcing, предварительного напряжения,
вращения основания, distributed damping или вязкости узла.

Материальная координата \(s\in[0,L_i]\) направлена от внешней заделки к
узлу. Канонический порядок независимых полей одного плеча:

\[
q=(u,w,v,\Phi,\psi,\theta,c),\qquad
U=(u,w,v)^T,\quad z=(\Phi,\psi,\theta)^T.
\]

Здесь \(u\) — продольное, \(w\) — внутриплоскостное и \(v\) —
внеплоскостное перемещения оси; \(\Phi\) — twist, \(\psi,\theta\) —
независимые повороты сечения; \(c\) — инженерное изменение поперечного
масштаба. Не накладываются `theta=w_s`, `psi=v_s`, `c=-nu*u_s` или
нерастяжимость. Динамическая депланация, бимомент и новая warping DOF
не вводятся.

## 2. Источники, входной draft и границы переноса

| Основание | Что используется | Что не переносится |
| --- | --- | --- |
| [Jang-type single rod](mindlin_herrmann_timoshenko_single_rod.md) и [production rigid joint](mindlin_herrmann_timoshenko_rigid_joint.md) | Принятый линейный reduced M-H/Timoshenko предел и common-DOF closure | Нелинейный закон ниже не приписывается Jang; численное source Jang kappa не восстановлено |
| [Перевод обозначений Ярцева](../anisotropic_rods/yartsev_ch2_notation_translation.md), [его joint note](../anisotropic_rods/yartsev_ch2_rigid_angular_joint.md) и existing helper | Изотропный Timoshenko/torsion предел, signs и existing generalized \(C_T\) | HMS/DX-209, monoclinic coupling, complex moduli, source book slope clamp |
| [Crespo da Silva 1988](../literature/source_index.md#crespo_da_silva_1988_flexural_torsional_extensional_formulation), DOI 10.1016/0020-7683(88)90087-X | Структура вариации, boundary work и различение amplitude order; PDF3/(1), PDF5–6/(8)–(14), PDF7–8/(15)–(19) | Авторская no-shear кинематика, \(u=O(\epsilon^2)\), \(EA=O(\epsilon^{-1})\), буквальное совпадение с (18) |
| [Supplied expansion](../../data/input/cubic_seven_field_expansion_supplied.md) | Отдельный аналитический comparison artifact | Перечисленные в draft проверки не считаются ранее выполненными проектными tests |

Пользовательский файл найден как
`C:/Users/Nikita/Downloads/cubic_seven_field_expansion.md`; его исходные
байты сохранены в `data/input/cubic_seven_field_expansion_supplied.md`.
SHA256:
`87d7ce3c0ba1f00e1f282fdf68f6392f6e73f1bc178fdec60aec614e1dad063f`.
Draft не перезаписывается generated appendix. Полное source qualification
Crespo da Silva находится в
[тематической карте](../literature/nonlinear_inplane_outofplane_sources.md).
Его полный авторский вывод независимо не воспроизводился.

## 3. Правый базис и перевод прежних соглашений

Используется постоянный исходный ортонормированный правый базис

\[
B=[t,n,k],\qquad t\times n=k,\qquad k=-E_Z.
\]

В действующем project contract \(t\) направлен от clamp к joint,

\[
t_1=E_X,\quad n_1=-E_Y,\quad
t_2=-\cos\beta\,t_1+\sin\beta\,n_1,\quad
n_2=-\sin\beta\,t_1-\cos\beta\,n_1.
\]

При beta=0 второе плечо ориентировано противоположно первому; при
beta=90 его tangent совпадает с \(n_1\). Эти формулы восстанавливаются
из unchanged `reddy_inplane_geometry`; новая angle convention не создаётся.

Прежний right-handed rectangular/Reddy basis имеет столбцы

\[
B_{old}=[t,+E_Z,n],\qquad
B=B_{old}Q_B,\qquad
Q_B=\begin{pmatrix}1&0&0\\0&0&-1\\0&1&0\end{pmatrix},
\quad Q_B^TQ_B=I,\quad\det Q_B=1.
\]

Таким образом, \(U_{old}=Q_BU=(u,-v,w)^T\): смена порядка осей и знак
внеплоскостного перемещения явны. Локальный rotation vector задаётся как

\[
P=\operatorname{diag}(1,-1,1),\qquad
a=Pz=(\Phi,-\psi,\theta)^T.
\]

Физические малые rotation components равны

\[
Ba=\Phi t-\psi n+\theta k.
\]

Положительный \(\theta\) остаётся project in-plane rotation about \(-E_Z\);
положительный \(\psi\) соответствует повороту about (-n), как в
book out-of-plane translation. Координатный \(P\) не является вращением
геометрии: его знак сохраняет различие именованных полей и physical axes.

## 4. Конечная кинематика и амплитудная степень

\[
r=r_0+BU,\qquad A=BR(a),\qquad R(a)=\exp([a]_\times),
\qquad [a]_\times b=a\times b,
\]
\[
\Gamma=R^T(e_1+U_s)-e_1,\qquad
[\chi]_\times=R^TR_s,\qquad[\Omega]_\times=R^TR_t.
\]

Текущие directors — столбцы \(A\). Поперечный mass director \(d_2\)
масштабируется на \(1+c\); \(c\) не Cartesian component и не является
наклоном оси. Все семь полей имеют первую амплитудную степень:

\[
U/\ell_{ref}=O(\epsilon_a),\qquad z=O(\epsilon_a),\qquad
c=O(\epsilon_a).
\]

Нормированные derivatives имеют ту же степень. Геометрия, material
coefficients, length scale и time scale фиксированы при изменении
\(\epsilon_a\). Этот параметр не равен slenderness, `mu` или `tau`;
продольное перемещение и contraction не понижаются до второй степени.

## 5. Исходное сечение, коэффициенты и размеры

Ширина (b) соответствует прежней \(+E_Z\) оси, толщина (h) — \(n\):

\[
A_0=bh,\quad I_\parallel=bh^3/12,\quad
I_\perp=hb^3/12,\quad I_p=I_\parallel+I_\perp.
\]
\[
m=\rho A_0,\quad j_\parallel=\rho I_\parallel,\quad
j_\perp=\rho I_\perp,\quad G=\frac{E}{2(1+\nu)},\quad\kappa=5/6,
\]
\[
C=\frac{EA_0}{1-\nu^2},\quad H=\kappa GI_\parallel,\quad
S=\kappa GA_0,\quad B_\parallel=EI_\parallel,\quad
B_\perp=EI_\perp,
\]
\[
D_\Gamma=\operatorname{diag}(C,S,S),\qquad
D_\chi=\operatorname{diag}(C_T,B_\perp,B_\parallel).
\]

Численное \(\kappa=5/6\) — сохранённое project decision, а не
восстановленная величина Jang. \(C,H,S,B_\parallel,B_\perp\) считаются
по исходному сечению и постоянны внутри плеча.

| Величина | SI dimension |
| --- | --- |
| (m) | kg/m |
| \(j_\parallel,j_\perp\) | kg m |
| (C,S) | N |
| \(H,B_\parallel,B_\perp,C_T\) | N m² |
| \(U\); \(z,c\); \(\chi,\Omega\) | m; dimensionless; 1/m, 1/s |

\(C_T\) берётся из existing generalized torsional section reduction
Ярцева, при нулевой isotropic coupling \(\overline S_{16}=0\). В
G20 control `Geometry(a=0.20,b=0.05)` переводит book \(I_y\) именно в
новое \(I_\perp\); буквы book `a,b` не равны новым rotation vector \(a\)
и rectangle `b,h`. `C_T` не заменяется на `G*Ip`. Section warping
остаётся сконденсированной в жёсткости без independent dynamic field.

## 6. Принятые неусечённые энергии

Все densities относятся к единице исходной длины:

\[
J(c)=\operatorname{diag}\bigl(j_\perp+(1+c)^2j_\parallel,
j_\perp,(1+c)^2j_\parallel\bigr),
\]
\[
T=\frac m2U_t^TU_t+\frac{j_\parallel}{2}c_t^2
+\frac12\Omega^TJ(c)\Omega,
\]
\[
V^{(0)}=\frac12\Gamma^TD_\Gamma\Gamma+\nu Cc\Gamma_1
+\frac C2c^2+\frac H2c_s^2+\frac12\chi^TD_\chi\chi.
\]

Инерционный закон соответствует аффинному mass field

\[
x=r+\eta(1+c)d_2+\zeta d_3
\]

с интегрированием по исходной массе. Изменяется распределение массы
относительно осей, а не масса; дополнительный общий множитель `1+c`
в \(T\) не вводится. \(\eta\) в этой локальной формуле — материальная
координата точки сечения, а не project thickness-mismatch parameter.

### Статус \(V^{(0)}\)

Это наше принятое нелинейное продолжение reduced law. Кинематика сама
по себе не определяет его однозначно. Его основания: заданный линейный
предел, объективность упругой энергии к общему жёсткому движению,
положительность исходной quadratic strain energy, отсутствие скрытого
prestress и новых подбираемых нелинейных коэффициентов.

Положительность normal block показывает точное completion of square:

\[
\frac C2(\Gamma_1^2+2\nu\Gamma_1c+c^2)
=\frac{EA_0}{2}\Gamma_1^2+\frac C2(c+\nu\Gamma_1)^2,
\qquad EA_0=C(1-\nu^2).
\]

Не включаются дополнительные elastic terms, например
\(c\chi_3^2,\Gamma_1\chi_1^2,c\Gamma_2^2\). Это конститутивное
допущение; не утверждается, что такие terms обязательно выше cubic
order. Зависимая от \(c\) инерция и постоянные исходные elastic
coefficients — две разные принятые редукции. Их сочетание не
объявляется единым полным выводом из 3D Hooke law и не доказывает
точность будущего nonlinear threshold реального соединения.

## 7. Вариация и неусечённые balances

При \(\delta R=R[\xi]_\times\), \(g=e_1+\Gamma\) независимые
вариации дают

\[
\delta\Gamma=R^T\delta U_s-\xi\times g,\qquad
\delta\chi=\xi_s+\chi\times\xi,\qquad
\delta\Omega=\xi_t+\Omega\times\xi.
\]

Определим

\[
N=D_\Gamma\Gamma+\nu Cc e_1,\quad M=D_\chi\chi,\quad
R_c=Hc_s,\quad \ell=J(c)\Omega.
\]

В \(\delta T\) вращательная часть равна
\(\ell\cdot\xi_t-(\Omega\times\ell)\cdot\xi\);
в \(\delta V^{(0)}\) coupling к rotation равен
\(-[g\times N+\chi\times M]\cdot\xi\), вместе с \(M\cdot\xi_s\).
После интегрирования по \(s,t\) внутренняя часть вариации действия
\(\int(T-V^{(0)})\,ds\,dt\) имеет знак
\(-E_U\cdot\delta U-E_{body}\cdot\xi-E_c\delta c\), где

\[
E_U=mU_{tt}-\partial_s(RN),
\]
\[
E_{body}=\partial_t\ell+\Omega\times\ell
-\partial_sM-\chi\times M-g\times N,
\]
\[
E_c=j_\parallel c_{tt}-\partial_s(Hc_s)+C(c+\nu\Gamma_1)
-\frac12\Omega^TJ_{,c}\Omega.
\]

Последний contribution равен
\(-j_\parallel(1+c)(\Omega_1^2+\Omega_3^2)\); его знак следует из
`T-V`, а не выбирается по ожидаемому пределу.

### Coordinate moment residual

Правый Jacobian экспоненциальных координат связывает вариации:

\[
\xi=J_r(a)\delta a=J_r(a)P\delta z,\qquad
J_r(a)=I-\frac12[a]_\times+\frac16[a]_\times^2+\cdots.
\]

Поэтому

\[
E_{body}\cdot\xi
=\bigl(P^TJ_r(a)^TE_{body}\bigr)\cdot\delta z,\qquad
E_z=P^TJ_r(a)^TE_{body}.
\]

Транспонирование определяется dual virtual work. Хотя \(P^T=P\),
материальный residual \(E_{body}\) не равен coordinate residual \(E_z\).
Euler–Lagrange equations по (z) сравниваются именно с последним.

## 8. Действие до четвёртой степени

Обозначим \(\Gamma=\gamma_1+\gamma_2+\gamma_3+O(\epsilon_a^4)\):

\[
\gamma_1=U_s-a\times e_1,
\]
\[
\gamma_2=-a\times U_s+\frac12a\times(a\times e_1),
\]
\[
\gamma_3=\frac12a\times(a\times U_s)
-\frac16a\times[a\times(a\times e_1)].
\]

Для \(d=s,t\), \(h_{d,1}=a_d\), \(h_{d,2}=-a\times a_d/2\),
\(h_{d,3}=a\times(a\times a_d)/6\). При \(d=s\) это \(\chi_n\),
при \(d=t\) — \(\Omega_n\). Положим

\[
J=J_0+J_1+J_2,
\quad J_0=\operatorname{diag}(j_\perp+j_\parallel,j_\perp,j_\parallel),
\]
\[
J_1=2c\operatorname{diag}(j_\parallel,0,j_\parallel),\qquad
J_2=c^2\operatorname{diag}(j_\parallel,0,j_\parallel).
\]

Тогда \(T_{\le4}=T_2+T_3+T_4\), где

\[
T_2=\tfrac m2|U_t|^2+\tfrac{j_\parallel}{2}c_t^2
+\tfrac12\Omega_1^TJ_0\Omega_1,
\]
\[
T_3=\Omega_1^TJ_0\Omega_2+\tfrac12\Omega_1^TJ_1\Omega_1,
\]
\[
T_4=\tfrac12\Omega_2^TJ_0\Omega_2+\Omega_1^TJ_0\Omega_3
+\Omega_1^TJ_1\Omega_2+\tfrac12\Omega_1^TJ_2\Omega_1.
\]

Аналогично \(V_{\le4}=V_2+V_3+V_4\):

\[
V_2=\tfrac12\gamma_1^TD_\Gamma\gamma_1+\nu Cc\gamma_{1,1}
+\tfrac C2c^2+\tfrac H2c_s^2+\tfrac12\chi_1^TD_\chi\chi_1,
\]
\[
V_3=\gamma_1^TD_\Gamma\gamma_2+\nu Cc\gamma_{2,1}
+\chi_1^TD_\chi\chi_2,
\]
\[
V_4=\tfrac12\gamma_2^TD_\Gamma\gamma_2+\gamma_1^TD_\Gamma\gamma_3
+\nu Cc\gamma_{3,1}+\tfrac12\chi_2^TD_\chi\chi_2+\chi_1^TD_\chi\chi_3.
\]

В quartic energy существенно присутствуют не только squares второго
порядка, но и произведения первого и третьего. Из
\(L_{\le4}=T_{\le4}-V_{\le4}\) получаются

\[
E_q=\partial_t L_{q_t}+\partial_s L_{q_s}-L_q
=E_q^{[1]}+E_q^{[2]}+E_q^{[3]}.
\]

Это точные coordinate equations полиномиального действия. Относительно
full reduced model отброшены residuals степени четыре и выше. Ускорения
остаются в mass form; усечённая mass matrix не инвертируется как
несогласованное по порядку приближение.

## 9. Два независимых пути и supplied comparison

**Путь A:** vector kinematic series предыдущего раздела, quartic action,
затем Euler–Lagrange differentiation. Все derivatives полей имеют
первую степень; total derivatives переводят jets согласованно.

**Путь B:** отдельный ряд matrix exponential \(R\), матричные
производные \(R^TR_s,R^TR_t\), независимое построение body balances
из §7 и преобразование \(P^TJ_r^T E_{body}\). Остатки и fluxes пути A
не являются входом пути B.

Installed SymPy отсутствует. Реализована точная sparse polynomial
арифметика `Fraction`: 42 jet symbols, затем 10 symbolic constant
coefficients; material coefficients имеют amplitude degree zero.
Canonical symbol/monomial ordering, рациональные коэффициенты и
coefficient-difference evidence сохраняются в audit artifacts. Lean proof
не выполнялся; результат — exact polynomial computation с независимыми
математическими путями, а не formal proof assistant certificate.

Все 7×3 homogeneous residual comparisons дали нулевые коэффициентные
разности. Отдельная строгая адаптация только 21 supplied polynomial sum
также дала `MATCH` во всех 21 блоках. Общий LaTeX parser не создавался.
Это фактически выполненная новая сверка; прежние заявления draft о
checks не переименованы в project tests.

Полная развёрнутая система вынесена в
[generated appendix](weakly_nonlinear_spatial_rod_expansion_generated.md),
получаемый тем же CLI, и не переписывается вручную в этой note.

## 10. Boundary work и смысл generalized moments

Из того же quartic potential определяются

\[
F=\frac{\partial V_{\le4}}{\partial U_s}=[RN]_{\le3},\quad
p_z=\frac{\partial V_{\le4}}{\partial z_s}
=[P^TJ_r^TM]_{\le3},\quad R_c=Hc_s.
\]

Эти равенства проверены как точные полиномиальные identities. Положительная
boundary work записывается

\[
[F\cdot\delta U+p_z\cdot\delta z+R_c\delta c]_0^L;
\]

в вариации `T-V` corresponding boundary term имеет противоположный
знак. Endpoint sign для nodal efforts: `sigma=-1` при \(s=0\),
`sigma=+1` при \(s=L\). Оба joint ends имеют `sigma=+1`.

\(p_z\) — covector, сопряжённый именованным rotation coordinates,
а не Cartesian physical moment. В full model при обратимом Jacobian

\[
M=J_r^{-T}Pp_z,\qquad m_{spatial}=BRM.
\]

При cubic truncation преобразование physical moment также должно
усекаться согласованно по степени. Нельзя применять inverse либо
переименовывать `p_z` без этого преобразования.

## 11. Жёсткий точечный узел и dual assembly

Точные compatibility conditions:

\[
r_1(L_1)=r_2(L_2),\quad
A_1(L_1)B_1^T=A_2(L_2)B_2^T,\quad c_1(L_1)=c_2(L_2).
\]

В общей малой окрестности единичного relative rotation, где rotation
log однозначен, они эквивалентны common coordinates

\[
B_1U_1=B_2U_2=d_J,\quad
B_1Pz_1=B_2Pz_2=\alpha_J,\quad c_1=c_2=c_J.
\]

Это 3+3+1=7 scalar kinematic conditions. Для ordering `(U,z,c)`

\[
q_{global}=T_{7,i}q_{local},\qquad
T_{7,i}=\operatorname{diag}(B_i,B_iP,1).
\]

Матрица ортогональна. Force/covector mapping не назначается отдельно:

\[
q_{local}=T_{7,i}^Tq_{global},\quad
f_{global}=T_{7,i}f_{local},\quad
f_{local}^T\delta q_{local}=f_{global}^T\delta q_{global}.
\]

Нулевые coefficients boundary variation по common nodal DOFs дают
ещё 7 scalar balances:

\[
\sum_i\sigma_i B_iF_i=0,\qquad
\sum_i\sigma_i B_iPp_{z,i}=0,\qquad
\sum_i\sigma_iR_{c,i}=0.
\]

Средняя сумма — generalized moment balance, сопряжённый общему
rotation vector \(\alpha_J\). В full chart он связан с physical moment
balance общим Jacobian; это не основание называть каждый `BP*p_z`
сырым Cartesian physical moment. Суммарная мощность узла сокращается
в силу этих dual balances, без выбора знаков по спектру.

**Qualification сохраняется:** `c1=c2` — принятое reduced
variational/common-DOF closure, supported by adopted published frame
assembly structure; это не прямой вывод из 3D elasticity конечной области
сварного углового стыка.

Внешняя заделка задаёт `U=0,z=0,c=0`, то есть все семь essential
conditions. Никакие дополнительные `u_s=w_s=v_s=0` не назначаются.
Section-rotation clamp отличается от source `book_slope_clamp` Ярцева;
nonzero-angle historical spectra с другой clamp не объявляются
same-BC reference для новой модели.

## 12. Линейный предел и изолированные подпространства

Степень один даёт требуемые семь equations:

\[
m u_{tt}=[C(u_s+\nu c)]_s,\qquad
j_\parallel c_{tt}=(Hc_s)_s-C(c+\nu u_s),
\]
\[
m w_{tt}=[S(w_s-\theta)]_s,\qquad
j_\parallel\theta_{tt}=(B_\parallel\theta_s)_s+S(w_s-\theta),
\]
\[
m v_{tt}=[S(v_s-\psi)]_s,\qquad
j_\perp\psi_{tt}=(B_\perp\psi_s)_s+S(v_s-\psi),
\]
\[
(j_\parallel+j_\perp)\Phi_{tt}=(C_T\Phi_s)_s.
\]

In-plane block совпадает с unchanged project Jang-type M-H/Timoshenko
после перестановки `(u,c,w,theta)`; out-of-plane block — с isotropic
Yartsev limit после явного axis/sign mapping. General beta совместимость
смешивает physical components через geometry, тогда как в local linear
arm blocks эти subsystems разделены.

### Чистое продольное движение

При `w=v=Phi=psi=theta=0` и нулевых derivatives этих полей остаются
только исходные линейные M-H equations по \(u,c\). Quadratic/cubic terms
исчезают точно. Утверждение относится к прямому плечу; произвольный
angular joint нельзя заменить этим ограничением без преобразования
physical displacement compatibility.

### Плоское подпространство

При `v=Phi=psi=0` соответствующие три residuals равны нулю. Независимая
full planar запись для остальных fields:

\[
\Gamma_1=(1+u_s)\cos\theta+w_s\sin\theta-1,\quad
\Gamma_2=-(1+u_s)\sin\theta+w_s\cos\theta,
\]
\[
N=C(\Gamma_1+\nu c),\qquad Q=S\Gamma_2,
\]
\[
m u_{tt}=[N\cos\theta-Q\sin\theta]_s,\quad
m w_{tt}=[N\sin\theta+Q\cos\theta]_s,
\]
\[
\partial_t[j_\parallel(1+c)^2\theta_t]
=(B_\parallel\theta_s)_s+(1+\Gamma_1)Q-\Gamma_2N,
\]
\[
j_\parallel c_{tt}=(Hc_s)_s-C(c+\nu\Gamma_1)
+j_\parallel(1+c)\theta_t^2.
\]

Её cubic expansion совпадает с restricted spatial derivation.
Reflection `(v,Phi,psi)->(-v,-Phi,-psi)` оставляет действие неизменным;
residuals и fluxes имеют соответствующую parity. Инвариантность
плоского подпространства не доказывает его устойчивость.

Вращение возбуждает \(c\) через inertia даже при \(\nu=0\).
Поэтому `nu=0` не разрешает автоматически удалить `c`; чистая
torsion motion при `c=0` также не объявляется invariant nonlinear
подпространством. No-shear/EB и inextensible limits не развиваются
в отдельное исследование и не реализуются произвольным огромным \(S\).

## 13. Объективность, масса и energy identity

Для full \(V^{(0)}\) постоянное жёсткое движение удовлетворяет

\[
U_s=R(a)e_1-e_1,\quad a_s=0,\quad c=0,
\]

так что \(\Gamma=\chi=0\) и упругие resultants исчезают точно. Общий
постоянный spatial rotation/translation не изменяет strain measures.
Для quartic/cubic модели проверяется сохранённая амплитудная степень;
точная объективность при произвольных больших rotation coordinates
от полинома не требуется.

Hessian \(T_{\le4,q_tq_t}\) симметричен точно; в нуле это

\[
\operatorname{diag}(m,m,m,j_\parallel+j_\perp,j_\perp,j_\parallel,j_\parallel).
\]

Положительность проверяется в нуле и объявленной малой окрестности;
sampled neighborhood check не является доказательством global
positivity усечённой энергии. У full rotational kinetic block
используется \(P^TJ_r^TJ(c)J_rP\), без matrix inverse.

Для quartic energy exact polynomial identity:

\[
\partial_t(T_{\le4}+V_{\le4})
-\partial_s(F\cdot U_t+p_z\cdot z_t+R_cc_t)
-\sum_{\alpha=1}^7 q_{\alpha,t}E_\alpha=0.
\]

Неподвижные essential clamps имеют нулевую мощность. Общие nodal
velocities и dual balances обнуляют сумму powers joint. Это energy
check, не modal energy classification.

Отрицательные контроли строятся в отдельных тестовых expressions:
удаление `Jr^T` из moment residual и удаление \(c\)-dependent inertia
term из \(E_c\) дают ненулевые коэффициентные разности. Правильные
production expressions ради этих контролей не меняются.

## 14. Реализация, вычислительная политика и пределы проверки

Единственный новый isolated helper:
[weakly_nonlinear_spatial_rod.py](../../scripts/lib/weakly_nonlinear_spatial_rod.py).
Он различает full energies/residual, quartic action, cubic coordinate
residual, boundary flux и mass Hessian. Exact model использует Rodrigues
и SO(3) right Jacobian, с analytic directional derivatives; near-zero
entire series предотвращает cancellation. Грубая finite difference для
full derivative не используется. Symbolic generation lazy, не запускается
при обычном импорте helper. Installed SymPy не требуется и не ставился.

[CLI](../../scripts/analysis/verify_weakly_nonlinear_spatial_rod.py) объединяет
source/provenance check, derivation, bounded audit, appendix generation
и cache reuse. [Config](../../data/input/weakly_nonlinear_spatial_rod.json)
фиксирует параметры и numerical policies до интерпретации результатов.
Точность symbolic checks — нулевой coefficient difference, не floating
tolerance. Старые solver tolerances сохраняются.

Manufactured seven-field profiles и analytic jets проверяют residual
`full-cubic` при фиксированных length/time scales и амплитудах
`0.04,0.02,0.01,0.005,0.0025`. Ожидается порядок \(O(\epsilon_a^4)\)
остатка equations и \(O(\epsilon_a^5)\) действия. Компоненты около
заранее заданного numerical floor помечаются отдельно; они не служат
доказательством порядка и не включаются в slope assertion.

Straight split checks сопоставляют direct rod и artificial interface
на двух фиксированных splits с переводом signs. Проверяются bounded
linear controls и nonlinear boundary-power/action-additivity identities;
полная nonlinear PDE trajectory не решается. General-angle checks
касаются geometry, work duality и rank; старые source slope-clamp spectra
не превращаются в same-BC nonlinear validation.

## 15. Фактические результаты и статусы

Исходный checkout: main, HEAD `7b3d667d5418cde247a6ca09fe9186945177537b`,
чистое дерево. Содержимое затрагиваемых документов, hashes 658 tracked files,
Git index и 104 файлов прежних reference bundles сохранены в OS temp до
редактирования. Старые physics/source variants/tolerances не менялись.

Canonical local bundle:
`results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/`.
Полная раскрытая запись остаётся в tracked
[generated appendix](weakly_nonlinear_spatial_rod_expansion_generated.md);
полиномы, коэффициентные разности, manufactured jets и profiles — в bundle.
Supplied MD сохранён отдельно побайтово. Все21 его parts дали `MATCH`.
Строгий adapter читает только эти21 выражения и отвергает неизвестные
symbols/syntax/степени; это не универсальный LaTeX parser.

| Поле | Degree1 A−B | Degree2 A−B | Degree3 A−B | Supplied parts1/2/3 |
| --- | --- | --- | --- | --- |
| u | 0, PASS | 0, PASS | 0, PASS | MATCH / MATCH / MATCH |
| w | 0, PASS | 0, PASS | 0, PASS | MATCH / MATCH / MATCH |
| v | 0, PASS | 0, PASS | 0, PASS | MATCH / MATCH / MATCH |
| Phi | 0, PASS | 0, PASS | 0, PASS | MATCH / MATCH / MATCH |
| psi | 0, PASS | 0, PASS | 0, PASS | MATCH / MATCH / MATCH |
| theta | 0, PASS | 0, PASS | 0, PASS | MATCH / MATCH / MATCH |
| c | 0, PASS | 0, PASS | 0, PASS | MATCH / MATCH / MATCH |

Все7 flux identities, linear/axial/independent planar/reflection,
quartic energy identity, mass symmetry, rigid-motion retained-order flux
и straight reversal/action identities имеют нулевые точные коэффициенты.
Проверены SI dimensions каждого одночлена. Два deliberate omission controls
дают ненулевые расхождения, как и намеренная правка coefficient supplied MD.
Усечённая mass matrix positive в 24 predeclared small-neighborhood samples:
min eigenvalue1.68575e-6. Глобальная положительность усечённой энергии не заявлена.

### Manufactured amplitude check

Все7 полей имеют первую степень. Три заранее заданных analytic sine/cosine
набора проверены в точках (s,t)=(.23,.17),(.61,.39); все jets получены
аналитически, а не конечными разностями. Full evaluator применяет Rodrigues,
точный SO(3) right Jacobian и analytic directional derivatives со stable
near-zero evaluation. Он не является cubic polynomial под другим именем.
Норма собирается по всем6 samples; fixed residual scales — EA/l для
translation и EA для rotations/c, здесь .01. Flux scales — EA и EA*l.

| epsilon_a | ||E_full−E_cubic|| fixed scale | Previous/current | p=log2(ratio) | Boundary difference |
| --- | --- | --- | --- | --- |
| .04 | 6.7675342e-8 | — | — | 5.5472734e-9 |
| .02 | 4.2297129e-9 | 15.9999847 | 3.9999986 | 3.4674371e-10 |
| .01 | 2.6435693e-10 | 16.0000076 | 4.0000007 | 2.1672578e-11 |
| .005 | 1.6522301e-11 | 16.0000073 | 4.0000007 | 1.3547013e-12 |
| .0025 | 1.0326434e-12 | 16.0000063 | 4.0000006 | 8.4585231e-14 |

Заранее выбранный aggregate acceptance interval3.7–4.3 не менялся.
Все component errors/orders сохранены. Последние Phi/psi/theta/c errors
7.26e-14/4.04e-14/3.68e-14/1.01e-14 ниже консервативного predeclared
reporting floor1e-13 и помечены: их последние slopes не используются как
evidence gate. Aggregate остаётся выше floor; предыдущие component values
согласуются с четвёртым порядком. Повышенная точность не потребовалась;
fields, sampling points и амплитуды после результата не подбирались.
Это порядок mass-form residual на prescribed jets, не ошибка nonlinear
частоты, периода или траектории.

### Bounded linear / split controls

G20 E=rho=1,nu=.3,b=.20,h=.05,L=1,κ=5/6;
C_T=2.7001245991303573e-6 из existing generalized rectangular torsion,
I_perp=3.333333333333334e-5, Sbar16=0. Jacobian нового действия сравнен
с independently embedded old M-H/Tim/Yartsev operators: max scaled
coefficient difference1.11e-16. New state shooting использует short-step
expm/positive-diagonal QR; profiles отдельно строятся через spatial
eigenvectors с bounded exponentials. Reference не строится из новой матрицы.

| Sorted position | Linear family | omega |
| --- | --- | --- |
| 1 | in-plane bending | .317474290788 |
| 2 | in-plane bending | .856560057752 |
| 3 | torsion | .867436981301 |
| 4 | out-of-plane bending | 1.038923588955 |
| 5 | in-plane bending | 1.634084000059 |
| 6 | torsion | 1.734873962604 |
| 7 guard | out-of-plane bending | 2.378101731158 |

Old saturated min-max/scalar count9/9 below omega2.7 certifies the required
prefix; sign scans alone are not the certificate. First3 M-H family omega
3.146029487067,6.290908679850,9.433461117452 additionally match old saved
eigenpairs (max relative1.64e-13), since this family is absent from prefix7.
Max direct prefix frequency difference1.09e-12, projected SVD ratio1.02e-12.
Ten family-specific direct profiles match saved old/bounded PDE basis,
old Yartsev state expm and exact torsion, max component L2 difference1.55e-14.

For splits.50/.35: max relative frequency difference1.71e-12;
first3 M-H on each split4.81e-13; selected profile L2 difference7.87e-13;
joint transmission residual1.31e-12, including c/R and physical moments.
Nonlinear prescribed-field action additivity errors≤3.98e-23; internal
covector balance/power cancel, fixed clamps have zero power. No nonlinear
trajectory was solved. Rank14 and virtual-work checks pass at beta0/45/90,
max relative work error2.95e-16.

Available old production in-plane first3 at45/90 match within8.97e-13.
Old out-of-plane angular **same-clamp spectrum reference is UNAVAILABLE**:
the historical Yartsev joint uses book_slope_clamp. No forced comparison,
new angular spectrum study or unsupported analytic-basis extension was run.
Its linear operator, coordinate maps/rank/duality are checked independently.

| Status | Result |
| --- | --- |
| NLSP_MODEL_SPECIFICATION | PASS |
| NLSP_VARIATIONAL_REDERIVATION | PASS |
| NLSP_SUPPLIED_EXPANSION_COMPARISON | MATCH |
| NLSP_BOUNDARY_AND_JOINT_WORK | PASS |
| NLSP_LINEAR_LIMIT | PASS, bounded controls; optional out-of-plane angular reference unavailable |
| NLSP_AXIAL_AND_PLANAR_LIMITS | PASS |
| NLSP_ENERGY_AND_SYMMETRY | PASS |
| NLSP_AMPLITUDE_TRUNCATION_ORDER | PASS |
| NLSP_STRAIGHT_SPLIT_CHECK | PASS |
| NLSP_CUBIC_MODEL_AUDIT | PASS, mathematical/computational audit only |

246 unique targeted/regression checks PASS:57 new,185 unchanged
M-H finite/source/beta0/general/Bishop/Yartsev checks and4 rectangular
Timoshenko checks. The combined238-test run took41.27s; final57 new checks
1.31s;4 rectangular checks.83s. One initial test referred to an absent
metadata key in the old count helper; its assertion was corrected to the
actual certificate schema. Physics and tolerances did not change. Independent
review also made frozen-reference artifact hashes part of cache identity and
guarded the narrow supplied adapter against silent high-degree truncation.
Earlier developmental bundles remain separate provenance, not current cache.
Final bounded compute:2.07s,715 boundary evaluations, one lazy derivation.
Matching compute/report reuse:0 roots,0 derivations. No old studies rerun.

## 16. Воспроизведение из корня проекта

```powershell
python scripts/analysis/verify_weakly_nonlinear_spatial_rod.py --check-sources
python scripts/analysis/verify_weakly_nonlinear_spatial_rod.py --compute
python scripts/analysis/verify_weakly_nonlinear_spatial_rod.py --report-only results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295
python -m pytest tests/test_weakly_nonlinear_spatial_rod.py -q
```

`python` здесь означает рабочий project interpreter. Fingerprint учитывает
config, supplied artifact, source/code hashes и версии среды. Повторный
audit/report переиспользует только совпадающий проверенный artifact;
generated appendix содержит provenance. Не выполняются старые parameter
studies или полный дорогостоящий test suite.

Проверенный interpreter: `D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe`,
Python3.12.4/NumPy2.1.3/SciPy1.15.2. SymPy не установлен; используется
точный sparse Fraction ring с canonical jet/coefficient/monomial order,
не floating samples вместо symbolic proof. Lean не использовался.
Новых пакетов не установлено. Saved single-rod references проверяются по
artifact hashes. При их отсутствии limited saved-profile validation
недоступна и CLI сообщает конкретный отказ; игнорируемые bundles не
объявляются автоматически частью свежего clone.

## 17. Текущая остановка и ограничения

Зафиксированы кинематика, принятое \(V^{(0)}\), boundary work и cubic
operator. Алгебраическая/численная проверка не повышает claim до exact
3D elasticity, физической валидации nonlinear joint или универсальной
пригодности M-H. `c1=c2` сохраняет reduced-frame qualification.

В этом этапе не выполняются nonlinear time integration, periodic orbits,
Galerkin/modal reduction, Floquet analysis, critical-amplitude search,
frequency/amplitude/geometry maps, branch tracking или energy-based modal
classification. Не вводятся forcing, rotational spring или Kelvin–Voigt:
они уже существуют в отдельной EB/RLB ветви и сюда не переносились.
FEM/3D elasticity не запускались, coefficients не подбирались.

Исторические решения сохраняются:
`LONGITUDINAL MODEL QUESTION CLOSED IN THE ADOPTED 1D SCOPE` и
`EB/RLB-KV: PAUSED_FOR_SUPERVISOR_DIRECTION`. Пользователь явно разрешил
отдельную текущую nonlinear spatial branch; из выполнения этого аудита
не следует автоматическое разрешение следующего исследования.
Нелинейное движение и критическая амплитуда пока не рассчитывались.
