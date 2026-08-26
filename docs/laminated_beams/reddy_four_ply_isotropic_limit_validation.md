# Изотропный четырёхслойный предел стержней Reddy

## 1. Область проверки

Этап RLB-1C-ISO проверяет изотропный предел модели двух симметрично
слоистых стержней Reddy с идеальным жёстким узлом. Каждый стержень содержит
ровно четыре слоя одинаковой толщины. Все слои изготовлены из одного
изотропного материала, а углы укладки проходят обычное преобразование
\(Q\rightarrow\overline Q\).

Основной результат Reddy сопоставляется с независимым замкнутым решением для
изотропной балки Тимошенко прямоугольного сечения. Независимый метод развивает
прежнюю circular-реализацию проекта, но принимает размерные свойства сечения.
Он не использует матрицы слоистого материала, редуцированные свойства Reddy,
матрицу состояния Reddy или матрицу узла RLB.

Оба спектральных метода основаны на независимых граничных определителях.
Метод Рэлея--Ритца в RLB-1C-ISO не применяется. Для восьми заранее
объявленных случаев построены и до сопоставления заморожены два coupled root
inventories. Проверены первые 12 корней и root 13 guard, локальные пространства
решений и выбранные физические формы.

## 2. Замороженный baseline

Работа начата от commit
051cfe049536640d7768db5263182b3bed608045 ветви main, subject
Version 0.4.2. Перед научной реализацией рабочее дерево было чистым.

Сохраняются следующие результаты предыдущих этапов:

    RLB-0-SINGLE-BEAM-MODEL: PASS
    RLB-0A-MATHEMATICAL-MODEL: PASS
    RLB-0A-INDEPENDENT-RITZ: PASS
    RLB-0-OVERALL-SCIENTIFIC: PASS_WITH_SOURCE_TABLE_LIMITATIONS
    RLB-1G-COORDINATES: PASS
    RLB-1J-JOINT-THEORY: PASS
    RLB-1A-BETA0-HOMOGENEOUS: PASS
    RLB-1B-BETA0-STEPPED: PASS
    RLB-1-BETA0-ROOT-INVENTORY: PASS

Исторический двухплечевой Ritz-аудит имеет отдельный статус

    RLB-1C0-BETA0-RITZ-BRIDGE: FAIL_CONVERGENCE_AT_N_LE_16

Этот статус означает, что разрешённое глобальное полиномиальное пространство
при \(N\leq16\) не обеспечило заданную сходимость первых 13 частот. Он не
изменяет governing equations, transfer formulation, coordinate contract или
теорию узла. Поэтому исторический Ritz-статус не входит в OVERALL этапа
RLB-1C-ISO.

Замороженными остаются файлы
scripts/lib/reddy_symmetric_laminated_beam.py,
scripts/lib/reddy_inplane_geometry.py и
scripts/lib/reddy_symmetric_coupled_beams.py. Их уравнения, порядок
состояния, \(K=5/6\), свойства ламината, соглашение по \(z\) и результаты
предыдущих этапов в настоящей проверке не изменяются.

## 3. Аудит прежнего изотропного способа

### 3.1. Канонический замкнутый базис Тимошенко

Каноническая regime-complete реализация прежней модели находится в
scripts/lib/variable_length_timoshenko.py. Круглое сечение определено в
строках 44--65, а прежний коэффициент сдвига -- в строке 117. Фабрика
сечения в строках 142--156 задаёт

\[
A=\pi r^2,\qquad I=\frac{\pi r^4}{4},\qquad
G=\frac{E}{2(1+\nu)}.
\]

Связь прежнего параметра \(\Lambda\) с угловой частотой записана в строке
163. Частота cutoff и соответствующее значение \(\Lambda\) вычисляются в
строках 167--172.

Пространственные корни и три режима базиса построены функцией timo_basis
в строках 214--276. Две изгибные функции, непрерывные в пределе cutoff,
заданы в строках 280--354. Третья функция является продольным решением
\(\sin(k_a x)\). В строках 357--374 из этих функций вычисляются

\[
Q_{\mathrm{old}}=KGA\,(w_{\mathrm{old}}'-\psi_{\mathrm{old}}),\qquad
M_{\mathrm{old}}=EI\,\psi_{\mathrm{old}}',\qquad
N_{\mathrm{old}}=EA\,u_{\mathrm{old}}'.
\]

Матрица соединения прежнего метода собирается в строках 377--410. Первое
плечо вычисляется при \(x_1=L_1\), а второе -- при \(x_{2,\mathrm{old}}=-L_2\),
что явно видно в строке 394. Поля второго плеча строятся на сетке
\(0\rightarrow-L_2\) в строке 520.

Такое направление подтверждает независимый reconstruction helper. В
scripts/lib/analytic_coupled_rods_shapes.py отрицательные аргументы второго
плеча заданы в строках 71--76. Endpoint quantities и шесть residuals
соединения вычисляются в строках 189--259. Helper предназначен для диагностики
полей и не заменяет замкнутый базис Тимошенко.

### 3.2. Геометрия отображения

В scripts/lib/in_plane_shape_geometry.py первый Timoshenko display basis
задан в строках 32--36, а второй -- в строках 39--47. Для второго плеча

\[
\boldsymbol t_{2,\mathrm{old}}=(\cos\beta,\sin\beta),\qquad
\boldsymbol n_{2,\mathrm{old}}=(\sin\beta,-\cos\beta).
\]

Эти векторы равны \(-\boldsymbol t_2\) и \(-\boldsymbol n_2\) frozen
RLB-контракта. Геометрическая координата display curve восстанавливается из
\(x_{2,\mathrm{old}}+L_2\) в строках 138--165. Следовательно, отрицательная
локальная координата не является знаком длины. Она согласует внешний конец
\(x_{2,\mathrm{old}}=0\) и узел \(x_{2,\mathrm{old}}=-L_2\) с display basis,
направленным от узла к внешней заделке.

### 3.3. Отделение Euler--Bernoulli модели

Файл src/my_project/analytic/formulas.py хранит прежний
Euler--Bernoulli determinant. Его матрица находится в строках 34--75, а
определитель -- в строке 78. Sign-scan этого определителя реализован в
src/my_project/analytic/solvers.py, строки 10--56.

Эти функции прочитаны для установления истории обозначений, но не входят в
RLB-1C-ISO. Euler--Bernoulli model не является точным компаратором FSDT при
конечной толщине. Её slender-limit проверка должна составлять отдельный этап.

### 3.4. Круглые и общие части прежней модели

К круглому сечению относятся \(A=\pi r^2\), \(I=\pi r^4/4\), прежний
коэффициент сдвига, параметризация радиуса через \(\epsilon\) и сохранение
массы через \(\tau_i\). Общими для балки Тимошенко являются продольное
решение, два изгибных пространственных корня, три режима этих корней,
clamped endpoint columns и шесть условий соединения.

Новая реализация находится в
scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py. Размерные
свойства сечения задаёт SectionProperties в строках 110--215. Фабрики
прямоугольного и контрольного круглого сечений находятся в строках 227--274.
Модуль не импортирует RLB, прежний circular solver, Yartsev или
ShellBuckling и не использует matrix exponential.

## 4. Четырёхслойный изотропный материал

Для каждого слоя принимаются

\[
E_1=E_2=E,\qquad \nu_{12}=\nu,\qquad
G_{12}=G_{13}=G_{23}=G=\frac{E}{2(1+\nu)},\qquad \rho=\mathrm{const}.
\]

Тогда \(\nu_{21}=\nu_{12}E_2/E_1=\nu\). В плоском напряжённом состоянии

\[
Q=
\begin{bmatrix}
\dfrac{E}{1-\nu^2} & \dfrac{\nu E}{1-\nu^2} & 0\\[2mm]
\dfrac{\nu E}{1-\nu^2} & \dfrac{E}{1-\nu^2} & 0\\
0&0&G
\end{bmatrix}.
\]

Для transverse shear ordering \((yz,xz)\)

\[
Q_s=\begin{bmatrix}G&0\\0&G\end{bmatrix}.
\]

Изотропия приводит к

\[
\overline Q(\theta)=Q,\qquad
\overline Q_s(\theta)=Q_s
\]

при любом угле \(\theta\). В расчёте это равенство не используется как
shortcut. Каждый слой проходит функции transformed_reduced_stiffness и
transformed_transverse_shear_stiffness frozen single-beam module, строки
375--425.

Основная стопка имеет вид

\[
[0/90/90/0].
\]

Все четыре слоя имеют толщину \(h/4\). Дополнительные algebraic controls
\([0/0/0/0]\) и \([17/-38/-38/17]\) также состоят ровно из четырёх слоёв.
Стопки задаются снизу вверх в project coordinate \(z_P=-z_R\) и проходят
обычную функцию integrate_laminate, строки 432--477. Однослойная модель в
численном contract этого этапа не используется.

## 5. Интегральные характеристики ламината

Для однородной изотропной толщины аналитически получим

\[
A=hQ,\qquad B=0,\qquad D=\frac{h^3}{12}Q,\qquad A_s=hQ_s.
\]

Массовые моменты на единицу ширины равны

\[
I_0=\rho h,\qquad I_1=0,\qquad I_2=\frac{\rho h^3}{12}.
\]

Free-edge reduction выполняется двумя эквивалентными способами: через
элемент compliance matrix и через Schur complement. Соответствующие функции
находятся в строках 505--592 frozen single-beam module. После умножения на
ширину \(b\)

\[
\begin{aligned}
EA_{\mathrm{axial}}&=Ebh,\\
EI_{\mathrm{bending}}&=\frac{Ebh^3}{12},\\
KGA_{\mathrm{shear}}&=KGbh,\\
m&=\rho bh,\\
J&=\frac{\rho bh^3}{12}.
\end{aligned}
\]

Здесь и далее обе сравниваемые модели используют одно принятое значение
\(K=5/6\). Проверка не устанавливает его универсальность для анизотропных
ламинатов.

## 6. Численный contract и preliminary algebraic gates

Используется набор

\[
E=1,\qquad \rho=1,\qquad \nu=0{,}3,\qquad
L_{\mathrm{total}}=2,\qquad L_{\mathrm{ref}}=1.
\]

Две геометрии имеют параметры

| Геометрия | \(b\) | \(h\) | \(b/h\) | \(L_{\mathrm{ref}}/h\) |
|---|---:|---:|---:|---:|
| G20 | 0,20 | 0,05 | 4 | 20 |
| G10 | 0,40 | 0,10 | 4 | 10 |

Предварительная проверка выполнена в памяти через обычный Qbar pipeline для
обеих геометрий и трёх четырёхслойных стопок. Для относительной матричной
разности использован максимальный компонент, нормированный на максимальный
компонент аналитической матрицы. Нулевые \(B\) и \(I_1\) нормированы на
\(\lVert A\rVert h\) и \(I_0h\) соответственно.

| Gate | Максимальный residual | Порог | Статус |
|---|---:|---:|---|
| In-plane \(\overline Q(\theta)\), \(\theta=0,17,45,90,-38^\circ\) | \(4{,}0412\cdot10^{-16}\) | \(10^{-12}\) | PASS |
| Transverse-shear \(\overline Q_s(\theta)\) | \(1{,}4433\cdot10^{-16}\) | \(10^{-12}\) | PASS |
| \(A,D\) относительно аналитических матриц | \(1{,}4799\cdot10^{-16}\) | \(10^{-12}\) | PASS |
| Масштабированный \(B\) | \(1{,}3004\cdot10^{-17}\) | \(10^{-12}\) | PASS |
| Shear matrix | 0 | \(10^{-12}\) | PASS |
| \(I_0,I_2\) | 0 | \(10^{-12}\) | PASS |
| Масштабированный \(I_1\) | 0 | \(10^{-12}\) | PASS |
| \(EA,EI,KGA,m,J\) | \(4{,}0658\cdot10^{-16}\) | \(10^{-11}\) | PASS |
| Равенство трёх четырёхслойных стопок | \(3{,}4694\cdot10^{-16}\) | \(10^{-11}\) | PASS |
| Compliance/Schur | \(2{,}7756\cdot10^{-16}\) | \(10^{-11}\) | PASS |

Для G20 получены аналитические section targets

\[
A_{\mathrm{sec}}=0{,}01,\quad
I_y=2{,}083333333333334\cdot10^{-6},\quad
KGA=0{,}003205128205128205.
\]

Здесь \(EA=m=0{,}01\) и \(EI=\rho I_y=2{,}083333333333334\cdot10^{-6}\).
Для G10

\[
A_{\mathrm{sec}}=0{,}04,\quad
I_y=3{,}333333333333334\cdot10^{-5},\quad
KGA=0{,}01282051282051282,
\]

причём \(EA=m=0{,}04\) и
\(EI=\rho I_y=3{,}333333333333334\cdot10^{-5}\).

Эти данные устанавливают preliminary gates
RLB-ISO-4PLY-CONSTITUTIVE: PASS и
RLB-ISO-SECTION-REDUCTION: PASS. Они не являются результатом сравнения
coupled spectra.

## 7. Ориентация прямоугольного сечения

Ширина \(b\) направлена вдоль локальной оси
\(\boldsymbol e_y=-\boldsymbol k=\boldsymbol E_Z\), то есть перпендикулярно
плоскости конструкции. Толщина \(h\) направлена вдоль
\(\boldsymbol e_z=\boldsymbol n_i\) в плоскости конструкции.

Изгиб происходит относительно оси \(\boldsymbol e_y\). Поэтому

\[
A_{\mathrm{sec}}=bh,\qquad I_y=\frac{bh^3}{12}.
\]

Величина \(hb^3/12\) и полярный момент инерции в этом этапе не используются.
Вращательная инерция на единицу длины равна \(\rho I_y\).

## 8. Независимая размерная модель Тимошенко

### 8.1. Уравнения одного плеча

Независимый comparator получает \(E,\nu,\rho,b,h,K,L_i,\beta\) без чтения
RLB properties. Он вычисляет

\[
G=\frac{E}{2(1+\nu)},\quad
EA,\quad EI_y,\quad KGA,\quad \rho A,\quad \rho I_y
\]

непосредственно из прямоугольного сечения. При гармоническом движении с
угловой частотой \(\omega>0\) прежнее знаковое соглашение имеет вид

\[
\begin{aligned}
N&=EA\,u', &N'&=-\rho A\,\omega^2u,\\
Q&=KGA\,(w'-\psi), &Q'&=-\rho A\,\omega^2w,\\
M&=EI_y\,\psi', &M'&=-Q-\rho I_y\,\omega^2\psi.
\end{aligned}
\]

Продольный волновой параметр равен

\[
k_a=\omega\sqrt{\frac{\rho A}{EA}},\qquad u(x)=C_a\sin(k_ax).
\]

Эта функция точно удовлетворяет условиям \(u(0)=0\). Она не использует
прежние скрытые значения \(E=1\), \(\rho=1\), radius или \(\epsilon\).

### 8.2. Пространственные корни изгиба

Обозначим \(z=r^2\). Пространственные корни определяются уравнением

\[
EI_yz^2+\omega^2\left(\frac{EI_y\rho A}{KGA}+\rho I_y\right)z
+\frac{\rho I_y\rho A}{KGA}\omega^4-\rho A\omega^2=0.
\]

Устойчивая формула корней реализована в
isotropic_rectangular_timoshenko_coupled_beams.py, строки 413--492.
Частота cutoff равна

\[
\omega_c=\sqrt{\frac{KGA}{\rho I_y}}.
\]

При \(\omega<\omega_c\) один корень \(z_a\) положителен, а второй \(z_b\)
отрицателен. Возникает mixed hyperbolic/trigonometric regime:

\[
a=\sqrt{z_a},\qquad b=\sqrt{-z_b},\qquad
\alpha=\frac{\rho A\omega^2}{KGA},\qquad
h_a=a+\frac{\alpha}{a},\qquad q=b-\frac{\alpha}{b}.
\]

При \(\omega=\omega_c\) выполняется \(z_a=0\); используется аналитический
cutoff limit без подстановки малой искусственной величины. При
\(\omega>\omega_c\) оба корня отрицательны. Тогда

\[
a=\sqrt{-z_a},\qquad b=\sqrt{-z_b},\qquad
h_a=a-\frac{\alpha}{a},\qquad q=b-\frac{\alpha}{b}.
\]

Этот случай образует two-trigonometric-pairs regime.

### 8.3. Clamped bending basis

В mixed regime введём

\[
d_0=a^2+b^2,\qquad r_h=\frac{h_a}{q},\qquad d_1=a-r_hb.
\]

Две канонические пары \((w,\psi)\) имеют вид

\[
\begin{aligned}
w_0&=\frac{\cosh ax-\cos bx}{d_0},&
\psi_0&=\frac{h_a\sinh ax+q\sin bx}{d_0},\\
w_1&=\frac{\sinh ax-r_h\sin bx}{d_1},&
\psi_1&=\frac{h_a(\cosh ax-\cos bx)}{d_1}.
\end{aligned}
\]

В two-trigonometric regime полагается

\[
d_0=b^2-a^2,\qquad r_h=\frac{h_a}{q},\qquad d_1=a-r_hb,
\]

и используются функции

\[
\begin{aligned}
w_0&=\frac{\cos ax-\cos bx}{d_0},&
\psi_0&=\frac{-h_a\sin ax+q\sin bx}{d_0},\\
w_1&=\frac{\sin ax-r_h\sin bx}{d_1},&
\psi_1&=\frac{h_a(\cos ax-\cos bx)}{d_1}.
\end{aligned}
\]

При cutoff аналитический предел записан как

\[
\begin{aligned}
w_0&=\frac{1-\cos bx}{b^2},&
\psi_0&=\frac{\alpha x+q\sin bx}{b^2},\\
w_1&=\frac{\sin bx}{b},&
\psi_1&=-\frac{q(1-\cos bx)}{b}.
\end{aligned}
\]

Все функции удовлетворяют \(w(0)=\psi(0)=0\). Для сохранения непрерывной
ориентации старых столбцов применяется заранее заданное обратимое
масштабирование, реализованное в строках 495--522. Оно не зависит от корней
coupled determinant. Производные, \(Q\), \(M\) и \(N\) вычисляются
аналитически в строках 525--669.

### 8.4. Независимый coupled determinant

Для каждого плеча endpoint rows обозначаются
\(\boldsymbol u_i,\boldsymbol w_i,\boldsymbol\psi_i,\boldsymbol N_i,
\boldsymbol Q_i,\boldsymbol M_i\). Каждый row содержит три столбца в порядке

\[
[\text{bending 0},\text{bending 1},\text{axial}].
\]

Второе плечо вычисляется при \(x_{2,\mathrm{old}}=-L_2\). При
\(c=\cos\beta\), \(s=\sin\beta\) прежняя матрица равна

\[
D_{\mathrm{old}}=
\begin{bmatrix}
\boldsymbol w_1&-c\boldsymbol w_2+s\boldsymbol u_2\\
\boldsymbol u_1&-s\boldsymbol w_2-c\boldsymbol u_2\\
\boldsymbol\psi_1&-\boldsymbol\psi_2\\
\boldsymbol M_1&-\boldsymbol M_2\\
-\boldsymbol Q_1&c\boldsymbol Q_2-s\boldsymbol N_2\\
\boldsymbol N_1&-s\boldsymbol Q_2-c\boldsymbol N_2
\end{bmatrix}.
\]

Матрица строится в строках 676--722 нового comparator. Условие существования
нетривиального решения имеет вид

\[
\det D_{\mathrm{old}}(\omega)=0.
\]

Raw и scaled matrices сохраняются раздельно. Scaling использует только
конечные положительные множители строк и столбцов, строки 770--865. Поэтому
он не изменяет нули определителя. Matrix exponential в comparator отсутствует.

## 9. Circular backward-compatibility gate

Контрольное круглое сечение строится только для сопоставления с замороженным
прежним способом:

\[
A=\pi r^2,\qquad I=\frac{\pi r^4}{4},\qquad
K_{\mathrm{circ}}=\frac{6+12\nu+6\nu^2}{7+12\nu+4\nu^2}.
\]

При \(\mu=\eta=0\) принимаются \(L_1=L_2=1\), \(r=2\epsilon\) и
\(\omega=\epsilon\Lambda^2\). Сопоставляются spatial roots, regime labels,
endpoint columns и полные coupled matrices при \(\beta=0^\circ\) и
\(30^\circ\).

Предварительная isolated-проверка использовала \(\epsilon=0{,}0025\) и
\(0{,}05\). Для mixed, exact-cutoff, двух точек в окрестности cutoff и
two-trigonometric regimes максимальные относительные различия spatial-root
data, endpoint columns и coupled matrices равны нулю в машинной арифметике.
AST-аудит подтвердил отсутствие импорта прежнего circular module в новом
comparator.

Дополнительная независимая проверка первых шести корней legacy circular case
дала максимальную относительную разность (2{,}207\cdot10^{-13}). Вместе с
below-cutoff, exact-cutoff и above-cutoff проверками это даёт статус
`RLB-ISO-LEGACY-RECTANGULAR-ADAPTER: PASS`. Прежний circular solver не
изменялся.

## 10. Каноническая система и преобразование состояний

### 10.1. Глобальный физический базис

Используется frozen RLB-1G contract

\[
\boldsymbol E_X=(1,0,0),\qquad
\boldsymbol E_Y=(0,1,0),\qquad
\boldsymbol E_Z=(0,0,1),\qquad
\boldsymbol k=-\boldsymbol E_Z.
\]

Для RLB

\[
\begin{aligned}
\boldsymbol t_1&=\boldsymbol E_X,&
\boldsymbol n_1&=-\boldsymbol E_Y,\\
\boldsymbol t_2&=-c\boldsymbol t_1+s\boldsymbol n_1,&
\boldsymbol n_2&=-s\boldsymbol t_1-c\boldsymbol n_1.
\end{aligned}
\]

Обе координаты \(x_i\) направлены от внешней заделки к узлу. Физические
перемещение и поворот равны

\[
\boldsymbol d_i=u_i\boldsymbol t_i+w_i\boldsymbol n_i,\qquad
\boldsymbol\vartheta_i=-\psi_i\boldsymbol k.
\]

### 10.2. Первое плечо

Для первого плеча legacy display basis совпадает с
\((\boldsymbol t_1,\boldsymbol n_1)\), а координаты совпадают:

\[
x_{1,\mathrm{old}}=x_1.
\]

Из равенства физических перемещений следует
\(u_1=u_{1,\mathrm{old}}\) и \(w_1=w_{1,\mathrm{old}}\). Старое сдвиговое
соотношение содержит \(w_{\mathrm{old}}'-\psi_{\mathrm{old}}\), тогда как
RLB использует \(w'+\psi\). Поэтому
\(\psi_1=-\psi_{1,\mathrm{old}}\).

Применение constitutive relations даёт полную карту

\[
\begin{bmatrix}
u_1\\w_1\\\psi_1\\N_1\\Q_1\\M_1
\end{bmatrix}_{\mathrm{RLB}}\!(x_1)
=
\operatorname{diag}(1,1,-1,1,1,-1)
\begin{bmatrix}
u_1\\w_1\\\psi_1\\N_1\\Q_1\\M_1
\end{bmatrix}_{\mathrm{old}}\!(x_1).
\]

В частности,

\[
Q_1=KGA(w_{1,\mathrm{old}}'-\psi_{1,\mathrm{old}})=Q_{1,\mathrm{old}},
\qquad M_1=-M_{1,\mathrm{old}}.
\]

### 10.3. Второе плечо

Для второго плеча legacy display basis противоположен RLB basis:

\[
\boldsymbol t_{2,\mathrm{old}}=-\boldsymbol t_2,\qquad
\boldsymbol n_{2,\mathrm{old}}=-\boldsymbol n_2.
\]

Координаты связаны равенством

\[
x_{2,\mathrm{old}}=-x_2.
\]

После приведения перемещений и поворота к одному физическому базису получим

\[
\begin{bmatrix}
u_2\\w_2\\\psi_2\\N_2\\Q_2\\M_2
\end{bmatrix}_{\mathrm{RLB}}\!(x_2)
=
\operatorname{diag}(-1,-1,-1,1,1,1)
\begin{bmatrix}
u_2\\w_2\\\psi_2\\N_2\\Q_2\\M_2
\end{bmatrix}_{\mathrm{old}}\!(-x_2).
\]

Действительно,

\[
\frac{d[-w_{\mathrm{old}}(-x_2)]}{dx_2}
=w_{\mathrm{old}}'(-x_2),
\]

и потому \(Q_{2,\mathrm{RLB}}=Q_{2,\mathrm{old}}\). Для момента две смены
знака, в \(\psi\) и в координате, компенсируются:
\(M_{2,\mathrm{RLB}}=M_{2,\mathrm{old}}\).

### 10.4. Сопоставление строк узла

Обозначим через \(r_j^{\mathrm{RLB}}\) шесть scalar residuals frozen матрицы
узла в порядке displacement \(u,w\), rotation, force \(N,Q\), moment. Через
\(r_j^{\mathrm{old}}\) обозначим строки \(D_{\mathrm{old}}\). Подстановка
полученных карт даёт

\[
\begin{aligned}
r_1^{\mathrm{old}}&=r_2^{\mathrm{RLB}},&
r_2^{\mathrm{old}}&=r_1^{\mathrm{RLB}},&
r_3^{\mathrm{old}}&=-r_3^{\mathrm{RLB}},\\
r_4^{\mathrm{old}}&=-r_6^{\mathrm{RLB}},&
r_5^{\mathrm{old}}&=-r_5^{\mathrm{RLB}},&
r_6^{\mathrm{old}}&=r_4^{\mathrm{RLB}}.
\end{aligned}
\]

Это signed row relation получено до сравнения спектров. Оно не является
подбором по корням. Для первого плеча локальная виртуальная работа сохраняет
знак. Для второго плеча она меняет знак вместе с ориентацией endpoint, что
согласуется с указанной ориентацией строк прежнего определителя.

## 11. Local solution-space gate

До coupled root search должен быть сопоставлен локальный clamped-arm solution
space. Для заданных \(E,\nu,\rho,b,h,K,L,\omega\) строятся две независимые
матрицы endpoint columns:

- RLB clamp map через одну frozen transfer matrix;
- legacy rectangular map через две замкнутые bending-функции и одну axial
  sine-функцию.

Legacy columns сначала переводятся в RLB state по картам раздела 10. После
единого физического nondimensional scaling сравниваются ранги, проекторы
трёхмерных column spaces и главные углы. Отдельно проверяются outer clamp,
дифференциальные уравнения, resultants и virtual-work pairing.

Заранее установлены пороги \(10^{-10}\) вне окрестности cutoff и \(10^{-8}\)
в объявленной окрестности cutoff. Проверка должна охватить mixed regime,
аналитический cutoff limit и two-trigonometric regime. Coupled roots нельзя
вычислять до прохождения этого gate.

Максимальный residual пространства clamped endpoint states равен
(1{,}9064\cdot10^{-13}). Проверка охватила mixed regime, аналитический
cutoff limit и two-trigonometric regime. Поэтому
`RLB-ISO-LOCAL-ARM-EQUIVALENCE: PASS`.

## 12. Политика независимых спектров

После algebraic, circular-backcompat и local-arm gates каждый метод должен
построить собственный inventory без данных другого метода. Рассматриваются
ровно восемь случаев:

| Case | Геометрия | \(L_1\) | \(L_2\) | \(\beta\) |
|---|---|---:|---:|---:|
| ISO-01 | G20 | 1,0 | 1,0 | \(0^\circ\) |
| ISO-02 | G20 | 0,7 | 1,3 | \(0^\circ\) |
| ISO-03 | G20 | 1,0 | 1,0 | \(30^\circ\) |
| ISO-04 | G20 | 1,0 | 1,0 | \(90^\circ\) |
| ISO-05 | G20 | 1,0 | 1,0 | \(-30^\circ\) |
| ISO-06 | G20 | 0,7 | 1,3 | \(30^\circ\) |
| ISO-07 | G20 | 1,3 | 0,7 | \(30^\circ\) |
| ISO-08 | G10 | 1,0 | 1,0 | \(30^\circ\) |

Каждый solver выполняет seed-free global scan, determinant sign brackets и
поиск локальных минимумов нормированного \(\sigma_{\min}\). Затем выполняются
refinement, duplicate handling, cluster/nullity audit и completeness check.
Ищутся первые 12 положительных корней и root 13 guard.

Inventories RLB и legacy rectangular сначала сохраняются и хешируются
раздельно. Корни одного метода не используются как seeds, интервалы или
локаторы другого. Для isolated roots установлен порог относительной разности
\(10^{-8}\). Для clusters сравниваются число корней, multiplicity, nullity и
cluster centre с тем же порогом.

Оба inventory содержат по 104 строки: восемь случаев по 12 основных корней и
одному guard root. Их semantic hashes равны соответственно

    RLB: 4E0DDE518A134F5C58A8A36F25B6C40A9750CFBE361288A0212F9664360EA837
    legacy: 6F2A06D39EC66D2E49DE8D519D201BFBF42150C1BA9D6E493EB2EFEB44F534D6

Все 104 межмодельных сопоставления проходят порог \(10^{-8}\). Максимальная
относительная разность равна \(4{,}56834\cdot10^{-9}\) для ISO-02, root 13.
Близкая пара ISO-04 при \(\Omega=309{,}048931\) и \(309{,}078378\) сохранена
как два различных accepted canonical candidates. Первый корень занимает slot
13, второй входит в post-guard candidate tail.

Старый direct-\(\beta=0\) fixed--fixed 3x3 diagnostic оказался чувствителен к
численной обусловленности mixed-regime closed-form basis. Его числовые строки
сохранены, но получили роль `AUXILIARY_DIAGNOSTIC`, статус
`SUPERSEDED_DIAGNOSTIC_ARTIFACT` и признак
`used_for_scientific_status=false`. Этот diagnostic не входит в основной
scientific gate: эквивалентность прямого coupled-разбиения и одного
fixed--fixed стержня уже доказана в завершённом RLB-1 baseline.
Зафиксированная причина классификации —
`CLOSED_FORM_MIXED_REGIME_BASIS_CONDITIONING`.

## 13. Политика сравнения форм

Формы RLB должны восстанавливаться из nullspace RLB determinant и
распространяться frozen transfer matrices. Формы comparator восстанавливаются
только из nullspace \(D_{\mathrm{old}}\) и замкнутых пространственных функций.
Ritz-реконструкция не используется.

Для ISO-03 заранее выбраны sorted positions 1, 2, 3 и 6. Для ISO-04 выбраны
positions 1, 2 и 3. Sorted position обозначает только место в упорядоченном
спектре и не определяет modal descendant.

После применения карт раздела 10 обе формы переводятся в поля

\[
[u_1,w_1,\psi_1,u_2,w_2,\psi_2]
\]

frozen RLB frame. Массовая норма имеет вид

\[
\sum_{i=1}^{2}\int_0^{L_i}
\left[\rho A(u_i^2+w_i^2)+\rho I_y\psi_i^2\right]\,dx_i=1.
\]

Для isolated roots вычисляется mass MAC. Gate требует
\(1-\mathrm{MAC}\leq10^{-6}\). При close или repeated roots сравниваются
mass-orthonormal subspaces, singular values overlap matrix и главные углы.
Индивидуальные SVD vectors внутри кратного пространства не сопоставляются.

Для каждого метода независимо проверены outer clamps, совместность перемещения
и поворота, глобальное равновесие сил и моментов, дифференциальные уравнения,
energy identity и сходимость сеток 401/801. Все 14 pointwise comparisons и 28
строк физических residuals имеют статус `PASS`. Максимум \(1-\mathrm{MAC}\)
равен \(4{,}441\cdot10^{-16}\), joint compatibility —
\(6{,}343\cdot10^{-12}\), joint equilibrium — \(3{,}727\cdot10^{-12}\),
energy identity — \(6{,}123\cdot10^{-10}\), а grid convergence —
\(2{,}690\cdot10^{-10}\).

## 14. Симметрии

После независимого cross-method comparison проверяются только две необходимые
симметрии:

- ISO-03 при \(\beta=30^\circ\) против ISO-05 при \(\beta=-30^\circ\);
- ISO-06 при \((L_1,L_2)=(0{,}7,1{,}3)\) против ISO-07 при
  \((L_1,L_2)=(1{,}3,0{,}7)\).

Сравниваются первые 12 корней, guard root, clusters и multiplicities. Порог
равен \(10^{-9}\). Эти симметрии не заменяют сопоставление двух методов и не
обобщаются на анизотропные ламинаты или \(B\ne0\).

RLB angle reflection и arm exchange, а также legacy angle reflection проходят
исходный порог. В legacy arm exchange единственное превышение относилось к
ISO-06/ISO-07, slot 12: \(1{,}7902\cdot10^{-9}\). Оба корня были повторно
уточнены независимо только внутри собственных frozen sign-changing brackets,
без root зеркального случая как seed. После этого разность равна
\(4{,}2224\cdot10^{-10}\). Статус auxiliary-проверки —
`PASS_AFTER_BOUNDED_LOCAL_REFINEMENT`; он характеризует точность локализации
корня, а не ошибку математической модели.

## 15. Итоговые статусы

Основные scientific gates:

    RLB-ISO-4PLY-CONSTITUTIVE: PASS
    RLB-ISO-SECTION-REDUCTION: PASS
    RLB-ISO-LEGACY-RECTANGULAR-ADAPTER: PASS
    RLB-ISO-LOCAL-ARM-EQUIVALENCE: PASS
    RLB-ISO-COUPLED-SPECTRUM: PASS
    RLB-ISO-MODE-SHAPES: PASS

Вспомогательные статусы:

    AUX-LEGACY-DIRECT-BETA0-FIXED-FIXED:
        INCONCLUSIVE_NUMERICAL_BASIS_CONDITIONING
    AUX-LEGACY-ARM-EXCHANGE:
        PASS_AFTER_BOUNDED_LOCAL_REFINEMENT

Исторический статус исходного перегруженного контракта сохраняется как
`ORIGINAL-CONTRACT-OVERALL: PARTIAL_PASS`. Главный научный статус:

    SCIENTIFIC_OVERALL:
        PASS_WITH_AUXILIARY_NUMERICAL_QUALIFICATIONS

Разрешённая формулировка результата:

    FOUR_PLY_ISOTROPIC_LIMIT_VALIDATED_FOR_DECLARED_FINITE_CASES_
    AGAINST_INDEPENDENT_RECTANGULAR_TIMOSHENKO_DETERMINANT

## 16. Ограничения и исключения

Этап представляет конечную validation-проверку, а не parameter study. Даже
полный PASS не устанавливает эквивалентность для любого угла, любой геометрии
или любой stacking sequence. Проверка не распространяется на анизотропное
coupling, \(B\ne0\), несимметричные ламинаты или произвольный коэффициент
сдвига.

Не выполняются broad beta sweep, beta refinement, length-ratio sweep,
thickness sweep, stacking-sequence study и \(K\)-sensitivity. Также исключены
кручение, затухание, complex frequencies, FEM, физическая зона узла, масса и
податливость узла, branch tracking, veering/localization study и подготовка
статьи.

Euler--Bernoulli comparison и Rayleigh--Ritz не входят в RLB-1C-ISO.
Численный workflow использует ровно четыре слоя. Исторические RLB outputs и
предыдущий Ritz FAIL не изменяются и не переинтерпретируются.
