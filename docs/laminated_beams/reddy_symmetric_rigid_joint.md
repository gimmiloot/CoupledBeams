# Идеальный жёсткий узел двух симметрично слоистых стержней Reddy

## Область рассмотрения

Документ фиксирует условия идеального точечного соединения двух стержней
Reddy в плоскости рисунка. Узел считается жёстким, безмассовым и не имеющим
собственной вращательной инерции, эксцентриситета или упругой податливости.
В узле отсутствуют сосредоточенные внешние силы и моменты.

Рассматривается состояние каждого плеча

\[
\boldsymbol y_i=
[u_i,w_i,\psi_i,N_i,Q_i,M_i]^{\mathsf T}.
\]

Здесь не вводятся кручение, затухание, несимметричные ламинаты или
пространственная модель зоны соединения. Спектральная проверка данного этапа
ограничена прямым пределом \(\beta=0\); алгебраические формулы узла записаны
для произвольного допустимого угла.

## Замороженный координатный контракт

Используется контракт RLB-1G из
`docs/laminated_beams/reddy_inplane_coordinate_contract.md`. Глобальный
правый базис и нормаль к плоскости рисунка заданы как

\[
\boldsymbol E_X=(1,0,0),\qquad
\boldsymbol E_Y=(0,1,0),\qquad
\boldsymbol E_Z=\boldsymbol E_X\mathbin\times\boldsymbol E_Y,
\qquad
\boldsymbol k=-\boldsymbol E_Z.
\]

Локальные оси обоих плеч направлены от внешней заделки к узлу. При
\(c=\cos\beta\) и \(s=\sin\beta\)

\[
\begin{aligned}
\boldsymbol t_1&=\boldsymbol E_X,
&\boldsymbol n_1&=\boldsymbol k\mathbin\times\boldsymbol t_1,\\
\boldsymbol t_2&=-c\,\boldsymbol t_1+s\,\boldsymbol n_1,
&\boldsymbol n_2&=-s\,\boldsymbol t_1-c\,\boldsymbol n_1.
\end{aligned}
\]

Физический базис Reddy для плеча \(i\) имеет вид

\[
\boldsymbol e_{x,i}=\boldsymbol t_i,\qquad
\boldsymbol e_{y,i}=-\boldsymbol k,\qquad
\boldsymbol e_{z,i}=\boldsymbol n_i.
\]

Программная реализация получает эти векторы и локально-глобальные
преобразования только из `scripts/lib/reddy_inplane_geometry.py`.
Координатные формулы в coupled-модуле не дублируются. На публичной границе
программного интерфейса угол обозначается `beta_rad`.

## Физические переменные

Перемещения, повороты, силы и моменты записываются в глобальном физическом
пространстве как

\[
\begin{aligned}
\boldsymbol d_i&=u_i\boldsymbol t_i+w_i\boldsymbol n_i,\\
\boldsymbol\vartheta_i&=\psi_i\boldsymbol e_{y,i}=-\psi_i\boldsymbol k,\\
\boldsymbol F_i&=N_i\boldsymbol t_i+Q_i\boldsymbol n_i,\\
\boldsymbol m_i&=M_i\boldsymbol e_{y,i}=-M_i\boldsymbol k.
\end{aligned}
\]

Эти определения сохраняют сопряжённость локальных координат и усилий:

\[
\boldsymbol F_i\mathbin\cdot\delta\boldsymbol d_i
=N_i\,\delta u_i+Q_i\,\delta w_i,
\qquad
\boldsymbol m_i\mathbin\cdot\delta\boldsymbol\vartheta_i
=M_i\,\delta\psi_i.
\]

## Инвариантные условия соединения

Непрерывность перемещения и поворота, а также равновесие сил и моментов
задаются непосредственно в глобальном пространстве:

\[
\boldsymbol d_1-\boldsymbol d_2=\boldsymbol0,
\qquad
\boldsymbol\vartheta_1-\boldsymbol\vartheta_2=\boldsymbol0,
\]

\[
\boldsymbol F_1+\boldsymbol F_2=\boldsymbol0,
\qquad
\boldsymbol m_1+\boldsymbol m_2=\boldsymbol0.
\]

Для получения одной канонической ориентации строк первые две векторные
группы проецируются на \((\boldsymbol t_1,\boldsymbol n_1)\), а поворот и
момент — на \(\boldsymbol e_{y,1}\). В частности, вторая и пятая строки
проецируются на \(\boldsymbol n_1=-\boldsymbol E_Y\), а не на
\(\boldsymbol E_Y\).

Если

\[
\boldsymbol q_i=[u_i,w_i,\psi_i]^{\mathsf T},
\qquad
\boldsymbol r_i=[N_i,Q_i,M_i]^{\mathsf T},
\]

то соответствующие физические карты в базис первого плеча равны

\[
G_1=I_3,
\qquad
G_2=
\begin{bmatrix}
-c&-s&0\\
s&-c&0\\
0&0&1
\end{bmatrix}.
\]

Поэтому условия узла имеют компактный вид

\[
G_1\boldsymbol q_1-G_2\boldsymbol q_2=\boldsymbol0,
\qquad
G_1\boldsymbol r_1+G_2\boldsymbol r_2=\boldsymbol0.
\]

## Скалярное разложение

Раскрытие инвариантных условий даёт

\[
\begin{aligned}
u_1+c u_2+s w_2&=0,\\
w_1-s u_2+c w_2&=0,\\
\psi_1-\psi_2&=0,\\
N_1-c N_2-s Q_2&=0,\\
Q_1+s N_2-c Q_2&=0,\\
M_1+M_2&=0.
\end{aligned}
\]

Знаки этих уравнений следуют из физического координатного контракта и не
выбираются по совпадению частот в прямом пределе.

## Полная матрица узла

Для объединённого порядка

\[
\boldsymbol Y_J=
[u_1,w_1,\psi_1,N_1,Q_1,M_1,
 u_2,w_2,\psi_2,N_2,Q_2,M_2]^{\mathsf T}
\]

условия записываются как \(J_{\mathrm{RLB}}(\beta)\boldsymbol Y_J=0\), где

\[
J_{\mathrm{RLB}}(\beta)=
\begin{bmatrix}
1&0&0&0&0&0& c&s&0&0&0&0\\
0&1&0&0&0&0&-s&c&0&0&0&0\\
0&0&1&0&0&0&0&0&-1&0&0&0\\
0&0&0&1&0&0&0&0&0&-c&-s&0\\
0&0&0&0&1&0&0&0&0&s&-c&0\\
0&0&0&0&0&1&0&0&0&0&0&1
\end{bmatrix}.
\]

Непосредственно следует

\[
J_{\mathrm{RLB}}J_{\mathrm{RLB}}^{\mathsf T}=2I_6,
\]

поэтому ранг матрицы равен шести для любого конечного допустимого угла.
В реализации эта матрица строится двумя независимыми способами: из
физических карт и по приведённой замкнутой формуле. Перестановки или смена
знака строк при сравнении не допускаются.

## Виртуальная работа

Для совместного виртуального движения узла существует
\(\delta\boldsymbol q_J\), для которого

\[
\delta\boldsymbol q_i=G_i^{\mathsf T}\delta\boldsymbol q_J.
\]

Тогда суммарная граничная виртуальная работа равна

\[
\begin{aligned}
\delta W_J
&=\sum_{i=1}^{2}\boldsymbol r_i^{\mathsf T}
  \delta\boldsymbol q_i\\
&=\left(G_1\boldsymbol r_1+G_2\boldsymbol r_2\right)^{\mathsf T}
  \delta\boldsymbol q_J=0.
\end{aligned}
\]

В локальных компонентах это тождество имеет вид

\[
\sum_{i=1}^{2}
\left(N_i\,\delta u_i+Q_i\,\delta w_i+M_i\,\delta\psi_i\right)=0.
\]

Следовательно, кинематические условия и равновесие образуют согласованную
пару относительно граничной виртуальной работы.

## Предел \(\beta=0\)

При \(\beta=0\)

\[
\boldsymbol t_2=-\boldsymbol t_1,
\qquad
\boldsymbol n_2=-\boldsymbol n_1,
\]

и условия узла дают

\[
u_2=-u_1,\quad w_2=-w_1,\quad\psi_2=\psi_1,
\qquad
N_2=N_1,\quad Q_2=Q_1,\quad M_2=-M_1.
\]

Иными словами,

\[
\boldsymbol y_{2,J}=R_0\boldsymbol y_{1,J},
\qquad
R_0=\operatorname{diag}(-1,-1,1,1,1,-1),
\]

и для любого состояния первого плеча

\[
J_{\mathrm{RLB}}(0)
\begin{bmatrix}\boldsymbol y_{1,J}\\R_0\boldsymbol y_{1,J}\end{bmatrix}
=\boldsymbol0.
\]

## Предел \(\beta=90^\circ\)

При \(\beta=90^\circ\) скалярные условия принимают вид

\[
u_1+w_2=0,
\qquad
w_1-u_2=0,
\qquad
\psi_1-\psi_2=0,
\]

\[
N_1-Q_2=0,
\qquad
Q_1+N_2=0,
\qquad
M_1+M_2=0.
\]

Этот предел используется только как алгебраическая проверка матрицы узла.
Спектральные расчёты при ненулевом угле в данном этапе не выполняются.

## Исключения и границы вывода

Полученная матрица относится только к идеальному безмассовому точечному
узлу двух симметрично слоистых стержней Reddy с замороженным порядком
физического состояния. Она не описывает конечную зону соединения,
податливость или массу узла, кручение, затухание, несимметричные ламинаты и
пространственные эффекты. Алгебраическая проверка при ненулевых углах не
заменяет спектральную проверку углового соединения.

## Статус

Для углов \(\beta=-90,-30,0,30,90^\circ\) и 100 углов с фиксированным
seed максимальное расхождение инвариантного и замкнутого построителей равно
\(2.2204460493\times10^{-16}\) при допуске \(10^{-14}\). Во всех 105
случаях ранг равен шести. Максимальная нормированная невязка 1003 проверок
виртуальной работы равна \(1.9735944830\times10^{-16}\) при допуске
\(10^{-12}\).

```text
RLB-1J-JOINT-THEORY: PASS
```
