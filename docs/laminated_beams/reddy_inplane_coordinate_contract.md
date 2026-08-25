# Координатный контракт RLB-1G для двух плеч в плоскости рисунка

```text
RLB-1G-COORDINATES: PASS
```

## 1. Область проверки

Контракт вводит физические базисы двух будущих плеч слоистой балки до
построения условий соединения. Проверены направления локальных координат,
трёхмерное вложение обозначений Reddy, преобразование локальных физических
величин и инвариантность виртуальной работы. Матрица соединения,
характеристический определитель и расчёт корней в RLB-1G не строятся.

Проверка выполнена на baseline
`ce1f9ba689bb96cbfb07065348fe3a2fd25542d4` (`main`, commit subject
`Version 0.4.0`). Исходное рабочее дерево было чистым.

## 2. Глобальный базис рисунка

Фиксируется правый ортонормированный базис

\[
\boldsymbol E_X=(1,0,0),\qquad
\boldsymbol E_Y=(0,1,0),\qquad
\boldsymbol E_Z=\boldsymbol E_X\mathbin\times\boldsymbol E_Y=(0,0,1).
\]

На рисунке \(\boldsymbol E_X\) направлен вправо, а \(\boldsymbol E_Y\) —
вверх. Наблюдатель расположен со стороны \(+\boldsymbol E_Z\) и смотрит
вдоль \(-\boldsymbol E_Z\). Вектор

\[
\boldsymbol k=-\boldsymbol E_Z=(0,0,-1)
\]

направлен от наблюдателя в плоскость рисунка и обозначается символом
\(\otimes\). Вектор \(-\boldsymbol k=\boldsymbol E_Z\) направлен к
наблюдателю и обозначается символом \(\odot\).

## 3. Локальные базисы плеч

Обе локальные координаты \(x_i\) направлены от внешней заделки к узлу
будущего соединения. Для первого плеча

\[
\boldsymbol t_1=\boldsymbol E_X,\qquad
\boldsymbol n_1=\boldsymbol k\mathbin\times\boldsymbol t_1=-\boldsymbol E_Y.
\]

При положительном \(\beta\) геометрический луч от узла к внешней заделке
второго плеча равен

\[
\boldsymbol g_2=\cos\beta\,\boldsymbol E_X+\sin\beta\,\boldsymbol E_Y.
\]

Вектор \(\boldsymbol g_2\) направлен от узла к внешней заделке, а
\(\boldsymbol t_2\) — от внешней заделки к узлу. Поэтому

\[
\begin{aligned}
\boldsymbol t_2
  &=-\boldsymbol g_2
   =-\cos\beta\,\boldsymbol t_1+\sin\beta\,\boldsymbol n_1,\\
\boldsymbol n_2
  &=\boldsymbol k\mathbin\times\boldsymbol t_2
   =-\sin\beta\,\boldsymbol t_1-\cos\beta\,\boldsymbol n_1.
\end{aligned}
\]

Явные контрольные случаи имеют вид

| \(\beta\) | \(\boldsymbol g_2\) | \(\boldsymbol t_2\) | \(\boldsymbol n_2\) |
|---:|---|---|---|
| \(0^\circ\) | \((1,0,0)\) | \((-1,0,0)=-\boldsymbol t_1\) | \((0,1,0)=-\boldsymbol n_1\) |
| \(30^\circ\) | \((\sqrt{3}/2,1/2,0)\) | \((-\sqrt{3}/2,-1/2,0)\) | \((-1/2,\sqrt{3}/2,0)\) |
| \(90^\circ\) | \((0,1,0)\) | \((0,-1,0)=\boldsymbol n_1\) | \((-1,0,0)=-\boldsymbol t_1\) |

Для обоих плеч выполнены

\[
\lVert\boldsymbol t_i\rVert=\lVert\boldsymbol n_i\rVert=1,\qquad
\boldsymbol t_i\mathbin\cdot\boldsymbol n_i=0,\qquad
\boldsymbol t_i\mathbin\times\boldsymbol n_i=\boldsymbol k.
\]

## 4. Физическое вложение координат Reddy

Для каждого плеча используется тройка

\[
\boldsymbol e_{x,i}=\boldsymbol t_i,\qquad
\boldsymbol e_{y,i}=-\boldsymbol k=\boldsymbol E_Z,\qquad
\boldsymbol e_{z,i}=\boldsymbol n_i.
\]

Она является правой:

\[
\boldsymbol e_{x,i}\mathbin\times\boldsymbol e_{y,i}
  =\boldsymbol e_{z,i},\qquad
\det[\boldsymbol e_{x,i},\boldsymbol e_{y,i},\boldsymbol e_{z,i}]=+1.
\]

Для первого плеча \(\boldsymbol e_{x,1}=+\boldsymbol E_X\),
\(\boldsymbol e_{z,1}=-\boldsymbol E_Y\), а
\(\boldsymbol e_{y,1}=+\boldsymbol E_Z\) направлен к наблюдателю. Это
согласуется с изображением source convention Reddy: \(x\) направлен вдоль
балки, \(y\) — по ширине, а \(z_R\) по толщине положителен вниз относительно
этого изображения. В проектном вложении каждого плеча положительному
направлению \(z_R\) поставлен в соответствие вектор
\(\boldsymbol e_{z,i}\).

Следует различать три объекта:

- \(z_R\) — физическая координата источника Reddy;
- \(z_P=-z_R\) — проектная координата хранения для интегрирования по слоям;
- \(\boldsymbol e_{z,i}\) — физический единичный вектор в глобальном
  пространстве.

Соглашение хранения RLB-0 не изменено. Его исходное определение остаётся в
[source contract](reddy_ch4_source_contract.md).

## 5. Физические переменные и виртуальная работа

Порядок локального состояния каждого плеча фиксируется как

\[
\boldsymbol y_i=[u_i,w_i,\psi_i,N_i,Q_i,M_i]^{\mathsf T}.
\]

Физические векторы определяются соотношениями

\[
\begin{aligned}
\boldsymbol d_i&=u_i\boldsymbol t_i+w_i\boldsymbol n_i,\\
\boldsymbol\vartheta_i&=\psi_i\boldsymbol e_{y,i}=-\psi_i\boldsymbol k,\\
\boldsymbol F_i&=N_i\boldsymbol t_i+Q_i\boldsymbol n_i,\\
\boldsymbol m_i&=M_i\boldsymbol e_{y,i}=-M_i\boldsymbol k.
\end{aligned}
\]

При фиксированной геометрии, то есть при \(\delta\beta=0\),
ортонормированность базиса непосредственно даёт энергетические пары

\[
\boldsymbol F_i\mathbin\cdot\delta\boldsymbol d_i
  =N_i\,\delta u_i+Q_i\,\delta w_i,
\qquad
\boldsymbol m_i\mathbin\cdot\delta\boldsymbol\vartheta_i
  =M_i\,\delta\psi_i.
\]

Эти равенства проверены отдельно для поступательной и вращательной частей,
а затем для их суммы.

## 6. Сопоставление с прежними соглашениями проекта

Старый модуль `scripts/lib/in_plane_shape_geometry.py` преобразует локальные
поля только в координаты отображения. Он не задаёт физический локальный
RLB-базис.

Для Timoshenko в старом модуле

\[
\begin{aligned}
\boldsymbol t_{1,\mathrm{old}}&=\boldsymbol t_1,&
\boldsymbol n_{1,\mathrm{old}}&=\boldsymbol n_1,\\
\boldsymbol t_{2,\mathrm{display}}&=-\boldsymbol t_2,&
\boldsymbol n_{2,\mathrm{display}}&=-\boldsymbol n_2.
\end{aligned}
\]

Следовательно, касательная старого display basis параметризует второе плечо
от узла к внешней заделке, то есть противоположно новому направлению
\(x_2\).

Старый EB helper использует другой знак поперечного поля:

\[
\boldsymbol n_{1,\mathrm{EB}}=-\boldsymbol n_1,\qquad
\boldsymbol t_{2,\mathrm{EB,display}}=-\boldsymbol t_2,\qquad
\boldsymbol n_{2,\mathrm{EB,display}}=\boldsymbol n_2.
\]

Поэтому ни `rod2_display_basis`, ни `eb_rod2_display_basis` не переносятся
в новый helper как физический базис.

Production FEM использует отдельное локально-глобальное преобразование

\[
\boldsymbol q_{\mathrm{global}}=\boldsymbol T\boldsymbol q_{\mathrm{local}},
\]

для которого трёхмерное продолжение правого продольного направления равно
\((\cos\beta,\sin\beta,0)=-\boldsymbol t_2\), а правого поперечного
направления — \((-\sin\beta,\cos\beta,0)=\boldsymbol n_2\). FEM-нумерация
второго плеча также идёт от узла к внешней заделке. Это преобразование
сохраняется без изменений и не отождествляется с парой
\((\boldsymbol t_2,\boldsymbol n_2)\).

Старые analytic fields используют знаковую координату
\(x_2:0\rightarrow-L_2\), после чего точки второго плеча разворачиваются
для отображения в порядке узел — внешняя заделка. Эти координаты и их
производные не переносятся в RLB-1G неявно.

Существующий `scripts/lib/yartsev_ch2_coupled_rods.py` также не используется
в RLB-1G. Его численный вектор `e_z=(0,0,+1)` задан в собственной системе,
где направление \(+\boldsymbol e_z\) названо положительным вниз. Модуль
содержит условия соединения, поэтому его `joint_basis` нельзя отождествлять
с новым глобальным базисом рисунка.

Аудит охватывал старый geometry helper и его тесты, production FEM transform
и его regression test, `docs/project_rules.md`, `scripts/README.md`,
`docs/theory/equations.tex` и существующие отчёты об исправлении отображения
Timoshenko. Старые файлы не изменены.

## 7. Реализация и проверка gate

Изолированная реализация находится в
`scripts/lib/reddy_inplane_geometry.py`. Она содержит только глобальные
векторы, базисы плеч и преобразования физических векторов. Модуль не
импортирует старый display helper. В реализации и тестах RLB-1G старый
helper используется только для comparison-проверки.

В файле `tests/test_reddy_inplane_geometry.py` выполнены следующие проверки:

- Проверяются \(\boldsymbol E_X\mathbin\times\boldsymbol E_Y=\boldsymbol E_Z\)
  и \(\boldsymbol k=-\boldsymbol E_Z\).

- Проверяются равенство
  \(\boldsymbol n_i=\boldsymbol k\mathbin\times\boldsymbol t_i\),
  ортонормированность и правая тройка Reddy.

- Проверяются пределы \(\beta=0^\circ\) и \(90^\circ\), а также явный
  случай \(30^\circ\).

- Проверяются точные отношения со старым Timoshenko display basis и отличие
  от знака поперечного поля EB.

- Выполняется реконструкция случайных векторов при фиксированном seed
  `20260825`.

- Проверяется равенство локальной и глобальной виртуальной работы.

Все девять целевых тестов прошли. Для геометрических сравнений задан
абсолютный допуск `5e-14`. Для суммы поступательной и вращательной частей
виртуальной работы задан допуск `2e-13`. Поэтому
координатный gate имеет статус `PASS`. Этот статус разрешает переход к
отдельному будущему этапу, но сам по себе не вводит матрицу соединения и не
запускает расчёт корней.
