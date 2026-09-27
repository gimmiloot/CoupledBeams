# KV: выбор симметрийно-редуцированного и полного solver-пути

2026-09-27, `main`, исходный HEAD `d1cd61fc7db816e81b6c77f7760d16a3966b7547`
(`Version 0.5.5.1`). Исходное дерево содержало 12 изменений D16/K17;
они сохранены. Дополнение — working-tree version `eb-kv-routing-v1`.
Разрешение — [D17](../memory/decisions.md#rlb-d17), результат —
[K18](../memory/knowledge.md#rlb-k18).

## 1. Причина изменения и область

[K17](inplane_kelvin_voigt_ac_diagnostics.md) связал отказы A/C с численной
согласованностью переноса и реакций. Его точная задача на одном плече
подтвердила прежние комплексные корни и физические формы. Теперь эта же
реализация используется по умолчанию для **точно одинаковых EB-плеч**.
Для неодинаковых плеч сохраняется полная задача. Это техническая регрессия
на известных решениях, без новых физических точек и без повторения
[внешней проверки](inplane_kelvin_voigt_literature_benchmarks.md).

Сохранены H/L/L/H, chi=.4, b=.20, h=.05, kappa=1 и reference K12:
A=.011, D=2.979166666666667e-6, m=.010000000000000002,
t_ref=69.2820323027551, k_theta=2.083333333333334e-6.
В EB 1/S=J=0; p=−alpha+i*omega_d, z=p*t_ref.
Секция получена штатной послойной редукцией.

## 2. Точные классы и выбор пути

При c_h=cos(beta/2), s_h=sin(beta/2) используются прежние условия K17:

| Класс | Перемещения | Усилия | Момент |
|---|---|---|---|
| eta=+1 | c_h*u−s_h*w=0 | s_h*N+c_h*Q=0 | M+2*(k_theta+c_theta*p)*psi=0 |
| eta=−1 | s_h*u+c_h*w=0 | c_h*N−s_h*Q=0 | M=0 |

Неизвестные — три реакции заделки одного плеча. Внешний результат содержит
обе физические формы: y2=eta*F*y1, F=diag(1,−1,−1,1,−1,−1), с общей
массовой нормировкой на 129 узлах каждого плеча. Класс eta=−1 не содержит
k_theta,c_theta: при заданном упругом eigenvalue выполняется проверка
подстановкой без вызова комплексного Newton. Возвращается
`EXACT_INACTIVE_BY_SYMMETRY`, без искусственного отрицательного Re(p).

Новый вход — [solver.solve_mode](../../scripts/lib/inplane_kelvin_voigt_solver.py):
`solver_path='auto'` выбирает `SYMMETRY_REDUCED` при точном равенстве
L,A,D,m и одного EB-закона обоих плеч. API фиксирует одинаковые внешние
заделки; другие опоры этим контрактом не поддерживаются. Явный mu!=0
исключает редукцию даже при округлившихся к одинаковым длинах.
Для остальных случаев выбирается `FULL_TWO_ARM`. Принудительный `'full'`
доступен всегда, `'reduced'` при неодинаковых плечах отклоняется.
Редуцированный вызов требует явного eta: sorted-номер его не заменяет.
При mu!=0 форма не проецируется, а дефект отражения не является gate.

## 3. Единственное изменение полного переноса

В историческом `kv.Provider.transfer(..., derivative=True)` оба результата
T,T_z получались из `expm_frechet`. Новый `FullProvider` вычисляет
**T=expm(X)** и только производную берёт как
`expm_frechet(X,X_z,compute_expm=False)`; X — прежняя масштабированная H*L.
Формула H_z, физический узел и знаки не изменены. SciPy 1.15.2 поддерживает
этот API. При запросе производной физический T теперь побайтно совпадает
с T без запроса производной в адресных тестах A/C.

Фиксированное row/column equilibration K17 задаётся один раз при predictor.
Если B_bal=diag(1/rows)*B_hat*diag(1/cols), обратное отображение есть
a=b/cols, r=diag(reaction_units)*a. Эта формула **не изменена**;
тесты проверяют её, реакции при x=0 и нормировку целой конструкции.
Full-формы по-прежнему восстанавливаются независимо в двух плечах
пошаговым expm K12. Относительные расхождения физических реакций full
с K17 после согласования общей фазы: A 2.78e-12, C 2.20e-12.

Исторический [inplane_kelvin_voigt.py](../../scripts/lib/inplane_kelvin_voigt.py)
оставлен побайтно неизменным для воспроизведения K12/K16/K17. Новый
production-вход оборачивает его корректор и диагностику. В
[K17 helper](../../scripts/lib/inplane_kelvin_voigt_symmetry_diagnostics.py)
выделен общий `AnalyticHalfProvider`; прежний `ClosedHalfProvider` сохранил
условный diagnostic API. Второго набора формул переноса не добавлено.
Теоретические основания — метки `eq:Hp`, `eq:Jmatrix`, `eq:scaledB`,
`eq:newton`, `eq:derivativeB`, `eq:modalidentity` в
[KV-заметке](inplane_kelvin_voigt_joint_theory.tex) и блоковая проверка K17.
Общего KV-блока в старом analytic code/equations.tex нет; иной моделью
он здесь не подменяется. Проверки обоих классов и производных сохранены.

## 4. Точность и контрольные результаты

До прохода зафиксировано regression agreement по комплексному p:
relative ≤1e-9, либо absolute ≤1e-10 при |p_reference|<1e-6.
Gates K12 неизменны: null/sigma/physical ≤1e-9, compatibility ≤1e-10,
energy ≤1e-6, MAC≥.95; raw second-singular separation ≥1e-8.
Это разные проверки, не единая оценка ошибки eigenvalue.

В таблице R — reduced, F — full; physical — максимум шести условий узла.
Для R это также невязка поднятой полной формы. Все корни взяты с полной
точностью из источников и проверены текущей матрицей. Корректор принял
начальные значения без обновлений. **Нулевая разность ниже отражает reuse
корня, а не независимую оценку его точности.**

| Контроль | d | Relative p difference к источнику, R/F | Physical R/F | Energy R/F | Статус R/F |
|---|---:|---:|---:|---:|---|
| K12 ACTIVE, 5° | .0005423772776686932 | 0 / 0 | 4.38e-16 / 1.35e-15 | 9.89e-10 / 9.89e-10 | принято / принято |
| K12 INACTIVE, 5° | то же | 7.58e-32 / — | 1.38e-15 / — | 2.40e-12 / — | exact reuse / — |
| B, 45°, sorted_02 | .001 | 0 / 0 | 5.36e-15 / 4.77e-15 | 1.85e-9 / 1.85e-9 | принято / принято |
| B, 45°, sorted_02 | .005 | 0 / 0 | 2.72e-14 / 1.98e-14 | 1.84e-9 / 1.84e-9 | принято / принято |
| A, 0°, sorted_05 | .001 | 0 / 0 | 6.39e-10 / 4.27e-10 | 1.03e-8 / 1.03e-8 | принято / принято |
| C, 75°, sorted_05 | .001 | 0 / 0 | 3.30e-11 / 2.51e-11 | 4.92e-8 / 4.92e-8 | принято / rank qualification |
| K11, 5°, mu=.01, p01 | 0 | — / 0 | — / 3.53e-15 | — / 9.90e-10 | — / принято |

Общий максимум null residual 9.71e-16, sigma ratio 2.80e-16,
conjugate residual 9.76e-16; минимум MAC к соответствующему упругому seed
0.999982608. У K11 используются L1=.99,L2=1.01 и две независимые формы;
это только упругий asymmetric control, не проверка вязкого продолжения
неодинаковых плеч.

## 5. A/C и оставшаяся квалификация

| Точка | z, сохранённый K16=K17=новый R=новый F | Physical K16 full | Новый full |
|---|---|---:|---:|
| A/.001 | −.17856506546575265+i*76.57755928391374 | 2.79e-7 | 4.27e-10 |
| C/.001 | −.006633333405021945+i*100.01655751370787 | 3.90e-8 | 2.51e-11 |

Оба новых full-восстановления проходят прежний physical gate 1e-9.
У C сохраняется `POSSIBLE_MULTIPLICITY`: raw full second-singular ratio
8.61e-9 ниже 1e-8. Поэтому CSV оставляет `accepted=false` этой строки.
Root equation и form recovery проходят; отдельно записан `rank_status=QUALIFIED`.
Это **квалификация raw rank gate, а не оставшийся отказ физической формы**.
Кратность здесь заново не исследовалась; объяснение двумя блоками остаётся
результатом [K17](inplane_kelvin_voigt_ac_diagnostics.md).

Итоги в пределах данного контрольного набора:
`SYMMETRY_REDUCED_PRODUCTION_PASS` (6/6) и
`FULL_TWO_ARM_PRODUCTION_PASS_WITH_HIGH_MODE_QUALIFICATION`
(5/6 безусловно принятых строк; C с указанным raw rank flag).
Они обозначают пригодность путей для проверенных условий, не общий PASS
всей KV-модели. Старые K16 rejected/NOT_RUN не переименованы.
После одного изменения full-пути дополнительное уточнение не выполнялось.

## 6. Данные, проверки и остановка

[CSV](../../results/laminated_beams/inplane_kelvin_voigt_solver_architecture/solver_regression.csv),
[diagnostics](../../results/laminated_beams/inplane_kelvin_voigt_solver_architecture/diagnostics.json),
[manifest](../../results/laminated_beams/inplane_kelvin_voigt_solver_architecture/run_manifest.json)
и `control_shapes.npz` лежат в отдельном локальном results-каталоге,
который по правилам репозитория исключён из Git. CSV содержит 12 контрольных
строк, provenance и раздельные статусы; JSON — все невязки и затраты.
Новый entry point — техническая регрессия двух solver-путей с другим
контрактом результата, а не preset физического sweep K16 или условного
диагностического триггера K17. Он переиспользует их helpers и сохраняет
старые runners для воспроизведения; произвольные параметры не принимает.
SHA-256 36 прежних научных файлов совпадают с исходными; старые отчёты,
теория и результаты не переписаны.

Использован `D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe`,
Python 3.12.4, NumPy 2.1.3, SciPy 1.15.2. Один проход: 0.63086 s,
B=83, B_z=47 (full 54/30, half 29/17), direct expm=108,
Frechet=60, shape expm=12, analytic transfers=803, восстановлений=12.
11 вызовов корректора дали 0 Newton updates; eta=−1 имеет 0 таких вызовов.
Новых физических roots/targets — 0; build-equivalents=1173, по каждой
точке соблюдён штатный защитный предел 2000. Тестовые вызовы учитываются
отдельно в manifest. Missing-only не выполняет root/matrix/form calls
и не меняет сохранённые файлы.

Пройдены 69 различных целевых тестов (1.70+0.69+0.60 s): 19 новых routing,
18 K17, 13 K16, 17 EB/общих K12 и два EB-контроля длин/массы K11.
После добавления отдельных root/form/rank полей и проверки eta повторены
только 19 routing-тестов (1.25 s); новых спектральных вычислений не было.
RLB-тесты не запускались. Команды из корня репозитория:

```powershell
$py = 'D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe'
& $py -B scripts/analysis/laminated_beams/check_inplane_kelvin_voigt_solver_architecture.py --compute
& $py -B -m pytest -q -s -p no:cacheprovider tests/test_inplane_kelvin_voigt_solver_routing.py tests/test_inplane_kelvin_voigt_symmetry_diagnostics.py tests/test_inplane_kelvin_voigt_targeted_weak_damping.py
& $py -B -m pytest -q -s -p no:cacheprovider tests/test_inplane_kelvin_voigt.py -k 'not RLB and not rlb and not unstable_trial'
& $py -B -m pytest -q -p no:cacheprovider 'tests/test_inplane_spring_robustness.py::test_actual_lengths_in_transfer_and_shapes[EB]' 'tests/test_inplane_spring_robustness.py::test_mass_contains_each_length_and_only_rlb_rotation[EB]'
```

**Место остановки:** техническая стабилизация завершена с локальной
квалификацией C full. A/.005 и C/.005 — `NOT_RUN`; новые d, beta,
RLB roots и физические исследования — 0. High precision, новые алгоритмы
и каскады уточнений не применялись. Выбор следующего физического
продолжения требует отдельного решения пользователя.
