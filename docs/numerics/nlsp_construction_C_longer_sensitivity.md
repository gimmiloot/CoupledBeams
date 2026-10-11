# Конструкция C: один период внеплоскостного изгиба и вычислительная стоимость

**Итог: `CONSTRUCTION_C_LONGER_EXPLORATORY_P12`.** Продолжение до одного
низшего внеплоскостного периода дало новый sampled maximum1.07026572
в исходном Sdisp3. Уточнение p16 не поместилось в оставшийся бюджет;
новый интервал имеет numerical qualification PARTIAL. Сильного физического
выхода из плоскости в рассчитанном p12 отклике не наблюдается.

Успешный bundle:
`results/nlsp_construction_C_longer_sensitivity/6f504c8b30b27df3/`.
Все старые trajectories и нормировка приt=0 сохранены, первые2Tref не
пересчитывались. Отдельно разрешённый повтор относится только к новому
planar interval после output-layer failure.

## Цель и неизменная постановка

Это отдельное ограниченное продолжение [модального опыта C](nlsp_construction_C_modal_sensitivity.md).
На прежнем интервале двух reference-периодов получено
`Gperp,disp^(3)=1.0530286991` при p16 и `1.0530316774` при p12.
Новый горизонт выбран по низшей сохранённой внеплоскостной частоте:
`omega_OUT1=.07773966990954037`, `t_end=80.82341119393537`,
то есть `4.031266620508157 Tref`. Здесь `Tref=20.049135619749023`.
Первоначальный интервал охватывал около половины этого внеплоскостного периода.
Период линейной формы не объявляется периодом нелинейного движения.

Остаются прежние beta30deg, L1=.8, L2=1.2, оба сечения b=.05/h=.10,
E=rho1, nu=.3, kappa5/6. Полная плоская собственная форма возбуждения
имеет max|w|=.005 и нулевые начальные скорости. Статической нагрузки нет.
Принятые V0/V4, нелинейная инерция, переменная масса, RHS/Jacobian,
Shen basis, value clamps и семь условий идеального узла, включая c1=c2,
не изменяются. Шесть физических kernels совпадают с историческими hashes.
Сохраняется [семиполевая теория](../theory/weakly_nonlinear_spatial_rod.md).

## Источники и настоящее продолжение

Основной parent:
`results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/`,
manifest SHA256 `b54b9271260b55841b260692d88f8802f430d780fe21a5d91fe049f91705f3dd`.
Его manifest, артефакты и родительские источники проверяются до новых ODE.
Исторический статус C_MODAL_SENSITIVITY_RESOLVED не переписывается.

Сохранены все координаты и скорости плоской базы при t0=2Tref,
конечная tangent matrix трёх направлений, исходный basis и фиксированная
reference-energy metric. Полное состояние p12 содержит322 компоненты,
плоская подсистема184, OUT response138; при p16 —434/248/186.
Отсутствует внутреннее adaptive state интегратора. Поэтому используется
новый integrator от сохранённого физического состояния при абсолютном t0,
а не побайтовое продолжение прежнего solver path.

Исходный basis Sdisp3 содержит displacement starts OUT forms1/2/5.
Tangent matrix при restart не ортонормируется и не заменяется исходным
basis. Gain по-прежнему относится к исходным состояниям при t=0 и той же
фиксированной E0 norm. Все OUT output DOFs сохраняются; ограничение
касается только трёх начальных направлений. Исторический prefix переносится
без интерполяции и без повторного решения первых двух периодов.

## Решения, зафиксированные до результатов

Основной новый расчёт — p12. Для p16 заранее задан computational trigger:
новый глобальный maximum превосходит прежний более чем на .005 либо
наблюдается выраженный растущий участок в последние .2Tref. Для второго
варианта требуется положительный least-squares slope, рост endpoint относительно
начала окна больше исторической spatial curve difference9.457687954e-4
и endpoint не ниже maximum этого окна более чем на ту же величину.
Это правила решения о вычислительных затратах, не физические пороги.
Историческая spatial difference не является error bound нового интервала.

Прогноз по сохранённым costs и pointwise acceleration OUT-оператора:
p12 base445.85s+tangent318.51s; p16 base723.27s+tangent635.58s.
С резервом15% четыре продолжения оценены в2441.69s при budget2700s.
Прогноз не является гарантией. После каждого completed prefix оставшийся
расчёт снова проходит resource gate. Допуски, componentwise atol и max_step
переиспользуются из actual records. Новые FEM, полный finite-perturbation
control, другие геометрии/амплитуды, Floquet и maps запрещены.

## Отказ архивирования и отдельное разрешение на повтор

Первый p12 continuation:
`results/nlsp_construction_C_longer_sensitivity/2bd3ecdc5caa2f28/`.
Он сохранил2033 полных sampled states на40.09827123949805...80.82341119393537.
Начальное состояние совпадает с parent endpoint; все состояния конечны,
OUT=0, energy drift нового интервала1.42122e-11, относительно исходного
initial energy2.85856e-11. Энергетика не заменяет temporal certification.

После вычисления trajectory старый save_dense_records отклонил native
cubic coefficients: его контракт требовал t_old[0]=0. Это техническая
ошибка orchestration preflight данного продолжения. Отказ физических
уравнений или Radau не установлен. Dense coefficients и native counter
archive не сохранены; sampled states ими не подменяются. Первый attempt
имеет статус BLOCKED_BY_DENSE_ARCHIVE_CONTRACT. Его manifest SHA256:
`73f4d225d2cb9daefc2b12b16667b355969baa11a601a17b3c75fe6b788dcf60`.

Writer исправлен добавлением явно заданного t_start при прежнем default0;
границы и непрерывность native intervals по-прежнему проверяются.
Пользователь отдельно разрешил ровно один технический повтор p12.
Он начинается с того же t0, без replay первых двух периодов. Первый
attempt сохранён, в общем budget ему консервативно начислены450s;
filesystem elapsed до trajectory archive составлял443.621s. Это измерение
включает расчёт sampled energy и запись NPZ; точный native timing record
первого attempt отсутствует. Native steps/RHS/LU counters первого запуска
не восстанавливаются по предположениям.

## Аудит вычислительной стоимости

| Historical 2Tref component | p12 | p16 |
|---|---:|---:|
| Planar base, seconds | 438.984 | 712.140 |
| Accepted planar steps | 38030 | 54275 |
| Planar RHS evaluations | 266776 | 382597 |
| Planar analytic Jacobian evaluations | 2 | 2 |
| Radau LU factorizations | 376 | 1824 |
| Matrix tangent, seconds | 415.753 | 962.528 |
| Accepted tangent steps | 6028 | 9286 |
| Tangent RHS evaluations | 86582 | 145346 |
| Full Jacobian evaluations in tangent operator cache | 79867 | 133734 |

В успешно архивированном новом p12 base на дополнительном интервале
получены39044 accepted steps,273790 RHS evaluations,2 analytic Jacobian
calls и326 LU factorizations. Wall time478.240s включает сохранение dense
polynomials. Actual RHS time327.589s (68.50%), отдельные safety checks
116.388s (24.34%); LU factorization .160s и solve12.568s (вместе2.66%).
Jacobian calls суммарно .0116s. Внутри RHS potential route занимает125.811s,
inertial terms73.640s, mass solve121.652s, из них rotational mass assembly
и factorization113.408s. Эти вложенные timings не суммируются с RHS;
сборка массы и её Cholesky не разделены как production wall fractions.
Остаток включает housekeeping, dense/sample output и прочую работу,
а не исключительно disk I/O. Во время base также выполнялись read-only
tests; timings являются наблюдаемыми wall costs этого запуска.

Writer remediation проверена на одинаковом continuation: times, полные
state arrays и energy arrays двух solves побитово совпали. У первого
не сохранены native adaptive states, поэтому их побитовое сравнение
не заявляется. Прежний0...2Tref prefix не пересчитывался ни в одном solve.

FFT сохранённых Radau polynomials на равномерных grids dt=.001/.0005
показывает основной изгиб около.0491 cycles/time, тогда как quarter-span
c velocity имеет dominant frequencies5.84/5.79 cycles/time у двух плеч.
В этих точках около97.3/98.2% Hann-window spectral power скоростей c
лежит выше1 cycle/time. Это pointwise field diagnostic, не новый spectrum
или eigenanalysis. Он подтверждает наличие быстрых внутренних компонент,
но не устанавливает оптимальные solver tolerances. При действующей tight
policy accepted base increments существенно меньше cap .00720939.

В плоском расчёте главная повторяющаяся работа связана с RHS, переменной
массой, инерционными членами и safety checks. В tangent route основной
кандидат для ускорения — многократная сборка Jacobian. DOP853 не имеет
Radau LU; его njev=0 не означает отсутствия analytic Jacobian внутри RHS.
Basis, quadrature и compiled expressions создаются заранее, а не выводятся
символически на каждом шаге. Сохранённые скорости — state variables.

Pointwise profiling выполнен без новых ODE на трёх сохранённых состояниях.
Вложенные cumulative timings не складываются как независимые wall fractions.
При p16 gradient+Hessian potential evaluation занимает около1.81ms,
gradient около.090ms, energy value около.019ms. Это отдельный benchmark,
не доли production runtime. В microbenchmark assembled mass cache
инвалидирован, но local arm caches могут оставаться warm. Их реальная
стоимость при меняющемся состоянии требует production measurements.

Radau linear systems отдельно оценены на сохранённом Jacobian и среднем
историческом шаге. При p16 median real/complex LU factorization составляет
.377/1.115ms, solve .0165/.0991ms. Эти показатели не определяют точную
долю LU в историческом runtime; times per solve и число Newton iterations
могут меняться. В продолжении добавлены observational timers настоящих
LU calls, без изменения аргументов и результатов.

## Единственная узкая оптимизация

Scoped OUT evaluator использует прежние local force_tangents, полное
assembled acceleration и тот же variable-mass solve. Он собирает только
нужные OUT input columns, сохраняя все global rows до mass solve.
Локальное ускорение из isolated-arm equations не подставляется.
Ни одна из трёх propagated columns не проецируется обратно на initial space.

Сравнение на шести сохранённых состояниях каждого p даёт различия лишь
на уровне float64 rounding. Для zero state проверена identity
`A.T E0+E0 A=0`, relative residual1.88e-14/2.96e-14 приp12/p16.
Median pointwise speedup1.326/1.538; это ускорение одного operator evaluation,
не обещание такой же скорости всего эксперимента. Физические kernels,
математические уравнения и tolerances не изменены. Серьёзный refactor
или reduced-order model в этом этапе не создаются.

На40 одинаковых saved states нового p12 интервала, с cold mass/potential
caches, median full/OUT-only operator times составили5.555/4.407ms,
speedup1.2605. Ни одной дополнительной ODE для этого сравнения нет.
В production tangent фактически391.612s из397.207s, то есть98.59%,
заняла OUT operator evaluation; выполнено81316 cache misses. Для этого
route ускорение operator важнее ускорения записи небольших sampled arrays.
Whole-IVP speedup на одинаковом интервале не измерен: исходный tangent
route заново не запускался. Разность costs старого и нового интервалов
не интерпретируется как чистый эффект оптимизации.

## Новый отклик и numerical qualification

| Quantity | Parent p12, 0...2Tref | Continued p12, 0...T_OUT1 |
|---|---:|---:|
| Sampled max Gdisp3 | 1.05303167743216 | 1.07026571974363 |
| Time of global max / Tref | .746 | 2.747 |
| Physical time of global max | 14.9566551723 | 55.0749755475 |
| Maximum only in added interval | — | 1.07026571974363 |
| Global maximum increase | — | .01723404231148 |
| Gain at final time | — | 1.04702635129771 |
| Max physical v, own fixed optimum | 3.1732215326e-5 | 3.2756578851e-5 |
| Rv = max abs v / .005 | .006346443065 | .006551315770 |
| Max abs Phi, rad | 3.7100578558e-5 | 4.0963109433e-5 |
| Max abs psi, rad | 5.1046276218e-5 | 5.4714200113e-5 |

Physical maxima — sampled maxima по общей201-point coordinate grid каждого
плеча. Две колонки используют собственные fixed optimal исходные комбинации,
а не одинаковый initial direction. Одно общее масштабирование задаёт
max{|v(0)|,Lref|Phi(0)|,Lref|psi(0)|}=5e-5, Lref1; все initial velocities0.
Новый max v составляет .65513% основной плоской amplitude. Это касательный
прогноз для малого конечного масштаба; новый full finite-perturbation solve
не выполнялся и данным интервалом не подтверждён.

Для нового fixed optimal input prefix/added physical maxima равны:
v3.18710e-5/3.27566e-5, Phi3.58273e-5/4.09631e-5,
psi5.11043e-5/5.47142e-5. Такое сравнение уже сохраняет один initial input
на обоих интервалах. Оно отделено от вновь оптимизируемой по времени gain curve.

Новая initial combination в неизменном p12 E0 basis:
(.9987011120,-.00794795245,.05032811197), прежняя:
(.9969334430,-.01352238393,.07707694454). Направление остаётся преимущественно
первой изгибной формой; reference-energy angle1.56884deg. Sign eigenvector
выбран по положительному scalar product с прежним optimum, без phase fitting.
At global peak83.2034% квадрата E0 norm приходится на скорости. Поэтому
7.03% norm gain не равен7.03% увеличения прогиба.

В последние .2Tref gain увеличился1.02725225...1.04702635, endpoint rise
.01977410, regression slope+.10347110 per Tref. Это растущий локальный
участок, а не новый global record у конца или доказательство secular growth.
Оба computational triggers для p16 выполнены. Однако1562.680s forecast
with reserve превышает1374.553s remaining budget. P16 не запускался.
Никакой другой пространственный level или temporal refinement не введён.

Новый maximum increase в18.22 раза больше historical sampled p12/p16
gain-curve difference9.45769e-4. Это только сопоставление с прежним indicator:
он не ограничивает ошибку на новом интервале. Direct p12/p16 comparison и
новый tight temporal control имеют NOT_RUN. Parent optimal temporal change
1.93e-9 и parent full consistency относятся лишь к0...2Tref. Поэтому
новый p12 gain не объявляется точным эталоном. Исторический prefix gain
воспроизведён с maximum difference0; restart sampled states непрерывны,
но independent restart-error solve не выполнен.

## Фактическое исполнение

| New production attempt | Duration, s | Accepted steps | RHS | Jacobian / LU |
|---|---:|---:|---:|---:|
| First p12 base; output archive failed | ~443.62; budget charge450 | NOT_AVAILABLE | NOT_AVAILABLE | NOT_AVAILABLE |
| Authorized p12 base repeat | 478.24018 | 39044 | 273790 | 2 / 326 |
| p12 matrix tangent continuation | 397.20729 | 6130 | 88154 | solver njev0 / LU0; OUT operator81316 |
| p16 base/tangent continuation | NOT_RUN | — | — | — |

Всего3 новых IVPs:2 planar (включая отдельно разрешённый output repeat)
и1 matrix tangent с тремя columns. Full seven-field perturbation/FEM/static
Newton/new spectra/maps0. Budget charge1325.447s из2700s, включая450s
первого attempt; successful production875.447s. Process peak working set
.641045GiB поOS observer при ceiling4GiB. Jobs последовательны, один поток.
Observed p12 gain/physical postprocessing1.864s; figures не имели отдельного
production timer. Actual native end обоих успешных jobs80.82341119393537.
Base steps .000764454... .001371786, tangent .004085551... .007209393.
Скорости — сохранённые state variables; native timestamps, dense polynomials
и actual final states находятся в соответствующих case folders.

Финальная проверка native/sample endpoint difference даёт5.29e-23 для
base и0 для tangent. Original basis иtangent restart matrix сохранены
точно; это проверка архивирования, а не independent restart-error bound.
200 targeted/related tests PASS (28.50s),27 новых navigation links проверены,
git diff --check PASS. Cached compute/report-only/plot-only с запрещёнными
scientific entry points выполняются с0calls и неизменным manifest. Чистое
повторное построение трёх figures из saved arrays занимает1.707s; это
отдельный postprocessing measurement без ODE. README/CHANGELOG, journal,
navigation иD23/K24 обновлены; старыеD/K сохранены побайтовым prefix.

Суммарный observed prefix+extension cost одной C: p12 base917.224s,
tangent812.960s, вместе1730.184s на4.031Tref без первого неудачного output
attempt. Это около7.15min наTref и28.84min на одну точку данного горизонта.
При p16 прежние2Tref стоили1674.668s (13.96min/Tref); actual continuation
cost отсутствует. Эти значения не гарантируют стоимость другой геометрии.
Сто p12 точек при таком последовательном rate потребовали бы около48h
только на две ODE-части, без новых controls и остальных стадий. Для регулярных
карт нужна предварительная оптимизация повторяемых вычислений.

Сохранены три PDF+PNG:
[gain](../../results/nlsp_construction_C_longer_sensitivity/6f504c8b30b27df3/C_longer_restricted_gain.pdf),
[physical OUT response](../../results/nlsp_construction_C_longer_sensitivity/6f504c8b30b27df3/C_longer_physical_OUT_response.pdf),
[numerical indicator and cost](../../results/nlsp_construction_C_longer_sensitivity/6f504c8b30b27df3/C_longer_numerical_cost.pdf).
Grey zero-base line — conservative theoretical reference, независимо
проверенный ранее; новая zero-base ODE на extended horizon не выполнялась.
Raw NPZ иCSV, including p12_gain.csv/computational_cost.csv, сохранены
рядом с [manifest](../../results/nlsp_construction_C_longer_sensitivity/6f504c8b30b27df3/manifest.json).

## Воспроизведение и ограничения

Для будущих повторяемых геометрических расчётов наиболее перспективно
вычисление Jacobian-vector products для всей matrix из трёх направлений,
без сборки ненужных входных столбцов. Оно должно сохранить M_q*a,
общую узловую кинематику и полный OUT выход. Текущий scoped evaluator
лишь исключает PL columns; matrix-free propagation не реализована.
Повторный расчёт unused potential gradient внутри local Hessian route
также заслуживает отдельного измерения; вторую оптимизацию здесь не вводим.

Reduced-order representation плоской базы могла бы сократить многократные
RHS evaluations. Однако простое исключение c или быстрых внутренних форм
не обосновано: поле c и velocities имеют незакрытую spatial qualification.
Будущая редукция потребует контроля не только w и energy, но и tangent
gain вдоль reduced base. Ни новая модель, ни новые amplitudes, ни карта
в настоящем продолжении не строятся. Ослабление tolerances не считается
безусловно безопасным ускорением после известного numerical-floor отказа.

Новый runnable script не нужен по Script Proliferation Control:
существующий CLI получает scoped preset --construction-c-longer,
переиспользуя прежние Radau/DOP853, assembly, SavedPlane и gain recovery.
Новые [config](../../data/input/nlsp_construction_C_longer_sensitivity.json)
и отдельно разрешённый [output repeat](../../data/input/nlsp_construction_C_longer_sensitivity_output_retry.json)
не ослабляют исторические failed-attempt guards.

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
& 'D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe' scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --construction-c-longer --compute --config data/input/nlsp_construction_C_longer_sensitivity_output_retry.json
```

Report/plot replay и cached compute читают сохранённые arrays и manifests,
не запускают ODE, static Newton, eigenanalysis, Gmsh или CalculiX.
Новый интервал сам по себе не получает independent temporal certification.
Историческая all-field spatial qualification, особенно c и velocities,
остаётся PARTIAL. Полная семиполевая finite-perturbation проверка здесь
запрещена; её прежний успех относится к0...2Tref.

DEFERRED_TWO_PARAMETER_ACTION_EXPANSION, физическая применимость ideal
point joint/c1=c2, permanent imperfections, полный amplification operator
и устойчивость периодических движений сохраняются как отложенные задачи.
Конечновременный restricted gain не доказывает асимптотической устойчивости
или неустойчивости. После этого C контроля дальнейшие расчёты при данной
геометрии и амплитуде останавливаются.
