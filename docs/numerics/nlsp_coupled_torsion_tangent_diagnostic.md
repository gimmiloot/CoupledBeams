# Адресная диагностика крутильного tangent/full контроля

**Итог: `NUMERICAL_FLOOR_IDENTIFIED`.** Историческая немонотонность вызвана
преимущественно численной погрешностью full nonlinear интегрирования при
фиксированной tolerance policy, усиливающейся после деления на delta.
Уточнение tangent изменяет отклик значительно меньше. При уточнении full
восстановлено наблюдаемое отношение ошибок около4 для delta→delta/2;
дополнительный temporal level подтверждает устойчивость малого остатка.
Исторический `BLOCKED_BY_TANGENT_CHECK` не переписан. Новый вывод ограничен
крутильным направлением A,p16,.1Tref; прежний пилот A/B не возобновлён.

## Объект и сохранённые источники

Проверяется один отказавший крутильный контроль конструкции A, p16,
на `0…0.1Tref`, где `Tref=35.40765544870329` и конечное время
`3.540765544870329`. Это отдельное продолжение
[пилота двух стержней](nlsp_coupled_planar_sensitivity_pilot.md).
Геометрия B, усиления, параметрические карты и устойчивость здесь не исследуются.

Исторический bundle `d024ee84989ce6b5` остаётся неизменным. Его manifest:
`7a149b76f64165a300e27d33f7cb1f2876a5f3094d5193d72f0c591c57564d80`.
Исторический статус `BLOCKED_BY_TANGENT_CHECK` и прежние thresholds сохраняются.
Новая диагностика находится в
`results/nlsp_coupled_torsion_tangent_diagnostic/a10f1110ffbce7ca/`.
Основание отдельного разрешения и границы вычислений:
[NLSP-D20](../memory/decisions.md#nlsp-d20).

Сохраняются quartic action, variable mass, nonlinear inertia, полные
семиполевые RHS/Jacobian, материал, геометрия, endpoint lifting, basis и
идеальный жёсткий узел. Новых статических задач или собственных форм нет.
Все начальные состояния загружаются непосредственно из parent bundle.

## Воспроизведение исторического отказа

Крутильное направление — column 2 сохранённой E0-orthonormal basis.
Начальные скорости равны нулю. Обе амплитуды заданы в глобальных coupled
coordinates: `X_delta(0)=X_base(0)+delta*z0`. Проверка воспроизводит этот
вектор побитово; добавки к плоским полям отсутствуют. Joint translation,
finite orientation и contraction constraints сохраняются.

`delta0=7.958509826438746e-7`, `delta0/2=3.979254913219373e-7`.
Сравниваются прямые solver state variables `(v,Phi,psi)` и их скорости.
Скорости не получаются дифференцированием восстановленных перемещений.
Матрица E0 norm совпадает с текущими K_perp,0 и M_perp,0 побитово.
Relative denominator `1.0001689730987355` — прежний фиксированный максимум tangent E0 norm,
с нижним ограничением 1; новые кривые по нему не перенормируются.

| Величина | delta0 | delta0/2 |
|---|---:|---:|
| Historical relative mismatch | `6.895227577955444e-8` | `8.871387187778917e-8` |
| Physical time at sampled maximum | `3.540765544870329` | `3.5053578894216253` |
| Coordinate E0 norm at total maximum | `3.42619e-8` | `7.53222e-9` |
| Velocity E0 norm at total maximum | `5.98510e-8` | `8.84086e-8` |
| Velocity fraction of squared norm | `75.3181%` | `99.2794%` |
| Error E0 norm before division by delta | `5.48850e-14` | `3.53075e-14` |

Уменьшение delta уменьшает абсолютную разность, но деление на меньшую
амплитуду увеличивает нормированную разность. Само это наблюдение не
устанавливает её численный источник.

В максимуме для меньшего delta вклад `psi_t` составляет `92.3634%` квадрата
полной error norm, `Phi_t` — `6.4625%`, `v_t` — `0.4535%`.
Следовательно, остаток преимущественно связан с изгибным поворотом,
хотя начальное направление преимущественно крутильное.
На первое/второе плечо приходится `51.3853%/48.6147%` квадрата нормы.

Проекция на шесть сохранённых low-mode initial states содержит лишь
`.3678%` error energy. Остальная часть не идентифицирована по eigenmodes:
полная линейная база не сохранена, нового eigenfrequency scan нет.
Rayleigh frequency velocity error равна `188.003`, существенно выше
частоты исходного torsional seed `1.54813`. Это характеристика спектрального
масштаба остатка, а не идентификация отдельных высокочастотных крутильных мод.
Десять крупнейших вычислительных velocity DOFs сохраняются в JSON;
ни один компонент из основной нормы не исключается.
Они содержат94.21% total squared error при delta/2. Крупнейшая whitened
velocity coordinate214 даёт70.40%; её связь с raw узловым alphaY не означает
локализацию ошибки в физическом узле. Mass whitening смешивает lifting и
внутренние функции; raw group contributions и восстановленные поля нужны
для физического толкования. Чисто крутильная кинетическая часть мала.

Для internal/joint и displacement/rotation groups сохранены diagonal
projected norms, cross terms и signed bilinear shares. Из-за coupling
последние могут быть отрицательными: это алгебраическое распределение полной
квадратичной формы, не независимые положительные энергии групп.
У joint share при меньшем delta значение `-.4431%`; internal share `100.4431%`.
Это не отрицательная энергия идеального узла.

## Время, dense output и точность исходных данных

Обе full trajectories имеют точно те же 101 comparison timestamps,
от 0 до `.1Tref`. Native accepted-step grids различаются.
Full state оценивается прежними кубическими Radau polynomials.
Повторное чтение dense output на этих же timestamps даёт normalized
разность `1.861e-17/4.247e-17`, значительно меньше historical mismatch.
Time/phase shifts, fitting или новый interpolation method не применяются.

| Full run | Accepted steps | Actual min step | Actual max step | RHS / Jac / LU |
|---|---:|---:|---:|---:|
| delta0 | 1829 | `.0011154516291114192` | `.0023282689967417947` | 12817 / 2 / 6 |
| delta0/2 | 1680 | `.0009002336743471062` | `.0026850014537608204` | 11774 / 2 / 6 |

Old tangent: 5104 accepted steps, 36954 RHS, 2 Jacobian, 886 LU.
Его accepted timestamps и dense polynomials не сохранены (`NOT_AVAILABLE`);
actual min/max step и полноценное old native-time сравнение не восстанавливаются
догадкой. Сохранены sampled states и точные settings. Output dense calls
на 101-точечной сетке — 100 по действующему алгоритму.

Old `rtol=1e-10`, `max_step=.0036046964938503193`. Полные componentwise
atol сохраняются в source JSON. Norm positive atol vector после деления
на delta равна `1.167e-6/2.334e-6`. Это диагностический масштаб policy,
не достигнутая ошибка, глобальная оценка или доказательство причины failure.

Pointwise symmetric RHS differences в четырёх сохранённых base states
при delta factors16,8,4,2,1,.5 показывают порядок около2. После исключения
ведущего delta² члена relative Jacobian difference не превышает `6.31e-14`.
Проверены OUT, velocity/acceleration и raw joint rows. Initial mapping,
normalization и ошибка повторного dense evaluation не объясняют historical
failure. Ошибка динамического интегрирования ещё требует независимого контроля.

Новая read-only diagnostic implementation первоначально использовала
неподходящую advanced indexing форму NumPy. Она исправлена до первого ODE;
original code/error/manifest сохранены в `read_only_implementation_error`.
Исторический metric от этой ошибки не зависел и воспроизводился точно.

## Замороженное временное уточнение

До новых integration results сохранён `frozen_refinement.json`:

- `rtol=1e-12`;
- full и tangent componentwise atol уменьшаются в100 раз;
- `max_step=.0018023482469251596`, вдвое меньше old;
- прежние initial states, geometry A, p16, horizon и relative denominator;
- прежняя сохранённая base trajectory, без нового equilibrium/base solve.

Вычисляется один torsional tangent column, с сохранением всех186 OUT DOFs,
и две full434component trajectories. Линейные column equations независимы;
выбор исходного направления не обрезает пространство отклика. Изменение
batch size влияет на численный error control и явно входит в диагностику.
Прогноз648.253s использует даже прежнюю стоимость всех шести tangent columns
и запас1.5; budget1200s,600s/run,maximum6 new ODE calls. Прогноз не гарантирует
стоимость. После настоящего solver failure автоматического retry нет.

## Независимое временное уточнение

Все три primary runs и один conditional run дошли до `3.540765544870329`.
Основная таблица использует точно прежние101 timestamps и один исторический
denominator. Для обоих delta primary settings одинаковы:

| Сочетание | delta0 | delta0/2 |
|---|---:|---:|
| old full / old tangent | `6.89522758e-8` | `8.87138719e-8` |
| old full / tight tangent | `6.92912345e-8` | `8.91000765e-8` |
| tight full / old tangent | `3.73873846e-8` | `9.34917832e-9` |
| tight full / tight tangent | `3.73876652e-8` | `9.34777630e-9` |
| D_full: normalized change old→tight | `5.81115727e-8` | `8.84864915e-8` |
| D_tangent: change old→tight | `4.24537180e-10` | `4.24537180e-10` |

Изменение tangent составляет менее.5% historical mismatch для меньшего delta.
Изменение full сопоставимо со всем старым mismatch. Одно лишь уточнение
tangent не устраняет немонотонность; уточнение full устраняет её даже при
старом tangent. Эти независимые сочетания определяют основной источник
failure, без предположения о шуме до получения новых trajectories.

Primary errors дают отношение `3.999632` и наблюдаемый порядок
`1.9998672631308596`. Это finite-interval, sampled comparison двух амплитуд.
Условное исключение ведущего delta² члена оставляет offset `1.569e-10`;
оно служит дополнительной диагностикой и не заменяет исходные данные.

После Stage2 выбран только вариант B. Наибольшая sensitivity была у
full trajectory delta0/2, поэтому для неё выполнен один extra level:
`rtol=1e-13`, componentwise atol=primary*0.1,
`max_step=.0009011741234625798`. Решение и стоимость заморожены в
`conditional_decision.json` до extra results. Прогноз этой части239.145s,
суммарно783.359s; negative-delta variant не выполнялся.

Extra level меняет нормированный full response на `1.51939605e-10`,
то есть примерно1.63% primary mismatch. Его mismatch с tight tangent —
`9.34625980e-9`. Малый остаток устойчив относительно этого дополнительного
temporal control. Extra point относится только к delta0/2; равномерной
extra-level pair для обеих амплитуд нет.

## Что установлено о причине отказа

В интервале старого sampled maximum доля atol в Radau error-control scale,
взвешенная вкладом velocity error, составляет `99.9558%/99.9757%`.
Использована фактическая формула установленного SciPy:
`atol+rtol*max(abs(y_start),abs(y_end))`. На последнем old increment конец
восстановлен dense polynomial; в остальных случаях сохранён как следующий
accepted start. Эта доля характеризует scale, не достигнутую ошибку.

Full accepted max steps меньше заданных max_step во всех old/tight/extra
runs. Поэтому step cap непосредственно не ограничивал эти шаги; их величины
определялись адаптивным tolerance control. Совместное уточнение rtol/atol
уменьшило actual steps и разрешило чрезвычайно малый OUT-отклик.
Влияние rtol и atol не разделено отдельными jobs. Нельзя приписать весь
эффект только изменению одного atol или только размера шага.

Уменьшение исходного возмущения при неизменном atol не обеспечивало
пропорционального уменьшения абсолютного numerical error. После деления
на delta этот остаток стал больше физического finite-amplitude remainder.
Большая доля psi_t и высокий Rayleigh scale объясняют, почему визуально
близкие координатные кривые не гарантировали прохождение полной state-нормы.
Оператор, IC mapping и metric mismatch как причина failure не подтверждаются.

Native-time проверка также не объясняет failure interpolation ошибкой.
Для old full / tight tangent sampled errors `6.929e-8/8.910e-8` на union
actual accepted timestamps становятся `7.042e-8/8.997e-8`: немонотонность
сохраняется. Для tight pair maxima на union совпадают с основной таблицей.
Union имеет19446/19297 точек для old full и23362/23620 для tight full.
Одна сторона в каждом union time оценивается native dense polynomial;
это не actual native output обеих сторон в каждой точке. Новые full
dense-vs-saved differences<=2.51e-17. Old tangent native polynomials
отсутствуют; полноценное old/old native comparison остаётся NOT_AVAILABLE.

Прежний floor1e-8, absolute gate1e-5 и trend/order policy не изменены.
Меньшая новая ошибка лежит ниже прежнего объявленного floor; её разрешимость
здесь подкреплена независимым temporal control. Это не ретроактивный PASS
старого контроля и не новое универсальное правило для сколь угодно малых delta.

## Фактические расходы и ограничения

| New run | rtol | Steps | RHS / Jac / LU | Actual min / max step | Время, s |
|---|---:|---:|---:|---:|---:|
| Tight tangent, one column | `1e-12` | 17617 | 124707 / 2 / 944 | `.0000955870 / .000305811` | 387.518 |
| Tight full delta0 | `1e-12` | 5745 | 40229 / 2 / 6 | `.000378557 / .000736263` | 76.981 |
| Tight full delta0/2 | `1e-12` | 6003 | 42047 / 2 / 12 | `.0000295925 / .000849072` | 79.715 |
| Extra full delta0/2 | `1e-13` | 10571 | 74023 / 2 / 12 | `.000133232 / .000477468` | 137.656 |

Суммарно4 новых ODE,681.869s<1200s: три nonlinear и одна tangent trajectory.
Новых base/static/eigen/FEM solves нет; jobs sequential, one thread.
Max relative energy drift full trajectories `1.04–1.13e-14`; safety calls
не выявили нарушения прежних bounds. Малый drift не используется вместо
temporal comparison.

Первоначальный frozen forecast использовал грубый tolerance exponent1/5
с запасом1.5. Последующая проверка установленного Radau predictor выявила
exponent1/4. Read-only пересчёт прогноза даёт816.102s, также меньше1200s.
Первоначальные648.253s и actual execution evidence не переписаны;
helper исправляет лишь будущий unfrozen forecast, не настройки расчётов.

Нет независимого spatial control p16, новой проверки base trajectory или
отдельной upper bound для floating-point effects. Один дополнительный full
level не является полной temporal certification. Полученная согласованность
относится к двум delta, одной геометрии и sampled finite interval.
Физическая ошибка reduced quartic model относительно3D этой диагностикой
не проверяется; прежние strict float64 и FEM qualifications сохранены.

Для ограниченного возврата к пилоту A/B имеется численное основание:
крутильный blocker объяснён и наблюдаемая согласованность восстановлена.
Mixed/sign, spatial и geometry controls исходного пилота ещё не выполнены;
их продолжение требует отдельного решения и нового cost preflight.
Карты чувствительности и Stage E/F автоматически не запускаются.

## Данные и рисунки

`historical_error_decomposition.npz`, `tight_error_decomposition.npz`
содержат raw coordinate/velocity errors, E0 matrix, delta, times и norms.
`read_only_diagnostics.json` хранит initial, pointwise, arm/field/internal-joint
и top-DOF checks. `primary_temporal_comparison.npz/json` сохраняет все четыре
сочетания; `conditional_comparison.npz/json` — extra control.
Полные новые states, dense polynomials, settings, counters и source snapshots
сохранены отдельно от parent, с attempt ledger и manifest hashes.

CSV: `error_time_histories.csv`, `temporal_pairings.csv`, `actual_execution.csv`.
Два основных PDF+PNG: `torsion_coordinate_velocity_error` и
`torsion_error_amplitude_refinement`. Нет smoothing, time/phase alignment,
изменения амплитуд или исключения высокочастотных компонент.

53 focused tests PASS,4.72s: diagnostic algebra/config/freeze/dispatch/cache
и адресные regressions прежней сборки/протокола; real IVPs внутри tests нет.
Final source/cache/navigation/Git checks сохраняются в `verification.json`.
Первый рисунок не показывает exact t=0 на log scale, где velocity error
равна нулю; все исходные данные, initial checks и comparison maxima включают0.

Финальные проверки: cached compute/report/plot имеют zero scientific/render
calls и неизменный manifest. Старый failed CLI compute также не запускает ODE.
Проверены20 новых navigation links; `git diff --check` и whitespace PASS.
Parent manifest неизменен. Original mathematical source hashes и старые
байты D/K/current/canonical pilot report сохранены; новый раздел добавлен отдельно.

Cwd/root `D:\PHD\CoupledBeams\CoupledBeams`, branch `main`,
HEAD `642758754250608f03b18e16ff1874bcce0d071b` сохранены. Исходные20 dirty
paths предыдущего этапа и пустой index защищены snapshot. Этот этап меняет
19 files: четыре новые (config,note,helper,test) и15 существовавших
(две scoped helper extensions и13 документационных/навигационных файлов).
Полный список и оба Git status — в verification.json; staging/commit/push нет.
Python `D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe`,3.12.4;
NumPy2.1.3,SciPy1.15.2,Matplotlib3.9.2,pytest8.3.4. README/CHANGELOG обновлены.
Evidence содержит98 артефактов; final manifest SHA256
`d889f8ef5e090f9b1de11a8ba9dbdaec60963e42e31f5a4750e6857a7249a2f9`.
Режим остаётся EXPLORATORY_NOT_CERTIFIED, admitted=False; выявление numerical
floor не снимает прежних строгих и физических qualifications.

## Воспроизводимость

~~~powershell
$py = 'D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe'
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
& $py scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --torsion-diagnostic --compute --through-stage read-only
& $py scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --torsion-diagnostic --compute --through-stage primary
& $py scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --torsion-diagnostic --compute --through-stage all
& $py scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --torsion-diagnostic --report-only results/nlsp_coupled_torsion_tangent_diagnostic/a10f1110ffbce7ca
& $py scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --torsion-diagnostic --plot-only results/nlsp_coupled_torsion_tangent_diagnostic/a10f1110ffbce7ca
~~~

Один existing CLI получил scoped diagnostic flag/config/helper, согласно SPC.
Переиспользуются tangent/nonlinear Radau и saved dense helpers; новый solver
не создан. Старый default и historical failed-attempt guard сохранены.

Отложены `DEFERRED_TWO_PARAMETER_ACTION_EXPANSION`, physical ideal-joint
validation, sensitivity maps, stability и 3D FEM angular construction.
Даже успешный diagnostic verdict не возобновляет автоматически Stage E/F.
