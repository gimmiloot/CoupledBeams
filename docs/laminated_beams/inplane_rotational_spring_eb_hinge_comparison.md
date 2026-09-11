# Ветви EB с внутренним шарниром и сравнение трёх жёсткостей

## Постановка и критерии до обработки

2026-09-11, исходный HEAD `2ab2ca1ab8ef12308f4ff6b0e3606390e59bc95d`,
`main`, чистое дерево. Дополнение — working-tree version.
Используется [прежняя EB-постановка](inplane_rotational_spring_joint.md#шесть-условий-и-матрица-узла)
и [метод форм и отражения](inplane_rotational_spring_eb_tracked.md#формы-и-зеркальная-симметрия).
Новых физических предположений нет: `L1=L2=l=1`, `b=.20`, `h=.05`, `E=rho=1`,
`Ag=b*h`, `Ig=b*h^3/12`, `A=E*Ag`, `D=E*Ig`, `m=rho*Ag`.
`k_theta=kappa_theta*D/l`, `Omega=omega*l^2*sqrt(m/D)`, `Lambda=sqrt(Omega)`.
Точный шарнир имеет `k_theta=0`, `M1=M2=0` с сохранением передачи сил
и поступательной совместности; это не две независимые консоли.

Локальная политика `frequency-map-v1 / fast_plot`, `tracked_branches`.
Базовые 201 угол: 0…10° через .1°, 10.5…30° через .5°, 31…90° через 1°.
Для kappa=0 сначала импортируются принятые корни с guard из
[sorted-карты](inplane_rotational_spring_eb_beta_maps.md).
Исходные каталоги sorted и tracked kappa=1,100 защищены SHA-256 всех файлов.
Старые частоты, формы, идентификаторы и шесть событий не пересчитываются.

Сохранены 129 узлов и Simpson-квадратура; массовое произведение включает
только u,w обоих плеч в постоянных материальных координатах. Угловое
назначение: полный массовый MAC, глобальная биекция, одинаковый класс
отражения; MAC >= .95 и отрыв >= .20. Эти же критерии сначала применяются
к прямому сопоставлению seed. При неоднозначности разрешён путь по kappa
только при beta=0, не более 12 дополнительных значений. До 50 новых углов
шарнира, двух восстановлений на событие, 6000 полных B и блоков суммарно
на точку; пул до десяти позиций с целой кратной группой.
Корневые gates сохранены: sigma `1e-9`, nullity `1e-12`, физическая невязка
`1e-9`, поступательная `1e-10`, нуль-вектор `1e-9`. Частотный ориентир `1e-6`
по Omega не является строгой апостериорной оценкой. Проверки массовой
ортогональности кратной пары имеют допуск `1e-10`, эквивалентности блоков — `1e-12`.

## Кратность при beta=90°

При `k_theta=0` и `beta=pi/2` блоки классов связаны знаками, а не равны
покомпонентно. Для `R=diag(1,-1,-1)`, `L=diag(1,1,-1)`:
`B_minus=L*B_plus*R`. Это следует из трёх условий классов в предыдущем
отчёте, `F=diag(1,-1,-1,1,-1,-1)`, коммутирования `F*H=H*F`
и согласованного преобразования реакций заделки. Поэтому два независимых
класса имеют одинаковую характеристическую задачу в этой граничной точке.
Ниже проверяются полная B, два проектора `(I+eta*S)/2`, массовая независимость
и одностороннее продолжение от beta<90. Частоты разных разрешённых корней
не смешиваются и не усредняются; совпадающие строки сохраняют разные shape_key.

## Соответствие начальных мод

Reference: первые шесть форм при `beta=0,kappa=1` задают comparison_mode_01…06.
После сопоставления seed таблица соответствия замораживается; далее используется
готовое продолжение каждого исходного branch_id по beta. Текущий sorted номер
не переименовывается. Независимость от любого пути в плоскости beta–kappa
не предполагается. Промежуточные углы каждой кривой принадлежат только её
собственной сетке; объединение сеток и интерполяция частот не используются.

## Фактический результат

Шарнир: **201/201 BASE-точка**, 1408 переиспользованных ROOT/GUARD-позиций
(1404 различных сочетания угол–частота), 1408 восстановленных форм.
Все корни взяты с полной точностью из локальных CSV/diagnostics, без уточнения.
Новых углов и поиска корней для kappa=0 — **ноль**. Получены 1206 записей
шести ветвей: 6 SEED_CONFIRMED и 1200 TRACKED, без пропусков и неоднозначных
назначений. На этом наборе текущие позиции остаются 1…6; это результат
продолжения форм, а не правило их идентификации. Guard при 90° сохранён целиком
в позициях 7/8. Два одинаковых значения частоты имеют разные ключи форм.

Минимальные угловые MAC `.9995424141`, margin `.9991383714`; максимальные
физический остаток `1.399e-10`, остаток полного нуль-вектора `3.879e-16`,
sigma ratio `8.743e-15`, дефект симметрии `2.249e-10`.
Ошибка единичной массовой нормы не превышает `6.662e-16`.
Это проверки согласованности восстановления с исходной характеристической
задачей; независимый механический метод здесь не вводился.

### Граничные вырождения и внутренние события

| Пара ветвей при kappa=0 | Omega при 90° | Lambda при 90° | Min MAC перехода 89° → 90° |
| --- | ---: | ---: | ---: |
| mode_01 / mode_02 | 15.370552823848 | 3.920529661136 | .9999999198 |
| mode_03 / mode_04 | 49.319329725728 | 7.022772225107 | .9999963027 |
| mode_05 / mode_06 | 94.738791770001 | 9.733385421836 | .9999724700 |
| guard 7/8 | 114.031537130972 | 10.678555011375 | не отслеживается как целевая ветвь |

Показанные числа округлены только в тексте. Во всех четырёх группах полная
nullity равна 2, проекторы отражения выделяют eta=+1 и eta=-1. Ошибка
массовой Gram-матрицы относительно I не более `4.476e-16`. Численная ошибка
знакового соотношения блоков равна 0 в использованных матрицах; относительная
ошибка эквивалентности шести физических строк и классов `6.064e-17`, ранг 6.
Обе формы каждой пары проходят физические gates, включая нулевые моменты.
Таким образом, три целевые пары имеют **ENDPOINT_DEGENERACY_SUPPORTED**;
аналогичная кратность guard записана отдельно. Односторонние метки задаются
продолжением форм, не произвольным выбором SVD-базиса.

На 201 узле не обнаружена внутренняя смена порядка ни одной из 15 целевых пар.
Trigger относительного сближения <= .02 срабатывает для пар 1/2 с 43°,
3/4 с 71°, 5/6 с 86°; все заканчиваются указанными граничными вырождениями.
Зазоры Omega при 89° равны соответственно `.0033292700`, `.0457187077`,
`.4327789834`, то есть разрешены на принятом частотном уровне.
Избегаемое пересечение этим набором не установлено. Отсутствие всех внутренних
пересечений между узлами сетки или во всём спектре не доказывается.
Угол за 90° и новые угловые окрестности не рассчитывались.
Старые шесть CROSSING_SUPPORTED при kappa=1,100 сохранены из предыдущего
этапа, их скобки не уточнялись и не считаются вновь проверенными.

### Проверенное соответствие seed

Прямой MAC reference→kappa=0 проходит критерии для всех шести форм.
Для kappa=100 пятая форма имеет прямой MAC `.9359222162`, margin `.9107935325`:
исходный статус SEED_MAPPING_AMBIGUOUS сохранён. Порог не понижался.
При beta=0 выполнен **один** дополнительный случай kappa=10; готовое kappa=.1
лежит вне интервала 1…100 и не заменяет этот шаг. Получены семь корней из
тех же симметрийных блоков и проверены по полной B, восстановлены семь форм.
Ни один старый корень kappa=1 или 100 не вычислялся заново.

| comparison mode | eta | branch при kappa=0 / 1 / 100 | MAC 1 → 0 | MAC 1 → 10 | MAC 10 → 100 |
| --- | ---: | --- | ---: | ---: | ---: |
| 01 | +1 | mode_01 / mode_01 / mode_01 | .99582769 | .99566779 | .99978747 |
| 02 | -1 | mode_02 / mode_02 / mode_02 | 1 | 1 | 1 |
| 03 | +1 | mode_03 / mode_03 / mode_03 | .99109136 | .97975947 | .99808934 |
| 04 | -1 | mode_04 / mode_04 / mode_04 | 1 | 1 | 1 |
| 05 | +1 | mode_05 / mode_05 / mode_05 | .99221488 | .96599469 | .99490727 |
| 06 | -1 | mode_06 / mode_06 / mode_06 | 1 | 1 | 1 |

Минимальные отрывы по шагам 1→0, 1→10, 10→100: `.98662151`, `.95167433`,
`.99356695`. Итоговые 18 соответствий подтверждены, перестановка оказалась
тождественной. Полная точность, конкуренты и исходный прямой отказ сохранены
в таблице и manifest. Для kappa=100 таблица использует оба проверенных шага
всех шести форм. Это продолжение по конкретному пути, не утверждение о любом пути.

Вспомогательный случай kappa=10 имеет max физический остаток `3.331e-11`,
max null residual `5.751e-17`, sigma ratio `9.555e-15`. Сохранён отказ detector
`BOUNDARY_MINIMUM` около Omega=108.7967118541, выше шестой Omega=104.2476964589.
Его интервал отделён от шестой позиции на `4.5490153952`; найденный guard
Omega=108.8279618541 отделён от правой границы поиска на `.35`.
Исходный helper status CONFIRMED не переписан; дополнительная явная запись
ROOT=CONFIRMED, GUARD=QUALIFIED_DETECTOR_WARNING сохраняет это предупреждение.
Для seed используются только шесть проверенных форм. Отказ не объявляется
устранённым; дальнейшее уточнение guard не запускалось.

### Что видно в сравнении

Моды 1,3,5 заметно зависят от жёсткости пружины; близость кривых в отдельных
угловых участках не означает их одинаковости на всём диапазоне.
Кривые мод 2,4,6 совпадают в разрешаемой точности: max относительное различие
Omega на 201 общем BASE-угле соответственно `3.008e-14`, `9.657e-14`,
`1.909e-12`. Это согласуется с независимым от k_theta блоком eta=-1
(см. три условия классов в предыдущем отчёте). Из одной картинки точное
равенство не выводится. Совпадение показывается без искусственного сдвига;
различимые линии и редкие маркеры обозначают три источника данных.

## Данные, рисунки и воспроизведение

Общая таблица содержит **3966 записей**: 1206 новых tracked-записей шарнира
и 2760 неизменённых записей прежних kappa=1,100. Все частоты и текущие
sorted позиции побуквенно сверены с исходными CSV. Каждая кривая использует
свою сетку: 201 угол kappa=0, по 230 углов kappa=1 и 100 с их прежними
локальными точками. Неподтверждённых вершин нет; интерполяция не применялась.

| Рисунок | PNG, 300 dpi | Векторный PDF |
| --- | --- | --- |
| Мода 1 | [PNG](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode01_kappa_comparison.png) | [PDF](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode01_kappa_comparison.pdf) |
| Мода 2 | [PNG](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode02_kappa_comparison.png) | [PDF](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode02_kappa_comparison.pdf) |
| Мода 3 | [PNG](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode03_kappa_comparison.png) | [PDF](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode03_kappa_comparison.pdf) |
| Мода 4 | [PNG](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode04_kappa_comparison.png) | [PDF](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode04_kappa_comparison.pdf) |
| Мода 5 | [PNG](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode05_kappa_comparison.png) | [PDF](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode05_kappa_comparison.pdf) |
| Мода 6 | [PNG](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode06_kappa_comparison.png) | [PDF](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/eb_tracked_mode06_kappa_comparison.pdf) |
| Шесть ветвей шарнира | [PNG](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_k0/eb_spring_tracked_k0.png) | [PDF](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_k0/eb_spring_tracked_k0.pdf) |

Все семь PNG просмотрены; легенды вынесены под оси, подписи не обрезаны.
Это семь содержательных рисунков в двух форматах. Старые обзоры не перерисованы.

Новый набор шарнира:
[корни](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_k0/verified_roots.csv),
[ветви](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_k0/tracked_branches.csv),
[формы](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_k0/shapes.npz),
[diagnostics](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_k0/tracking_diagnostics.json),
[события](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_k0/crossing_events.csv),
[manifest](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_k0/run_manifest.json).
В shape_key входит позиция, а не только частота; 1408 ключей дают 4224 массива
states/reactions/vector. Старые два набора лежат совместно в
[tracked-каталоге kappa=1,100](../../results/laminated_beams/inplane_rotational_spring_eb_tracked/tracked_branches.csv).

Сравнение:
[соответствие seed](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/seed_mode_mapping.csv),
[общая таблица](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/comparison_branches.csv),
[manifest с исходными версиями, семью новыми корнями kappa=10 и диагностикой пути](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/comparison_manifest.json),
[семь форм вспомогательного seed](../../results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison/seed_shapes.npz).
Результаты локальные и игнорируемые Git; прочитаны фактически. Содержательный
итог и код версионируются. Все 29 файлов обоих исходных каталогов, включая
частоты, формы, manifests, события и изображения, побайтно сохранены.
Старый source HEAD `ee40f33…` у kappa=1,100 не заменён новым HEAD дополнения.

Параметризован [прежний entry point](../../scripts/analysis/laminated_beams/track_inplane_rotational_spring_eb.py):
`import_saved`, `track`, `hinge_endpoint_checks`, `prepare_seed_mapping`,
`continue_seed_mapping`, `prepare_comparison_table`, `render_comparison`.
[Helper](../../scripts/lib/inplane_rotational_spring_eb_modes.py) дополнен
`saved_root_classes`; production physics не изменена. Нового runner нет.

```powershell
$beamPython = 'D:\python\Pycharm\pythonProject\.venv\Scripts\python.exe'
$beamTrack = 'scripts/analysis/laminated_beams/track_inplane_rotational_spring_eb.py'
& $beamPython -B $beamTrack --mode compute --kappas 0
& $beamPython -B $beamTrack --mode resume --kappas 0
& $beamPython -B $beamTrack --mode seed-map --kappas 0
& $beamPython -B $beamTrack --mode seed-continue --kappas 0
& $beamPython -B $beamTrack --mode comparison-table --kappas 0
& $beamPython -B $beamTrack --mode plot-only --kappas 0
& $beamPython -B $beamTrack --mode comparison-plot-only --kappas 0
& $beamPython -B -m pytest tests/test_inplane_rotational_spring_eb_modes.py -q -p no:cacheprovider
```

Все перечисленные режимы действительно запускались. `--output` задаёт каталог
углового tracking; `--hinge-output`, `--legacy-tracked`, `--comparison-output`
задают источники и каталог сравнения. Значения по умолчанию — ссылки выше.
`seed-continue` нужен только при сохранённой неоднозначности прямого назначения;
здесь он использует один зафиксированный путь 1→10→100 при beta=0.
Готовый результат переиспользуется. Для изменения подписей, цветов или линий
достаточно соответствующего plot-only: он не запускает seed matching,
tracking, формы или solver и не меняет научные CSV/NPZ.

## Затраты, проверки и ограничения

Python 3.12.4 из `D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe`;
numpy 2.1.3, scipy 1.15.2, matplotlib 3.9.2, pytest 8.3.4. Путь без разделителя
перед `.venv` недоступен; установки и изменения PATH не выполнялись.

Шарнир: 1404 полных B, 8 блоков и одна проверка физической матрицы условий —
1413 построений, максимум 13 на точку; 1408 expm и 2832 аналитических
восстановления плеч. Обработка с атомарными checkpoint/NPZ заняла 77.259 с,
собственно восстановление форм — .749 с внутри этого времени.
Проверочный resume занял 1.760 с и дал **ноль** новых матриц/форм/корней/tracking.
Вспомогательный kappa=10: 878 построений (7 B + 871 блок), 872 expm,
14 аналитических восстановлений плеч, .348 с на локальный путь с проверками.
Первичное прямое сопоставление сохранённых seed заняло .114 с.
Итого 2291 построение с проверкой матрицы условий, 2280 expm,
1415 новых восстановленных форм конструкции и 2846 вычислений форм плеч;
новые частоты — только семь позиций вспомогательного seed, новых углов нет.
Бюджеты не увеличивались, дополнительные восстановительные попытки не нужны.
Вызовы в целевых тестах отделены от этих счётчиков workflow.

Отрисовка обзора — .783 с; финальных шести сравнений — 2.323 с.
Предыдущая отрисовка сравнений 2.292 с сохранена в render history: после
просмотра добавлены только редкие маркеры совпадающих линий.
В plot-only matrix/root/shape/tracking = **0/0/0/0**, seed matching также 0.
Пройдены **34 целевых теста**: прежние 21 и дополнения кратности, независимых
классов, ключей форм, seed-биекции, неизменности исходных данных, своих сеток,
постоянства частот при mapping, root/guard qualification и отсутствия
вычислений при отрисовке/повторном импорте. Полный pytest не запускался.

Сопоставление конкретного пути подтверждено; независимость от произвольного
двумерного пути, полнота всего спектра и отсутствие всех скрытых пересечений
не заявляются. Старые source/Ritz qualifications и непринятые EB-случаи
kappa=10000/legacy RIGID guard остаются в силе. Вязкость, RLB и новые
геометрии не рассчитывались. Следующее исследование автоматически не начато.
