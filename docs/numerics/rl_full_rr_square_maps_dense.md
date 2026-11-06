# Плотные карты квадратных стержней: Λ(β), Λ(μ), Λ(τ)

2026-11-03. **Все пять карт PASS.** Сохранены 735 точек сеток
(731 различная геометрия), 6615 значений рангов 1–9.
Выполнено 544 новых решения, переиспользована 191 точка прежних сеток.
Неразрешённых случаев нет. Научная интерпретация кривых не проводилась.

[Все рисунки: общий вид и панели 1–4 / 5–8, PNG и PDF](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/README.md).
[Общий CSV](../../results/rl_full_rr_maps/7b99bfdcc7975b2d/spectrum.csv),
[расчётный manifest](../../results/rl_full_rr_maps/7b99bfdcc7975b2d/manifest.json),
[manifest рисунков](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/manifest.json).

| Карта | Фиксировано | Сетка, включая концы | Точек | Новых / прежних | min Lᵢ/hᵢ |
| --- | --- | --- | ---: | ---: | ---: |
| BETA | μ=.1, τ=0 | 0…90°, шаг 1° | 91 | 0 / 91 | 51.961524 |
| MU | β=5°, τ=0 | 0… .8, шаг .005 | 161 | 144 / 17 | 11.547005 |
| MU | β=30°, τ=0 | 0… .8, шаг .005 | 161 | 144 / 17 | 11.547005 |
| TAU | β=5°, μ=.1 | −.8… .8, шаг .01 | 161 | 128 / 33 | 35.118846 |
| TAU | β=30°, μ=.1 | −.8… .8, шаг .01 | 161 | 128 / 33 | 35.118846 |

Геометрические отношения Lᵢ/hᵢ и Lᵢ/bᵢ сохранены в каждой строке CSV;
для квадрата они совпадают. Ни одна точка не нарушает предварительный
ориентир Lᵢ/hᵢ≥10. Это не оценка физической точности относительно 3D.

| Карта | Две панели: PNG / векторный PDF | Отдельная таблица / manifest |
| --- | --- | --- |
| BETA | [PNG](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_beta_mu01_tau0_split.png) / [PDF](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_beta_mu01_tau0_split.pdf) | [CSV](../../results/rl_full_rr_maps/a3c8b2c6b1f2f2fa/spectrum.csv) / [manifest](../../results/rl_full_rr_maps/a3c8b2c6b1f2f2fa/manifest.json), прежние данные |
| MU, 5° | [PNG](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_mu_beta5_tau0_dense_split.png) / [PDF](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_mu_beta5_tau0_dense_split.pdf) | [CSV](../../results/rl_full_rr_maps/8ef758bc176d1149/spectrum.csv) / [manifest](../../results/rl_full_rr_maps/8ef758bc176d1149/manifest.json) |
| MU, 30° | [PNG](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_mu_beta30_tau0_dense_split.png) / [PDF](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_mu_beta30_tau0_dense_split.pdf) | [CSV](../../results/rl_full_rr_maps/1897644eb33cd196/spectrum.csv) / [manifest](../../results/rl_full_rr_maps/1897644eb33cd196/manifest.json) |
| TAU, 5° | [PNG](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_tau_beta5_mu01_dense_split.png) / [PDF](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_tau_beta5_mu01_dense_split.pdf) | [CSV](../../results/rl_full_rr_maps/508a363ed17ba31a/spectrum.csv) / [manifest](../../results/rl_full_rr_maps/508a363ed17ba31a/manifest.json) |
| TAU, 30° | [PNG](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_tau_beta30_mu01_dense_split.png) / [PDF](../../results/rl_full_rr_maps/figures/7b99bfdcc7975b2d_1f3615bd8ca5/square_tau_beta30_mu01_dense_split.pdf) | [CSV](../../results/rl_full_rr_maps/f8f6841cfc780106/spectrum.csv) / [manifest](../../results/rl_full_rr_maps/f8f6841cfc780106/manifest.json) |

Общие изображения с суффиксом `_all` доступны в галерее. Всего созданы
10 PNG и 10 PDF. Формат 14×9 дюймов, штатные маркеры 4 pt, PNG 220 dpi;
на каждом узле оставлен маркер. Отрезки соединяют рассчитанные значения,
сглаживания нет. Изображены только sorted-позиции 1–8, guard9 — в таблицах.
Фиксированные параметры указаны в заголовках. PDF не содержат растровых
изображений; цвета, подписи и компоновка проверены по PNG.

**Неизменный контракт и воспроизводимость.** Использованы прежние
[сценарий](../../scripts/analysis/map_rl_full_rr.py),
[квадратный адаптер](../../scripts/lib/rl_full_rr_aspect.py),
[K8-интерфейс](../../scripts/lib/rl_full_rr_sections_prefix.py),
[допуски](../../data/input/rl_full_rr_geometry_c.json) и
[проверенный контроль сечения](../../results/rl_full_rr_square/8c3fc90e3fb0fbfd/manifest.json).
Физическая постановка и нормировка соответствуют
[отчёту исходных квадратных карт](rl_full_rr_square_maps.md):
RL_FULL/Тимошенко, квадрат a=1, RIGID/RIGID, ε=.005,
E=ρ=l=1, ν=.3, k_s=5/6, Λ⁴=ρA₀ω²l⁴/(EI₀), ω=.005Λ².
Применены `frequency-map-v1`, `fast_plot`, `sorted_positions`, K_plot=8,
K_guard=9. Исходники, физические модули и старая конфигурация не изменены.

Создана только новая [конфигурация плотных сеток](../../data/input/rl_full_rr_square_maps_dense.json).
Значения μ сформированы как i/200, i=0…160; τ — как i/100,
i=−80…80, через Decimal с однократным преобразованием в числа JSON.
Проверены точные десятичные шаги, границы, строгий порядок и уникальность.
В [предварительной проверке](../../results/rl_full_rr_maps/_audit_dense_20261103/preflight.json)
подтверждены совпадение контрактов, хешей исходников и 191 доступная старая
точка без вызова решателя. Новый fingerprint зависит от новых сеток;
переиспользование отдельных точек определяется полным прежним контрактом.

**Численные проверки.** Все 731 уникальные геометрии имеют нулевой счёт
ниже первой частоты и счёт 9 на контрольной границе guard9; все сохранённые
значения простые. Не возникло неразрешённых интервалов или кратных кластеров.
Действующие процедуры их обработки и допуски сохранены.
Максимумы по всем девяти рангам:

| Проверка | Максимальная ошибка | Прежний допуск |
| --- | ---: | ---: |
| Граничные условия и узел | 2.674·10⁻¹⁰ | 10⁻⁹ |
| Уравнения движения | 1.394·10⁻¹² | 10⁻⁹ |
| Энергетический баланс | 4.556·10⁻¹³ | 5·10⁻⁸ |
| Массовый Gram | 7.386·10⁻¹² | 5·10⁻⁷ |

Max |M−1|=2.23·10⁻¹⁶; ошибка сохранения объёма ≤7.78·10⁻¹⁶.
Новая часть потребовала 6053 count-запроса (11–17 на точку), без повторов
запросов и расширений правой границы. В переиспользованных старых записях
сохраняются десять прежних `INSUFFICIENT_RIGHT_COUNT`; это история тех
расчётов, не новые отказы и не неразрешённые точки.

Контроль соседних значений применён на всех пяти заданных сетках. Триггер
изменения ω более .2 сработал только при β=2°,3°; обе проверки переиспользованы
из прежнего архива, статус `REPRODUCED`. **Новых адресных повторов — 0.**
Дополнительных узлов, Ritz-расчётов, certified_audit и tracking не было.
Счёт, BC/ODE, энергетические и массовые проверки пройдены в каждой новой
точке. [Сводка диагностики](../../results/rl_full_rr_maps/_audit_dense_20261103/numerical_summary.json),
[журнал расчёта](../../results/rl_full_rr_maps/_audit_dense_20261103/compute.log).
Время расчётного прохода — 147.08 с; построения 20 файлов рисунков — 6.87 с.

**Работа без новых решений.** Отдельные таблицы карт собраны из общего
комплекта; их manifests закономерно показывают ноль новых решений.
Статистика 544/191 относится к первому общему расчётному проходу.
При сборке отдельных таблиц, `tables-only` и построении рисунков вызовы
solve/count/catalogues/Ritz/физического базиса были программно запрещены:
[проверка завершения](../../results/rl_full_rr_maps/_audit_dense_20261103/finalization.json).

```powershell
python scripts/analysis/map_rl_full_rr.py --config data/input/rl_full_rr_square_maps_dense.json --mode ALL
python scripts/analysis/map_rl_full_rr.py --config data/input/rl_full_rr_square_maps_dense.json --plot-only results/rl_full_rr_maps/7b99bfdcc7975b2d
python scripts/analysis/map_rl_full_rr.py --tables-only results/rl_full_rr_maps/7b99bfdcc7975b2d
```

Первая команда переиспользует завершённый комплект. Для выбора карты доступны
`--mode MU` / `--mode TAU` и `--name`, например
`square_mu_beta5_tau0_dense`. Диапазоны изменяются в `maps[].grid`,
фиксированные параметры — в `maps[].fixed`, оформление — в `plot`
или отдельном `--style`. Изменения оформления требуют только `plot-only`.

**Сохранность.** До запуска зафиксированы Git-статус и SHA256 существующих
исходников, конфигураций, закрытых документов и затрагиваемых старых архивов.
Проверены хеши новых manifest/CSV, точечных записей и NPZ форм.
Все 1719 переиспользованных Λ совпали с исходным комплектом точно.
Итоговая [проверка сохранности](../../results/rl_full_rr_maps/_audit_dense_20261103/preservation.json)
подтверждает неизменность старых архивов и защищённых материалов;
из прежних документов дополнены только README, CHANGELOG и индекс результатов.
Исторические отчёты, D/K, outline/claims/materials не редактировались.
Commit/push, FEM и новые исследования не выполнялись.
