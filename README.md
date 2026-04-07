# CoupledBeams

CoupledBeams — исследовательский репозиторий по задаче о собственных колебаниях двух упругих стержней, сопряжённых под углом. В репозитории собраны теория, литература, аналитические программы и вычислительные результаты.

Базовые обозначения для theory-facing формул согласованы с `docs/literature/pdf/Статья-Дорофеев-2025.pdf`: `\Lambda` — безразмерная частота, `\beta` — угол сопряжения, `\varepsilon` — параметр толщины. Локальный параметр `\mu` используется как параметр несимметрии длины в постановках с разными длинами стержней.

При сравнении с `docs/literature/pdf/2003JSVb.pdf` нужно учитывать известную проблему со знаками в determinant-like записи матрицы `T` формулы `(39)`: печатный sign pattern нельзя переносить в теорию или код без ручной сверки с условиями непрерывности. В локальной теории часть знаков уже поглощена аргументами вида `-\Lambda(...)`.

## Структура проекта

- Теория: `docs/theory/`
- Литература и заметки по источникам: `docs/literature/`
- Журнал проекта: `docs/project_log/journal.md`
- Аналитический код: `src/my_project/analytic/`
- FEM-код: `src/my_project/fem/`
- Запускные скрипты: `scripts/`
- Вычислительные результаты: `results/`

## Аналитический код

После структурного рефакторинга общий математический слой аналитических программ вынесен в:

- `src/my_project/analytic/formulas.py` — сборка общих параметров, матрицы и determinant без изменения формул.
- `src/my_project/analytic/solvers.py` — общая solver-логика: поиск корней, трекинг ветвей, уточнение близких пар и служебные численные функции.

Сценарии верхнего уровня по-прежнему запускаются отдельно:

- `src/my_project/analytic/FreqFromAngle.py` — сценарий по углу `\beta`.
- `src/my_project/analytic/FreqFromMu.py` — сценарий по параметру несимметрии `\mu`.

## Запуск

Из корня репозитория:

```bash
python src/my_project/analytic/FreqFromAngle.py
python src/my_project/analytic/FreqFromMu.py
```

Baseline FEM-сценарий запускается отдельно:

```bash
python src/my_project/fem/python_fem.py
```

Для FEM-скрипта нужны `numpy` и `scipy`. Скрипт не требует входных файлов и сохраняет таблицу спектра в `results/fem_spectrum.csv`.

Для быстрой проверки аналитического слоя есть smoke test:

```bash
python -m unittest discover -s tests -p "test_analytic_smoke.py"
```

Smoke test находится в `tests/test_analytic_smoke.py` и проверяет, что после рефакторинга determinant, контрольные корни и базовые численные значения не изменились.
