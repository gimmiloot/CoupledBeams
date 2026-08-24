# Симметрично слоистая балка Reddy

Раздел содержит первую автономную модель одного прямого симметрично
слоистого стержня. Реализованы изгибная FSDT-подсистема, независимая
продольная подсистема и их блочно-диагональное объединение.

Текущий результат:

```text
RLB-0-SOURCE: PASS_WITH_DOCUMENTED_SOURCE_RECONSTRUCTION
RLB-0A-BENDING: PARTIAL_PASS
RLB-0B-AXIAL: PASS
RLB-0C-COMBINED-UNION: PASS
OVERALL: PARTIAL_PASS
```

`PARTIAL_PASS` для RLB-0A связан с 22 расхождениями с Table 4.3.3 за
пределами точности трёх напечатанных десятичных знаков. Аналитические
изгибные gates, Eq. (4.3.51), продольная модель и объединение спектров
прошли. Параметры, граничные условия и допуски не подбирались по отдельным
случаям.

## Соглашения

У Reddy координата по толщине положительна вниз, а слои нумеруются сверху
вниз. Python API хранит стопку снизу вверх в проектной координате
`z_project = -z_Reddy`. Различие явно записано в source contract. Оно не
меняет текущие палиндромные симметричные benchmark-стопки, но существенно
для будущей модели с \(B\ne0\).

Коэффициент сдвига `K` является обязательным входом. Benchmark Table 4.3.3
задаёт `K=5/6` с provenance
`INFERRED_FROM_INTERNAL_NUMERICAL_CONSISTENCY`; это значение не считается
универсальным свойством произвольного ламината.

## Файлы

- `reddy_ch4_source_contract.md` — формулы, страницы, нормировка и политика
  реконструкции источника;
- `reddy_symmetric_single_beam_validation.md` — численные результаты и
  ограничения;
- `scripts/lib/reddy_symmetric_laminated_beam.py` — узкий вычислительный API;
- `scripts/analysis/laminated_beams/validate_reddy_symmetric_single_beam.py`
  — воспроизводимый CLI;
- `tests/test_reddy_symmetric_laminated_beam.py` — целевые регрессии;
- `tests/data/reddy_ch4_table_4_3_3.json` — машинная транскрипция источника.

Generated CSV, JSON, report и две диагностические фигуры находятся в
`results/laminated_beams/reddy_symmetric_single_beam/`. Каталог игнорируется
Git.

## Запуск

Из корня репозитория:

```powershell
python scripts/analysis/laminated_beams/validate_reddy_symmetric_single_beam.py --source-check-only
python -m pytest -q -p no:cacheprovider tests/test_reddy_symmetric_laminated_beam.py
python scripts/analysis/laminated_beams/validate_reddy_symmetric_single_beam.py --full-validation
python scripts/analysis/laminated_beams/validate_reddy_symmetric_single_beam.py --plot-only
```

`--source-check-only` не импортирует научный solver. `--plot-only` читает
только сохранённые CSV/JSON и не выполняет matrix exponential или поиск
корней. Режим `--smoke` намеренно возвращает для неисполненного RLB-0C статус
`PARTIAL_PASS`, а внутренний статус — `NOT_RUN`.

## Границы этапа

Продольная подсистема является project derivation from symmetric CLT и
стандартной динамики стержня, а не результатом, напечатанным в §4.3.4.
Combined state не содержит искусственной продольно-изгибной связи. На этом
этапе не реализованы coupled rods, angular joint, \(B\ne0\), torsion,
damping, FEM и complex roots.
