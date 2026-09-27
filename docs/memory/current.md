# Текущий контекст

## 1. Текущее научное состояние

2026-09-27: [D20](decisions.md#rlb-d20) завершено,
[K21](knowledge.md#rlb-k21) — **RLB_KV_PRODUCTION_PASS** в ограниченном
регрессионном наборе. Общий production API поддерживает EB и проектную
RLB: точно одинаковые плечи → reduced eta±; неодинаковые → full 6×6.
Eta=−1 возвращает упругое решение без complex Newton.

Воспроизведены сохранённые RLB-контроли K12, проверены full/reduced,
точный invS=J=0 предел и один упругий несимметричный контроль K11.
EB-регрессии проходят. Численных блокеров в новом наборе нет;
старые source/Ritz qualifications и raw full-флаг C не пересмотрены.
Это техническая валидация, не новый результат о влиянии RLB на демпфирование.

## 2. Основания следующей работы

- [K21 и отчёт](../laminated_beams/inplane_kelvin_voigt_rlb_solver_architecture.md) — общий API, контрольный набор, точный предел и команды.
- [K12](knowledge.md#rlb-k12) / [K14](knowledge.md#rlb-k14) — KV-теория, внутренний пилот и внешняя проверка с source-print qualifications.
- [K15](knowledge.md#rlb-k15) / [K19](knowledge.md#rlb-k19) — завершённые EB screening и сравнение трёх ACTIVE-состояний.
- [K18](knowledge.md#rlb-k18) — происхождение production routing и прежние EB-квалификации.
- [K20](knowledge.md#rlb-k20) — локальная чётность по вязкости; не требует одинаковых плеч.

## 3. Открытое решение

Возможный следующий этап — sparse elastic RLB screening и сравнение
EB/RLB по Omega, Delta_psi, s_joint и zeta_slope_pred.
Он требует отдельного решения пользователя. Углы, состояния и дальнейшая
вязкость этим техническим этапом не выбираются.

## 4. Действующие границы и остановка

Production RLB validation завершена; новый screening не начат.
Не разрешены новые физические d/beta sweeps, asymmetric positive-d roots,
crossing/veering, strong damping, FRF, FEM или уточнение точности.
Обычный повтор завершённой регрессии — missing-only без solver/matrix/form
вызовов; он не восстанавливает отсутствующие старые результаты автоматически.

## 5. Среда

Проверено на этом этапе: D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe,
Python 3.12.4, NumPy 2.1.3, SciPy 1.15.2; matplotlib и pytest доступны.
52 уникальных целевых теста пройдено; полный pytest не запускался.
Подробные затраты, исходные HEAD и хеши — в отчёте и manifest K21.
