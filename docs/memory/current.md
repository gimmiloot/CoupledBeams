# Текущий контекст

## Circular EB reviewer diagnostic — завершён, 2026-10-05

[D25](decisions.md#rlb-d25) / [K26](knowledge.md#rlb-k26):
[единая научная заметка](../laminated_beams/circular_eb_rotational_spring_rigid_limit.md).
Одна геометрия круглых изотропных EB-стержней: mu=.30, beta=15°,
r=.005 м, l=1 м. Это отдельное продолжение темы поворотного узла,
без переноса результата на RLB, анизотропию, Bishop или вязкость.

Exact RIGID согласован с baseline для первых шести корней и guard 7
(max|Delta Lambda|=3.016e-10). Generic solver разрешил событие около
204.287 Гц как один простой корень позиции 5. Sorted spectrum растёт
1→10→100→RIGID, отличие от RIGID уменьшается для всех шести позиций.

28 основных форм восстановлены с M=1. Ограниченное продолжение по
kappa подтвердило все шесть локальных seeds до RIGID:
1→3.162278→5.623413→10→17.782794→31.622777→100→RIGID.
Min accepted MAC=.968142, margin=.936732; пороги неизменны,
endpoint conflicts и смены sorted position отсутствуют.

## Сохранённые qualifications

- Полная 12-root эквивалентность не установлена: отдельная численная
  SVD-nullity classification на позициях 11–12 около Lambda=18.139529.
- Direct 1→10 для 02/03/05/06 остаётся UNRESOLVED. Успешное
  continuation — отдельный результат, а не переименование direct.
- |Delta psi| и s уменьшаются вдоль проверенного пути у 5/6.
  Seed 06 от 1 до 10 сначала увеличивает их на 2.95634% и 4.53756%.
  В exact RIGID относительный поворот проходит physical gate; s не задана.
- Старые [EB-JOINT-K01](knowledge.md#eb-joint-k01) /
  [K02](knowledge.md#eb-joint-k02) сохраняют остановки прежнего workflow,
  без доказательства физической кратности.
- Локальные raw results игнорируются Git; смысл и ограничения сохраняет
  научная заметка. Ни реальный k_theta, ни C, ни finite-r/l ошибка не найдены.

## Где остановились и следующий вопрос

1D spring/RIGID reviewer diagnostic завершён. Следующий выбранный вопрос —
теоретическое обоснование transmission/equilibrium conditions и границ
thin-joint asymptotics. В заметке отдельно записаны литературные
асимптотические условия, баланс на области узла, условный compact-joint
scaling и конечный 1D численный результат.

Локальная FEM-калибровка реального узла отложена пользователем.
Новые корни, kappa, геометрии, r/l sweep и исправление high-spectrum
classification автоматически не разрешены. Текущий этап — документационная
синхронизация без расчётов или изменения моделей.

## Независимый статус EB/RLB Kelvin–Voigt

Пауза `PAUSED_FOR_SUPERVISOR_DIRECTION` от 2026-09-27 относится к этой
отдельной линии и сохраняется: [D22/K23](../laminated_beams/inplane_kelvin_voigt_eb_rlb_complex_confirmation.md),
[сводка для руководителя](../laminated_beams/inplane_kelvin_voigt_research_status_for_supervisor.md).
Заключительное сравнение пяти EB/RLB-пар подтверждает прогноз отношения
затуханий в своих исходных границах; новых расчётов по нему не требуется.
Подготовлены [Figure 3](../laminated_beams/figure03_crossing_veering.md),
[Figure 4 revised с сохранённой v1](../laminated_beams/figure04_eb_rlb_shapes.md),
[Figure 5 — огибающие](../laminated_beams/figure05_damping_envelopes.md);
[прежняя столбчатая версия](../laminated_beams/figure05_frequency_damping_comparison.md)
остаётся историей. Новая физическая программа KV не выбрана.
