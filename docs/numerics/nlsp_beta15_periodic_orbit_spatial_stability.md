# β15: уточнённая плоская периодическая ветвь и спектральный OUT-кандидат

**OUT_GROWTH_CANDIDATE_RESOLVED**, 2026-10-11, для одной точки: **A=0.005**,
β=15°, μ=0.1, p12, основная ветвь типа II с Ψ=0.
Уточнённая орбита даёт σ⊥=0.003411410067 в спектральной постановке Хилла.
Показатель сохраняется при включении всех 69 OUT-мод p12 и изменении N=7/9.
Другие три амплитуды не прошли строгий HB gate при H≤11; их устойчивость не оценивалась.
Это кандидат дискретной модели около контролируемой приближённой орбиты,
а не доказанная неустойчивость точной физической конструкции.

## Источники и постановка

Использованы [принятое действие](../theory/weakly_nonlinear_spatial_rod.md),
[сохранённые полиномы](../theory/weakly_nonlinear_spatial_rod_expansion_generated.md),
[канонический аудит](nlsp_out_of_plane_instability_feasibility.md),
[остановленный одномодовый PL](nlsp_beta15_nonlinear_planar_resonance.md) и
[двухмодовая ветвь](nlsp_beta15_planar_1to3_internal_resonance.md).
Исторический bundle `results/nlsp_beta15_planar_1to3_internal_resonance/b59944c2ef43d0b7/`
и его manifest SHA256 `dd1f5d5300783fc67f0386081715a0f2a64865d1dabe8e2f7ce0362f06157df6`
не изменены. Новый bundle:
`results/nlsp_beta15_periodic_orbit_spatial_stability/6b93bc782316296d/`.

L₁=.9,L₂=1.1; оба b=.05,h=.10; E=ρ=1,ν=.3,κ=5/6.
Все physical kernels, S≤4, cubic mass-form equations, nonlinear inertia,
variable mass, Galerkin space, local/global maps, внешние заделки и
идеальный узел и c₁=c₂ сохранены. Совпадение действия и Jacobian проверено.
Новых ODE, tangent propagators, FEM, геометрий, p16 или иных ветвей нет.
Частоты сохранённых p12 форм: ω₂=.4348207304902722,ω₄=1.3063293949106531.
Текущие sorted positions используются только в этой фиксированной геометрии.

## Уточнение полной плоской системы

Все 92 плоские координаты (u,w,θ,c) представлены косинусными гармониками0…H.
Их coefficients и Ω решаются совместно Newton/AFT, без пассивного Picard-исключения.
H=9 имеет921 unknowns; H=11 —1105. Исходные Fourier coefficients каждой точки
взяты из прежней Ψ=0 ветви; x=U₁,₂ и y=U₃,₄ остаются положительными.
Косинусная постановка следует из обратимости принятого действия и выбранного
состояния поворота с нулевыми скоростями. Фаза фиксирована этой постановкой.
A — точный максимум пространственного полинома w(s,0), с его стационарными точками.
Никакого amplitude/phase fitting нет; форма4 не задаётся дополнительно.

Решается исходная mass-form невязка R=M(q)q̈+inertia(q,q̇)+∇V(q).
Её Jacobian использует прежние analytic force tangents с фактическим HB-ускорением.
Независимая проверка по направлению дала относительную разность9.46e−10.
На независимых53 фазах со сдвигом.173 проверяется полный R, не только его
проекция на retained harmonics. Нормировка R — прежний инерционный масштаб
Ω²√(x²+81y²) в resting-mass coordinates.

До спектральных результатов зафиксированы: full relative residual≤1e−7,
energy variation≤1e−8, |Ω₁₁−Ω₉|≤2e−7, amplitude error≤2e−12,
relative active-coordinate change≤.02. Operator/Hamiltonian gate2e−9,
Hill eigen-residual1e−9, edge weight1e−6, spectral-control gate5e−6,
numerical floor1e−8; growth must exceed5×observed uncertainty.
Эти значения не менялись после просмотра знака показателя.

| A | Ω, полный HB11 | Independent full residual | Re λ OUT | Re λ PL, ненейтральные | Статус точки |
|---:|---:|---:|---:|---:|---|
| 0.005 | 0.435115586723 | 8.855e-09 | 0.003411410 | 4.05e-15 | OUT_GROWTH_CANDIDATE_RESOLVED |
| 0.006 | 0.435238797100 | 1.914e-07 | — | — | NOT_EVALUATED_HB_GATE |
| 0.008 | 0.435471178146 | 2.407e-05 | — | — | NOT_EVALUATED_HB_GATE |
| 0.010 | 0.435686011681 | 2.260e-04 | — | — | NOT_EVALUATED_HB_GATE |

Ω изменились относительно исходного двухмодового приближения на
+9.47e−12,+6.06e−11,+4.84e−10,−3.22e−9.
Разности9/11 гармоник:5.06e−13,9.37e−12,1.85e−10,4.01e−9.
Малая разность частот сама по себе не заменяет full residual gate.
При больших A основная невязка находится в отброшенной12-й гармонике;
её коэффициентная норма растёт до3.99e−8 приA=.010.
Поэтому точные algebraic roots внутри H11 не делают эти орбиты пригодными
для текущей проверки устойчивости. H не повышалось сверх объявленного11.

Для принятой точки energy variation=4.260e-11,
minimum mass eigenvalue=0.985729796, joint compatibility≤1.12e−16.
Скорости приt=0 нулевые; cosine closure за2πτ точное в представленном усечении.
Scaled Newton Jacobian condition=5730.8;
для всех четырёх точек4288…7019. Последние две retained гармоники принятой
орбиты имеют relative energy-coordinate coefficient norm1.10e−8.
Доля формы4 в начальной квадратичной энергии активных координат остаётся.8995%.
Это не отдельно сохраняемая энергия формы4.

## Касательная система и канонические переменные

Плоский и OUT блоки полного Jacobian разделены: запрещённые блоки точно нулевые.
Вторая вариация сохранённого действия совпадает coefficient-exactly с
линейной по OUT частью mass-form equations; нечётных OUT-мономов нет.
Для каждого блока собраны M,C,Q,Ṁ,Ċ из прежних T≤4,V≤4.
Π=Mη̇+Cη, S=[[I,0],[C,M]], и канонический оператор учитывает
ṠS⁻¹+SAS⁻¹. Член Ṡ не отброшен.

Action/полный velocity Jacobian agreement≤1.48e−14;
канонический Hamiltonian defect≤8.31e−17;
transform identity≤2.58e−16.
При использовании Ṡ именно HB-орбиты и полного physical Jacobian дополнительная
разность≤1.49e−14. Это также контролирует off-orbit discrepancy.
Left/right duality≤1.24e−13. Для спектральной проекции применяются
canonical energy coordinates (ωQ,P), с mass-normalized physical eigenvectors.
Положительные/отрицательные частотные ветви сохранены.

96 фаз дают Fourier coefficients m=−16…16. Реальность F₋m=conj(Fm)
сохранена без исправляющей симметризации. На12 других фазах reconstruction
error=8.721e-15.
Нормы OUT operators в этих одинаково нормированных координатах:
||F₁||F=3.75291,||F₂||F=.103614,||F₃||F=.227210 в1/model time.
Третья гармоника присутствует; она не заменяется single-harmonic ansatz.
Эти полные нормы не являются отдельно выделенными резонансными g.

## Метод Хилла и численная разрешённость

Подстановка z=exp(λt)Σzₙexp(inΩt) даёт
Hₙₘ=Fₙ₋ₘ−inΩδₙₘI. Знак подтверждён constant-oscillator limit.
λ иλ+ikΩ — replicas. Показатели приводятся в strip[−Ω/2,Ω/2),
с одновременным сдвигом Fourier-вектора; accepted counts не выдаются за
число различных физических мод. Требуется малый вес крайних coefficients
и исходного, и сдвинутого вектора. Неконтролируемые boundary replicas исключены.

| OUT-моды p12 | N | σ⊥,1/model time |
|---:|---:|---:|
| 8 | 7 | .0033618895 |
| 16 | 7 | .0033650410 |
| 24 | 7 | .0033663575 |
| 48 | 7 | .0034099308 |
| 69 | 7 | .003411410066593 |
| 69 | 9 | .003411410066595 |

После первичного положительного результата объявлен ограниченный full-space
control:48/69 мод иN7/9, без изменения thresholds. Dense OUT matrices занимают
65.4/104.9MiB; полная PL N7 —116.2MiB. Стоимость оценена до их сборки.
Переход24→48 заметнее, чем16→24; поэтому первый restricted результат не принимался
за полную p12 характеристику. Последний результат включает все69OUTмод.

Последняя modal difference=1.479e-06;
temporal difference=1.854e-15;
HB9/11 orbit sensitivity в24 модах=4.456e-12.
Это observed indicators, не строгие апостериорные error bounds.
Их максимум с floor1e−8 даёт uncertainty indicator=1.479e-06.
Положительный показатель превышает его более чем в2300 раз.
Все69 пространственных мод входят в оператор; σ здесь — максимум по
разрешённым Hill eigenpairs. Дальние высокочастотные показатели, для которых
перенос в fundamental strip требует Fourier support внеN, не квалифицированы.
В последней Hill eigenpair residual относительно matrix norm1.49e−17,
относительно собственного λ и action6.56e−15; independent convolution совпадает
с dense matrix action. Quartet defect1.07e−14;10accepted replicas критической
семьи имеют real-part spread4.57e−12 и относятся к одному механизму.
Zero-base σ≤4.44e−15; его canonical operator совпал с линейным до3.63e−14.

## Плоский контроль и форма опасного возмущения

После8/16PLмод выполнено сравнение полного пространства92 мод,N5/7.
Ненейтральный σ∥=4.052e-15; положительного ненейтрального
роста в разрешённом спектре не обнаружено. Близкие нулевые показатели,
включая real splitting±2.73e−9, находятся ниже numerical floor и не
трактуются как неустойчивость. Полный canonical time-shift direction
проверен convolution приN11: relative residual=9.221e-12.
Предварительный16-modal нейтральный residual.00230 был projection error;
он не принимался за характеристику полной системы. Neutral band5e−6 отделён
от остальных показателей. Неполнота finite temporal spectrum не позволяет
утверждать глобальную плоскую устойчивость.

OUT critical vector состоит преимущественно из sortedOUT2:
93.39% energy-coordinate Fourier weight; OUT4 —5.77%, OUT5 —.785%.
Это modal composition, не branch tracking. В phase-zero physical complex
state reference-energy bending/torsion fractions равны97.30%/2.70%.
Канонические импульсы обратно переведены в физические скорости при помощи
текущих M,C. Профили v,Φ,ψ сохранены; фаза и величина eigenvector произвольны,
он не является finite-amplitude spatial trajectory.

Условный tₑ=293.133918model time, Nₑ=20.299757 периодов.
Это уже разумный масштаб для одного будущего адресного контроляA=.005,
но такой запуск требует отдельного задания и здесь не выполняется.

## Интерпретация и ограничения

Прежний harmonic M0 не давал роста приA=.005; уточнённая Type-II орбита даёт
разрешённый OUT2-dominated spectral candidate. Исторический интервал
.00555… .00988 не подтверждён и не подставлен вместо новых результатов.
Менялись одновременно нелинейная частота, средняя конфигурация, форма4 и
гармонический состав. Причинный вклад внутреннего1:3 отдельно не установлен.
Адресное сравнение Fourier contributions могло бы разделить влияние;
орбита с искусственно удалённой формой4 не становится автоматически свободной орбитой.
Три остальные A остаются NOT_EVALUATED, а не stable/unstable.

GateA PASS; GateB PARTIAL (1 из4); GateC/D/E PASS только дляA=.005.
Общий статус OUT_GROWTH_CANDIDATE_RESOLVED относится к наличию одного
кандидата; вся четырёхамплитудная проверка не объявляется завершённой без ограничений.
Остались existence/accuracy exact periodic orbit, spatialp-refinement,
динамический контроль и независимая3Dvalidation. Type-I initial excitation
не заменяется Type-II periodic branch. DEFERRED_TWO_PARAMETER_ACTION_EXPANSION
и physical idealjoint/c₁=c₂ остаются отложенными.

## Воспроизводимость и технические проверки

Расширен существующий scoped helper; новый runnable solver не создан:

```powershell
python scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --beta15-periodic-spatial-stability --compute
python scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --beta15-periodic-spatial-stability --report-only results/nlsp_beta15_periodic_orbit_spatial_stability/6b93bc782316296d
python scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py --beta15-periodic-spatial-stability --plot-only results/nlsp_beta15_periodic_orbit_spatial_stability/6b93bc782316296d
```

Completed replay использует только cache. Tests compute/report/plot работают
на temporary copies; historical PDF/manifest/source hashes сохранены.
Первый qualification writer ошибочно применял refined-orbit absolute gate
также к HB9 comparator. Применение исправлено; все thresholds прежние,
HB solutions не повторялись. Исправлен выбор replicas у края fundamentalstrip
и проверка одновременно original/shifted edges. Неразрешённые coarse spectra
не использовались для классификации. Evidence исходных проверок сохранена
в qualification_writer_correction; никакой physical kernel не исправлялся.

Научные вычисления и conservative failed-attempt charge около121s из900s;
учёт с regression runtime приведён в numerical_costs/verification.
48Hill eigensolves, включая28 техническихповторов;20retained problems.
5 недорогих повторных OUT diagonalizations одного неизменногоK/M (cache reuse
исправлен),2600production HB phase evaluations;686+4canonical phase evaluations,
14pointwise Jacobian evaluations. Префлайт и regression probes учитываются отдельно.
Новые nonlinear/tangent ODE/FEM calls —0.
58 адресных tests PASS (19 новых и39 прежних); дополнительный cache recheck3PASS.
Консервативный общий учёт с tests и короткими probes — 146.629s из900s.
Postprocessing двух рисунков — 1.146s.
Тесты, git diff --check и final source verification сохранены в verification.json.

Рисунки: [HB orbit controls](../../results/nlsp_beta15_periodic_orbit_spatial_stability/6b93bc782316296d/HB_orbit_controls.pdf)
и [Hill growth controls](../../results/nlsp_beta15_periodic_orbit_spatial_stability/6b93bc782316296d/Hill_growth_controls.pdf).
Имеются PNG-копии. Указаны единицы нормированной модели; кривые не сглажены.
После этапа STOP: новые spatial dynamics, геометрии и амплитуды не запускаются.
