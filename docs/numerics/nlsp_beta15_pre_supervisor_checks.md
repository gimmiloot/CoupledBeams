# β15: завершающие проверки перед обсуждением с руководителем

Дата: 2026-10-11–12. Общий статус:
`PRE_SUPERVISOR_CHECKS_COMPLETE_WITH_PL_PARTIAL`.
Исходный положительный OUT-результат сохранён; причины непрохождения PL
neutral gate установлены, но сам строгий gate не закрыт. Обычное начальное
OUT2-отклонение с нулевыми скоростями возбуждает растущее подпространство.
За 5T его физический прогиб остаётся малым.

Bundle: `results/nlsp_beta15_pre_supervisor_checks/d3f6375a27dba883/`.
[Самостоятельный научный отчёт](../reports/nlsp_beta15_parametric_instability/report.tex),
[PDF, 5 страниц](../reports/nlsp_beta15_parametric_instability/report.pdf).
Технический текст не подменяет научный отчёт для руководителя.

## Постановка и происхождение

Фиксированы p14, β=15°, μ=0.1, L₁=0.9, L₂=1.1, b=0.05, h=0.10,
E=ρ=1, ν=0.3, κ=5/6. Тип II, Ψ=0, max|w(s,0)|=0.005,
внешние заделки и прежний идеальный жёсткий узел, включая c₁=c₂.
Никаких новых constitutive assumptions, физического действия, инерционных
членов, RHS/Jacobian, массы, узла или локальных преобразований не менялось.

Использованы сохранённые источники:

| Bundle | SHA256 исторического manifest |
|---|---|
| p12 HB/Hill `6b93bc782316296d` | `0ec30aad088ed9c3d85464bdcacfe024ec88ff9f6b8f901a43d34f7637721d3b` |
| p12 dynamic `49dc29161a5dda78` | `437088a5f20df823d4ac78f8fe30d0d059f9f09aea1e045c850bab8363f91a72` |
| p14 spatial `2f77d2952671122b` | `3a2e35815c950fa21d9fcfbc6ba1efd001633e5c3e29057d097540a55ab6f1f4` |

Source hashes проверены до работы; все файлы этих трёх исторических bundles,
включая PDF и manifests, проверяются относительно начального snapshot.
Исторические результаты не пересчитывались и их статусы не переписаны.
Git main, HEAD `c66743a33cb8539ee3bc7d82b37103a5f0e366c0`, staging пустой.
До работы имелись 39 dirty paths; snapshot и hashes пользовательских файлов
сохранены в `results/_smoke/nlsp_beta15_pre_supervisor_checks/initial_state.json`.

## A. Источник PL time-shift residual

Для X=(q,v) приближённой HB-орбиты вычислены r=Ẋ−f(X), её производная
и независимое действие полного аналитического Jacobian. Проверена identity

\[
\ddot X_{\rm HB}-Df(X_{\rm HB})\dot X_{\rm HB}=\dot r.
\]

В mass-form R=M(q)a+I(q,v)+∇V(q), a=q̈HB:
Ṙ=QHB v+VHB a+M q⃛HB, ṙv=M⁻¹(Ṙ−Ṁrv).
QHB,VHB взяты из прежних аналитических mass-form tangents. Для канонического
Z=S(t)z использованы оба члена ṠS⁻¹+SAS⁻¹, включая переменную массу.
Отдельно сохранена малая off-orbit разность между оператором второй вариации
действия и преобразованным Jacobian около не вполне точной орбиты.

Масштаб относительного остатка: Ω·RMS||Ztime||=2.74136926305e−5
в прежних канонических энергетических координатах.

| Контроль | Результат |
|---|---:|
| Исторический N13 raw residual | 7.58599272745e−8 |
| Pointwise канонический остаток | 7.59294257572e−8 |
| Преобразованная производная r | 7.59294223777e−8 |
| Off-orbit canonical difference | 4.50016e−13 |
| Identity по полному Jacobian, max | 9.49403e−13 |
| Canonical decomposition identity, max | 6.05122e−13 |
| Независимый central FD, hphase=0.001 | 7.00997e−10 |
| Независимый central FD, hphase=0.0005 | 1.65002e−9 |
| MAC сохранённого near-zero PL eigenvector / time shift | 0.999999997643 |
| Сохранённый ненейтральный PL max Re λ | 6.22885e−15 |

FD-проверки проходят прежний 2e−9 identity gate. При уменьшении шага
остаток возрастает из-за вычитания малых величин; это не утверждение о
сходимости второго порядка FD на данном численном полу.

| N | Raw Hill neutral residual | После учёта известной производной r |
|---:|---:|---:|
| 7 | 4.57541e−7 | 4.57540e−7 |
| 11 | 4.26463e−10 | 1.64625e−11 |
| 13 | 7.58599e−8 | 1.64629e−11 |
| 17 | 7.59294e−8 | 1.64635e−11 |
| 25 | 7.59294e−8 | 1.64639e−11 |
| 33 | 7.59294e−8 | 1.64646e−11 |

Доминирует 12-я гармоника производной орбитальной невязки. N=11 её
не видит: малый проецированный остаток при таком усечении не закрывает
физический neutral gate. N≥13 воспроизводит pointwise остаток, а проверка
баланса с ṙ остаётся малой. N7 имеет ещё и заметную ошибку усечения самого
нейтрального направления. Fourier aliasing или ошибка Ṡ не объясняют N13
остаток; он обусловлен конечной точностью исходного HB-представления.

Порог **2e−9 сохранён**. Статус `PL_STABILITY_REMAINS_PARTIAL`.
Коррекция диагностического баланса не является исправлением исходной орбиты,
новой нулевой модой или изменением физических коэффициентов. В ранее
разрешённом ненейтральном PL спектре роста не обнаружено; near-zero splitting
±2.62e−8 согласуется с нейтральным направлением и чувствительной нейтральной
структурой приближённой консервативной орбиты. Полная PL-устойчивость не
утверждается. Новые PL eigenvalue solves не выполнялись.

## B. Displacement-only OUT2: нормировка и доступность

Все поля v,Φ,ψ взяты из сохранённой p14 линейной OUT2. В физическом состоянии
все OUT-скорости нулевые; PL-начальное состояние совпадает с прежней HB13.
Lref=1, max(|v|,Lref|Φ|,Lref|ψ|)=5e−5. Отсутствуют искусственно добавленные
Floquet-скорости или подбор другого начального возмущения.

Пусть ru,rs --- вещественные растущий и обратный Floquet-векторы в физических
velocity coordinates, нормированные ruᵀE₀ru=rsᵀE₀rs=1,
E₀=diag(K₀|OUT,M₀|OUT), ||x||E₀=√(xᵀE₀x). S переводит их в сопряжённые импульсы, J --- прежняя
симплектическая матрица. Использована dual projection

\[
c_u=\frac{(Sr_s)^T J(Sx_0)}{(Sr_s)^T J(Sr_u)},\qquad
c_s=-\frac{(Sr_u)^T J(Sx_0)}{(Sr_s)^T J(Sr_u)}.
\]

При ином масштабе/фазе eigenvector коэффициент меняется, но физический
projector cu·ru инвариантен. Обычное евклидово произведение правых векторов
не использовалось. Получено cu=−1.87794806382405e−7, cs≈cu; разность со
сохранённым результатом −5.29e−23. Диагностический indicator uncertainty
1.12166e−10 подтверждает ненулевую проекцию, но не является строгой границей.

| Начальная составляющая | max|v| | max|Φ| | max|ψ| |
|---|---:|---:|---:|
| Полное OUT2-отклонение | 1.3285261e−5 | 1.4296469e−5 | 5.0000000e−5 |
| Растущая составляющая | 6.5047060e−6 | 1.4277470e−5 | 2.1640881e−5 |
| Обратная составляющая | 6.5047060e−6 | 1.4277470e−5 | 2.1640881e−5 |
| Остальные составляющие | 1.3029616e−6 | 3.1927977e−5 | 9.6309686e−6 |

Φ,ψ --- углы; их сравнение с перемещением применяется только через Lref.
Растущий общий физический масштаб составляет **43.2818%** полного начального,
а max|v| --- **48.9618%** полного начального max|v|. Это не доли энергии:
компоненты не ортогональны, поля взаимно компенсируются. Их E₀-norms:
1.926794e−7 (total), 1.877948e−7 (unstable), 1.877948e−7 (stable),
6.770891e−8 (remaining); последняя составляет 35.1407% total norm.
Сумма этих норм не обязана равняться total norm.

Для растущей составляющей в той же фазе через целое число периодов:

| Интервал | exp(σnT) | max|v| растущей составляющей |
|---:|---:|---:|
| 3T | 1.159258951 | 7.540639e−6 |
| 5T | 1.279285359 | 8.321375e−6 |
| 10T | 1.636571030 | 1.064541e−5 |
| 20T | 2.678364736 | 1.742198e−5 |

Множитель отрицательный: знак меняется каждый T. Таблица --- только
стробоскопический прогноз растущей составляющей, не полная OUT-траектория и
не максимум за период. 10T/20T не интегрировались.

## C. Один короткий tangent control

В p14 архиве сохранены только критический и обратный eigenvectors,
а векторы остальных нейтральных Floquet-компонент отсутствуют. Поэтому
полную временную реакцию выбранного displacement-only начального состояния
нельзя восстановить из критической пары. Выполнен один разрешённый IVP
до 3T и native-checkpoint продолжение той же траектории до 5T.

Использовано всё пространство 81 OUT-моды, а не прямая двухмодовая модель.
Сохранённый канонический Fourier operator m=−16…16 включает массу,
инерцию и Ṡ. На трёх новых независимых фазах его сравнение с полным
аналитическим Jacobian через L(Ṡ+S Jphysical)S⁻¹R даёт relative error
4.18963e−14. Вещественное вычисление F₀+2ReFm cos−2ImFm sin согласуется
с комплексным Fourier representation до 2.82e−16. Физические состояния
восстановлены S(t)⁻¹R; временно меняющиеся координаты не трактуются как
неизменная physical-energy metric.

DOP853, rtol=1e−11, componentwise atol от прежнего norm-scaled rule,
исходное направление и коэффициент масштаба сохранены при продолжении.
Сохранены native accepted grid и endpoint, dense sampled observables T/256.
Нет retries, других initial states или новой PL base integration.

| Величина | 3T | 5T |
|---|---:|---:|
| Усиление растущей dual-проекции, tangent | 1.159258910 | 1.279285429 |
| Спектральный прогноз | 1.159258951 | 1.279285359 |
| Общая E₀-norm amplification в конце | 1.103175867 | 1.217985617 |
| Максимум max|v|/A за интервал | 0.2752784% | 0.2863587% |

Ошибка прогноза projection gain ≤6.95e−8. Finite-window lograte
0.003411316428 согласуется с saved σ=0.003411315508. Это не измерение
показателя Флоке полной нелинейной динамики. Максима по времени sampled,
не строгие continuous-time supremum bounds. Физические spatial maxima
для каждого состояния находятся по прежним полиномиальным профилям.

За 5T max|v|=1.4317936e−5, max|Φ|=4.7841438e−5,
max|ψ|=6.2088013e−5. Начальный max|v|/A=0.2657052%; увеличение
максимального прогиба невелико. Рост компонента 27.9%, total norm 21.8%
и рост физического максимума --- разные показатели. Растущая составляющая
доступна и разрешима уже на 3–5T, но значительного выхода из плоскости нет.
На 10–20T спектральный рост мог бы быть заметнее на фоне остальных
колебаний; для будущего эксперимента это осмысленный масштаб, а не
гарантия большого прогиба или разрешение нового расчёта.

## Проверки, стоимость и результаты

Science production/preflight/diagnostics: **61.39633 s**, из них IVP 3T
25.49805 s и продолжение 3→5T 18.51458 s. Всего один initial-state control,
два вызова tangent ODE с native restart; full nonlinear ODE=0,
новые HB/Hill eigensolves=0, FEM=0. Никаких p16, геометрий, амплитуд,
тип-I, длинной динамики или нормальной формы не выполнялось.

Адресные regressions проверяют исторические ссылки/хэши, OUT2 amplitude и
нулевые velocities, symplectic duality и projector invariance, исходную
time-shift identity на дополнительной независимой фазе, неослабленный
neutral gate, Fourier signs/realness, negative/reciprocal multiplier,
native continuation и физические observables. Existing p14/dynamic tests
также проверяют Hamiltonian/Floquet conventions, mass/joint, eigen-residual,
spatial/temporal qualification и canonical Ṡ. Cached compute/report/plot
выполняются на временных копиях с scientific functions запрещёнными.
Один первоначальный test harness failure (неверное имя config field)
исправлен без изменения gate. В итоге **58 различных tests PASS**:
16 новых и 42 существующих p14/dynamic regressions. Новый набор после
исправления проходит полностью: 16 PASS, 6.85 s pytest / 7.38022 s wall.
С учётом первоначального прогона, smoke checks и конечной проверки сохранности
общий консервативный расход **92.77216 s / 600 s**. Точный budget ledger
сохранён в `verification.json`. `git diff --check` и whitespace checks
новых текстовых файлов PASS. Все 483 исторических файла, шесть физических
kernels и 24 не затронутых задачей начальных dirty files сохранили SHA256.
HEAD и пустой staging не изменились; final Git status содержит 48 file paths
при `--untracked-files=all` против начальных 39. Commit/push не выполнялись.

Existing CLI, без нового runnable script:

```powershell
python -m scripts.analysis.pilot_nlsp_coupled_planar_sensitivity --pre-supervisor-checks --compute
python -m scripts.analysis.pilot_nlsp_coupled_planar_sensitivity --pre-supervisor-checks --report-only <temporary-copy>
python -m scripts.analysis.pilot_nlsp_coupled_planar_sensitivity --pre-supervisor-checks --plot-only <temporary-copy>
```

Завершённый compute/report читает cache без scientific calls; plot строит
три рисунка из готового `report_data.npz`. CLI проверяет source manifests.
Научный LaTeX собирается отдельно штатным XeLaTeX, без package installation;
исходник и рисунки редактируемы. Все пять страниц финального PDF отрендерены
при 110 dpi и просмотрены; нет clipping, collisions, missing glyphs,
overfull boxes или unresolved references. Недостающие optional polyglossia,
caption/booktabs и hyperref dependencies заменены уже установленными babel
и стандартными LaTeX средствами, без изменения научного содержания.
Code source snapshots различают исходную
science implementation и завершающие renderer/report additions;
physical kernel hashes совпадают на всех этапах.

## Рисунки, границы и остановка

В bundle три научных рисунка (PNG/PDF): geometry, physical modes,
growth confirmation. Геометрия следует s от внешних заделок к узлу:
t₁=EX, n₁=−EY, t₂=(−cosβ,−sinβ), n₂=(−sinβ,cosβ), k=−EZ.
β не является малым внутренним углом V между плечами. Физические профили
получены из сохранённых коэффициентов, а не заменены синусоидами.
Growth figure разделяет старый Floquet-selected p12 initial state и
новый displacement-only p14 initial state. Нет карты неустойчивости.

Уровень научного заключения прежний: численное подтверждение в принятой
стержневой модели. OUT-рост p12/Hill/tangent/full и p14/Hill согласуются;
здесь добавлена доступность и короткая observability-проверка естественного
OUT2-отклонения. Нет теоремы полной p-сходимости или независимой 3D/physical
validation. PL spectrum qualification остаётся отдельным PARTIAL.

README/CHANGELOG, navigation, journal и D32/K33 обновлены по полученному
результату. `DEFERRED_TWO_PARAMETER_ACTION_EXPANSION`, physical ideal joint
с c₁=c₂, independent 3D validation, further p-refinement, Type-I dynamics,
damping/experimental observability и large-OUT saturation остаются отложенными.
**STOP: дальнейшие исследования только после обсуждения с руководителем.**
