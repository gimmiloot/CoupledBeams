# Mindlin–Herrmann + Timoshenko: локальные источники и границы переноса

Дата: 2026-10-05. Литературно-документационный этап. Здесь зарегистрированы
постановки и проверки **авторов**, без воспроизведения расчётов, нового
вывода модели или выбора условий углового узла CoupledBeams.
Метаданные, версии, пути, SHA256 и объём чтения — в
[source index](source_index.md); ключи — в [bibliography.bib](bibliography.bib).

## Статус направления

- [Standalone Bishop reproduction](../theory/bishop_literature_reproduction.md)
  сохранён вместе с кодом, тестами, benchmarks и раздельными статусами
  численной проверки и совпадения с печатью. Bishop остаётся reference /
  diagnostic longitudinal model.
- [Аудит Timoshenko + Bishop](../theory/timoshenko_bishop_single_rod.md)
  сохраняет статус `COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED`: исчезновение
  cross terms у центрированного прямоугольника не обеспечило одновременного
  восстановления обеих неизменённых исходных теорий из буквальной общей
  3D-кинематики. Произвольное hybrid closure не принято.
- Bishop больше не является предпочтительным кандидатом основной
  продольной части будущей combined in-plane model. Текущий **кандидат** —
  **Mindlin–Herrmann axial + Timoshenko bending**, с независимой переменной
  поперечного сокращения. Реализация и валидация этой комбинации в проекте
  **не выполнены**; наличие литературных formulations не меняет этот статус.
- Новые источники показывают конкретные способы редуцированного сочетания
  теорий. Они не отменяют предыдущий Bishop audit и не подтверждают модель
  для нашей геометрии. Условия нашего углового узла не выбраны.

## Сравнение источников

Ссылки ведут к canonical records. «Нет локальной связи» относится только
к соответствующему прямому элементу; связь внутри M-H между axial motion
и lateral contraction при этом сохраняется.

| Source | Axial theory | Bending theory | Independent lateral-contraction DOF? | Local axial-bending coupling? | How combined? | Cross-section/material | Problem type | Validation у авторов | Correction-factor status | Relevance to our project | Transfer limitations |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| [Rucka 2010](source_index.md#rucka_2010_l_joint_guided_waves) | M-H | Timoshenko | Да, source ψ | Нет в локальных блоках; есть mode conversion в L-joint | Четыре поля; блочные D, E, μ; time-domain SEM, поворот и сборка | Сталь, квадрат 6×6 mm | Высокочастотные волновые пакеты, обнаружение надреза | Straight rod, затем intact/notched L-joint; эксперимент и SEM | K₁ᴹᴴ=1.1, K₂ᴹᴴ=2.1, K₁ᵀⁱᵐ=.95 fitted; K₂ᵀⁱᵐ=12K₁ᵀⁱᵐ/π² по cutoff | Прямой опубликованный M-H-Tim frame precedent | Калибровка 100/120 kHz; не low-frequency defaults и не наши условия узла |
| [Jang–Park–Lee 2014](source_index.md#jang_2014_timoshenko_composite_patch_guided_waves) | M-H по толщине | Timoshenko | Да, source ψ_b | Bare isotropic base: раздельные axial/bending contributions; patched system: coupled | Displacement fields + **разные reduced constitutive blocks**; frequency-domain SEM | Прямоугольная металлическая основа, односторонняя composite patch | Частоты, FRF, guided waves, debonding | Собственный 1D FEM и ANSYS 2D plane-stress; численная проверка | κ_b, κ_c — source shear corrections; численный default здесь не установлен | Основной источник для понимания reduced closure и четырёх полей | Patch coupling не доказывает обязательной связи bare beam; не единый неизменённый полный 3D law |
| [Banerjee–Ananthapuvirajah 2019](source_index.md#banerjee_2019_rayleigh_love_timoshenko) | Rayleigh–Love | Timoshenko | Нет самостоятельного M-H DOF | Приняты uncoupled, выводятся независимо | Отдельные DSM, simple superposition, 6×6 plane-frame DSM | Изотропные стержни; рама задана жёсткостями и инерциями | Свободные колебания rods/frame | Аналитический rod control, опубликованный stepped-bar пример, сравнение теорий для рамы | В frame example k=2/3; source input, не M-H fit | Published precedent for modular assembly of refined 1D theories | Не M-H equations и не Bishop + Timoshenko validation |
| [Liu et al. 2021](source_index.md#liu_2021_multibody_beams_rigid_bodies) | В demonstration: classical / Rayleigh–Love | В demonstration: Euler–Bernoulli / Timoshenko | Нет в показанном элементе | Uncoupled local demonstration; global coupling при сборке | Beam DS, coordinate transformations, перенос к центрам rigid bodies, assembly | Изотропные beams, включая круглые, и rigid bodies | Свободные колебания multibody structures | Conventional DSM, опубликованные результаты, ANSYS FEM | Например k=1 в §3.1; заданный source input | Различие local block structure и structural coupling | Bishop только упомянут в обзоре; нет прямой M-H или Bishop + Timoshenko validation |

## Rucka 2010: что именно объединено

Прочитаны §2.1–2.3, pp. 1761–1764 (PDF 2–5); experimental setup §3,
p. 1764; §4, pp. 1766–1768 (PDF 7–9); §5, pp. 1768–1775
(PDF 9–16); §6, pp. 1776–1778 (PDF 17–19). Обозначения и матрицы
сверены с изображениями PDF 3–5, straight-rod context — PDF 7.

В (1) M-H displacement field имеет вид `u_bar≈u(x)`, `v_bar≈y ψ(x)`.
**В самом PDF стоит ψ, не латинская c**; извлечение текста теряет этот
символ. Это независимое поле поперечного сокращения (смысл `c(x,t)` из
формулировки задачи). В (5) Timoshenko использует `u_bar≈−y φ(x)`,
`v_bar≈v(x)`. Source φ — поворот, source ψ — contraction; эти обозначения
не переносятся на проектный ψ-поворот.

В §2.3 порядок полей `q=(u,ψ,v,φ)ᵀ`. В (8)–(13):

```
ε_MH = (u_x, ψ, ψ_x)ᵀ,       ε_Tim = (v_x − φ, φ_x)ᵀ,
D = diag(D_MH, D_Tim),        E = diag(E_MH, E_Tim),
μ = diag(μ_MH, μ_Tim),
E_MH = [[EA/(1−ν²), νEA/(1−ν²), 0],
        [νEA/(1−ν²), EA/(1−ν²), 0],
        [0, 0, K₁_MH G I]],
μ_MH = diag(ρA, K₂_MH ρI),
E_Tim = diag(K₁_Tim G A, E I), μ_Tim = diag(ρA, K₂_Tim ρI).
```

Это транскрипция структуры **источника**, не утверждённые формулы проекта.
D имеет размер 5×4, E — 5×5, μ — 4×4. Source I — момент площади в данной
плоскости, его нельзя автоматически заменять полярным Ip кругового Bishop.
Ненулевая normal coupling внутри E_MH связывает u и ψ, а не M-H с изгибом.

Метод — **time-domain** spectral elements: полиномы Лагранжа и GLL
quadrature, локальная диагональная mass matrix (25). Это отличается от
frequency-domain SEM Jang. В (26)–(28), p. 1764, элемент переводится в
глобальные координаты и агрегируется. Напечатано `K_global=Tᵀ K_local T`
в конвенции источника; направление T и знаки проекта здесь не меняются.

Эксперимент: straight rod L=1000 mm и сварной L-joint с длиной каждого
участка 995 mm между осями, квадрат 6×6 mm; сталь, ρ=7556 kg/m³ и
E=200.11 GPa измерены, ν=.33 задано. Свободный образец, целое состояние и
надрез; источник возбуждает продольный пакет 120 kHz и изгибный 100 kHz.
В §4 longitudinal и flexural waves сначала сопоставляются с экспериментом
отдельно на прямом стержне, Figs. 5–6. Затем §5 применяет ту же frame
formulation к L-joint. §5.2, p. 1770 (PDF 11), описывает conversion
longitudinal→flexural; §5.3, p. 1774 (PDF 15), flexural→longitudinal.
Глобальная конверсия не противоречит локальным раздельным блокам.

### Поправки Rucka — SOURCE-SPECIFIC / NOT A PROJECT DEFAULT

Точные основания выбора приведены в §2.1–2.2, p. 1762 (PDF 3):

| Параметр | Значение в статье | Как получен | Область |
| --- | --- | --- | --- |
| K₁ᴹᴴ | 1.1 | Совместно с K₂ᴹᴴ: least squares по измеренным скоростям axial wave | Измерения при 100 и 120 kHz для образца статьи |
| K₂ᴹᴴ | 2.1 | Та же экспериментальная идентификация | Тот же диапазон |
| K₁ᵀⁱᵐ | 0.95 | Least squares по измеренным скоростям flexural wave | Измерения при 100 и 120 kHz |
| K₂ᵀⁱᵐ | 12K₁ᵀⁱᵐ/π² | Согласование cutoff с Lamb modes; не отдельный least-squares fit | Выбор в formulation автора |

Figs. 4 и последующие графики не превращают двухчастотную калибровку в
проверку всего диапазона или наших низких собственных частот. Дополнительно
§4, p. 1766, подбирает mass-proportional damping по соотношению амплитуд
отражений: 1000 s⁻¹ для longitudinal, 2000 s⁻¹ для flexural waves.
Сопоставление сигналов поэтому не является проверкой при всех независимо
известных параметрах. Эти damping и correction factors проект не принимает.

## Jang–Park–Lee 2014: кинематика и разные constitutive reductions

Прочитаны §2.1–2.2, pp. 249–250 (PDF 2–3), Appendix A целиком,
pp. 258–259 (PDF 11–12); §5.1, pp. 254–255 (PDF 7–8), описание
dispersion comparison §5.2, p. 255, и Conclusions p. 258.
PDF 2, 11–12 сверены визуально, включая (1), (5), (A1)–(A10).

Для isotropic base beam (1), p. 249:

```
u_b(x,z,t) = u_b0(x,t) − z θ_b(x,t),
w_b(x,z,t) = w_b0(x,t) + z ψ_b(x,t).
```

u_b0 и w_b0 — axial/transverse mid-plane displacements; θ_b — rotation
около y; ψ_b — **независимая** lateral contraction по толщине z.
Ширина b направлена вдоль y. Это planar thickness formulation; она не
равна круговой радиальной кинематике Bishop, где contraction связана с u′.

Appendix A.2, (A5)–(A6), делит полные strains следующим образом (штрих — x):

| Компонента | Axial / M-H contribution (a) | Bending / Timoshenko contribution (b) |
| --- | --- | --- |
| ε_xx | u_b0′ | −z θ_b′ |
| ε_zz | ψ_b | Нет отдельного вклада в (A5) |
| γ_xz | z ψ_b′ | w_b0′−θ_b |

Однако общая кинематическая запись **не означает** применения одного
неизменённого полного 3D-закона к сумме всех strains. В Appendix A.1:

| Блок | Зануляемые 3D stresses | Сохраняемые stresses | Reduced coefficients |
| --- | --- | --- | --- |
| M-H, (A1)–(A2) | σ_yy=τ_xy=τ_yz=0 | σ_xx, σ_zz, τ_xz | C₁₁*=E/(1−ν²), C₁₃*=Eν/(1−ν²), C₅₅*=G; в (A1) второй normal diagonal также C₁₁* |
| Timoshenko, (A3)–(A4) | Те же, **дополнительно σ_zz=0** | σ_xx, τ_xz | Q₁₁*=E, Q₅₅*=G |

Appendix A.3 (A7) прямо записывает `V_b=V_b^(a)+V_b^(b)` как два
stress-work contributions. В основной (5), p. 249, bending rigidity
содержит `I_0b Q₁₁*`, axial normal term — `A_b C₁₁*`; есть M-H internal
term `2 A_b C₁₃* u_b0′ ψ_b`. Для bare centered isotropic base отсутствуют
mixed axial-bending terms в этой записи и в kinetic energy (6), p. 250.
Это характеристика **принятой авторами редуцированной formulation**.

Именно дополнительная stress reduction для bending возвращает Q₁₁*=E,
а не жёсткость буквальной 3D-подстановки из предыдущего Bishop audit.
Отсюда следует полезность источника для выбора будущего reduced closure;
не следует, что проблема прежнего аудита исчезла без изменения предпосылок.
Здесь не утверждается единственность или строгая полная 3D-совместимость
новой combination и не выводятся уравнения CoupledBeams.

Shear correction κ_b появляется в (5) и в resultants (A8)–(A9), p. 259:
`Q_b=κ_b Q₅₅* A_b(w_b0′−θ_b)` и `R_b=κ_b C₅₅* I_0b ψ_b′`.
В **этом** источнике correction применяется к обоим shear contributions;
это не разрешение умножать Bishop H на Timoshenko κ. Численное значение
κ_b/κ_c и процедура его выбора в проверенных местах и поиске по тексту
не установлены; значение 5/6 не приписывается статье и не вводится в проект.

### Предупреждения к печати Appendix A

Изображения pp. 258–259 подтверждают следующие места, которые нельзя
переносить в будущий код буквально:

- В (A3) первая компонента правого вектора напечатана `σ_xx^(b)`, а не
  strain `ε_xx^(b)`. Название stress–strain relation, размерность Q₁₁*=E
  и strains в (A5)–(A7) указывают на опечатку символа. Печатная запись
  сохранена здесь; программная реконструкция не выполнялась.
- В (A10) напечатано `A_b=b∫[-h_b/2,h_b/2] z dz=b h_b`. Этот интеграл
  нечётной функции равен нулю, поэтому две части равенства несовместимы.
  По определению площади множитель z избыточен. PDF не исправлен.
- В (A7) κ_b явно отсутствует, тогда как (5), (A8), (A9) его содержат.
  Место введения correction зарегистрировано; формальный переход от
  uncorrected stress-work к corrected energy требует отдельного аудита
  перед реализацией. Здесь эти записи не объявляются буквально тождественными.

Это ограниченные source warnings, не полный независимый аудит всех формул.

### Patched example и пределы validation

После perfect-bond constraints (3) статья получает шесть coupled fields,
§2.2 (9): u_b0,w_b0,θ_b,ψ_b,θ_c,ψ_c. Односторонняя patch, её материал
и interface constraints создают coupling; её нельзя приписывать обязательной
связи однородной симметричной основы без patch.

§5.1 сравнивает SEM с собственным 1D FEM (Appendix C) и ANSYS **2D
eight-node quadratic plane-stress** model. Table 1 и Fig. 6 относятся к
частотам/FRF cantilever с металлической основой (E=69 GPa, ν=.33,
ρ=2700 kg/m³) и graphite/epoxy patch; Fig. 8 — к guided-wave response
при 100 kHz. Это численная проверка постановки авторов, не выполненный
нами расчёт и не экспериментальная валидация нашей балки. Детальный вывод
composite Appendix B и debonding study не принимаются в область проекта.

## Banerjee–Ananthapuvirajah 2019: modular assembly precedent

Ключ `banerjee_2019_rayleigh_love_timoshenko` сохранён. Дополнительно
прочитаны конец §2, рукопись p. 5 (PDF 6); §2.1, p. 6 (PDF 7), (1),(2);
§2.2, pp. 9–10 (PDF 10–11), (22)–(28); §2.3, p. 15 (PDF 16), (65),(66);
§4.3, pp. 24–25 (PDF 25–26). Это страницы accepted manuscript,
не финальные журнальные pp. 337–347.

Авторы прямо принимают Rayleigh–Love axial и Timoshenko bending как
uncoupled, treated independently. Отдельно получают dynamic stiffness
matrices и в §2.3 объединяют их посредством simple superposition в 6×6
plane-frame DSM. §4.3 применяет её к свободным колебаниям рамы; параметры
заданы через EI, EA, kAG, ρA, ρIp, ν и k=2/3. Это опубликованный пример
модульного сочетания refined axial/bending 1D theories. Он не содержит
M-H equations, независимой contraction DOF или доказательства общей
Bishop–Timoshenko кинематики. Его результаты в этом этапе не пересчитывались.

## Liu et al. 2021: обзор доступных DS и реально показанный случай

Прочитаны title/abstract, §2.1 (PDF 8–9), описание преобразования
beam/rigid-body system в §2.2–2.3 (PDF 9–16), §3.1 representative
validations (PDF 19–21), conclusion (PDF 25–26), Appendix (PDF 26–30).
Локальный файл — accepted manuscript с репозиторной обложкой. Для
указателей используются **PDF pages**, не журнальные страницы; структура
uncoupled matrix (2) сверена по изображению PDF 9.

В §2.1, PDF 8, обзор перечисляет classical, Rayleigh–Love, Rayleigh–Bishop
axial DS и Euler–Bernoulli, Timoshenko, higher-order bending DS. Но текст
перед/после (2), PDF 9, ограничивает **demonstration** uncoupled local case:
axial classical/Rayleigh–Love и bending Euler–Bernoulli/Timoshenko.
Именно эти четыре сочетания обсуждаются в §3.1, PDF 21, и Appendix:
(27) uncoupled DS, (28)–(33) axial, (34)–(49) bending, PDF 26–30.
Обзорное упоминание Bishop не является его демонстрацией с Timoshenko.

Общая DS форма (1) допускает локальную связь, но выбранный пример (2)
её не имеет. Преобразования координат, относительные положения узлов и
центров rigid bodies (§2.2), затем assembly (§2.3) дают structural coupling.
Это полезное разграничение локальных теорий и глобальной структуры, а не
готовые дополнительные условия нашего углового узла.

§3.1 сравнивает собственный метод с conventional DSM для двух beams и
эксцентричного rigid body, затем с опубликованным TMM и ANSYS FEM;
в примерах есть круглые изотропные beams и структуры, заданные EI/EA/ρA.
Заданный k=1 (PDF 19, 21) — параметр тех примеров. Независимая проверка
точности, сравнений быстродействия и всех формул accepted manuscript
здесь не выполнялась. Это не источник M-H equations и не validation
Rayleigh–Bishop + Timoshenko.

## Доступность первоисточников и provenance чтения

Все четыре основных PDF найдены; Banerjee уже был зарегистрирован,
остальные три — новые файлы пользователя. Проверены 76 PDF локального
каталога, включая untracked: имена, извлечённые титульные страницы,
содержание подходящих работ и SHA256. Новых дубликатов не обнаружено.
Ранее существующая пара `A2-4.pdf` / `работа2.pdf` с одинаковым SHA256
не относится к этой регистрации и сохранена. Имена и байты PDF не менялись.

| Cited source | Что доступно | Статус |
| --- | --- | --- |
| Mindlin, R. D.; Herrmann, G. — *A one dimensional theory of compressional waves in an elastic rod* | Rucka, p. 1779, ref. [35]: Proceedings of First US National Congress of Applied Mechanics, 1950, pp. 187–191 | **cited / full text unavailable**; год и выходные данные здесь переданы как напечатанная ссылка Rucka, оригинал не прочитан |
| Martin, M.; Gopalakrishnan, S.; Doyle, J. F. — *Wave propagation in multiply connected deep waveguides* | Rucka, p. 1779, ref. [34]: JSV 174 (1994), 521–538 | **cited / full text unavailable**; содержание оригинала не восстанавливалось по пересказам |

Новых canonical keys/BibTeX-записей по этим двум вторичным упоминаниям
не создано. Уже зарегистрированные локальные
[Liu 2022](source_index.md#liu_2022_longitudinal_dynamic_stiffness) и
[Shatalov et al. 2011](source_index.md#shatalov_2011_longitudinal_rod_theories)
содержат M-H в своей тематике; прежний объём чтения этих записей не расширен
и не заменяет отсутствующие оригиналы. Liu 2022 и Liu 2021 — разные работы.

Метаданные взяты из PDF; адресные внешние проверки ограничены выпуском
Rucka ([издатель](https://doi.org/10.1016/j.jsv.2009.12.004): 329(10)),
названием журнала/томом без выпуска Jang
([издатель](https://www.sciencedirect.com/science/article/abs/pii/S1359836813007804))
и online date Liu ([репозиторий](https://openaccess.city.ac.uk/id/eprint/26468/):
2020-09-16, при citation year 2021). Online dates Rucka 2009-12-30 и Jang
2014-01-03 напечатаны на первых страницах; Banerjee сохраняет ранее
проверенную дату 2018-10-10. Широкого поиска литературы и OCR не было.
Текст и изображения для чтения создавались только во временной папке ОС;
в `results/` новые данные не создавались. Численные расчёты, тесты моделей,
FEM и реализация M-H на этом этапе не запускались.
