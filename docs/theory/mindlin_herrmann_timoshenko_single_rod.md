# Mindlin–Herrmann + Timoshenko: source audit одного прямого стержня

Последующий [reduced joint и beta0 transparency gate](mindlin_herrmann_timoshenko_rigid_joint.md)
отдельно принимает published common-DOF closure и подтверждает artificial
interface для того же G20 rod. Source audit и single-rod solver ниже не
менялись; nonzero-angle spectrum/3D joint closure этим не валидированы.

**Текущий отдельный production decision (2026-10-06):**
`PRODUCTION_MHTIM_FORMULATION_SELECTED=PASS`, formulation
`JANG_BARE_ISOTROPIC_REDUCED_MH_TIMOSHENKO`;
`PRODUCTION_MHTIM_KAPPA_RESOLVED` (project rectangular κ=5/6).
`MHTIM_SINGLE_ROD_FINITE_SPECTRUM=PASS`,
`HIERARCHY_SINGLE_ROD=PARTIAL_PASS`: numerical comparison passes, но
пониженные axial models не разрешают отдельную contraction constraint
полной четырёхполевой заделки. [Новый этап, §14–20](#14-production-formulation-decision).
Прежние source/preset audits ниже относятся к своим historical questions
и не переписаны как один общий PASS.

Дата: 2026-10-06. Diagnostic-only, один однородный изотропный прямоугольный
стержень, малая линейная кинематика и одна плоскость изгиба. Литературные
проверки здесь — dispersion/group velocity; конечный собственный спектр
с новыми концевыми условиями не вводится.

**Итог: `MHTIM_VARIANT_DEPENDENT`.** Семейства энергий Rucka и Jang имеют
точное отображение коэффициентов, но опубликованные поправки Rucka не
равны варианту Jang. Статус сравнения этих вариантов:
`MH_SOURCE_VARIANTS_NOT_EQUIVALENT`. Коэффициенты будущей основной модели:
**`PRODUCTION_MH_COEFFICIENTS_UNRESOLVED`**. M-H + Timoshenko остаётся
кандидатом, не новой production baseline CoupledBeams.

Дополнение 2026-10-06: проверка новых локальных PDF завершилась
**`RECTANGULAR_MH_PRESET_MAPPING_UNRESOLVED`**. Ng прямо подтверждает роли
двух обсуждаемых factors, но его normal block отличается от текущего;
ожидаемый Fernandes et al. не найден, второй PDF — circular-theory SSRN
preprint. **`PRODUCTION_MH_COEFFICIENTS_UNRESOLVED` сохранён**, preset не
добавлен. Полная транскрипция и hard gate — в [§13](#13-production-rectangular-m-h-correction-prescription).

## 1. Вопрос и сохранённый Bishop gate

Проверяется опубликованный способ сочетания refined axial и planar bending
для bare isotropic rectangle: поля, энергия, поправки, локальное разделение,
дисперсия и низкочастотный предел. Новый общий 3D-закон или гибридная
кинематика не предлагаются. [Предыдущий Bishop audit](timoshenko_bishop_single_rod.md)
закрыт со статусом `COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED`; он не
переинтерпретируется как доказательство физической связи. [Standalone
Bishop reproduction](bishop_literature_reproduction.md), его код, fixtures
и tests сохранены как diagnostic/reference theory.

M-H вводит самостоятельную contraction coordinate. Jang дополнительно
показывает явное разделение reduced constitutive laws, которым обычный
Timoshenko self term сохраняется. Это литературно определённое замыкание
1D energy, а не утверждение о точном едином 3D-поле для всех эффектов.

## 2. Источники и проверка проекта

Начальное состояние: cwd/git root `D:/PHD/CoupledBeams/CoupledBeams`, main,
HEAD `7ef59a5d8340d46735363d789ec557cca5144422`. Checkout содержит прежнюю
литературную регистрацию и пользовательские изменения другого reviewer
этапа, включая staged theory note, journal/navigation, memory и PDF.
Начальные diff и hashes сохранены во временной папке ОС. Эти изменения
сохранены; memory другой ветви не обновлялась.

| Локальный источник | Реально использованные места | Назначение |
| --- | --- | --- |
| Rucka 2010, `rucka_2010_l_joint_guided_waves`, [PDF](../literature/pdf/j.jsv.2009.12.004.pdf), DOI 10.1016/j.jsv.2009.12.004 | PDF 2–5, pp. 1761–1764, (1)–(13), §3; PDF 7, p. 1766, §4/Fig. 4; (26)–(28) только для понимания assembly | M-H, source factors, Timoshenko и straight-rod dispersion |
| Jang–Park–Lee 2014, `jang_2014_timoshenko_composite_patch_guided_waves`, [PDF](../literature/pdf/j.compositesb.2013.12.050.pdf), DOI 10.1016/j.compositesb.2013.12.050 | PDF 2–3, pp. 249–250, (1),(5)–(7); Appendix A, PDF 11–12, pp. 258–259; §5.1–5.2, PDF 7–9, Fig. 9(a) | Bare base energy и separate constitutive reductions |
| Banerjee–Ananthapuvirajah 2019, [PDF](<../literature/pdf/Rayleigh_Love accepted version.pdf>) | PDF 6, конец §2, manuscript p. 5; PDF 16, §2.3, (65),(66) | Modular Rayleigh–Love/Timoshenko precedent; не M-H validation |
| Liu et al. 2021, [PDF](<../literature/pdf/Journal Paper 122-accepted version.pdf>) | PDF 8–9, §2.1, (1),(2); прежняя регистрация Appendix | Available DS theories и отличие local blocks от global assembly |

Формулы Rucka PDF 3–4 и Jang PDF 2–3,11–12, графики PDF 7/9 сверены по
изображениям без OCR. SHA256, исходные числовые строки, страницы и правила
поправок сохранены в [source fixtures](../../data/input/mindlin_herrmann_timoshenko_sources.json).
Снимки страниц: `results/mindlin_herrmann_timoshenko_literature/source_audit/`.
[Литературная карта](../literature/mindlin_herrmann_timoshenko_sources.md)
и существующие canonical keys не изменены. Нового поиска литературы не было.

Frozen `equations.tex`, раздел «Вид решения», и
`src/my_project/analytic/formulas.py` содержат elementary axial/EB baseline;
общего M-H блока там нет. Он добавлен отдельно. Новый Timoshenko block
берёт EA/EI/KGA/rhoA/rhoI из существующего
[rectangular helper](../../scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py).
Его знаки, API и матрицы не изменяются. Bishop/Rayleigh–Love находятся в
неизменённом `bishop_longitudinal.py`; новая задача не является его preset:
M-H имеет независимую coordinate и двухветвевую дисперсию.

## 3. Координаты и source-faithful поля

В проектном локальном обозначении b вдоль y, h вдоль z, изгиб в (x,z):
`A=bh`, `I=I_y=bh³/12`, `I_z=hb³/12`, `I_p=I_y+I_z`.
**В источниках M-H используется I_y, не I_p.** Source Rucka transverse y
переименована здесь в z; её ψ обозначается c, её v — w, φ — θ.
Jang u_b0,w_b0,θ_b,ψ_b соответствуют u,w,θ,c. Локальный порядок
`q=(u,c,w,θ)`; он не меняет порядок unknowns baseline модели.

Rucka (1): Ux=u, transverse displacement=z c; c независимо от u′.
Её Timoshenko (5): Ux=−zθ, transverse displacement=w.
Jang (1) объединяет их в planar displacement field:

```
Ux = u − z θ,       Uz = w + z c.
```

u,w имеют m, c,θ безразмерны. Положительная растягивающая u′ в
квазистатическом acoustic пределе даёт c=−νu′, то есть contraction.
Это planar thickness formulation со suppressed width stresses, а не
осесимметрическая теория или полный учёт contraction в обеих transverse
directions. Не вводится утверждение Uy=0 совместно с полным 3D Hooke law.

Прямое дифференцирование Jang (A5),(A6) даёт

| Strain | Axial/M-H (a) | Bending/Timoshenko (b) | Total |
| --- | --- | --- | --- |
| ε_xx | u′ | −zθ′ | u′−zθ′ |
| ε_zz | c | Отдельный компонент исключён reduced law | c в (A5) |
| γ_xz | zc′ | w′−θ | zc′+w′−θ |

В Rucka generalized strain order (9): `(u′,c,c′,w′−θ,θ′)`; z² уже
учтено в интегрированном I. При чистом axial поле γ_xz=zc′. В отличие от
Bishop c не является заранее заданной функцией u′.

## 4. Constitutive reductions Jang

Appendix A.1, (A1)–(A4):

| Contribution | Suppressed stresses | Retained stresses | Coefficients |
| --- | --- | --- | --- |
| Axial/M-H | σ_yy=τ_xy=τ_yz=0 | σ_xx,σ_zz,τ_xz | C₁₁*=C₃₃*=E/(1−ν²), C₁₃*=νE/(1−ν²), C₅₅*=G |
| Bending/Timoshenko | Те же, дополнительно σ_zz=0 | σ_xx,τ_xz | Q₁₁*=E, Q₅₅*=G |

Проверка через compliance: при σ_yy=0 normal compliance в (x,z) равна
`E⁻¹[[1,−ν],[−ν,1]]`; её обратная — первая строка таблицы.
При σ_zz=0 остаётся ε_xx=σ_xx/E. Эти переходы проверены точной Fraction
арифметикой. Shear γ — engineering strain, τ_xz=Gγ_xz.

Источник (A7) определяет V_b как **сумму двух stress-work contributions**,
а (5) и (A8),(A9) вводят κ_b в shear energy/resultants. Поэтому bending
rigidity — `Q₁₁* I=EI`, shear rigidity — `κ_b Q₅₅* A=κ_b GA`.
Поперечная strain, eliminated при stress reduction bending, не становится
новой независимой bending DOF. Это источник reduced closure; оно не
равносильно применению одного неизменённого 3D law к total strains таблицы.

Предупреждения прежней регистрации подтверждены: в (A3) RHS печатает σ
вместо ε; в (A10) лишний z в интеграле площади; в (A7) нет явного κ.
Реализована именно **corrected energy (5)**, согласованная с (A9), а не
неисправленная буквальная (A3)/(A10) или uncorrected (A7). Печатные записи
сохранены в fixtures/source index, PDF не исправлены. Это однозначно
названный source variant, без новых composite equations.

## 5. Correction factors и точное отображение

| Коэффициент | Rucka §2.1–2.2, p.1762 | Jang bare energy (5),(6) | Единицы / значение |
| --- | --- | --- | --- |
| K₁ᴹᴴ | 1.1 | κ_b | Dimensionless; multiplies GI contraction-gradient stiffness |
| K₂ᴹᴴ | 2.1 | 1 | Dimensionless; multiplies ρI contraction inertia |
| K₁ᵀⁱᵐ | .95 | κ_b | Dimensionless; multiplies GA bending shear stiffness |
| K₂ᵀⁱᵐ | 12·.95/π² | 1 | Dimensionless; multiplies ρI bending rotary inertia |

Rucka K₁ᴹᴴ,K₂ᴹᴴ совместно fitted методом least squares по axial velocities
при 100 и 120 kHz; K₁ᵀⁱᵐ fitted по flexural velocities при тех же двух
частотах. K₂ᵀⁱᵐ выбран по Lamb cutoff, не отдельный fit. Все четыре:
**SOURCE-SPECIFIC / NOT A PROJECT DEFAULT**. Это empirical corrections к
приближённым fields; в частности, literal velocity integral сам по себе
дал бы K₂=1. Исправленная kinetic energy источника принимается явно.

Jang использует один κ_b в двух shear terms и единичные inertia factors.
Численное κ_b **не установлено** по полному локальному тексту; оно не
восстанавливается из внешнего textbook или формы графика. Для условного
контроля явно задано `--jang-kappa 5/6`. Это документированный test input,
не заявленное значение авторов и не выбор production prescription.

Точное отображение семейства Rucka → Jang при одинаковых E,ρ,ν,A,I:

```
K1_MH = kappa_b, K2_MH = 1,
K1_Tim = kappa_b, K2_Tim = 1,
(u,psi_R,v,phi_R) -> (u_b0,psi_b,w_b0,theta_b).
```

Матрицы при этом равны точно (проверено побайтно для одинаковых floating
inputs). Опубликованный fitted вариант Rucka не удовлетворяет этому
mapping: K₂ᴹᴴ=2.1≠1 уже исключает эквивалентность при любом κ_b.
Масштабирование c не устраняет различие: нормальный c² coefficient и νC
cross coefficient одновременно фиксируют масштаб. Подбора коэффициентов
или density/geometry не выполнялось.

Project Timoshenko имеет r=ρI и Q=κGA(w′−θ). Его source reproduction
Rucka отдельно использует r=K₂ᵀⁱᵐρI и source κ=.95. Коррекции Rucka не
записываются в project helper; production parameter set не создан.

## 6. Энергии, размерность и локальная независимость

Определим локально `C=EA/(1−ν²)`, `H=K₁ᴹᴴ GI`, `m=ρA`, `j=K₂ᴹᴴρI`,
`B=EI`, `S=K₁ᵀⁱᵐGA`, `r=K₂ᵀⁱᵐρI`. Эти H,j — M-H coefficients,
не H,J standalone Bishop. Единицы: C,S — N; H,B — N·m²;
m — kg/m; j,r — kg·m. Все energy densities ниже имеют N=J/m.

```
T = 1/2 ∫ [m*u_t² + j*c_t² + m*w_t² + r*theta_t²] dx,
V = 1/2 ∫ [C*(u_x² + 2*nu*u_x*c + c²) + H*c_x²
             + B*theta_x² + S*(w_x-theta)²] dx.
```

Из Rucka (11)–(13) это следует из заданных D/E/μ blocks. У Jang
potential decomposition следует из (A7)/(5), kinetic — из (6) и
центрированного сечения. Разделение здесь использует **и source closure,
и centroidal symmetry**, а не одну декларацию о симметрии.

| Mixed term между axial и bending | До интегрирования | После интегрирования | Причина отсутствия |
| --- | --- | --- | --- |
| u_t θ_t | −ρz из (u_t−zθ_t)²/2 | −ρQ_z, Q_z=∫z dA | Q_z=0 для centered homogeneous rectangle |
| w_t c_t | +ρz из (w_t+zc_t)²/2 | +ρQ_z | Q_z=0 |
| u_t w_t, u_t c_t, c_t θ_t | Нет произведений разных Cartesian velocity components | 0 | Кинетическая энергия суммы квадратов component velocities |
| u′ θ′, c θ′ | Не включены как axial/bending cross stress-work в (A7) | 0 в source energy | Separate reduced constitutive contributions |
| c′(w′−θ) | Не включён в corrected source decomposition (5) | 0 | Separate shear self contributions; их κ явно задано источником |

∫y dA=∫z dA=∫yz dA=0; ∫z²dA=bh³/12, ∫y²dA=hb³/12
проверены точно. При формальном сдвиге reference axis Q_z≠0 кинетические
cross terms появляются — algebraic negative control, не новая geometry.
Исчезновение axial/bending kinetic coupling требует первого момента
относительно bending coordinate; principal axes нужны для разделения двух
bending directions, которых этот этап не рассматривает.

Normal quadratic form положительно определена при физическом isotropic
ν и C>0; H,S,B>0, mass matrix положительна. Вариационные mixed blocks
между (u,c) и (w,θ) равны нулю. Это свойство принятой **source 1D energy**,
не новая полная 3D-теория.

## 7. Уравнения, концевые переменные и оператор

Из Hamilton variation source energy:

```
m*u_tt = C*u_xx + nu*C*c_x,
j*c_tt = H*c_xx - C*c - nu*C*u_x,
m*w_tt = S*(w_xx-theta_x),
r*theta_tt = B*theta_xx + S*(w_x-theta).

N = C*(u_x+nu*c), R = H*c_x,
Q = S*(w_x-theta), M = B*theta_x.
```

Первое уравнение Rucka (2) использует `2GA/(1−ν)=EA/(1−ν²)=C`;
знак второго соответствует variation c, не Bishop Γ. Концевой член
variation V — `[N δu + R δc + Q δw + M δθ]`; наружные quantities
слева отрицательны, справа положительны. N,Q имеют N, R,M — N·m.
Возможны отдельные essential/natural pairs u/N, c/R, w/Q, θ/M, но
**условия нашего angular joint не выбираются**, включая contraction DOF.

Тот же bare limit читается непосредственно из Jang (9)–(11), PDF 3:
при исключении patch contribution остаются EA_u=EA_ψ=C, GI_ψ=H,
EI_θ=B, κGA=GA_b=K₇=S, K₃=νC, m₁₁=m₂₂=m, m₃₃=m₄₄=ρI.
Остальные patch coupling coefficients отсутствуют. После указанного
notation mapping четыре base equations совпадают с системой выше;
шестиполевой composite operator не реализовывался.

При exp(i(kx−ωt)) mass diag(m,j,m,r). Из общей D* E D матрицы:

```
K_MH = [[C*k², -i*nu*C*k], [i*nu*C*k, C+H*k²]],
K_Tim = [[S*k², i*S*k], [-i*S*k, S+B*k²]].
```

Обе Hermitian. После данного вывода допускается block diagonal operator.
Нули full generalized eigenproblem совпадают с обоими blocks; determinant
factorization проверяется при одной общей mass/spectral normalization.
Не утверждается буквальное равенство с определителями в иных базисах.

## 8. Дисперсия и устойчивый group velocity

Для s=k², λ=ω² characteristic polynomials:

```
MH: C*H*s² + [C²*(1-nu²)-(C*j+m*H)*lambda]*s
             + m*j*lambda²-m*C*lambda = 0,
Tim: S*B*s² -(S*r+m*B)*lambda*s + m*r*lambda²-m*S*lambda = 0.
```

Используется stable quadratic product для малого корня. Нет sinh/cosh,
transfer matrix products, root scans или iterative retries. Обе spatial
roots сохраняются как propagating/evanescent/zero; attenuation положительна
для затухания вдоль +x. В Jang plot imaginary axis показана снизу как
−attenuation·h, что является только display convention.

Group velocity получена аналитическим дифференцированием λ(k), а не
конечной разностью. Независимый контроль — Hellmann–Feynman:
`vg=(v* K′ v)/(2ω v* M v)`. Source labels acoustic/optical продолжены от
k=0 как axial/contraction и bending/shear. При ν=.33 gap M-H ненулевой;
labels не перескакивают между roots и не означают pure polarization при
любом k. Descendant tracker других направлений не применяется.

## 9. Литературные проверки и один общий случай

Rucka: L=1 m, b=h=.006 m, E=200.11 GPa, ρ=7556 kg/m³, ν=.33;
source factors таблицы §5. Fig. 4: 0–500 kHz. L сохранено в конфигурации,
но infinite-wave dispersion от него не зависит. Damping для этих curves
не нужен, time-signal/amplitude reproduction и L-joint не выполняются.

| Рассчитанная величина | Rucka source inputs |
| --- | --- |
| Classical acoustic axial speed | 5146.220866 m/s |
| Contraction cutoff | 345680.639263 Hz |
| Shear cutoff | 262945.871880 Hz |
| Axial vg при 100 / 120 kHz | 5075.467967 / 5039.080004 m/s |
| Bending vg при 100 / 120 kHz | 3012.059713 / 3094.707878 m/s |

На всём [100,120] kHz источник §4 утверждает один longitudinal и один
flexural propagating mode. Оба cutoff выше 120 kHz; analytic roots
подтверждают counts на всём интервале, не только в двух sampled points.
`SOURCE_STATEMENTS=PASS`. Остальные числовые строки таблицы выше — наши
расчёты, не напечатанные измерения или digital author data.
Округлённые source inputs не приобретают дополнительной физической точности
из-за показанных шести десятичных знаков: они нужны для воспроизводимости.
Вид Fig. 4 согласуется визуально; статус `QUALITATIVE_REPRODUCTION` без pixel tolerance.

Jang: bare metallic base b=.020 m, h=.002 m, E=69 GPa, ρ=2700 kg/m³,
ν=.33; Fig. 9(a) 0–2 MHz. κ_b=5/6 задано **только для conditional control**.
Акoustic speed=5055.250296 m/s, contraction cutoff=1476250.985206 Hz —
не зависят от κ_b и определяются source parameters. Shear cutoff при
заданном κ_b=779995.291328 Hz. Восстановлены real/imag kh и group curves
всех четырёх branches. Порядок onset и общая форма согласуются с Fig. 9(a);
статус `CONDITIONAL_QUALITATIVE_REPRODUCTION`. Значение κ авторов не
выведено из этого сходства; exact numeric reproduction всех curves остаётся
неустановленным. Composite patch, Table 1 и FEM не воспроизводились.

Один controlled comparison использует **геометрию и материал Rucka**:

| Величина | Rucka fitted variant | Jang conditional κ_b=5/6 на той же geometry |
| --- | --- | --- |
| Acoustic axial speed, m/s | 5146.220866 | 5146.220866 |
| Contraction cutoff, Hz | 345680.639263 | 500938.837742 |
| Shear cutoff, Hz | 262945.871880 | 264677.171157 |
| Formal M-H high-k speeds, m/s | 2283.675004, 5451.615273 | 2880.427709, 5451.615273 |
| Formal Timo high-k speeds, m/s | 3075.455204, 4788.349792 | 2880.427709, 5146.220866 |

График общего случая ограничен 0–500 kHz: Jang contraction cutoff
чуть выше этого интервала, её propagating branch в нём отсутствует;
evanescent root и cutoff сохранены. `comparison.csv` содержит также
фиксированные k=0,1,10,100,500,1000 m⁻¹ и phase/group velocities.
Эти operator controls и формальный k→∞ не являются high-frequency
экспериментальной validation за пределами source plots. Отличия при общих
E,ρ,ν,A,I вызваны corrections, не constitutive law, geometry или Hz/ω.
Никакая formulation не ранжируется как физически лучшая по fitted example.

## 10. Низкочастотный предел

Acoustic M-H branch:

```
c = -nu*u_x + higher-order terms,
omega² = (E/rho)*k² + nu²*(H-j*E/rho)/m*k⁴ + O(k⁶).
```

Contraction branch имеет nonzero cutoff `sqrt(C/j)/(2π)` и vg→0 при k→0.
Для Timoshenko acoustic bending `omega=sqrt(B/m)*k²+O(k⁴)`,
`vg=2sqrt(B/m)*k+O(k³)`. Shear cutoff `sqrt(S/r)/(2π)`, vg→0.
Оба acoustic branches имеют ω→0; это не две новые finite-rod mode lists.
У источников planar I, поэтому higher-order axial asymptotics не
отождествляются с standalone circular Bishop/Rayleigh–Love через I_p.

Все low-k controls kh=.001,.0001,.00001 подтверждают эти leading limits.
Positive corrected energies обеспечивают реальную неотрицательную λ(k)
в данной undamped source theory. Никаких выводов о coupled spectrum нет.

## 11. Численная проверка и воспроизводимость

Допуски записаны **до вычисления** в fixtures: spectral-scale eigenvalue
error 3e−10, mass-normalized equation residual 1e−10, HF group relative
error 3e−10, spatial roundtrip 5e−10, polynomial residual 1e−10,
existing project spatial roots 3e−11. Полиграфические графики numeric
tolerance не получают.

Независимая матрица строится из D* E D; analytic polynomial roots не
получаются вызовом того же eigensolver. На фиксированных k controls:
max equation residual ≤1.82e−16, eigenvalue spectral-scale error ≤3.52e−16,
group/HF relative error ≤1.10e−11, spatial roundtrip ≤4.16e−10.
Raw relative errors малых acoustic eigenvalues и spectral_scale/λ
сохранены: нормировка на spectral scale не заявляет высокой relative
точности плохо обусловленного малого eigenvalue у generic eigensolver.
Project Timoshenko spatial roots проверены ниже, на и выше cutoff.

Есть точные rectangle/compliance/cross-term checks, positive mass,
Hermitian/full block operator, determinant factorization, неизменность
bending при изменении только M-H factors, source counts и cache tests.
Полный тяжёлый набор исследований не запускается. SymPy/Lean недоступны;
использованы прямой вывод, Fraction и независимый numerical eigenproblem.
Начальный тестовый прогон выявил strict float equality на одном ulp и
eager formatting заголовка Rucka; исправлены тест с **прежним** fixed
numerical contract и renderer. Коэффициенты/численные допуски не менялись.

Финальный прогон: **74 tests PASS** (26 новых + 48 прежних Bishop literature/
kinematics), дополнительно **3 existing rectangular Timoshenko regressions
PASS**, итого 77. Source-check, compute, plot-only и matching-cache reuse
прошли; повторный compute выполнил zero root evaluations. `git diff --check`
без ошибок. Рабочее окружение: Python 3.12.4, NumPy 2.1.3, SciPy 1.15.2,
Matplotlib 3.9.2; фактический Python/version manifest имеет приоритет.
Проверенный interpreter: `D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe`;
пакеты не устанавливались, рабочее окружение не изменялось.

CLI: [один entry point](../../scripts/analysis/reproduce_mindlin_herrmann_timoshenko_literature.py),
[один helper](../../scripts/lib/mindlin_herrmann_longitudinal.py),
[targeted tests](../../tests/test_mindlin_herrmann_timoshenko_literature.py).
Команды из корня, в установленном проектном NumPy/SciPy окружении:

```powershell
python scripts/analysis/reproduce_mindlin_herrmann_timoshenko_literature.py --check-sources
python scripts/analysis/reproduce_mindlin_herrmann_timoshenko_literature.py --compute --case rucka
python scripts/analysis/reproduce_mindlin_herrmann_timoshenko_literature.py --compute --case all --jang-kappa 5/6
python scripts/analysis/reproduce_mindlin_herrmann_timoshenko_literature.py --plot-only --case all --jang-kappa 5/6
python -m pytest tests/test_mindlin_herrmann_timoshenko_literature.py tests/test_bishop_literature.py tests/test_timoshenko_bishop_single_rod.py -q
python -m pytest tests/test_rectangular_isotropic_models_vs_beta.py -k "old_timoshenko_and_four_ply_rlb_reduce_to_one_section_contract or frozen_old_vs_rlb_roots_match" -q
```

`--variant rucka_2010` / `--variant jang_2014_bare_isotropic` выбирает source
case. Jang без `--jang-kappa` сохраняет unresolved record с причиной и
возвращает exit 2, не выбирая число. Output root задаётся `--output-dir`.
Datasets: `results/mindlin_herrmann_timoshenko_literature/<fingerprint>/<case>/`;
`current.json` указывает актуальный запуск. Там manifest, exact parameters,
hashes/versions/Git/commands, roots/attenuation/group/phase CSV, reference
elementary/EB CSV, cutoffs, comparisons, numerical diagnostics и PNG/PDF.
Plots только Rucka Fig. 4, conditional Jang Fig. 9(a), один common comparison.
Графики читают saved CSV; tests запрещают compute в plot-only и reuse.

Fingerprint включает fixture, equations/helper/CLI versions, source hashes,
Python/dependencies, Git HEAD и явный κ input. Изменённые входы дают новый
namespace; cache требует совпадения identity и artifact hashes. SI parameters
и factors отдельно сверяются с printed strings: случайно изменённая geometry
не сохраняет source-reproduction label. Старые
bundles после технических исправлений сохранены как provenance, но не
читаются как актуальные. Поиск не повторяется: fixed analytic roots,
точные in-range cutoff points, zero iterative searches/retries. Исключение
сохраняет `failure.json` с case/frequency/причиной и прекращает этот случай.

## 12. Решение о будущей базе и ограничения

Опубликованная reduced linear energy математически определена, positive
и локально блочна; её параметризуемая diagnostic реализация пригодна для
последующего явного выбора coefficients и BC. **Production promotion
на этом этапе не выполнен.** Таблица §5 содержит разные prescriptions;
κ Jang не установлен численно, fitted 100/120 kHz Rucka не являются
универсальным low-frequency набором. Выбор coefficients — отдельное
научное решение, а не результат подгонки или успешного benchmark.

Banerjee подтверждает modular architecture Rayleigh–Love + Timoshenko,
но не M-H accuracy. У Rayleigh–Love три plane-frame DOF на конце и 6×6
DSM; у source M-H/Timoshenko четыре DOF на конце (включая c), то есть
возможный двухконцевой DS имел бы размер 8×8 до отдельно обоснованного
исключения переменных. Такая DS/condensation здесь не реализована.
Liu различает local independence и global structural
coupling; его фактические demonstration axial theories — classical/RL,
не Bishop или M-H. Точный источник первой M-H theory этим этапом не
заменяется. Source planar approximation не является 3D truth для любой
ширины/толщины rectangle.

**Angular-joint M-H conditions NOT yet selected. No coupled-beam
calculations yet.** Два стержня, L-joint production model, nonlinear
equations, другие плоскости, torsion, anisotropy, damping и 3D FEM не
реализованы/не запускались. Прежние scientific solvers, Bishop benchmarks,
article results, viscosity/Yartsev/nonlinear branches не изменены.

## 13. Production rectangular M-H correction prescription

### 13.1. Фактические источники и scope нового gate

Начальный main/HEAD `0f961acb227b504439408d3484b13a72a7acf3be`, tracked
diff пуст; untracked `hdl_85788.pdf` и `ssrn-5985611.pdf` сохранены.
Ng 2014 [accepted PDF](../literature/pdf/hdl_85788.pdf), 41 страница,
соответствует заданию. Второй файл — [Elishakoff–Tharu SSRN preprint](../literature/pdf/ssrn-5985611.pdf),
100 страниц, не Fernandes–Machado–Dutkiewicz 2022. Поиск по локальным
PDF проекта (92, включая ignored/untracked, title/DOI/content matching)
не нашёл DOI `10.3390/en15207725` или ожидаемое название Energies.
Этот источник не считается прочитанным или зарегистрированным по PDF.
Metadata, SHA256, объём чтения и предупреждения — в [source index](../literature/source_index.md#rectangular-m-h-prescription-audit-2026-10-06).

Ng: PDF 9, (1), `u_j≈bar(u_j)(x,t)`, `v_j≈y*bar(phi_j)(x,t)`.
Source bar(phi_j) соответствует независимой c, несмотря на словесное
rotational angle: это множитель поперечной координаты, не Timoshenko
rotation. На PDF 10 в (2) и непосредственно следующем абзаце напечатано:

```
(2*mu_j+lambda_j)*A_j*u_j,xx + lambda_j*A_j*phi_j,x = rho_j*A_j*u_j,tt
mu_j*I_j*S1*phi_j,xx - (2*mu_j+lambda_j)*A_j*phi_j
                     - lambda_j*A_j*u_j,x = rho_j*I_j*S2,j*phi_j,tt
mu_j=E_j/(2*(1+nu_j)), lambda_j=nu_j*E_j/((1+nu_j)*(1-2*nu_j))
A_j=b_j*h_j, I_j=b_j*h_j^3/12
S1=12/pi^2
S2,j=S1*((1+nu_j)/(0.87+1.12*nu_j))^2
```

Визуально проверены PDF 1–2,9–11,33. В (6), PDF 11, contraction-gradient
entry напечатан без k_j², в отличие от второй производной (2),(4).
Новый characteristic solver Ng не строится; mapping опирается на (2).
Численные эксперименты и Bayesian identification этой статьи не повторены.

### 13.2. Exact mapping factors и отдельное constitutive различие

| Source | Symbol / source expression | Section / place in equation | Project energy coefficient | Mapping / status |
| --- | --- | --- | --- | --- |
| Ng, PDF 10 (2) + following paragraph | S1=12/π² | Rectangle, I=bh³/12; μIS1 c_xx | H=K_MH1 GI | K_MH1=S1 directly because μ=G; **factor-role confirmed** |
| Ng, same page | S2=S1[(1+ν)/(.87+1.12ν)]² | Rectangle; ρIS2 c_tt | j=K_MH2 ρI | K_MH2=S2 directly; **factor-role confirmed** |
| Ng, (2) | D=(2μ+λ)A; F=λA | Normal diagonal / internal u_x c coupling | Current D=C, F=νC | **Different normal block** at fixed physical E,ν |
| Elishakoff–Tharu, PDF 20 (37),(38) | κ², κ₁² in circular stress-displacement relations | Circular axisymmetric r,z,a; u radial, w axial | No direct rectangular H/j mapping established | **Not a second confirmation**; do not identify squared symbols with S1/S2 |
| Expected Fernandes et al. 2022 | Not inspected | Local full text absent | Not established | **Unavailable**, no source substitution |

Таким образом, у Ng S1 и S2 не являются sqrt(K), K², 1/K или только
нормированными frequency factors. Их места в PDE устанавливают
соответствующие H,j в энергии. Но для всей normal energy из (2) следует:

```
T_Ng = 1/2 integral [rho*A*u_t^2 + rho*I*S2*c_t^2] dx
V_Ng = 1/2 integral [D*(u_x^2+c^2) + 2*F*u_x*c + G*I*S1*c_x^2] dx
D=(2*mu+lambda)*A, F=lambda*A
```

Это восстановление из напечатанных PDE, не цитата напечатанной энергии.
`D/EA=(1−ν)/[(1+ν)(1−2ν)]`, `F/EA=ν/[(1+ν)(1−2ν)]`.
Ratio `F/D=ν/(1−ν)` отличается от project `ν`. Нормированный mixed
coefficient `F/sqrt(D*D)` инвариантен при отдельном масштабировании u,c,
поэтому простое переименование/масштабирование DOF не устраняет различие.
При stationary contraction `c=−F*u_x/D`:

```
Ng:      D-F^2/D = EA/(1-nu^2)
project: C-(nu*C)^2/C = EA
```

Следовательно, source equation (2) как напечатана даёт acoustic speed
`sqrt(E/[rho*(1−nu²)])`, а текущая reduced model — `sqrt(E/rho)`.
Это algebraic convention audit, не обвинение статьи в ошибочности и не
новый dispersion benchmark. Различие исчезает при ν=0, но не при общем ν.
Точные Fraction checks подтверждают ratios и stationary reduction.
Lamé λ источника не заменён молча reduced plane-stress коэффициентом;
physical E,ν не перенормированы. Frozen equations.tex/analytic baseline,
current energy core и signs/unknown ordering сохранены.

### 13.3. Происхождение, статус и нерешённое условие переноса

Ng задаёт S1,S2 как published rectangular model formula, до specimen
damage identification; fit этих factors не описан. S1 не содержит b/h,
S2 содержит только ν. Это отличается от fitted Rucka factors. Два
напечатанных decimal constants .87 и 1.12 не объявляются точными
рациональными физическими константами. Prescription не названа уникально
правильной или оптимальной. Повторное употребление именно этой formula
двумя новыми доступными rectangular publications пока не установлено.

Ng PDF 33 ref. **37** цитирует Doyle 1997, 2nd ed. Elishakoff–Tharu PDF 95
ref. **10** цитирует Graff 1976 (как напечатано), ref. **15** — Doyle 1997.
Книги отдельно не проверены; circular review не заменяет ни оригинал книги,
ни отсутствующий Fernandes PDF. Его §2.2 использует другие геометрию и
систему squared corrections; все результаты препринта не валидированы.

**Итог нового этапа: `RECTANGULAR_MH_PRESET_MAPPING_UNRESOLVED`;
`PRODUCTION_MH_COEFFICIENTS_UNRESOLVED`.**
Candidate name `rectangular_literature_default` записан только в fixtures
как **unadopted candidate**, `adopted=false`, `code_preset=null`. Helper/
production constructor не добавлен. Требуются второй непосредственно
проверенный rectangular full source и явно обоснованный перенос к
текущему normal reduced block; совпадение H/j roles отдельно не доказывает
всей source-consistent energy. Никакой новой constitutive теории не выбрано.

Печатные формулы candidate при ν=.3 дают вычисляемые arithmetic sanity
values; они сохраняются в local audit, не становятся model defaults или
physical precision claims. Production cutoff/low-k sanity run **не
выполнен**, поскольку preset не принят. Существующие low-k и Timoshenko
regressions проверяются targeted tests без нового parameter study.
При ν=.3 arithmetic output: S1=1.2158542037080533,
S2=1.4127769143961029, вычислены из formula, не hard-coded в model.

M-H corrections остаются отделены от project Timoshenko κ и rotary
factor=1. Source Rucka/Jang cases, Jang `kappa_b_numeric=null` и тест
`test_jang_kappa_is_never_guessed` сохранены. Никакой coefficient fit,
two-beam assembly, contraction angular-joint BC, nonlinear derivation,
FEM или возврат к Bishop hybrid не выполнялись. Старые results не
пересчитываются; source-check проверяет также hashes новых records.

Новое локальное evidence находится в
`results/mindlin_herrmann_timoshenko_literature/rectangular_preset_audit/`:
page snapshots, source inventory/manifest и audit/sanity/verification JSON.
Это source audit, не ещё один numerical benchmark. README и script guides
не требуют нового workflow: CLI и helper API не изменены.

Проверки нового gate: **82 tests PASS** (31 M-H source/audit, 48 прежних
Bishop literature/kinematics, 3 существующих rectangular Timoshenko).
Fixed-point Rucka/Jang temporal/spatial outputs и model coefficients
совпали точно в JSON representation; прежние source records, cases и
numeric contract не изменены. Первый технический snapshot comparison
сравнивал Python tuple units с JSON list; после одинаковой сериализации
различий нет, numerical tolerances не вводились. Source hash-check passed;
защищённые PDF/baseline/CLI файлы побайтово сохранены. Git HEAD/staging
не изменены. Нет нового numerical benchmark или production sanity model.

## 14. Production formulation decision

Пользователь явно выбрал published Jang bare isotropic reduced closure
как основу последующей внутриплоскостной модели. Named preset
`project_jang_reduced_rectangular` отделён от `jang_2014_bare_isotropic`.
Это не установление численного source κ Jang, не fit и не утверждение
единственности/оптимальности теории. Rucka/Jang source cases и standalone
Bishop остаются прежними. `MHTIM_VARIANT_DEPENDENT` и
`MH_SOURCE_VARIANTS_NOT_EQUIVALENT` остаются истинными для source comparison.

Новый [Fernandes full text](../literature/source_index.md#fernandes_2022_spectral_tower_cable),
DOI 10.3390/en15207725, PDF/p.6–7 (19)–(23), повторяет Ng factors 12/π²
и ν-formula, но использует `(2mu+lambda)A,lambda*A`. Он поддерживает
provenance альтернативной ветви, не выбранный project normal block.
В этих местах I назван cross-section inertia, без rectangular b,h
prescription для source specimen; tower case не является нашим rectangle.
Ref.32 p.25 цитирует Doyle1997, книга не проверена. Historical
`RECTANGULAR_MH_PRESET_MAPPING_UNRESOLVED` сохраняется для переноса всего
Ng-type block. Unadopted candidate в старых fixtures не становится новым
production preset; общий production status прежнего этапа был historical.
Направление Ng/Fernandes/Liu-type сохраняется как отдельная literature
alternative. Здесь Liu означает зарегистрированный Liu2022 по longitudinal
DS (rectangular M-H вынесен в Appendix), не Liu2021 multibody framework.
Его Appendix не используется для выбора нового production closure и не
объявляется заново независимо воспроизведённым; в этом этапе прямое
подтверждение correction expressions получено по Ng/Fernandes.

## 15. Project kappa provenance и вариационные уравнения

Frozen [G20 contract](../../tests/data/reddy_four_ply_isotropic_limit_cases.json),
material.K=5/6, используется existing rectangular Timoshenko comparator.
Accepted provenance зафиксирована в
[RLB-2B note §2–3](../laminated_beams/rectangular_isotropic_models_vs_beta_note.md)
и проверяется `test_canonical_contract_and_four_equal_isotropic_plies` в
`test_rectangular_isotropic_models_vs_beta.py` (точное K=5/6).
Библиографический precedent: Kramer–Gfrerer2024, local `A2-1.pdf`, p.2
§2.1 после (4): κ=5/6 для rectangle, (5) γ=w′−Θ. Страница дополнительно
сверена визуально. Project helper сохраняет `Q=KGA*(w′−psi)` и
`KGA=K*G*A`; его API/rotary inertia не меняются. κ установлен **проектным
contract**, не догадкой о Jang. Один κ=5/6 входит в S=κGA и H=κGI;
`j=r=ρI` точно, `C=EA/(1−ν²)`. Section с иным K не принимается preset.

Для q=(u,c,w,θ), A=bh, I=I_y=bh³/12, production energy — §6 с
`K_MH1=K_Tim1=κ_project`, `K_MH2=K_Tim2=1`. Independent exact
quadratic differentiation (Fraction tests) даёт:

| Variation | Potential boundary coefficient | Volume coefficient in δV |
| --- | --- | --- |
| δu | N=C(u′+νc) | −N′ |
| δc | R=Hc′ | C(c+νu′)−R′ |
| δw | Q=S(w′−θ) | −Q′ |
| δθ | M=Bθ′ | −Q−M′ |

δT после интегрирования по времени даёт −m u_tt, −j c_tt, −m w_tt,
−r θ_tt. Из Hamilton variation получены:

```
m*u_tt = C*u_xx + nu*C*c_x
j*c_tt = H*c_xx - C*c - nu*C*u_x
m*w_tt = S*(w_xx-theta_x)
r*theta_tt = B*theta_xx + S*(w_x-theta)
N=C*(u_x+nu*c), R=H*c_x, Q=S*(w_x-theta), M=B*theta_x
```

Энд-член δV: `[Nδu+Rδc+Qδw+Mδθ]`; outward sign слева минус, справа плюс.
Mass/potential mixed blocks между (u,c) и (w,θ) нулевые из выбранной
source energy и centroid symmetry (§6), до spectrum union. Quasistatic
`c=−νu′` возвращает EA, acoustic speed √(E/ρ). Frozen equations.tex и
analytic baseline не содержат M-H и не меняются.

## 16. Finite single-rod boundary problem и hierarchy caveat

Один existing normalized G20 input: E=ρ=1, ν=.3, b=.20, h=.05,
L=L_ref=1; A=.01, I=.000002083333333333334. Это mathematical benchmark,
не измеренный материал и не HMS-DX209 specimen. L_total=2/two-arm assembly
не используется. f=ω/(2π) относится к formal model time units contract;
dimensionless `f*=f*L/√(E/ρ)`, здесь f*=f численно. Geometry не подбиралась.

На обоих концах **u=c=w=θ=0**. Source (1), PDF2 p.249:
`Ux=u−zθ`, `Uz=w+zc`; vanishing displacement при двух различных z
однозначно даёт все четыре conditions. Это full clamp **в planar
four-field approximation**, не утверждение о полном 3D Dirichlet model.
Jang (12), PDF3 p.250, явно включает essential δψ_b=0 как alternative
natural R_b pair; (28), PDF5 p.252, задаёт ψ_b1/ψ_b2 как end DOF.
Fig.5 PDF7 p.254 содержит clamped cantilever precedent. Новый bare CC
control не является воспроизведением composite Table1 или source bare CC
таблицы. Contraction constraint resolved: независимая c подавлена торцом;
это **не u′=0** в M-H. Natural free pairs — N=R=Q=M=0, без нового free run.

В hierarchy elementary/Rayleigh–Love берут те же section/material и
Timoshenko block. Planar Rayleigh–Love получается из thickness contraction
c≈−νu′: gradient energy MH опущена, lateral inertia сохранена;
`J_RL=ν²ρI_y`, `H=0`. Это **planar Rayleigh–Love approximation**, не
full polar/two-transverse contraction circular reference. Elementary
H=J=0; H не заменён epsilon. RL equation/effort:
`m u_tt−J_RL u_xxtt−EAu_xx=0`, `N_RL=(EA−J_RLω²)u′`.
Axial endpoint condition reduced models — только u=0 на каждом конце.

**Hierarchy physical qualification:** reduced theories не имеют independent
c, поэтому suppression c=0 full-clamp boundary layer не разрешается.
Перенос c=−νu′ к торцу дал бы u=u′=0, overconstrained для второго порядка;
он не выполняется. Numerical hierarchy использует model-specific
representations fixed-end fixture с теми же centroidal/bending clamps.
Exact одинаковая four-field/3D clamp kinematics A/B/C не доказана:
`HIERARCHY_SINGLE_ROD=PARTIAL_PASS`, хотя numerical gates проходят.
Это не мешает PASS finite M-H problem и не выбирает angular-joint BC.

## 17. State equations и exact bounded solution

При exp(iωt), y_MH=(u,c,N,R), y_T=(w,θ,Q,M):

```
u′ = N/C - nu*c              w′ = theta + Q/S
c′ = R/H                    theta′ = M/B
N′ = -m*omega²*u            Q′ = -m*omega²*w
R′ = (EA-j*omega²)*c+nu*N   M′ = -r*omega²*theta-Q
```

R′ следует из C(c+νu′)−jω²c, затем u′ definition; EA=C(1−ν²).
Общий order (u,c,w,θ,N,R,Q,M) даёт два blocks после energy derivation.
Baseline unknown ordering не меняется.

Primary: exact spatial polynomial §8, independent PDE amplitudes,
acoustic cos/sin и **bounded anchored exponentials** `exp(−αx)`,
`exp(−α(L−x))`. Positive length/row/column scaling — essential matrix4×4.
Нет неограниченных sinh/cosh или больших transfer products. Finite solver
намеренно ниже optical cutoff; ν=0 имеет отдельные decoupled columns.
Independent: first-order state из resultants, initial essential q=0 и
два независимых efforts; short-step `scipy.linalg.expm`, impedance scaling
и positive-diagonal QR продвигают двумерное пространство. Step exponent≤1,
budget512steps, projected boundary matrix2×2. Два представления одной 1D
theory, не две обёртки общей boundary matrix или experimental validation.

## 18. Root completeness и numerical verification

[Config](../../data/input/mindlin_herrmann_timoshenko_single_rod.json) задан
до расчёта: bounded interval/block, 400 scan intervals, root xtol1e−11 /
rtol1e−12; independent rel2e−8, scaled BC/ODE1e−9, energy5e−8,
mass orthogonality5e−7, nonzero singular condition≤1e8. Один doubled scan
только при count failure, без range/tolerance expansion; не потребовался.

Полнота не выводится из sign scan. Young inequality при η=.18>ν²:

```
C(u′²+2nu*u′c+c²)+Hc′²
 >= C(1-eta)u′² + C(1-nu²/eta)c² + Hc′².
```

Правый lower form имеет две exact scalar Dirichlet spectra. Min-max
upper eigenvalue count до search ceiling=7 (contraction lower count0);
verified MH roots=7, включая guard. Для Timoshenko ослабление θ essential
до natural M=0 даёт exact simply-supported lower spectrum: k=nπ/L,
две dispersion branches и uniform θ shear mode. Upper count=11, verified
CC roots=11 с guard. Search starts ниже lower first roots (count0).
Насыщение upper counts исключает пропущенные roots в заявленных ranges.
Это небольшой model-specific min-max gate, не WW/general framework.
Guard coverage подтверждает combined first12 prefix.

Max MH / Timo scaled BC residual: 1.55e−13 /1.75e−12;
ODE: 1.49e−15 /1.07e−14; energy relative: 2.00e−13 /2.18e−12;
mass Gram: 2.77e−13 /1.06e−12; primary/independent frequency difference:
1.06e−13 /1.71e−12. Nonzero singular condition≤1.990 /1.681.
Это actual diagnostics, не ужесточённые после расчёта tolerances.
Exact reductions проверены Fraction; SymPy/Lean отсутствуют, без installs.

## 19. Hierarchy comparison и spectral family inventory

| Family n | Elementary axial f* | Planar Rayleigh–Love f* | M-H f* | Common Timoshenko bending f* |
| --- | --- | --- | --- | --- |
| 1 | .500000000 | .499953743 | .500706144 | .050527603 |
| 2 | 1.000000000 | .999630095 | 1.001229213 | .136325767 |
| 3 | 1.500000000 | 1.498752436 | 1.501381967 | .260072546 |
| 4 | 2.000000000 | 1.997045678 | 2.000968673 | .416258846 |
| 5 | 2.500000000 | 2.494237017 | 2.499780436 | .599831008 |
| 6 | 3.000000000 | 2.990056680 | 2.997590005 | .806045922 |

Bending вычисляется один раз и переиспользуется во всех A/B/C; self terms
и existing-basis roots checked. Combined boundary determinant в общих
bounded bases имеет block factorization; zeros равны union **после**
energy decomposition. Равенство с иными нормировками determinants не заявлено.

First12 MH combined positions: B1,B2,B3,B4,A1,B5,B6,A2,B7,B8,A3,B9;
верхний f*=1.522253495. Family index отделён от sorted position.
Contraction **wave cutoff** f*=11.558994422, shear wave cutoff6.242570465,
оба выше bounded low inventory. Independent propagating contraction/shear
optical branches в нём отсутствуют. MH shapes содержат evanescent
contraction boundary contribution: c не ноль внутри rod. Wave cutoff
не называется первой finite CC optical eigenfrequency; upper finite
spectrum не рассчитывался. Contraction lower count scale8.173443340
тоже выше MH search ceiling. Более полная theory не обязана понижать
все frequencies; resolved clamp effect виден в первых MH roots, без fit.

## 20. Reproduction и remaining joint question

[Один finite CLI](../../scripts/analysis/verify_mindlin_herrmann_timoshenko_single_rod.py)
имеет новый finite-boundary/count contract, поэтому отделён от source
CLI. Reuses existing MH helper, rectangular section, exact Bishop H=0
boundary machinery и atomic writers. Source Jang требует explicit κ.

```powershell
python scripts/analysis/verify_mindlin_herrmann_timoshenko_single_rod.py --check-sources
python scripts/analysis/verify_mindlin_herrmann_timoshenko_single_rod.py --compute
python -m pytest tests/test_mindlin_herrmann_timoshenko_finite_rod.py tests/test_mindlin_herrmann_timoshenko_literature.py tests/test_bishop_literature.py tests/test_timoshenko_bishop_single_rod.py -q
```

Проверенный interpreter D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe,
Python3.12.4, existing NumPy/SciPy environment. Result root
`results/mindlin_herrmann_timoshenko_single_rod/<fingerprint>/`: manifest,
exact config/hashes/versions/Git, full result/count/diagnostics JSON,
hierarchy и mass-normalized profiles CSV. current.json указатель;
reuse проверяет identity/artifact hashes, zero root evaluations.
Page snapshots отдельно source_audit/. Новых plots/maps нет.

Финальные targeted checks: 19 новых finite +31 source MH +48 прежних
Bishop literature/kinematics =98 PASS; ещё4 existing rectangular Timo
regressions (включая exact G20 κ contract) PASS, итого102. Smoke source/
finite check/compute и matching-cache reuse прошли. Rucka/Jang fixed-point
outputs и прежние cases/source records/numeric contract сохранены точно;
13 protected initial files, включая три пользовательских PDF, baseline
helpers и source CLI, побайтово сохранены. Initial user diff сохранён
во временной папке ОС; прежние незакоммиченные source audit edits не удалены.

Production closure/project κ выбраны для дальнейшей работы; finite MH
gate passes, hierarchy qualified PARTIAL_PASS по clamp interpretation.
**Angular-joint conditions для c/R не выбраны.** Two-beam M-H assembly,
β, nonlinear equations, out-of-plane/torsion, coefficient fitting и
3D FEM не выполнялись. Ng/Fernandes variant, Rucka/Jang semantics и
закрытый Bishop gate сохранены. Автоматического coupled этапа нет.
