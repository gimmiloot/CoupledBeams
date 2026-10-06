# Mindlin–Herrmann + Timoshenko: source audit одного прямого стержня

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
