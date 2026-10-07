# M-H–Timoshenko: reduced rigid joint и beta=0 transparency

Последующий general-angle этап (2026-10-06), §§9–16:
**MHTIM_GENERAL_BETA_JOINT_GATE=PASS**. Одна общая assembly использует
действующий project β contract; structural gates, frozen beta0 regression
и fixed spectral pilot5°/45°/90° прошли. Разделы1–8 ниже сохраняют
исторический beta0 gate: их statements о невычисленном β≠0 относятся
к тому этапу. Новая проверка не является 3D joint proof или applicability
study; оговорка о c-continuity сохраняется без изменений.

Дата: 2026-10-06. Все семь раздельных gates ниже — **PASS**. Это один
однородный прямой fixed–fixed стержень G20, искусственно разбитый на два
плеча. **Спектр при β≠0 не рассчитывался и не объявлен валидированным.**

## 1. Scope и научная квалификация

По явному заданию принят published reduced frame closure: общий вектор
перемещения центра, общий поворот сечения и общая contraction coordinate.
Дополнительная совместность c принята из вариационной/сборочной структуры
опубликованной Mindlin–Herrmann–Timoshenko frame formulation. Это замыкание
редуцированной одномерной модели; оно **не заявляется прямым выводом из 3D
упругости конечной области реального сварного узла**. Не утверждается и
обратное, будто это условие ошибочно. При разных направлениях плеч c
относится к разным локальным поперечным осям; их общая скалярная DOF —
именно выбранное замыкание, а не доказанная локальная 3D кинематика узла.

Теория плеча, preset, коэффициенты, external clamps и direct reference
из [single-rod note](mindlin_herrmann_timoshenko_single_rod.md#14-production-formulation-decision)
сохранены. Ng/Fernandes 12/π² не принят; Rucka/Jang source prescriptions,
исторические variant-dependent/mapping статусы и закрытый Bishop audit
не изменены. Frozen equations.tex и src/my_project/analytic не содержат
этого M-H узла; расширение локализовано здесь, baseline не переписан.

## 2. Проверенные source equations и assembly inference

Локальные издательские PDF; metadata/citation keys не менялись:

| Source | Места и что прямо напечатано | Что из этого выведено здесь |
| --- | --- | --- |
| `rucka_2010_l_joint_guided_waves`, [PDF](../literature/pdf/j.jsv.2009.12.004.pdf) | PDF4/p.1763 (10)–(15): weak form, block energies/masses, node q=(u,ψ,v,φ); PDF5/p.1764 (26)–(28): nodal T, transformed stiffness/mass/load и standard aggregation; §5 pp.1768–1775 применяет тот же элемент к L-joint | ψ→c, φ→θ. Общие assembled scalar nodal DOFs дают c₁=c₂ и θ₁=θ₂; их сопряжённые nodal forces суммируются. **Скалярная формула c₁=c₂ отдельно не напечатана.** |
| `jang_2014_timoshenko_composite_patch_guided_waves`, [PDF](../literature/pdf/j.compositesb.2013.12.050.pdf) | PDF2/p.249 (1),(5), PDF3/p.250 (6),(12): fields, energies, natural/essential ψ_b pair; PDF5/p.252 (28)–(30),(36): end DOFs и weak boundary work; Appendix A pp.258–259 resultants | Bare base retains (u,c,w,θ) and pair (ψ_b,R_b). Jang поддерживает физику плеча/концевые пары; frame aggregation Rucka не приписывается Jang. |

Изображения PDF4–5 Rucka и PDF3,5 Jang повторно проверены вручную и
сохранены в `results/mindlin_herrmann_timoshenko_beta0_joint/source_audit/`.
Rucka SHA256 `d15786c6aa2bf71731c63b16bf594739a626e6c45461e63689d90a36d472ce5f`;
Jang SHA256 `8d945ae828f5a794090df67b968ca84901ff7e3441e6a2d5a4435c025b5c1929`.
Новых источников/web поиска/OCR нет. Jang patch terms не перенесены.
Rucka §5.1 также моделирует сварную область пониженным E; §5.2 применяет
frame element к mode conversion. Эта source-specific welded-zone модель,
fitted corrections и damping здесь не воспроизводятся: gate относится к
homogeneous artificial interface, а публикация даёт assembly precedent.

Транскрипция Rucka (26), в порядке (u,ψ,v,φ):

```text
T_i = [[cos α, 0, sin α, 0],
       [0,     1, 0,     0],
       [-sin α,0, cos α, 0],
       [0,     0, 0,     1]]
Kbar=Tᵀ K T, Mbar=Tᵀ M T, fbar=Tᵀ f.       (27)
K=aggregate(Kbar), M=aggregate(Mbar), f=aggregate(fbar). (28)
```

По (27) source T означает q_local=T q_global. В проектной реализации
ниже используется явно объявленный **local-to-global** map; он двойствен
force map. Translation components поворачиваются, contraction и rotation
скаляры остаются неизменными. Значения fitted corrections Rucka не
переносятся; source Jang numeric κ по-прежнему не угадывается.

## 3. Production arm и physical coordinate contract

`project_jang_reduced_rectangular`, κ=5/6, I=I_y=bh³/12:

```text
y=(u,c,w,theta,N,R,Q,M)
C=EA/(1-nu²), H=kappa*GI, m=rho*A, j=r=rho*I,
B=EI, S=kappa*GA
T=1/2 ∫[m*u_t²+j*c_t²+m*w_t²+r*theta_t²]dx
V=1/2 ∫[C*(u_x²+2nu*u_x*c+c²)+H*c_x²
        +B*theta_x²+S*(w_x-theta)²]dx
N=C*(u_x+nu*c), R=H*c_x, Q=S*(w_x-theta), M=B*theta_x
```

В drawing basis EX вправо, EY вверх, EZ к наблюдателю, k=−EZ.
b вдоль local y, h вдоль local z=n; x вдоль t. Positive w первого плеча
направлено вниз, как в действующей rectangular Timoshenko display convention.
Из Jang `Ux=u−zθ`, `Uz=w+zc` следует positive physical rotation θ*k:
k×t=n, k×n=−t. **Знак θ и M в PDE/state не меняется.** Поворот вокруг
общей оси k не меняет скаляр θ при повороте local frame. Не переносится
иной physical interpretation Reddy rotation автоматически.

За основу направления взят existing explicit physical geometry contract
[RLB coordinate note](../laminated_beams/reddy_inplane_coordinate_contract.md),
с независимым выводом Jang rotation/efforts из displacement/energy. Старый
rectangular comparator хранит x₂∈[0,−L₂], а display helper имеет
t₂=(cosβ,sinβ), n₂=(sinβ,−cosβ). Эти helpers не изменены. Здесь обе
координаты **положительны от outer clamp к joint**. Это явный выбор
параметризации, не подбор знаков по корням.

| Arm, β=0 | x interval | Clamp | Joint | Tangent t | Transverse n | Global single coordinate X |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | [0,L₁] | x₁=0 | x₁=L₁ | (1,0) | (0,−1) | X=x₁ |
| 2 | [0,L₂] | x₂=0 | x₂=L₂ | (−1,0) | (0,1) | X=L−x₂ |

На втором плече conversion к old negative-coordinate storage: x₂new=−x₂old,
u_new=−u_old, w_new=−w_old, θ_new=θ_old; новые c/R добавляются из собственной
M-H variation. При conversion к глобальному single-coordinate state:

```text
(u,c,w,theta,N,R,Q,M)_global = (-u,c,-w,theta,N,-R,Q,-M)_arm2.
```

Это следует также непосредственно из resultants и d/dX=−d/dx₂; R и M
меняют знак, N и Q — нет. c не является Cartesian translation.

## 4. Boundary variation и восемь условий узла

Точное дифференцирование квадратичной энергии даёт boundary δV:

```text
[N*delta u + R*delta c + Q*delta w + M*delta theta]_0^Li.
```

Следовательно nodal endpoint sign σ=−1 в x=0, σ=+1 в x=L_i. Оба joint
ends здесь right ends: σ₁=σ₂=+1, оба outer clamps left ends: −1.
Не следует смешивать internal N,R,Q,M с outward nodal efforts.

Инвариантные условия без предварительного выбора sin/cos:

| Pair | Kinematics | Variational equilibrium | Scalar count |
| --- | --- | --- | --- |
| Translations | d₁=d₂, d_i=u_i t_i+w_i n_i | Σ σ_i(N_i t_i+Q_i n_i)=0 | 2+2 |
| Rotation | θ₁=θ₂ | Σ σ_i M_i=0 | 1+1 |
| Contraction | c₁=c₂ | Σ σ_i R_i=0 | 1+1 |

Итого **8 scalar conditions**, rank8. При независимых common variations
δd_J,δc_J,δθ_J: δu_i=t_i·δd_J, δw_i=n_i·δd_J, δc_i=δc_J,
δθ_i=δθ_J. Подстановка в сумму endpoint work даёт

```text
delta W_J = [Σ σ_i(N_i*t_i+Q_i*n_i)]·delta d_J
          + [Σ σ_i R_i]*delta c_J + [Σ σ_i M_i]*delta theta_J.
```

Независимость этих variations даёт ровно указанные balances. Force map
двойствен displacement map; для orthogonal local-to-global T оба maps
используют T (local variations — Tᵀδq_global). Fraction checks проверяют
integration-by-parts endpoint signs и rational work identity точно;
48 arbitrary-state float checks имеют max scaled error6.67e−16.
Lean/SymPy не доступны; используются прямой вывод и exact rational algebra.

В actual row order (dX,dY,c,θ,FX,FY,Rnode,Mnode), при β=0:

```text
u1+u2=0;  -w1-w2=0;  c1-c2=0;  theta1-theta2=0;
N1-N2=0;  -Q1+Q2=0;  R1+R2=0;  M1+M2=0.
```

В частности R₁+R₂=0 **для этой** coordinate/endpoint convention, но в
single-coordinate representation это обычная непрерывность R. Без
объявленных signs нельзя произвольно писать [R]=0 или R₁+R₂=0.

## 5. Matrix form, beta0 proof и completeness

Local state order q-first из single-rod solver сохранён. Пусть T_i переводит
(u,c,w,θ) в (dX,c,dY,θ); P меняет порядок на (dX,dY,c,θ). Тогда

```text
J = [[P*T1, 0,       -P*T2, 0       ],
     [0,    P*σ1*T1, 0,     P*σ2*T2]]  (8x16 on (q1,p1,q2,p2)).
```

Primary matrix16×16 состоит из двух outer clamps по4 строки и J на joint
states. Каждое плечо использует **неизменённый** `finite_state_basis`:
cos/sin и bounded anchored exponentials. Unknowns arm1(MH4,Tim4),
arm2(MH4,Tim4). Положительная row scaling не двигает zeros. Отдельная
permutation даёт exact MH8×8 ⊕ Tim8×8; все mixed entries точно нулевые.
В одних этих bases determinant factorization имеет regular scaling и
permutation sign; равенство raw determinants с direct4×4 не заявляется.

**Прозрачность как exact form/domain result.** После указанного reflection
map right fields и energy densities совпадают с global single-rod fields.
Joint compatibility склеивает H¹ поля u,c,w,θ; два интеграла энергий
складываются в integral[0,L]. Никаких joint mass/stiffness/point energies нет.
Outer essential conditions совпадают. Force, R и moment balances точно
отменяют internal boundary work и являются natural transmission conditions
того же quadratic form. Restriction любой global solution даёт two-arm
solution и обратно. Это изометрия mass form и совпадение self-adjoint
energy domain; следовательно artificial split не меняет спектр/формы.
Это доказательство для **прямого однородного идеального 1D** interface,
не 3D и не validation при β≠0.

Поэтому прежние min-max upper counts direct problem применимы к assembly:
MH Young lower form η=.18 — upper count7, Tim relaxed rotation/SS lower
spectrum — upper count11. Отдельные seed-free bounded scans каждого
split насыщают оба counts; below-start upper count0. Это исключает
пропуски в заявленных диапазонах, а не только показывает 18 найденных roots.
Всего6 searches, по400 intervals, MH443/Tim467 determinant evaluations
на split, один allowed doubled scan только при count failure; retries и
failed intervals отсутствуют. Guard13 combined покрыт обеими family inventories.

## 6. Reference, numerical policy и результаты

Case тот же G20: E=ρ=1,ν=.3,b=.2,h=.05,total L=L_ref=1,κ=5/6.
Материал/геометрия не подбирались. Splits .50,.35,.65 — три algebraic
controls **одного** rod, не parameter map. Physical time units contract
формальны; ниже f*=f L/√(E/ρ), здесь численно равно f, не measured Hz.

Immutable reference `results/mindlin_herrmann_timoshenko_single_rod/342ce44bff81c36f/`
переиспользован после identity/input/code/dependency/artifact SHA checks;
0 reference root evaluations. Direct shapes восстановлены тем же unchanged
single-rod mode function при saved roots. Primary two-arm roots вычислены
независимо по matrix8×8, без seeded brackets от reference.

Третий уровень: самостоятельно assembled harmonic state, positive global
segments, short-step expm/positive-diagonal QR. Проверяется одинаковое
продвижение essential initial subspace сквозь два segments и весь L;
unsafe product H(L₂)H(L₁) full ill-conditioned matrices не выполняется.
Это ещё одна численная проверка той же 1D theory, не independent experiment.

До расчёта: inherited root xtol1e−11/rtol1e−12, 400 intervals и те же
bounded ranges; frequency relative2e−8, MAC loss1e−8, component L2 relative1e−6,
scaled clamp/interface1e−9, interface vs reference1e−6, short-step transfer1e−9,
nonzero singular condition≤1e8; Gauss200 на **каждом** subinterval.
Near-cluster gap threshold1e−6 потребовал бы subspace audit; clusters нет.
Tolerances после результатов не менялись.

| Position | Family/index | Direct f* | Two-arm 50/50 f* | Max absolute difference over3 splits |
| ---: | --- | ---: | ---: | ---: |
| 1 | bending1 | .050527602683514 | .050527602683507 | 7.64e−15 |
| 2 | bending2 | .136325767246192 | .136325767246196 | 6.22e−14 |
| 3 | bending3 | .260072546037018 | .260072546037019 | 8.11e−15 |
| 4 | bending4 | .416258846084940 | .416258846084941 | 1.12e−15 |
| 5 | axial1 | .500706143979489 | .500706143979483 | 6.00e−15 |
| 6 | bending5 | .599831008472830 | .599831008472830 | 2.23e−16 |
| 7 | bending6 | .806045921931286 | .806045921931286 | 1.12e−16 |
| 8 | axial2 | 1.001229212937640 | 1.001229212937643 | 4.22e−15 |
| 9 | bending7 | 1.030751272103510 | 1.030751272103510 | 0 |
| 10 | bending8 | 1.270446989448601 | 1.270446989448602 | 4.45e−16 |
| 11 | axial3 | 1.501381967307509 | 1.501381967307509 | 0 |
| 12 | bending9 | 1.522253494567150 | 1.522253494567607 | 4.57e−13 |
| 13 guard | bending10 | 1.783833631821124 | 1.783833631821124 | 0 |

Family inventories включают7 axial acoustic и11 bending modes; combined
sorted position не является family index/descendant branch. Optical wave
cutoffs contraction11.558994422/shear6.242570465 выше этих inventories;
их не называют finite eigenfrequencies. c(x) внутри acoustic MH modes
содержит evanescent contribution и проверяется, даже без optical roots.

## 7. Shapes, interface residuals и arm swap

Mass inner product ∫[m*u_a*u_b+j*c_a*c_b+m*w_a*w_b+r*θ_a*θ_b]dX;
coefficients mass-normalized, single overall sign alignment по overlap.
Все54 family modes проверены. MAC min≥.9999999999999995; max component
relative L2 error1.62e−12, **c error4.81e−13**. Inactive fields точно нулевые
по доказанному block decomposition; их относительная ошибка не делится на0.
Clamps max scaled residual2.15e−12; interface states vs direct2.16e−12.
Interface residual каждого компонента делится на max амплитуду соответствующего
global state component на quadrature, не на его (возможно нулевое) interface
значение. Raw residuals тоже сохранены; в inactive block они точно0.

| Joint row | Max scaled residual over54 modes |
| --- | ---: |
| dX | 1.86e−14 |
| dY | 1.27e−12 |
| c | 1.02e−15 |
| θ | 1.78e−12 |
| FX | 4.89e−14 |
| FY | 1.41e−12 |
| Rnode | 1.04e−15 |
| Mnode | 1.04e−12 |

Nonzero singular condition≤3.507; primary/reference frequency relative
error≤4.56e−13; segmented QR/direct projection difference≤4.60e−14.
Arm swap .35↔.65 с X→L−X и state reflection выше: max frequency difference
2.23e−16 relative, component L2≤2.06e−15, MAC loss≤machine precision.
Это отдельная проверка ориентации, не sign fitting.

## 8. Reproduction, статусы и ограничения

Один новый diagnostic CLI нужен для **joint/interface contract**; source и
direct single-rod CLIs не получают зависимости от joint module. Все arm/source
libraries, source fixtures, PDFs/BibTeX/index и прежние tests побайтово сохранены.

```powershell
python scripts/analysis/verify_mindlin_herrmann_timoshenko_beta0_joint.py --check-sources
python scripts/analysis/verify_mindlin_herrmann_timoshenko_beta0_joint.py --compute
python -m pytest tests/test_mindlin_herrmann_timoshenko_beta0_joint.py tests/test_mindlin_herrmann_timoshenko_finite_rod.py tests/test_mindlin_herrmann_timoshenko_literature.py tests/test_bishop_literature.py tests/test_timoshenko_bishop_single_rod.py -q
```

Если direct reference отсутствует, сначала выполнить unchanged single-rod
`--compute`; joint `--compute` также может создать reference тем же CLI.
`--check-sources` сам roots не считает. Фактический interpreter:
D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe; Python3.12.4,
NumPy2.1.3/SciPy1.15.2, новых зависимостей нет.

Bundle `results/mindlin_herrmann_timoshenko_beta0_joint/<fingerprint>/`:
manifest с working-state/Git/source/code/artifact hashes и reference provenance,
parameters, result с family inventories/brackets/counts/matrix diagnostics,
frequencies.csv (12+guard по3 splits), mode_profiles.csv (u,c,w,θ,N,R,Q,M
и direct fields). Matching-cache reuse проверяет все hashes, roots не считает.
Изменённые входы меняют identity; повреждённые artifacts не используются.
Новый plot не нужен: numerical profile comparison и saved data достаточны.

| Gate | Status |
| --- | --- |
| MHTIM_JOINT_SOURCE_CLOSURE | PASS, с reduced common-DOF qualification |
| MHTIM_JOINT_VARIATIONAL_FORM | PASS |
| MHTIM_BETA0_HOMOGENEOUS_SPECTRUM | PASS |
| MHTIM_BETA0_SPLIT_INVARIANCE | PASS |
| MHTIM_BETA0_MODE_SHAPES | PASS |
| MHTIM_BETA0_ARM_SWAP | PASS |
| MHTIM_BETA0_JOINT_GATE | PASS |

31 new joint tests +98 unchanged single/source/Bishop tests=129 PASS;
ещё4 relevant existing rectangular Timo regressions PASS: **133 total**.
Import/smoke/hash-validated reuse и diff whitespace checks проходят.
Previous HIERARCHY_SINGLE_ROD=PARTIAL_PASS сохраняет свою clamp qualification.
Статус β0 не переносится на угловую геометрию: вопросы finite 3D joint region,
неоднородного interface и nonzero-angle validation остаются за рамками.
β≠0 spectrum, sweeps, coefficient fit, nonlinear derivation, FEM truth,
article/Yartsev/viscosity и thematic memory другой ветви не затрагивались.
После этого gate работа остановлена.

## 9. General-angle geometry: действующий project contract

Расширение запрошено отдельно после beta0 PASS. Физика плеча, source
variants, single-rod solver и frozen bundle3059d70b1b50ea2e не изменяются.
General geometry берётся из `reddy_inplane_geometry.py` и
[coordinate contract §3](../laminated_beams/reddy_inplane_coordinate_contract.md#3-локальные-базисы-плеч).
Импортируются только t,n; constitutive/state/rotation physics Reddy не
переносятся. Frozen equations.tex/src analytic не содержат этого M-H block;
новая реализация согласуется с source energy и local state §§3–4 и
выделена в прежнем joint-helper, без ревизии baseline formulas.

Joint origin=(0,0); EX вправо, EY вверх, k=−EZ. β — signed отклонение луча
joint→right outer clamp от +EX против часовой стрелки. Left outer clamp
(−L₁,0), right outer clamp=(L₂cosβ,L₂sinβ). Обе положительные local x идут
**от заделки к узлу**. Поэтому, после восстановления project geometry:

```text
t1=(1,0), n1=(0,-1),
g2=(cos(beta),sin(beta)),
t2=-g2=(-cos(beta),-sin(beta)),
n2=k cross t2=(-sin(beta),cos(beta)).
```

В частности β=0 означает прямой rod, β=90° — прямой угол; opening angle
между outer rays при 0≤β≤90° равен180°−β. |t|=|n|=1, t·n=0,
t×n=k следуют из cos²β+sin²β=1. Относительные rotations proper;
translation matrix в drawing EX/EY имеет det−1 из-за local transverse
downward convention. Полный beam triad (t,EZ,n) right handed.
При mirror geometry и том же view axis k знаки rotation/M рассмотрены в §13.

## 10. Local/global transformation и dual forces

Сначала сохраняются invariant conditions §4. Пусть
A_i=[t_i,n_i], q_i=(u_i,c_i,w_i,θ_i), g_J=(dX,cJ,dY,θJ),
G_i — orthogonal matrix с A_i на translation indices(0,2) и единицами на
scalar indices(1,3). Тогда

```text
g_i=G_i*q_i, q_i=G_iᵀ*g_J;
outward f_i=σ_i*(N_i,R_i,Q_i,M_i);
f_global,i=G_i*f_i,
f_iᵀ*δq_i=(G_i*f_i)ᵀ*δg_J.
```

Force map **получен из virtual work**, а не задан отдельно. Оба joint
ends остаются local right ends, σ₁=σ₂=+1; external ends σ=−1.
c scalar reduced DOF, не Cartesian translation. θ invariant при proper
in-plane rotation вокруг общего k. General matrix — ровно J §5 с
G_i(beta); восемь conditions и их order сохранены. c/R и θ/M не получают
trig factors. Нормальная сила и transverse shear смешиваются геометрией
в global translation/force rows; это structural coupling, при неизменённых
local MH⊕Tim operators. При β0 exact permutation вновь даёт MH8⊕Tim8.

По orthogonality G_i: J Jᵀ=2I₈ при обоих joint signs++; rank8 и
singular values√2 для **любого** proper frame. Проверены β0/5/45/90 и
128 deterministic arbitrary-state work identities: max scaled error6.67e−16.
Exact rational duality checks тоже проходят; Lean/SymPy недоступны.

Одна production implementation `frame_boundary_matrix(...,beta_deg=...)`
собирает16×16 из прежнего bounded arm basis. `boundary_matrix` оставлена
как frozen beta0 guard/wrapper и **делегирует этой же функции**. Старый
CLI/test, запрещающий nonzero-angle call в старом API, сохранён; physics
fork отсутствует. General CLI использует общий API непосредственно.

## 11. Exact beta0 regression и beta→0 continuity

Frozen reference — immutable `results/mindlin_herrmann_timoshenko_beta0_joint/3059d70b1b50ea2e/`.
Проверены artifact hashes и unchanged arm/source/orchestration identities;
изменённый joint-helper — именно audit subject. Frozen CLI переиспользован
in memory, без записи в его старые bundles. Все3 splits .5/.35/.65:
54 family frequencies, full shapes/c/R residuals, split and arm-swap checks
проходят. Matrix max absolute difference=0; frequency relative=0 и
component L2 difference=0 в данном floating run. Общая assembly при β0
совпала с frozen physical matrices; сравнение raw determinants между
разными нормировками не используется как требование.

Отдельный zero-limit control: same equal arms и первый low spectrum до
f*=.7. Фиксированные малые углы — не scientific map/sensitivity study.

| β, deg | max relative difference first3 roots vs β0 | ||Jβ−J0||₂ |
| ---: | ---: | ---: |
| 1e−6 | 6.67e−15 | 1.745329252e−8 |
| 1e−4 | 3.76e−11 | 1.745329252e−6 |
| 1e−2 | 3.76e−7 | 1.745329250e−4 |

До расчёта allowances1e−9/1e−7/1e−3 соответственно; все соблюдены.
Это diagnostic convergence allowances, не пределы applicability и не
оценки sensitivity. Аналитически ||Jβ−J0||₂=2|sin(β/2)|≤|β| в radians.
Boundary matrix состоит из Jβ и неизменённых analytic arm columns и
поэтому непрерывна. Low β0 roots простые/разделённые; local implicit-root
continuity согласуется с расчётом. MAC между разными β не вычислялся.

## 12. Beta90 geometry diagnostic

t₂→(0,−1), n₂→(−1,0). Entries transforms становятся0/±1 с ошибкой
6.13e−17; c/θ diagonal entries остаются1. Из invariant maps, а не из
спектрального fit:

```text
d1=(u1,-w1), d2=(-w2,-u2);
Fnode1=(N1,-Q1), Fnode2=(-Q2,-N2).
```

Это понятные axis permutations; θ/c compatibility и M/R balances те же.
Отдельный right-angle virtual-work check входит в128 trials. Спектр
не использовался для назначения знаков.

## 13. Arm swap и reflection

Label swap при β45: переставляются **frames вместе с arm labels**, а не
только lengths. С column block swap joint operator меняет знак только
первых4 compatibility rows; balances остаются теми же. Matrix-map error0.
При возвращении в canonical drawing basis proper global rotation
−(180°+β) переводит swapped geometry к canonical −β; scalar c/θ не
меняются. Никакой signed fitting нет.

Зеркало EX: S=diag(1,−1). Чтобы сохранить t*×n*=k,
t*=S t, n*=−S n. Local state transformation:

```text
(u,c,w,theta,N,R,Q,M)*=(u,c,-w,-theta,N,R,-Q,-M).
```

Таким образом c/R scalars, θ/M signed planar pseudoscalars при improper
reflection. Global translation/force отражаются S; local shear/rotation
оба меняют знак, энергии/инерция сохраняются. Reflection operator identity
ошибка0; reflected frames приβ45 совпадают с canonical frames−45.

Для обязательного symmetry eigenpair gate **до full pilot** вычислены
6 roots до f*=.7 при canonical45, label-swap45 и mirror−45. Это одна
физическая geometry после isometric remapping, не extra angle study.
Max swap frequency relative difference2.23e−16, component L2≤3.72e−14;
mirror differences0; MAC min≥.9999999999999995. Force/M/R residual gates
проходят во всех remapped cases. Roots простые; близкий cluster потребовал
бы subspace comparison и остановил бы текущий индивидуальный-vector gate.
Эти same-geometry overlaps не являются across-beta tracking.

## 14. Exact energy count и bounded root inventory

Подсчёт direct-rod min-max не переносится на angular assembly. Вместо
этого выводится небольшой exact energy/Schur count для **этого** fixed
equal-arm below-cutoff problem, без universal spectrum framework/FEM.

Пусть aλ=V−λT, λ=ω² (quadratic forms без1/2). V₀ — поле с q=0 на
обоих концах каждого arm, включая joint. Its negative index J₀(ω) равен
сумме fixed–fixed arm eigenvalue counts. Off those Dirichlet poles существует
unique harmonic extension Eλ g общей four-coordinate joint DOF.
Integration by parts даёт aλ(Eλ g,z)=0 для z∈V₀. Поэтому

```text
aλ(z+Eλg,z+Eλg)=aλ(z,z)+gᵀ Sλ g,
Sλ=Σ G_i K_arm,i(ω) G_iᵀ,
N_full(<ω)=J0(<ω)+negative_inertia(Sλ).
```

K_arm — exact Dirichlet-to-Neumann map: решить essential endpoint
system q(0)=0,q(L_i)=given в прежнем analytic basis, получить p(L_i).
Реакция на joint имеет+ sign. Это independent4×4 energy matrix,
не determinant16×16 под другим именем. Symmetry следует из reciprocal
energy; raw skew residual проверяется **до** roundoff symmetric part.
Fixed positive congruence balances S и сохраняет inertia.

Для одного L_i=.5 catalog4 MH и9 Tim poles независимо проверен через
unchanged state-expm/QR. Saturated min-max counts certifies catalog до
4.5/4.568391355 f*, выше общего pilot ceiling3.75. J₀ удваивается для
одинаковых arms. Count samples исключают known pole neighborhoods;
в predefined symmetric exclusion window1e−7*max(1,ω) обе side counts
должны совпасть, иначе unresolved root/pole case и stop. Вне window
inertia margin≥4e−14, boundary condition≤1e10, raw Schur skew≤1e−10.
Ни одного unresolved pole/multiplicity case в этом run нет.

Primary roots — full16×16 determinant со positive row normalization.
Fixed400 intervals от .001πc₀/L до7.5πc₀/L; exact counts at boundaries
сертифицируют число roots в каждом interval. Count jump>1 вызывает единственный
предусмотренный bounded subdivision mechanism (depth≤20,budget2000),
затем brentq с unchanged xtol1e−11/rtol1e−12. Ни расширения диапазона,
ни tolerance relaxation, ни duplicate-root deletion по близости нет.
Unresolved multiplicity/conditioning означала бы stop/PARTIAL или FAIL.

## 15. Fixed spectral pilot, normalization и residuals

Только после PASS всех8 structural gates выполненыβ5/45/90. Same G20
E=ρ=1,ν=.3,b=.2,h=.05,total L=1, L₁=L₂=.5,κ=5/6, j=r=ρI.
External u=c=w=θ=0 неизменны. Ни μ, ни thickness mismatch, ни hierarchy
не вводятся. Ни одна position при одномβ не объявляется той же модой
при другомβ; таблица содержит independently sorted positions.

| Position | β5 f* | β45 f* | β90 f* |
| ---: | ---: | ---: | ---: |
| 1 | .055033967071 | .136051160059 | .134723376618 |
| 2 | .136322717201 | .155336542940 | .185632402245 |
| 3 | .260614727245 | .307572760372 | .383164055806 |
| 4 | .416184307973 | .409759368849 | .396591527265 |
| 5 | .500711404773 | .501134260130 | .502413042506 |
| 6 | .599648770566 | .586829718439 | .559866141189 |
| 7 | .806071027692 | .808241880224 | .817337287827 |
| 8 | .996894538416 | .935209839680 | .905590601426 |
| 9 | 1.035936564806 | 1.153232145595 | 1.241571900481 |
| 10 | 1.270389825139 | 1.265318983025 | 1.270814792498 |
| 11 | 1.501284203395 | 1.493128521996 | 1.464849754533 |
| 12 | 1.522217658810 | 1.519478633942 | 1.512533708799 |
| 13 guard | 1.783887293235 | 1.788406840272 | 1.805278176664 |

f*=f L/√(E/ρ), здесь численно f; это normalized contract, не эксперимент
в Hz. Оптические cutoffs не называются eigenfrequencies конечной системы.
В полном fixed search range найдено/сертифицировано23/24/24 positive roots;
ниже lower boundary count0. Это покрывает12+guard13. Subdivisions0/0/2,
failed intervals0; при90 два подразделения разрешили близкие simple roots.
Shapes/resultants сохранены для всех71 computed roots, не только prefix.

General mass normalization: Σ∫[m u²+j c²+m w²+rθ²]dx=1, Gauss200 per arm.
Mode checks из прежних PDE/basis дают ODE, energy ω² и mass Gram diagnostics.
Сдвиг angular frames не меняет local mass coefficients илиκ. Component
norms не превращаются в новую family-classification theory.

До расчёта зафиксированы scaled clamp/joint/ODE1e−9, energy5e−8,
Gram5e−7, nonzero singular condition≤1e8. Для mixed angular state использована
одна dimensionally consistent work metric: q=(u,Lc,w,Lθ),
p=(N,R/L,Q,M/L), conjugate generalized pairs. Joint displacement rows
(dX,dY,Lδc,Lδθ) делятся на max over-arm q amplitude; force rows
(FX,FY,δR/L,δM/L) — на max over-arm p amplitude. Векторы rotations/forces
таким образом сравниваются в общей физической системе. Raw residuals
и дополнительный own-amplitude c/R diagnostic тоже сохраняются.

| Max over full computed inventory | β5 | β45 | β90 |
| --- | ---: | ---: | ---: |
| Clamp scaled residual | 1.13e−12 | 8.03e−14 | 2.11e−13 |
| Joint scaled residual | 1.37e−11 | 4.37e−12 | 2.41e−12 |
| c compatibility (work metric) | 4.14e−15 | 1.45e−15 | 3.85e−14 |
| R balance (work metric) | 4.30e−16 | 6.39e−17 | 4.74e−16 |
| M balance (work metric) | 5.85e−14 | 1.39e−15 | 7.53e−15 |
| c compatibility / own c amplitude | 1.68e−12 | 1.59e−13 | 4.18e−13 |
| R balance / own R amplitude | 2.68e−12 | 4.84e−14 | 1.94e−13 |
| ODE relative residual | 1.50e−14 | 4.13e−15 | 3.81e−15 |
| Energy relative error | 6.55e−13 | 3.79e−13 | 7.00e−13 |
| Mass Gram error | 3.64e−11 | 1.88e−12 | 4.27e−12 |
| Nonzero singular condition | 47.14 | 92.68 | 461.22 |

Собственные формы сравнивались по MAC только для идентичной физической
geometry (frozen regression/swap/reflection), не междуpilotβ. Никаких
applicability/veering/safe-prefix/10% выводов из таблицы не делается.

## 16. General-stage reproduction, statuses и remaining limitations

```powershell
python scripts/analysis/verify_mindlin_herrmann_timoshenko_general_beta_joint.py --check-sources
python scripts/analysis/verify_mindlin_herrmann_timoshenko_general_beta_joint.py --compute
python -m pytest tests/test_mindlin_herrmann_timoshenko_general_beta_joint.py tests/test_mindlin_herrmann_timoshenko_beta0_joint.py tests/test_mindlin_herrmann_timoshenko_finite_rod.py tests/test_mindlin_herrmann_timoshenko_literature.py tests/test_bishop_literature.py tests/test_timoshenko_bishop_single_rod.py -q
```

Новый diagnostic entry point имеет другой contract: structural hard gates,
coupled exact count, symmetry and fixed3-angle pilot. Он переиспользует
existing joint-helper/arm solver, frozen beta0 orchestration и atomic writers;
расширять frozen beta0 CLI было бы небезопасно. Нет command для arbitrary
angle map, hierarchy или article figures. Existing Python3.12.4/NumPy2.1.3/
SciPy1.15.2 environment, interpreter D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe;
новых dependencies нет.

Bundle `results/mindlin_herrmann_timoshenko_general_beta_joint/<fingerprint>/`:
manifest/Git/source/code/frozen hashes, parameters, full gates/transforms/ranks,
zero-limit/symmetry/poles/count samples/brackets/SVD/mode diagnostics,
frequencies.csv и local full8-state mode_profiles.csv. Matching cache checks
hashes, roots не считает. Старые beta0/direct bundles не записывались;
source index/BibTeX/fixtures/source variants/single physics/old tests сохранены.
Generated conclusions продублированы в этом tracked документе; новых plots нет.

| MHTIM_GENERAL_BETA category | Status |
| --- | --- |
| COORDINATES | PASS |
| VIRTUAL_WORK_DUALITY | PASS |
| BETA0_REGRESSION | PASS |
| ZERO_LIMIT | PASS |
| JOINT_RANK | PASS |
| RIGHT_ANGLE_GEOMETRY | PASS |
| ARM_SWAP | PASS |
| REFLECTION | PASS |
| SPECTRAL_PILOT | PASS |
| JOINT_GATE | PASS |

37 new general tests +129 unchanged MH single/source/beta0/Bishop tests =166
PASS; ещё4 relevant rectangular Timo regressions PASS, **170 unique tests**.
Source/import/smoke/cache и diff checks проходят. Old beta07 gates и
historical source statuses сохранены. Теперь подтверждены general-angle
transformations и ограниченный finite pilot этой reduced model; это не
validation реальной сварной finite3D joint области или всех parameter cases.
**c₁=c₂ остаётся variational reduced-frame closure, не прямым следствием 3D
elasticity.** κ/coefficients не подбирались, FEM не использовался, nonlinear
equations/elementary–RL hierarchy/applicability maps/tracking не выводились.
На этом этапе работа останавливается; следующего parameter study нет.
