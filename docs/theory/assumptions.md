# Assumptions

For project-wide rules on source-of-truth priority, branch identity,
diagnostic workflow, and model-extension checks, see `../project_rules.md`.

Здесь фиксируются активные допущения и рабочие гипотезы проекта.

## Working Notes

## Out-of-plane EB + Saint-Venant torsion

- The out-of-plane subsystem is treated as a separate linear model for a
  planar two-beam system: Euler--Bernoulli out-of-plane bending plus
  Saint-Venant torsion. It is independent of the existing in-plane
  bending/axial subsystem in the ideal linear planar setting. Its coordinate
  convention is fixed in `docs/theory/out_of_plane_eb_torsion.md`: `e_z`
  points downward, `t_i` runs from clamp to joint, `n_i=e_z x t_i`, and
  tangent vectors must be written as `t_i` rather than `tau_i`.

## Diagnostic thickness mismatch

- The mass-preserving thickness-mismatch model with
  `eta=(r_2-r_1)/(r_1+r_2)` is diagnostic-only. It keeps the total length
  `2l`, the equal-radius base radius `r_0` at `eta=0`, and the total mass fixed
  through `tau_1=(1-eta)/sqrt(1+2 mu eta+eta^2)` and
  `tau_2=(1+eta)/sqrt(1+2 mu eta+eta^2)`. It does not replace the baseline
  equal-radius determinant, article model, or FEM model. See
  `docs/thickness_mismatch/README.md`.

- Базовые обозначения для theory-facing формул: `l_1`, `l_2` — длины плеч, `l=(l_1+l_2)/2` — базовая длина, `r` — радиус круглого сечения, `\Lambda` — безразмерная частота, `\beta` — угол сопряжения, `\varepsilon` — параметр толщины, как в `docs/literature/pdf/Статья-Дорофеев-2025.pdf`.
- Локальный параметр `\mu` используется как параметр несимметрии длины в текущем расширении задачи: `\mu=(l_2-l_1)/(l_1+l_2)`, `l_1=l(1-\mu)`, `l_2=l(1+\mu)`. Эквивалентно, в формулах с базовой длиной `L`: `L_1/L = 1-\mu`, `L_2/L = 1+\mu`. В опубликованной статье Дорофеева 2025 этого параметра нет, потому что там рассматриваются одинаковые стержни длины `L`.
- Случай `\mu=0` фиксирует симметричный эталон для рабочей нумерации ветвей. Нумерация относится к отслеживаемой ветви, а не к текущему месту этой ветви в отсортированном спектре.
- Пунктирные reference-линии на графиках `\Lambda(\mu)` соответствуют одиночным стержням с закреплением заделка--шарнир: `\Lambda_n^{(1)}(\mu)=\alpha_n/(1-\mu)` и `\Lambda_n^{(2)}(\mu)=\alpha_n/(1+\mu)`, где `\alpha_n` — корни `\tan\alpha=\tanh\alpha`. Эти линии используются только как ориентиры для интерпретации, а не как замена анализу форм колебаний.
- При использовании `docs/literature/pdf/2003JSVb.pdf` знаки в матрице `T` формулы `(39)` нужно сверять вручную с уравнениями непрерывности `(10)--(11)`; печатную запись этого определителя не стоит переносить в теорию или код без дополнительной проверки.
