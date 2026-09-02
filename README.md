# Low_Rank_IgA_Projection

Low-rank assembly of the system tensors arising in isogeometric analysis (IgA),
based on univariate projections and the tensor-train (TT) format.

This is the reference implementation accompanying the manuscript
*Projection-based low-rank assembly in IgA* by Tom-Christian Riemer and
Martin Stoll.

For a three-dimensional geometry map $G$, the mass and
stiffness operators are

```math
\mathbf{M}_{(\hat{\mathbf{i}},\hat{\mathbf{j}})} = \int_{\left[0,1 \right]^3} \hat{\beta}_{\hat{\mathbf{i}}} \left(  \hat{x}  \right)  \, \hat{\beta}_{\hat{\mathbf{j}}} \left(  \hat{x}  \right)  \, \omega  \left(  \hat{x}  \right)  \, \mathrm{d} \hat{x}, \qquad \mathbf{K}_{(\hat{\mathbf{i}},\hat{\mathbf{j}})} = \int_{\left[0,1 \right]^3}  \left(  Q  \left(  \hat{x}  \right)  \nabla \hat{\beta}_{\hat{\mathbf{i}}} \left(  \hat{x}  \right)  \right)  \cdot \nabla \hat{\beta}_{\hat{\mathbf{j}}} \left(  \hat{x}  \right)  \, \mathrm{d} \hat{x},
```
where $\omega \left( \hat{x} \right)  = \det \left( \nabla G \left( \hat{x} \right)  \right)$ and $Q \left( \hat{x} \right)  = \omega \left( \hat{x} \right)   \left( \nabla G \left( \hat{x} \right)  \right)^{-1}  \left( \nabla G \left( \hat{x} \right)  \right)^{-\top}$

The geometry-dependent weight functions $\omega$ and $Q$ are the only obstacle to a tensor-product
structure. This repository provides a method that recovers that structure by
projecting those factors mode by mode onto an auxiliary spline space, so that
the operator is assembled from univariate contractions and stored in TT format,
never as a full array.


## Method

For the **mass tensor**, the weight function $\omega = \det(\nabla G)$ is itself a
piecewise polynomial and therefore lies in a reduced spline product space. Its
coefficient tensor is obtained from the coordinate control-point tensors
$\mathbf{C}^{(a)}$, $a = 1,2,3$, by univariate coefficient transfer operators.
With exact quadrature and without TT truncation, the resulting representation is
an exact reformulation of the standard Galerkin mass tensor.

For the **stiffness tensor**, the matrix-valued weight function $Q$ is rational.
Its entries are split into polynomial numerators $N_{kl}$, which are treated in the
same way as $\omega$, and the reciprocal determinant $\rho = 1/\omega$, which is in
general not a spline function and is approximated by an $L^2$-projection onto a
tensor-product spline space $\mathbb{S}_{\rho}$. The projection system is solved by
the alternating minimal energy (AMEn) method.

Everything is carried out in the TT format. The coefficient tensors, of order nine
and twelve before grouping, are never formed explicitly, and the multivariate
integrals reduce to univariate integrals and contracted products.

The implementation also contains an adapted version of the interpolation-based
low-rank assembly of Bünger, Dolgov and Stoll, used as a comparison, and calls
GeoPDEs for full-format reference matrices.


## Scope and limitations

- The construction requires an **orientation-preserving tensor-product B-spline**
  geometry map. NURBS parameterizations are **not** covered. A determinant of
  constant negative sign is admissible up to a sign change in a single TT core.
- Efficiency depends on the coordinate control-point tensors $\mathbf{C}^{(a)}$
  having moderate TT ranks, and on a moderate number of univariate basis functions
  per direction. The univariate operators scale linearly in $n^{(d)}$ for fixed
  degree but do not disappear.
- The accuracy of the stiffness assembly is governed by the richness of the
  projection space $\mathbb{S}_{\rho}$. Refining it enlarges the univariate operators
  and the projection system, which dominates the cost of the stiffness assembly.


## Requirements

| Dependency | License | Purpose |
|---|---|---|
| MATLAB (R2020a or newer) | proprietary | `exportgraphics` and `string` are used in the plotting and driver code |
| [GeoPDEs](https://github.com/rafavzqz/geopdes) | GPL-3.0 | IgA framework: spline spaces, meshes, reference assembly (`op_u_v_tp`, `op_gradu_gradv_tp`) |
| NURBS toolbox | GPL-2.0-or-later | geometry handling (`geo_load`, `nrbmak`), shipped together with the GeoPDEs releases |
| [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) | MIT-style | `tt_tensor` / `tt_matrix` classes, TT rounding, cross approximation, AMEn solver |
| Linux | — | peak-memory measurement only; reads `/proc/self/status`, needs kernel 4.0 or newer built with `CONFIG_PROC_PAGE_MONITOR` |
| [matlab2tikz](https://github.com/matlab2tikz/matlab2tikz) | BSD-2-Clause | optional, used to export the manuscript figures |

Add GeoPDEs, the NURBS toolbox, the TT-Toolbox and this repository to the
MATLAB path before running anything.




## Reproducing the experiments in the paper


The results reported in the manuscript were produced on an Ubuntu Linux server
(kernel 5.15.0-41-generic, x86_64) with two AMD EPYC 9534 processors
(2 x 64 cores, 256 logical CPUs) at approximately 3.72 GHz and 1.5 TiB RAM.

Peak-memory figures are measured **per process**, by spawning one fresh
`matlab -batch` child process per configuration and reading its high-water mark
from `/proc`. In-process measurement is not comparable across methods, because
MATLAB pools memory and short-lived transients are missed. This part of the
harness is Linux-only.


## Citation

If you use this code, please cite the accompanying paper.

```bibtex
@misc{riemer2026projection,
  title         = {Projection-based low-rank assembly in {IgA}},
  author        = {Riemer, Tom-Christian and Stoll, Martin},
  year          = {2026},
  eprint        = {2609.01218},
  archivePrefix = {arXiv},
  primaryClass  = {math.NA},
  doi           = {10.48550/arXiv.2609.01218}
}
```


Please also cite [GeoPDEs](https://github.com/rafavzqz/geopdes) and the
[TT-Toolbox](https://github.com/oseledets/TT-Toolbox), on which this code depends.


## License

This project is licensed under the GNU General Public License v3.0. See
[`LICENSE`](LICENSE) for the full text.

GPL-3.0 was chosen because this code depends on GeoPDEs, which is itself
distributed under GPL-3.0.

