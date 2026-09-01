# Low_Rank_IgA_Projection

Low-rank assembly of the system tensors arising in isogeometric analysis (IgA),
based on univariate projections and the tensor-train (TT) format.

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


## Requirements

| Dependency | Purpose |
|---|---|
| MATLAB (R2020a or newer) | `exportgraphics` and `string` are used in the plotting and driver code |
| [GeoPDEs](https://github.com/rafavzqz/geopdes) | IgA framework: spline spaces, meshes, reference assembly (`op_u_v_tp`, `op_gradu_gradv_tp`) |
| NURBS toolbox | geometry handling (`geo_load`, `nrbmak`), shipped together with the GeoPDEs releases |
| [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) | `tt_tensor` / `tt_matrix` classes, TT rounding, cross approximation, AMEn solver |
| Linux | peak-memory measurement only; reads `/proc/self/status`, needs kernel 4.0 or newer built with `CONFIG_PROC_PAGE_MONITOR` |
| [matlab2tikz](https://github.com/matlab2tikz/matlab2tikz) | optional, used to export the manuscript figures |

Add GeoPDEs, the NURBS toolbox, the TT-Toolbox and this repository to the
MATLAB path before running anything.


