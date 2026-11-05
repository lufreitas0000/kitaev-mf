# Discovery and Pivot Log: kitaev-mf

| Date | Module/Context | Obstruction | Pivot/Resolution | Status |
| :--- | :--- | :--- | :--- | :---: |
| 2026-11-04 | `services/solver.py` (M1.3) | High-precision BZ integral for pure-Kitaev correlator $K\text{-}006$ required large meshes and showed boundary quadrature sensitivity | Introduced Romberg/tanh-sinh adaptive quadrature scheme in `scripts/k006_quadrature.py` with validated asymptotic error scaling | `Resolved` |
| 2026-11-05 | `services/solver.py` (M1.6) | Plain Picard iteration exhibits oscillatory divergence near intermediate-field phase boundaries | Pivoted solver hierarchy to joint DIIS and Grassmannian Newton-Krylov with linesearch damping | `Active` |
