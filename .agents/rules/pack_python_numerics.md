<!-- @CANONICAL-RULE: DO NOT EDIT DIRECTLY. MANAGED BY project-governance -->
# Pack: Python Quantitative Metrology & Bayesian Model Selection

This pack governs quantitative numerical physics, mean-field solvers, Monte Carlo, and spectral simulations implemented in Python (`numpy`, `scipy`, `torch`).

---

## 1. Conditioning, Residuals & Solver Invariants

Every reported numerical quantity and converged state must quantify algebraic conditioning and solver residuals:
1. **Matrix Conditioning ($\kappa$):**
   - Report 2-norm condition numbers $\kappa(\mathbf{A}) = \|\mathbf{A}\|_2 \|\mathbf{A}^{-1}\|_2$. If $\kappa(\mathbf{A}) > 10^8$, use iterative refinement, regularized pseudoinverses, or high-precision solvers, documenting precision bounds.
2. **Stationarity & Constraints:**
   - In constrained optimization (e.g. self-consistent Hartree-Fock-Bogoliubov, Augmented Lagrangian Multipliers), report KKT stationarity residuals $\|\mathbf{r}_{\mathrm{KKT}}\|_\infty$ and gauge constraint violations $\epsilon_{\mathrm{gauge}} = \max_i |c_i(\Omega)|$.
3. **Idempotency & Spectral Projectors:**
   - For one-body density matrices / projector algebras, verify purification idempotency: $\|\mathbf{P}^2 - \mathbf{P}\|_F / \sqrt{N} < \epsilon_{\mathrm{tol}}$.

---

## 2. Multi-Basin Sampling & Robust Confidence Intervals

Non-linear mean-field and self-consistent equations admit multiple local minima/basins of attraction:
1. **Multi-Basin Exploration:**
   - Require multi-seed initialization sweeps ($N_{\mathrm{seeds}} \ge 20$) spanning diverse configuration families (e.g. uniform, topological vortex/defect seeded, random Gaussian).
2. **Wilson Score Confidence Intervals (Discrete / Basin Survival):**
   - For discrete fraction or basin survival probability $\hat{p} = k/n$, report **95% Wilson score confidence intervals**:
     $$\mathrm{CI}_{95\%} = \frac{\hat{p} + \frac{z^2}{2n} \pm z \sqrt{\frac{\hat{p}(1-\hat{p})}{n} + \frac{z^2}{4n^2}}}{1 + \frac{z^2}{n}}, \quad z \approx 1.96$$
3. **Non-Parametric Bootstrap Confidence Intervals (Continuous Observables):**
   - For continuous physical observables (order parameters, gaps, string tensions), report non-parametric bootstrap 95% CIs with $N_{\mathrm{samples}} \ge 50$ and resample count $B \ge 1,000$ using the Bias-Corrected and Accelerated ($\mathrm{BC}_a$) or percentile method.

---

## 3. Bayesian Model Selection ($\Delta\mathrm{BIC}$ & $\Delta\mathrm{AICc}$)

When distinguishing between competing functional scaling laws (e.g. exponential localized decay vs. algebraic power-law, or linear string tension vs. logarithmic screening), perform formal Bayesian model comparison:
1. **Bayesian Information Criterion (BIC):**
   $$\mathrm{BIC}(\mathcal{M}_k) = n_{\mathrm{eff}} \ln\!\left(\frac{\mathrm{RSS}_k}{n_{\mathrm{eff}}}\right) + k \ln n_{\mathrm{eff}}$$
   $$\Delta\mathrm{BIC}_{01} = \mathrm{BIC}(\mathcal{M}_0) - \mathrm{BIC}(\mathcal{M}_1)$$
   - Strong empirical support for $\mathcal{M}_1$ over null/alternative $\mathcal{M}_0$ requires $\Delta\mathrm{BIC}_{01} > 10$ (approximate Bayes factor $B_{10} = \exp(\Delta\mathrm{BIC}_{01}/2) > 148$).
2. **Small-Sample Corrected Akaike Information Criterion ($\mathrm{AICc}$):**
   $$\mathrm{AICc}(\mathcal{M}_k) = n_{\mathrm{eff}} \ln\!\left(\frac{\mathrm{RSS}_k}{n_{\mathrm{eff}}}\right) + 2k + \frac{2k(k+1)}{n_{\mathrm{eff}} - k - 1}$$
   - Always report $\mathrm{AICc}$ alongside $\mathrm{BIC}$ when $n_{\mathrm{eff}}/k < 40$.
3. **Autocorrelation & Effective Sample Size ($n_{\mathrm{eff}}$):**
   - When data are measured along spatial lattices or time-series with spatial/temporal autocorrelation, never use raw point count $n$. Compute effective degrees of freedom:
     $$n_{\mathrm{eff}} = \frac{n}{1 + 2 \sum_{j=1}^{\infty} \rho(j)}$$
     where $\rho(j)$ is the empirical autocorrelation at lag $j$.

---

## 4. Finite-Size Scaling (FSS) & Two-Sided Stability

Extrapolations to the thermodynamic limit ($L \to \infty$) must pass two-sided holdout validation:
1. **Leave-$L_{\min}$-Out (Sub-Asymptotic Check):**
   - Verify that excluding the smallest cluster sizes does not significantly alter extracted critical exponents or scales, proving the fit is not dominated by pre-asymptotic lattice corrections.
2. **Leave-$L_{\max}$-Out (Thermodynamic Stability Check):**
   - Verify that the physical scales ($\xi, \sigma_{\mathrm{str}}, \Delta$) inferred when holding out $L_{\max}$ fall within the 95% confidence interval of the full-dataset extrapolation.
3. **Exact Diagonalization (ED) Oracle Calibration:**
   - Validate mean-field or variational predictions against exact diagonalization or DMRG benchmarks on small clusters before scaling up.
