$p = 'inst\rmarkdown\draft2.Rmd'
$t = [IO.File]::ReadAllText($p)
$nl = if ($t -match "`r`n") { "`r`n" } else { "`n" }
function Swap-Text($old, $new) {
  if (-not $script:t.Contains($old)) { Write-Host "MISSING: $old"; return }
  $script:t = $script:t.Replace($old, $new)
}
Swap-Text '| Battery replacement cycle | `t_repl_max` | 10 years | 8 years | Tunnel states ($T$) | Real-world LATITUDE / Registry |' ('| Maximum device life (tunnel states, $T$) | `t_repl_max` | 15 years | 15 years | Forced replacement at $T$ | Model assumption |' + $nl + '| Median battery life | `repl_median` | 10 years (SE 0.5) | 8.5 years (SE 0.41) | Normal; S-ICD 8.5 years from EMBLEM, 95% CI 7.9--9.5 | Real-world LATITUDE / Registry; EMBLEM |' + $nl + '| Weibull shape of battery life | `repl_shape` | 5 (SE 0.5) | 5 (SE 0.5) | Gamma | Assumption (increasing hazard of battery depletion) |')
Swap-Text 'For a device with a battery lifespan of $T = t_{\text{repl}}$ years ($T = 10$ for TV-ICD; $T = 8$ for S-ICD), the state space consists of $S = 2T + 4$ states:' 'For a device with a maximum battery lifespan of $T = t_{\text{repl\_max}}$ years ($T = 15$ for both devices), the state space consists of $S = 2T + 4$ states:'
Swap-Text 'Upon reaching the end of battery life at state $T$, surviving non-shocked patients transition back to State 1 in the following cycle, triggering a generator replacement cost (`c_icd_repl`).' 'In each tunnel state $t$, patients who remain alive are replaced with probability $h_t$ (see the next section) and transition back to State 1 in the following cycle, triggering a generator replacement cost (`c_icd_repl`); otherwise they advance to state $t+1$. Replacement is certain at state $T$ ($h_T = 1$).'
Swap-Text 'Surviving shocked patients progress to ICD tunnel state $t+1$ in the next cycle (or back to State 1 with device replacement if $t = T$).' 'Surviving shocked patients are replaced with probability $h_t$ (back to State 1) or progress to ICD tunnel state $t+1$ in the next cycle.'
Swap-Text 'Note that setting $T = 1$ and $\text{c\_icd\_repl} = 0$' 'Note that setting $T = 1$ (`t_repl_max = 1`) and $\text{c\_icd\_repl} = 0$'
$new = @'
## Survival-Curve Model of Device Battery Longevity

An earlier version of the model replaced every device deterministically after exactly 10 (TV-ICD) or 8 (S-ICD) years.
This forces an artificial step-function into the Markov trace, with all replacement costs, replacement-related complications and disutilities falling in a single cycle.
In practice battery depletion varies between devices.
We therefore modelled time to generator replacement with a Weibull survival function $S(t) = \exp\{-(t/\lambda)^{k}\}$, with shape $k > 1$ so that the hazard of battery depletion increases with device age.
The scale was parameterised from the median battery life $m$ as $\lambda = m / (\log 2)^{1/k}$, so that $S(m) = 0.5$.
The median was set to 10 years for TV-ICD and 8.5 years for S-ICD (EMBLEM; 95% CI 7.9 to 9.5), close to the fixed lifespans used previously, and the shape was set to $k = 5$.

The annual probability that a device in tunnel state $t$ (time since implant, in years) is replaced is the discrete-time hazard
$$h_t = 1 - \frac{S(t)}{S(t-1)}, \quad t = 1, \dots, T,$$
with $h_T = 1$ so that all devices are replaced by the maximum lifespan $T = 15$ years (the Weibull survival probability at 15 years is below 1% for both devices).
In each cycle a patient in ICD tunnel state $t$ who does not die and is not shocked is replaced with probability $h_t$ (returning to State 1) and otherwise moves to state $t+1$; patients in shock tunnel states follow the same rule after surviving the shock cycle.
Replacement costs and disutilities are therefore incurred in proportion to the occupancy of State 1 after the first cycle, rather than in a single cycle.
The deterministic fixed-lifespan model is recovered as the special case $h_t = 0$ for $t < T$.

In the probabilistic sensitivity analysis the replacement probabilities are treated in the same way as the MCMC posterior draws of the other transition probabilities.
For each simulation we draw the median battery life from a Normal distribution (truncated above zero) and the shape from a Gamma distribution, using the standard errors in Table 1, compute the vector $(h_1, \dots, h_T)$, and pass this as a simulation-specific row of the transition probability array.
The standard error of the S-ICD median (0.41 years) was derived from the reported 95% confidence interval; that of the TV-ICD median was assumed to be 5% of the median, and the shape standard error was assumed to be 10% of the shape.
The Weibull curve is assumed to apply from implantation; no peri-operative mortality is attached to replacement, and replacement carries the complication risk `p_repl_compl` and disutility `u_repl` as before.

'@
Swap-Text '## Bayesian Inference of Transition Probabilities' ($new.Replace("`r`n", "`n").Replace("`n", $nl) + '## Bayesian Inference of Transition Probabilities')
[IO.File]::WriteAllText($p, $t, (New-Object System.Text.UTF8Encoding($false)))
