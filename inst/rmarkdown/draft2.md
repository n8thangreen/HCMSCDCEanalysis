---
title: "Long-term cost-effectiveness of subcutaneous versus transvenous implantable cardioverter-defibrillators across sudden cardiac death risk algorithms in hypertrophic cardiomyopathy"
author:
  - Nathan Green, UCL
  - Yang Chen, UCL
date: "2026-10-22"
output:
  pdf_document: 
    keep_tex: yes
  bookdown::pdf_document2:
    extra_dependencies: ["amsmath", "arydshln"]
    toc: false
  word_document: default
  html_document:
    df_print: paged
params:
  data: ipd_risk
bibliography: bibliography.bib
editor_options:
  markdown:
    wrap: sentence
---







# Introduction

Hypertrophic cardiomyopathy (HCM) is a common inherited heart muscle disorder and a leading cause of sudden cardiac death (SCD) in young adults.
Patients at high risk of SCD need to be identified so they can be offered lifesaving prophylactic treatment with an implantable cardioverter-defibrillator (ICD).
Contemporary guidelines recommend that sudden death risk is assessed using quantitative clinical risk prediction models [@OMahony2014; @Ommen2020], yet real-world ICD implantation rates vary substantially across international healthcare systems [@Nauffal2021].

In our previous cost-effectiveness analysis [@Green2024], we evaluated the European Society of Cardiology (ESC) 5-year HCM Risk-SCD prediction algorithm [@OMahony2014] for primary prevention transvenous ICD (TV-ICD) implantation over a 12-year time horizon without device replacement.
However, HCM patients are typically implanted at a young mean age (~48 years) and face decades of living with a device, requiring multiple generator replacements over their lifetime.
Furthermore, the clinical landscape has shifted with the widespread adoption of the subcutaneous ICD (S-ICD) [@Lambiase2022; @Knops2020].
By placing the electrode and generator entirely extrathoracically, the S-ICD eliminates intravascular lead placement and substantially reduces long-term lead-related complications and systemic device infections, as demonstrated in the landmark ATLAS [@Healey2022] and 8-year PRAETORIAN-XL [@Knops2024] trials.
Conversely, S-ICD systems carry a higher upfront hardware acquisition cost and a shorter median battery longevity (~8 to 8.7 years) compared with contemporary single-chamber TV-ICDs (~10 years), necessitating more frequent generator replacements over an extended time horizon.

This paper addresses this critical evidence gap by extending our Bayesian Markov cost-effectiveness framework to a 30-year time horizon using tunnel states to explicitly track device battery longevity and generator replacement cycles.
We evaluate the joint clinical and economic trade-offs between device modality (S-ICD vs. TV-ICD) and risk stratification thresholds (ESC 5-year SCD risk $>4\%$ and $>6\%$, individual risk factor counts, and observed clinical practice) from a UK National Health Service (NHS) perspective.

All code and data are publicly available on GitHub at `https://github.com/n8thangreen/HCMSCDCEanalysis`.

# Data and Model Parameters

The study cohort has been described in detail previously [@OMahony2014; @Green2024].
In brief, the dataset comprises $N =$ 3672 individuals with HCM followed up across six European centres: Athens (474), Bologna (456), Coruna (590), London (1592), Murcia (404), and Naples (156).
The mean age (SD) at baseline was 48.3 (16.5) years.

Health state utilities, resource use, and unit costs were obtained from contemporary clinical trials, UK NHS tariffs, and St Bartholomew's Hospital (Barts Health NHS Trust) procurement data (Table 1).
Under the UK NHS Payment Scheme and Specialised Services Devices Programme (SSDP), implantable defibrillator hardware is a High-Cost Tariff-Excluded Device reimbursed in addition to the Healthcare Resource Group (HRG) implantation procedure tariff (`EY02B`, £4,666).
To ensure a like-for-like comparison between modalities, our base-case analysis uses the combined NHS procedure tariff and Barts Health NHS Trust device hardware cost (£15,317 for TV-ICD vs. £17,434 for S-ICD), with tariff-only and broader European cost bounds (£4,666 to £22,000) evaluated in sensitivity analyses.
Costs per shock episode (£1,766) were taken from the UK CareLink healthcare resource utilization study [@Llewellyn2023], and device-specific complication rates and battery longevity were parameterized using PRAETORIAN-XL [@Knops2024], ATLAS [@Healey2022], and EFFORTLESS [@Lambiase2022].

| Parameter | Variable | TV-ICD (`tvicd`) | S-ICD (`sicd`) | Sensitivity Range / Distribution | Source |
|:---|:---|:---|:---|:---|:---|
| **Health Utilities (QALYs)** | | | | | |
| HCM without ICD | `q_hcm` | 0.88 / yr | 0.88 / yr | 0.60 -- 1.00 | @Sanders2005 |
| Living with ICD (multiplier) | `u_icd` | 0.90 | 0.90 | 0.80 -- 1.00 (0.95 scen.) | @Magnusson2020; @Holbrook2020; @Knops2020 |
| Shock state (multiplier) | `u_shock` | 0.875 | 0.875 | 0.50 -- 1.00 | @Buxton2006 |
| Implant procedure decrement | `u_implant` | -0.048 | -0.048 | -0.096 -- 0.00 | @Holbrook2020 |
| Implant complication decrement | `u_compl` | -0.096 | -0.096 | Point (-0.096) | @Holbrook2020 |
| **Costs (£)** | | | | | |
| Initial implant + device | `c_icd` | £15,317 | £17,434 | TV: £4,666--£20,000; S: £13,788--£22,000 | NHS `EY02B` + Barts Trust SSDP |
| Generator replacement | `c_icd_repl` | £15,317 | £17,434 | TV: £4,666--£20,000; S: £13,788--£22,000 | NHS `EY02B` + Barts Trust SSDP |
| Annual ICD follow-up | `c_appt` | £145 | £145 | Point (£145) | NHS `WF02A` [@NHSEngland2021] |
| Shock episode (hospitalisation) | `c_shock` | £1,766 | £1,766 | £883 -- £6,098 | @Llewellyn2023 |
| Device / lead infection | `c_inf` | £37,116 | £37,116 | Point (£37,116) | @Thijssen2014; NHS `EY17A` |
| Lead dislodgement / revision | `c_dis` | £6,146 | £6,146 | Point (£6,146) | @Thijssen2014; NHS `EY17A` |
| **Device & Clinical Parameters** | | | | | |
| Battery replacement cycle | `t_repl` | 10 years | 8 years | Tunnel states ($T$) | Real-world LATITUDE / Registry |
| Initial implant complication | `p_compl` | 0.043 | 0.015 | Point | @Cunningham2012; @Knops2024; @Healey2022 |
| Initial infection probability | `p_inf_init` | 0.02277 | 0.007 | Point | @Thijssen2014; @Knops2024 |
| Initial dislodgement probability | `p_dis_init` | 0.00828 | 0.008 | Point | @Thijssen2014; @Knops2024 |
| Time horizon | `J` | 30 years | 30 years | 12 -- 30 years | Long-term replacement horizon |
| Discount rate (costs & QALYs) | `delta` | 3.5% | 3.5% | 0% -- 6.0% | NICE Reference Case |

: Base-case model parameters, sensitivity ranges, and data sources for Transvenous ICD (TV-ICD) and Subcutaneous ICD (S-ICD).

# Methods

## Tunnel-State Markov Model with Battery Replacement Cycles

To account for time-since-implant and periodic generator replacements over a $J = 30$-year horizon, we expanded our discrete-time annual Markov model [@Green2024] using tunnel states (Figure \@ref(fig:model)).
For a device with a battery lifespan of $T = t_{\text{repl}}$ years ($T = 10$ for TV-ICD; $T = 8$ for S-ICD), the state space consists of $S = 2T + 4$ states:

1. **HCM with ICD Tunnel States ($1, \dots, T$)**: Patients implanted with an ICD enter State 1 at cycle $j = 1$ and advance from state $t$ to $t+1$ each year they remain free of shock and death. Upon reaching the end of battery life at state $T$, surviving non-shocked patients transition back to State 1 in the following cycle, triggering a generator replacement cost (`c_icd_repl`).
2. **ICD Shock Tunnel States ($T+1, \dots, 2T$)**: If a patient in ICD tunnel state $t$ experiences an appropriate or inappropriate shock, they enter the corresponding shock tunnel state $T+t$ for one cycle (incurring shock cost `c_shock` and utility multiplier `u_shock`). Surviving shocked patients progress to ICD tunnel state $t+1$ in the next cycle (or back to State 1 with device replacement if $t = T$).
3. **All-Cause Mortality with ICD ($2T+1$)**: Absorbing sink state for non-arrhythmic mortality in the ICD arm.
4. **HCM without ICD ($2T+2$)**: Starting state for patients not selected for ICD implantation under a given risk threshold.
5. **Sudden Cardiac Death ($2T+3$)**: Absorbing sink state for fatal arrhythmic events in the non-ICD arm.
6. **All-Cause Mortality without ICD ($2T+4$)**: Absorbing sink state for non-SCD mortality in the non-ICD arm.

Note that setting $T = 1$ and $\text{c\_icd\_repl} = 0$ recovers the original 6-state memoryless model without battery replacement [@Green2024] as an exact special case.

<div class="figure" style="text-align: center">
<img src="../../images/model_diagram_tunnel.png" alt="Tunnel-state Markov model diagram incorporating device battery lifespan and replacement cycles." width="85%" />
<p class="caption">Tunnel-state Markov model diagram incorporating device battery lifespan and replacement cycles.</p>
</div>

## Bayesian Inference of Transition Probabilities

For each risk stratification rule (Observed practice, ESC 5-year SCD risk $>6\%$ and $>4\%$, $\ge 1, 2, 3$ risk factors, and No ICD), individual patient follow-up counts were aggregated into annual transition matrices and fitted in JAGS/WinBUGS [@Lunn2000; @R2021] to obtain 1,000 posterior MCMC draws of the transition probability matrices.

# Results

## Base-Case Cost-Effectiveness: TV-ICD and S-ICD (30-Year Tunnel Model)

Table 2 presents the 30-year expected total costs, QALYs, incremental outcomes, incremental cost-effectiveness ratios (ICERs), and incremental net monetary benefits (INB at a willingness-to-pay threshold of $k = £25,000$/QALY) per patient for both TV-ICD and S-ICD relative to the observed clinical practice strategy within each device modality.


Table: 30-year tunnel-state cost-effectiveness results per patient by device modality (reference: Observed strategy within each device).

|Device |Model         |Strategy  |Mean_Cost |Incremental_Cost |Mean_QALY |Incremental_QALY |ICER                     |INB (£25k) |
|:------|:-------------|:---------|:---------|:----------------|:---------|:----------------|:------------------------|:----------|
|TVICD  |tunnel (10yr) |Risk > 4% |£10,200   |£ 4,767          |13.03     |0.115            |£41,353/QALY             |£-1,885    |
|TVICD  |tunnel (10yr) |Risk > 6% |£ 4,966   |£  -467          |13.02     |0.108            |Dominates (£-4,332/QALY) |£ 3,165    |
|TVICD  |tunnel (10yr) |Observed  |£ 5,433   |Ref              |12.92     |Ref              |Reference                |Ref        |
|TVICD  |tunnel (10yr) |No ICD    |£     0   |£-5,433          |12.76     |-0.155           |£34,946/QALY             |£ 1,546    |
|SICD   |tunnel (8yr)  |Risk > 4% |£13,457   |£ 6,217          |13.03     |0.116            |£53,748/QALY             |£-3,325    |
|SICD   |tunnel (8yr)  |Risk > 6% |£ 6,537   |£  -703          |13.02     |0.108            |Dominates (£-6,519/QALY) |£ 3,400    |
|SICD   |tunnel (8yr)  |Observed  |£ 7,240   |Ref              |12.92     |Ref              |Reference                |Ref        |
|SICD   |tunnel (8yr)  |No ICD    |£     0   |£-7,240          |12.76     |-0.156           |£46,450/QALY             |£ 3,344    |

Across a 30-year time horizon with battery replacements:

* **Within each device modality, the $>6\%$ 5-year SCD risk threshold dominates observed clinical practice**:
  * For **TV-ICD** (10-year replacement cycle), the $>6\%$ risk threshold reduces mean lifetime cost by £467 per patient (£4,966 vs. £5,433) while increasing QALYs by $+0.108$ (13.02 vs. 12.92 QALYs), yielding a positive INB of **£3,165** per patient at $k = £25,000$/QALY.
  * For **S-ICD** (8-year replacement cycle), the $>6\%$ risk threshold saves £703 per patient (£6,537 vs. £7,240) for the same $+0.108$ QALY gain, yielding an INB of **£3,400** per patient relative to S-ICD under observed selection.
* **The $>4\%$ 5-year SCD risk threshold implants more than double the number of patients (1,103 vs. 542)**, which substantially increases initial implantation and repeat battery replacement costs over 30 years (£10,200/patient for TV-ICD; £13,457/patient for S-ICD) for only a marginal additional QALY gain over the $>6\%$ threshold ($+0.008$ QALYs), resulting in ICERs above the standard NICE threshold (£41,353/QALY for TV-ICD and £53,748/QALY for S-ICD vs. Observed).

<div class="figure" style="text-align: center">
<img src="../../images/state_pop_over_time_tvicd_tunnel.png" alt="Cohort state occupancy proportions over the 30-year time horizon under the tunnel-state model." width="95%" />
<p class="caption">Cohort state occupancy proportions over the 30-year time horizon under the tunnel-state model.</p>
</div>

## Joint Comparison: S-ICD vs. TV-ICD Across Risk Stratification Rules

When comparing S-ICD and TV-ICD Strategies jointly against the benchmark of current clinical practice (**TV-ICD: Observed**, Figure \@ref(fig:cecombined)):

* Switching from **TV-ICD: Observed** (£5,529, 13.83 QALYs in the PSA) to **S-ICD combined with the $>6\%$ Cox risk algorithm** (`S-ICD: Cox > 6%`, £6,574, 13.91 QALYs) is **cost-effective at the NICE threshold**, with an **ICER of £13,548/QALY**, an expected incremental benefit (EIB) of **+£883**, and a **56.3% probability of cost-effectiveness** at $k = £25,000$/QALY.
* However, within any *fixed* risk stratification rule (Table 3 and Figure \@ref(fig:inmbforest)), S-ICD incurs higher 30-year cumulative costs than TV-ICD ($+£1,556$ per patient under Risk $>6\%$; $+£1,823$ under Observed; $+£3,244$ under Risk $>4\%$). Although S-ICD reduces acute implant and lead complications (`1.5%` vs. `4.3%`), the combination of a £2,117 higher device cost per implant/replacement and an 8-year vs. 10-year battery replacement cycle (3 replacements vs. 2 replacements over 30 years) outweighs the acute complication savings when long-term living-with-ICD utility (`u_icd`) is assumed equal between devices.


Table: Head-to-head comparison of S-ICD versus TV-ICD within each of the 8 risk stratification rules over 30 years (WTP = £25,000/QALY).

|      |Strategy  |TV-ICD Cost |S-ICD Cost |Delta Cost (S - TV) |Delta QALY (S - TV) |INMB at £25k (95% CrI)     |
|:-----|:---------|:-----------|:----------|:-------------------|:-------------------|:--------------------------|
|2.5%  |Risk > 6  |£ 4,939     |£ 6,495    |+£1,556             |+0.0004             |£-1,546 (£-1,625, £-1,452) |
|2.5%1 |Risk > 4  |£10,177     |£13,421    |+£3,244             |+0.0008             |£-3,224 (£-3,343, £-3,096) |
|2.5%2 |any rf    |£16,453     |£21,739    |+£5,286             |+0.0013             |£-5,253 (£-5,397, £-5,099) |
|2.5%3 |over 1 rf |£ 4,449     |£ 5,863    |+£1,414             |+0.0004             |£-1,406 (£-1,490, £-1,323) |
|2.5%4 |over 2 rf |£   780     |£ 1,036    |+£  256             |+0.0001             |£  -255 (£  -275, £  -222) |
|2.5%5 |over 3 rf |£    12     |£    13    |+£    2             |+0.0000             |£    -2 (£    -5, £    -1) |
|2.5%6 |icd       |£ 5,459     |£ 7,283    |+£1,823             |+0.0004             |£-1,813 (£-1,875, £-1,735) |
|2.5%7 |no icd    |£     0     |£     0    |+£    0             |+0.0000             |£     0 (£     0, £     0) |

<div class="figure" style="text-align: center">
<img src="../../images/ce_plane_combined.png" alt="Combined Probabilistic Sensitivity Analysis (30-year tunnel model) comparing TV-ICD and S-ICD strategies against TV-ICD Observed practice: (left) Cost-Effectiveness Plane and (right) Cost-Effectiveness Acceptability Curves (CEAC)." width="48%" /><img src="../../images/ceac_combined.png" alt="Combined Probabilistic Sensitivity Analysis (30-year tunnel model) comparing TV-ICD and S-ICD strategies against TV-ICD Observed practice: (left) Cost-Effectiveness Plane and (right) Cost-Effectiveness Acceptability Curves (CEAC)." width="48%" />
<p class="caption">Combined Probabilistic Sensitivity Analysis (30-year tunnel model) comparing TV-ICD and S-ICD strategies against TV-ICD Observed practice: (left) Cost-Effectiveness Plane and (right) Cost-Effectiveness Acceptability Curves (CEAC).</p>
</div>

<div class="figure" style="text-align: center">
<img src="../../images/inmb_strategies_forest.png" alt="Incremental Net Monetary Benefit (INMB at WTP £25,000/QALY) of S-ICD vs. TV-ICD across all 8 risk stratification scenarios." width="75%" />
<p class="caption">Incremental Net Monetary Benefit (INMB at WTP £25,000/QALY) of S-ICD vs. TV-ICD across all 8 risk stratification scenarios.</p>
</div>

## Sensitivity Analyses

Figure \@ref(fig:tornado) displays the one-way deterministic sensitivity analysis (DSA) tornado plots for Expected Incremental Benefit (EIB at $k = £25,000$/QALY) for TV-ICD and S-ICD.
Consistent with our 12-year analysis [@Green2024], the dominant driver of cost-effectiveness across a 30-year horizon is the health-related quality of life multiplier associated with living with an ICD (`u_icd`), followed by baseline HCM utility (`q_hcm`) and device acquisition/replacement cost (`c_icd`).
In particular, if avoiding transvenous leads with an S-ICD yields even a modest chronic quality-of-life advantage (e.g., improving `u_icd` from `0.90` to `0.92`--`0.95` due to reduced lead anxiety or avoidance of long-term intravascular extraction morbidity), S-ICD becomes directly cost-effective over TV-ICD within the $>6\%$ risk group.

<div class="figure" style="text-align: center">
<img src="../../images/tornado_plots_EIB.png" alt="One-way deterministic sensitivity analysis tornado plots of Expected Incremental Benefit (EIB, £ at WTP £25,000) for (top) TV-ICD and (bottom) S-ICD under the &gt;4% and &gt;6% 5-year SCD risk thresholds." width="90%" />
<p class="caption">One-way deterministic sensitivity analysis tornado plots of Expected Incremental Benefit (EIB, £ at WTP £25,000) for (top) TV-ICD and (bottom) S-ICD under the >4% and >6% 5-year SCD risk thresholds.</p>
</div><div class="figure" style="text-align: center">
<img src="../../images/tornado_plots_EIB_sicd.png" alt="One-way deterministic sensitivity analysis tornado plots of Expected Incremental Benefit (EIB, £ at WTP £25,000) for (top) TV-ICD and (bottom) S-ICD under the &gt;4% and &gt;6% 5-year SCD risk thresholds." width="90%" />
<p class="caption">One-way deterministic sensitivity analysis tornado plots of Expected Incremental Benefit (EIB, £ at WTP £25,000) for (top) TV-ICD and (bottom) S-ICD under the >4% and >6% 5-year SCD risk thresholds.</p>
</div>

# Discussion

By expanding our Bayesian Markov model to a 30-year horizon with tunnel states for generator replacement cycles, this study provides two key insights for HCM clinical practice and health technology assessment:

1. **Risk Stratification Efficiency Can Fund Technological Innovation**:
   Adopting the ESC $>6\%$ 5-year SCD risk threshold rather than historic/unstructured clinical selection reduces unnecessary implantations in low-risk patients who never experience a life-saving shock.
   Over 30 years, this improved targeting makes **S-ICD under the $>6\%$ threshold cost-effective compared with current observed TV-ICD practice** (ICER £13,548/QALY), because the savings from avoiding low-yield implants partially offset the higher hardware and 8-year replacement costs of S-ICD.
2. **Battery Longevity and Chronic HRQoL Drive the S-ICD vs. TV-ICD Trade-off**:
   In a head-to-head comparison within the same risk stratum, S-ICD's lower perioperative complication and infection rates [@Knops2024; @Healey2022] are outweighed financially over 30 years by its £2,117 higher device cost and shorter battery life (8 vs. 10 years, requiring an extra generator replacement over three decades).
   Consequently, from a pure NHS payer perspective, S-ICD is most economically attractive when targeted to patients with higher infection/lead-failure risk, younger patients in whom transvenous lead failure accumulates over decades, or when S-ICD confers a chronic quality-of-life benefit (`u_icd`) over intravascular leads.

\newpage
# References

<div id="refs"></div>
