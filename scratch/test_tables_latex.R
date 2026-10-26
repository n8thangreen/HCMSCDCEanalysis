library(knitr)
library(kableExtra)
library(dplyr)

tex_doc <- '
\\documentclass[11pt]{article}
\\usepackage[margin=1in]{geometry}
\\usepackage{amsmath}
\\usepackage{graphicx}
\\usepackage{booktabs}
\\usepackage{float}
\\usepackage{array}

\\begin{document}

\\section{Table 4 Test}
'

tab_clin <- data.frame(
  Strategy = c("No ICD (Do Nothing)", "ESC Risk > 6%", "Observed Practice", "ESC Risk > 4%"),
  Implants = c("0 (0.0%)", "538 (14.7%)", "559 (15.2%)", "1106 (30.1%)"),
  `30-Yr Surv` = c("49.6%", "55.2%", "54.5%", "57.9%"),
  `Life-Years` = c("14.54 yrs", "15.11 yrs", "15.01 yrs", "15.38 yrs"),
  `SCD Deaths` = c("653", "397", "421", "275"),
  `SCDs Saved` = c("—", "255", "232", "378"),
  `NNT (SCD)` = c("—", "2.1", "2.4", "2.9"),
  `NNT (All)` = c("—", "2.6", "3.1", "3.6"),
  `LYs / Imp` = c("—", "+3.93", "+3.11", "+1.91"),
  check.names = FALSE
)

tbl4 <- kable(tab_clin, format = "latex", booktabs = TRUE, caption = "Clinical Outcomes") %>%
  kable_styling(font_size = 8, latex_options = "HOLD_position") %>%
  column_spec(1, width = "2.4cm")

tbl4_mod <- sub("\\\\begin\\{tabular\\}", "\\\\setlength{\\\\tabcolsep}{3pt}\\\\begin{tabular}", as.character(tbl4))

# Table 7
threshold_df <- data.frame(
  `Decision Context` = c(
    "Head-to-Head Parity (Risk > 6%)",
    "Head-to-Head Parity (Risk > 4%)",
    "Head-to-Head Parity (Observed)",
    "Practice Transition (S-ICD > 6% vs. TV-ICD Obs)"
  ),
  `Break-Even (£25k)` = c("£12,761 (-26.8%)", "£12,753 (-26.8%)", "£12,620 (-27.6%)", "£22,230 (+27.5%)"),
  `Break-Even (£35k)` = c("£12,773 (-26.7%)", "£12,765 (-26.8%)", "£12,631 (-27.5%)", "£25,469 (+46.1%)"),
  `Current ICER` = c("Dominated", "Dominated", "Dominated", "£10,198 / QALY"),
  `Economic Rationale` = c(
    "Must be discounted ~£2,550 below TV-ICD list price (£15,317) to compensate for 8-yr vs 10-yr battery replacement penalty.",
    "Must offset higher 30-year replacement frequency across 1,103 implanted patients.",
    "Must offset replacement burden under unguided historic selection.",
    "Risk > 6% clinical efficiency (avoiding low-risk implants) cross-subsidizes S-ICD hardware premium."
  ),
  check.names = FALSE
)

tbl7 <- kable(threshold_df, format = "latex", booktabs = TRUE, caption = "Threshold Analysis") %>%
  kable_styling(font_size = 8.5, latex_options = "HOLD_position") %>%
  column_spec(1, width = "3.0cm") %>%
  column_spec(5, width = "4.2cm")

tbl7_mod <- sub("\\\\begin\\{tabular\\}", "\\\\setlength{\\\\tabcolsep}{4pt}\\\\begin{tabular}", as.character(tbl7))

full_tex <- paste0(tex_doc, "\n", tbl4_mod, "\n\\newpage\n\\section{Table 7 Test}\n", tbl7_mod, "\n\\end{document}")
writeLines(full_tex, "scratch/test_tables.tex")
