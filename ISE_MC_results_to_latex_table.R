library(tidyverse) # Includes dplyr and stringr used below

# 1. Set path and read the CSV
path <- getwd()
df <- read.csv(file.path(path, "ISE_MC_results.csv"))

# 2. Process the data
df_processed <- df %>%
  mutate(
    # Calculate n from k (n = k(k+1)/2)
    n = k * (k + 1) / 2,
    # Format target variable
    Target = paste0("$m_", j, "$"),
    # Multiply metrics by 10^6 and round to integer
    Mean = round(mean_ISE_MC * 1e6),
    SD = round(sd_ISE_MC * 1e6),
    Median = round(median_ISE_MC * 1e6),
    IQR = round(IQR_ISE_MC * 1e6),
    # Format method to match template
    method_str = paste0("D-", method)
  ) %>%
  # Ensure custom ordering: GM, NW, LL
  mutate(method = factor(method, levels = c("GM", "NW", "LL"))) %>%
  arrange(j, n, method)

# 3. Create a helper to pad numbers with ~ for proper LaTeX alignment
pad_tilde <- function(x, width) {
  str_pad(as.character(x), width, side = "left", pad = "~")
}

# Calculate maximum string widths dynamically based on the dataset
w_mean <- max(nchar(as.character(df_processed$Mean)))
w_sd <- max(nchar(as.character(df_processed$SD)))
w_med <- max(nchar(as.character(df_processed$Median)))
w_iqr <- max(nchar(as.character(df_processed$IQR)))

df_processed <- df_processed %>%
  mutate(
    Mean_pad = pad_tilde(Mean, w_mean),
    SD_pad = pad_tilde(SD, w_sd),
    Median_pad = pad_tilde(Median, w_med),
    IQR_pad = pad_tilde(IQR, w_iqr)
  )

# 4. Construct LaTeX header (Note: backslashes are escaped in R strings)
latex_header <- "\\begin{table}[htbp]
\\caption{Comparison of the D-GM, D-NW, and D-LL methods based on the mean, median, standard deviation (SD), and interquartile range (IQR) of 100 $\\smash{\\widetilde{\\mathrm{ISE}}}$ values, multiplied by $10^6$, for regression functions $m_1$ through $m_6$ and sample sizes $n\\in \\{28, 55, 105\\}$. The integrals in the definition of the D-GM estimator are computed using the $\\mathsf{R}$ command \\texttt{adaptIntegrate} with a relative tolerance of $10^{-3}$.\\label{tab:1}}

\\bigskip
\\renewcommand{\\arraystretch}{0.82} % Adjust the value as needed
\\setlength{\\tabcolsep}{10pt} % Adjust the value as needed
\\centering
{\\footnotesize
\\begin{tabular}{ccccccc}
Target & $n$ & Method & Mean & SD & Median & IQR \\\\"

# 5. Construct LaTeX body
body_lines <- c()
targets <- unique(df_processed$j)

for (tj in targets) {
  body_lines <- c(body_lines, "\\toprule")
  df_target <- df_processed %>% filter(j == tj)
  ns <- unique(df_target$n)
  
  for (idx_n in seq_along(ns)) {
    # Add midrule between different 'n' blocks of the same target
    if (idx_n > 1) {
      body_lines <- c(body_lines, "\\midrule")
    }
    tn <- ns[idx_n]
    df_n <- df_target %>% filter(n == tn)
    
    for (idx_m in seq_len(nrow(df_n))) {
      row <- df_n[idx_m, ]
      if (idx_m == 1) {
        # First row for a specific n gets the multirow commands
        line <- sprintf("\\multirow{3}{*}{%s} & \\multirow{3}{*}{%d}\n   & %s & %s & %s & %s & %s \\\\", 
                        row$Target, row$n, row$method_str, 
                        row$Mean_pad, row$SD_pad, row$Median_pad, row$IQR_pad)
      } else {
        # Subsequent rows get empty slots for Target and n
        line <- sprintf("   &    & %s & %s & %s & %s & %s \\\\", 
                        row$method_str, 
                        row$Mean_pad, row$SD_pad, row$Median_pad, row$IQR_pad)
      }
      body_lines <- c(body_lines, line)
    }
  }
}

body_lines <- c(body_lines, "\\bottomrule")

# 6. Construct LaTeX footer
latex_footer <- "\\end{tabular}}\n\\end{table}"

# 7. Combine and write to file
latex_full <- paste(
  latex_header,
  paste(body_lines, collapse = "\n"),
  latex_footer,
  sep = "\n"
)

out_file <- file.path(path, "ISE_MC_results_table.tex")
writeLines(latex_full, out_file)
cat("LaTeX table successfully written to:", out_file, "\n")