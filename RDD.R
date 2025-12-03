library(rdrobust)
library(dplyr)
library(ggplot2)

set.seed(42)

# Data Generation Process
N <- 2000
cutoff_c <- 0
tau_RD <- 0.05  # 5% price jump

# Running Variable: Standardized Market Cap Rank
# Simulate ranking centered at c=0
X <- runif(N, min = -150, max = 150)
X <- X[order(X)] # Sort the running variable

# Treatment Assignment: Sharp RDD
# D_i = 1 if the stock is ranked >= 0 (e.g., included in R2000), 0 otherwise.
D <- ifelse(X >= cutoff_c, 1, 0)

# Potential Outcome Function (f(X_i)): Smooth, continuous trend (linear + quadratic)
# This represents the continuous effect of market capitalization on CAR absent the index shock.
f_X <- 0.001 * X - 0.00005 * X^2

# Noise: Standard error term
epsilon <- rnorm(N, mean = 0, sd = 0.02)

# Outcome Variable: Cumulative Abnormal Return (CAR)
# Y_i = f(X_i) + tau_RD * D_i + epsilon_i
Y <- f_X + tau_RD * D + epsilon

df_rdd <- data.frame(CAR = Y, Rank = X, Treated = D)

# Summary Statistics Table

# Filter data around a reasonable window for presentation
df_summary <- df_rdd %>% 
  filter(abs(Rank) <= 100) # Only consider firms within 100 ranks of the cutoff

# Calculate descriptive statistics by assignment status
summary_stats <- df_summary %>%
  group_by(Treated) %>%
  summarise(
    N_Obs = n(),
    Mean_Rank = mean(Rank),
    SD_Rank = sd(Rank),
    Mean_CAR = mean(CAR),
    SD_CAR = sd(CAR))

summary_stats <- summary_stats %>% 
  mutate(
    Group = case_when(
      Treated == 0 ~ "Control (R1000, X < 0)",
      Treated == 1 ~ "Treatment (R2000, X >= 0)"
    )
  ) %>%
  select(Group, N_Obs, Mean_Rank, Mean_CAR, SD_CAR)

cat("### Summary Statistics Table (Simulated Data, |Rank| <= 100)\n")
print(summary_stats)
cat("\n")

# RDD Estimation

# Run rdrobust to estimate LATE
# - c=0: Cutoff (Standardized Rank 1000)
# - p=1: Local Linear Regression is the preferred specification [1, 2]
# - kernel="triangular": Standard kernel for RDD [3]
# - h: Optimal bandwidth selection (CCT MSE-optimal bandwidth [4]) is automatic
rdd_estimate <- rdrobust(
  y = df_rdd$CAR,
  x = df_rdd$Rank,
  c = cutoff_c,
  p = 1,
  kernel = "triangular")

tau_rd     <- rdd_estimate$coef["Robust", "Coeff"]
se_robust  <- rdd_estimate$se["Robust", "Std. Err."]
z_value    <- rdd_estimate$z["Robust", "z"]
# Compute robust p-value
p_value    <- 2 * (1 - pnorm(abs(z_value)))

# Robust CI
ci_lower <- rdd_estimate$ci["Robust", "CI Lower"]
ci_upper <- rdd_estimate$ci["Robust", "CI Upper"]
bw_robust <- rdd_estimate$bws["b", "left"]

cat("     Robust RDD Treatment Effect     \n")
cat("τ_RD (LATE estimate): ", format(tau_rd, digits=4), "\n")
cat("Robust Std. Error:   ", format(se_robust, digits=4), "\n")
cat("z-value (Robust):    ", format(z_value, digits=4), "\n")
cat("P-value (Robust):    ", format(p_value, digits=4), "\n")
cat("95% CI (Robust):     [",
    format(ci_lower, digits=4), ", ",
    format(ci_upper, digits=4), "]\n", sep="")
cat("------------------------------------\n")
cat("Bandwidth (Robust):  ", format(bw_robust, digits=4), "\n")
cat("Effective N:         ", rdd_estimate$N_h[1] + rdd_estimate$N_h[2], "\n")

# RDD Plot
cat("### RDD Plot (Visualizing the Discontinuity at Cutoff)\n")
rdd_plot_output <- rdplot(
  y = df_rdd$CAR, 
  x = df_rdd$Rank, 
  c = cutoff_c,
  title = "Figure 1: Cumulative Abnormal Returns at Russell 1000/2000 Cutoff",
  x.label = "Standardized Market Cap Rank (X_i, Cutoff = 0)",
  y.label = "Cumulative Abnormal Returns (Y_i)",
  p = 1, # Local Linear Regression fit
  kernel = "triangular")

cat("\n----------------------------------------------------------------------\n")

# Robustness Check A: Pre-Treatment Covariate Balance Test
# Concept: Check for discontinuity in a predetermined covariate (e.g., T-1 CAR).
# A valid RDD requires the LATE on this covariate to be zero.

# Simulated Covariate Z should have a smooth trend f_Z(X) and NO DISCONTINUITY (tau_Z = 0)
f_Z <- 0.005 * X + 0.0001 * X^2 # Smooth trend
Z_T1_CAR <- f_Z + rnorm(N, mean = 0, sd = 0.015) # T-1 CAR with noise

extract_rd_results <- function(rdobj) {
  tau  <- rdobj$coef["Robust", "Coeff"]
  se   <- rdobj$se["Robust", "Std. Err."]
  z    <- rdobj$z["Robust", "z"]
  pval <- 2 * (1 - pnorm(abs(z)))
  ci_l <- rdobj$ci["Robust", "CI Lower"]
  ci_u <- rdobj$ci["Robust", "CI Upper"]
  h    <- rdobj$bws["b", "left"]
  
  return(list(
    tau = tau, se = se, z = z, p = pval,
    ci_l = ci_l, ci_u = ci_u, h = h
  ))}

# Balance test RDD
rdd_balance_check <- rdrobust(
  y = Z_T1_CAR,
  x = X,
  c = cutoff_c,
  p = 1,
  kernel = "triangular"
)

bal <- extract_rd_results(rdd_balance_check)

cat("\n### Robustness Check A: Covariate Balance (T-1 CAR)\n")
balance_table <- data.frame(
  Metric = c("LATE on T-1 CAR", "Robust SE", "z-value", "p-value", "95% CI", "h"),
  Value  = c(
    format(bal$tau, digits=4),
    format(bal$se, digits=4),
    format(bal$z, digits=4),
    format(bal$p, digits=4),
    paste0("[", format(bal$ci_l, digits=4), ", ", format(bal$ci_u, digits=4), "]"),
    format(bal$h, digits=4)))
print(balance_table)

placebo_cutoff <- 80

rdd_placebo <- rdrobust(
  y = df_rdd$CAR,
  x = X,
  c = placebo_cutoff,
  p = 1,
  kernel = "triangular"
)

pl <- extract_rd_results(rdd_placebo)

cat("\n### Robustness Check B: Placebo Cutoff at 80\n")
placebo_table <- data.frame(
  Metric = c("Placebo τ", "Robust SE", "z-value", "p-value", "95% CI", "h"),
  Value  = c(
    format(pl$tau, digits=4),
    format(pl$se, digits=4),
    format(pl$z, digits=4),
    format(pl$p, digits=4),
    paste0("[", format(pl$ci_l, digits=4), ", ", format(pl$ci_u, digits=4), "]"),
    format(pl$h, digits=4)))
print(placebo_table)

library(rddensity)
# McCrary Density Test at cutoff
mcc <- rddensity(X, c = cutoff_c)
summary(mcc)