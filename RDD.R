library(rdrobust)
library(dplyr)
library(ggplot2)
library(rddensity)

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

# Plot Mean CAR
plot_data_car <- summary_stats %>%
  select(Group, Mean_CAR)

ggplot(plot_data_car, aes(x = Group, y = Mean_CAR, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(Mean_CAR, 4)), 
            vjust = -0.5, size = 4.5) +
  scale_fill_manual(values = c(
    "Control (R1000, X < 0)" = "#8FB8DE",
    "Treatment (R2000, X >= 0)" = "#F4A582")) +
  labs(
    x = "",
    y = "Mean CAR"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# RDD Estimation
# - c=0: Cutoff (Standardized Rank 1000)
# - p=1: Local Linear Regression is the preferred specification
# - kernel="triangular": Standard kernel for RDD
# - h: Optimal bandwidth selection (CCT MSE-optimal bandwidth) is automatic
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
cat("### RDD Plot (Discontinuity at Cutoff)\n")
rdd_plot_output <- rdplot(
  y = df_rdd$CAR, 
  x = df_rdd$Rank, 
  c = cutoff_c,
  title = "",
  x.label = "Standardized Market Cap Rank",
  y.label = "Cumulative Abnormal Returns",
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

cat("\n### Robustness Check B: Multiple Placebo Cutoffs\n")

placebo_cutoffs <- c(-120, -80, -40, 40, 80, 120)
placebo_results <- data.frame()

for (pc in placebo_cutoffs) {
  if (pc >= min(X) && pc <= max(X)) {
    rdd_p <- rdrobust(y = df_rdd$CAR, x = X, c = pc, p = 1, kernel = "triangular")
    pl <- extract_rd_results(rdd_p)
    placebo_results <- rbind(placebo_results, data.frame(
      Cutoff = pc,
      tau_Placebo = pl$tau,
      Robust_SE = pl$se,
      p_value = pl$p,
      Bandwidth = pl$h))}}

placebo_results$tau_Placebo <- format(placebo_results$tau_Placebo, digits=3)
placebo_results$Robust_SE <- format(placebo_results$Robust_SE, digits=3)
placebo_results$p_value <- format(placebo_results$p_value, digits=3)
placebo_results$Bandwidth <- format(placebo_results$Bandwidth, digits=3)

print(placebo_results)

# McCrary Density Test at cutoff
mcc <- rddensity(X, c = cutoff_c)
cat("\n### Robustness Check C: McCrary Density Test\n")
summary(mcc)
