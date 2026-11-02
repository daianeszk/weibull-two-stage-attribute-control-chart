## ============================================================
## npˆs2_wb comparison based on LCL
## ============================================================

options(digits = 16)

## ============================================================
## Step 00 - Global parameters
## ============================================================

alpha_global_h0 <- 0.0027
ARL0_target <- 1 / alpha_global_h0

ASS_min <- 3
ASS_max <- 4

# Exact ASS0 values being investigated
ASS0_allowed <- c(3.0, 3.10, 3.20)

tol_ASS0 <- 1e-10

n_a_vec <- 1:ASS_max

n_b_max <- 20
n_b_vec <- 1:n_b_max

alpha_a_vec <- c(0.1)

## Table 6 shape parameters
k0_vec <- c(0.5, 3, 5, 10, 15, 20)

## Scale-invariant reference mean
mu0 <- 1

## ============================================================
## Table for decreases in mean
## ============================================================

d_vec <- c(
  -0.40,
  -0.30,
  -0.20,
  -0.10,
  -0.05,
  -0.01,
  0
)


cat("\n============================================================\n")
cat("GLOBAL SETTINGS - LCL COMPARISON\n")
cat("============================================================\n")
cat("alpha_global_h0 =", alpha_global_h0, "\n")
cat("ARL0_target     =", ARL0_target, "\n")
cat("ASS_min         =", ASS_min, "\n")
cat("ASS_max         =", ASS_max, "\n")
cat("ASS0 allowed    =", paste(ASS0_allowed, collapse = ", "), "\n")
cat("n_b_max         =", n_b_max, "\n")
cat("mu0             =", mu0, "\n")
cat("k0_vec          =", paste(k0_vec, collapse = ", "), "\n")
cat("d_vec           =", paste(d_vec, collapse = ", "), "\n")
cat("Table 6 benchmark n values: 5, 10 and 30\n")
cat("Chart direction: LOWER / LCL\n")


## ============================================================
## Step 01 - Helper functions
## ============================================================

.validate_inputs <- function(Y, n_sample, alpha) {
  
  if (
    length(n_sample) != 1 ||
    !is.finite(n_sample) ||
    n_sample <= 0 ||
    n_sample != as.integer(n_sample)
  ) {
    stop("n_sample must be a positive integer.")
  }
  
  if (
    length(alpha) != 1 ||
    !is.finite(alpha) ||
    alpha <= 0 ||
    alpha >= 1
  ) {
    stop("alpha must be in (0,1).")
  }
  
  if (
    any(!is.finite(Y)) ||
    any(Y != as.integer(Y)) ||
    any(Y < 1) ||
    any(Y > n_sample)
  ) {
    stop("Y must contain integers in the interval [1, n_sample].")
  }
  
  invisible(TRUE)
}


find_p1 <- function(Y, n_sample, alpha) {
  
  .validate_inputs(
    Y        = Y,
    n_sample = n_sample,
    alpha    = alpha
  )
  
  qbeta(
    alpha,
    shape1 = Y,
    shape2 = n_sample - Y + 1
  )
}


make_stage_design <- function(stage, n_sample, alpha, k0, lambda0) {
  
  Y_vec <- seq_len(n_sample)
  
  p1_vec <- find_p1(
    Y        = Y_vec,
    n_sample = n_sample,
    alpha    = alpha
  )
  
  ## Lower critical time - used for Table 6
  tcrit_low_vec <- qweibull(
    p1_vec,
    shape = k0,
    scale = lambda0
  )
  
  ## Kept only for reference/output
  tcrit_high_vec <- qweibull(
    1 - p1_vec,
    shape = k0,
    scale = lambda0
  )
  
  data.frame(
    stage       = stage,
    n_sample    = n_sample,
    alpha_stage = alpha,
    Y           = Y_vec,
    p1          = p1_vec,
    tcrit_low   = tcrit_low_vec,
    tcrit_high  = tcrit_high_vec,
    stringsAsFactors = FALSE
  )
}


make_design_pairs <- function(design_a, design_b) {
  
  design_a_perf <- data.frame(
    Y_a          = design_a$Y,
    p1_a         = design_a$p1,
    tcrit_low_a  = design_a$tcrit_low,
    tcrit_high_a = design_a$tcrit_high,
    stringsAsFactors = FALSE
  )
  
  design_b_perf <- data.frame(
    Y_b          = design_b$Y,
    p1_b         = design_b$p1,
    tcrit_low_b  = design_b$tcrit_low,
    tcrit_high_b = design_b$tcrit_high,
    stringsAsFactors = FALSE
  )
  
  idx <- expand.grid(
    row_a = seq_len(nrow(design_a_perf)),
    row_b = seq_len(nrow(design_b_perf))
  )
  
  out <- cbind(
    design_a_perf[idx$row_a, , drop = FALSE],
    design_b_perf[idx$row_b, , drop = FALSE]
  )
  
  rownames(out) <- NULL
  
  out
}


## ============================================================
## Step 02 - Out-of-control scenarios
## Decreases in mean
## ============================================================

scenario_grid <- data.frame(
  d = d_vec
)

scenario_grid$scenario_id <- paste0(
  "d",
  gsub("-", "m", gsub("\\.", "p", as.character(scenario_grid$d)))
)


## ============================================================
## Step 03 - TABLE 6 BENCHMARK
##
## TRUE ARL values from Table 6
## n = 5, 10 and 30
##
## DECREASES IN MEAN = LCL
## ============================================================

table6_benchmark <- data.frame(
  
  k0 = c(
    rep(0.5, 7),
    rep(3,   7),
    rep(5,   7),
    rep(10,  7),
    rep(15,  7),
    rep(20,  7)
  ),
  
  d = rep(d_vec, times = 6),
  
  
  ## ----------------------------------------------------------
  ## TRUE ARL - n = 5
  ## ----------------------------------------------------------
  
  table6_true_ARL_n5 = c(
    
    ## delta = 0.5
    126.582, 178.571, 210.970, 310.559,
    326.797, 362.318, 370.370,
    
    ## delta = 3
    2.458, 7.074, 25.667, 89.766,
    153.374, 251.256, 370.370,
    
    ## delta = 5
    1.042, 1.941, 8.007, 52.083,
    144.508, 299.401, 370.370,
    
    ## delta = 10
    1.000, 1.000, 1.351, 12.431,
    72.568, 280.898, 370.370,
    
    ## delta = 15
    1.000, 1.000, 1.000, 3.606,
    28.620, 214.592, 370.370,
    
    ## delta = 20
    1.000, 1.000, 1.000, 1.555,
    13.488, 201.612, 370.370
  ),
  
  
  ## ----------------------------------------------------------
  ## TRUE ARL - n = 10
  ## ----------------------------------------------------------
  
  table6_true_ARL_n10 = c(
    
    ## delta = 0.5
    66.666, 109.409, 179.856, 255.102,
    349.650, 367.647, 370.370,
    
    ## delta = 3
    1.110, 2.174, 8.647, 52.137,
    129.199, 284.090, 370.370,
    
    ## delta = 5
    1.000, 1.034, 2.279, 19.305,
    82.508, 268.817, 370.370,
    
    ## delta = 10
    1.000, 1.000, 1.000, 2.900,
    21.376, 183.823, 370.370,
    
    ## delta = 15
    1.000, 1.000, 1.000, 1.260,
    9.281, 163.934, 370.370,
    
    ## delta = 20
    1.000, 1.000, 1.000, 1.009,
    4.013, 151.057, 370.370
  ),
  
  
  ## ----------------------------------------------------------
  ## TRUE ARL - n = 30
  ## ----------------------------------------------------------
  
  table6_true_ARL_n30 = c(
    
    ## delta = 0.5
    19.076, 43.554, 91.911, 198.412,
    236.966, 344.827, 370.370,
    
    ## delta = 3
    1.000, 1.006, 1.611, 12.764,
    64.350, 274.725, 370.370,
    
    ## delta = 5
    1.000, 1.000, 1.007, 3.395,
    25.188, 201.612, 370.370,
    
    ## delta = 10
    1.000, 1.000, 1.000, 1.033,
    4.179, 127.226, 370.370,
    
    ## delta = 15
    1.000, 1.000, 1.000, 1.000,
    1.513, 77.881, 370.370,
    
    ## delta = 20
    1.000, 1.000, 1.000, 1.000,
    1.057, 43.706, 370.370
  ),
  
  
  ## ==========================================================
  ## TRUE LCL from benchmark
  ##
  ## These are the standardized LCL values reported by
  ## Fernandes et al.
  ## ==========================================================
  
  table6_true_LCL_n5 = c(
    rep(-0.973, 7),   # delta = 0.5
    rep(-2.610, 7),   # delta = 3
    rep(-2.884, 7),   # delta = 5
    rep(-3.179, 7),   # delta = 10
    rep(-3.218, 7),   # delta = 15
    rep(-3.233, 7)    # delta = 20
  ),
  
  table6_true_LCL_n10 = c(
    rep(-1.277, 7),   # delta = 0.5
    rep(-2.690, 7),   # delta = 3
    rep(-2.866, 7),   # delta = 5
    rep(-2.953, 7),   # delta = 10
    rep(-3.119, 7),   # delta = 15
    rep(-3.158, 7)    # delta = 20
  ),
  
  table6_true_LCL_n30 = c(
    rep(-1.745, 7),   # delta = 0.5
    rep(-2.770, 7),   # delta = 3
    rep(-2.874, 7),   # delta = 5
    rep(-3.124, 7),   # delta = 10
    rep(-2.945, 7),   # delta = 15
    rep(-2.926, 7)    # delta = 20
  ),
  
  
  ## Asymptotic LCL from Table 6
  table6_false_LCL_n5  = rep(-2.780, 42),
  table6_false_LCL_n10 = rep(-2.780, 42),
  table6_false_LCL_n30 = rep(-2.780, 42),
  
  stringsAsFactors = FALSE
)


write.csv(
  table6_benchmark,
  "table6_benchmark_n5_n10_n30.csv",
  row.names = FALSE
)


## ============================================================
## Step 04 - Double-sampling parameter grid
## ============================================================

double_param_grid <- expand.grid(
  k0      = k0_vec,
  n_a     = n_a_vec,
  n_b     = n_b_vec,
  alpha_a = alpha_a_vec
)

double_param_grid$lambda0 <- (
  mu0 / gamma(1 + 1 / double_param_grid$k0)
)

double_param_grid$mean_weibull_0 <- mu0

double_param_grid$var_weibull_0 <- double_param_grid$lambda0^2 * (
  gamma(1 + 2 / double_param_grid$k0) -
    gamma(1 + 1 / double_param_grid$k0)^2
)

double_param_grid$alpha_b <- (
  alpha_global_h0 / double_param_grid$alpha_a
)

double_param_grid <- double_param_grid[
  is.finite(double_param_grid$alpha_b) &
    double_param_grid$alpha_b > 0 &
    double_param_grid$alpha_b < 1,
]

double_param_grid$Avg_Sample0_per_cycle <- round(
  double_param_grid$n_a +
    double_param_grid$n_b * double_param_grid$alpha_a,
  10
)

double_param_grid <- double_param_grid[
  double_param_grid$Avg_Sample0_per_cycle >= ASS_min &
    double_param_grid$Avg_Sample0_per_cycle <= ASS_max,
]

double_param_grid <- double_param_grid[
  sapply(
    double_param_grid$Avg_Sample0_per_cycle,
    function(x) any(abs(x - ASS0_allowed) < tol_ASS0)
  ),
]

double_param_grid$ARL0_global <- (
  1 / (double_param_grid$alpha_a * double_param_grid$alpha_b)
)

if (nrow(double_param_grid) == 0) {
  stop("No double-sampling configurations found with allowed ASS0 values.")
}


## ============================================================
## Step 05 - Double-sampling sweep
##
## IMPORTANT:
## Table 6 from benchmark is LOWER-SIDED only.
## Therefore ALWAYS use tcrit_low.
## ============================================================

double_performance_list <- list()

counter <- 1

for (g in seq_len(nrow(double_param_grid))) {
  
  k0_g      <- double_param_grid$k0[g]
  lambda0_g <- double_param_grid$lambda0[g]
  
  weibull_mean0_g <- double_param_grid$mean_weibull_0[g]
  weibull_var0_g  <- double_param_grid$var_weibull_0[g]
  
  n_a_g     <- double_param_grid$n_a[g]
  n_b_g     <- double_param_grid$n_b[g]
  alpha_a_g <- double_param_grid$alpha_a[g]
  alpha_b_g <- double_param_grid$alpha_b[g]
  
  Avg_Sample0_g <- double_param_grid$Avg_Sample0_per_cycle[g]
  
  
  design_a_g <- make_stage_design(
    stage    = "A_first_sample",
    n_sample = n_a_g,
    alpha    = alpha_a_g,
    k0       = k0_g,
    lambda0  = lambda0_g
  )
  
  design_b_g <- make_stage_design(
    stage    = "B_confirmation_sample",
    n_sample = n_b_g,
    alpha    = alpha_b_g,
    k0       = k0_g,
    lambda0  = lambda0_g
  )
  
  design_pairs_g <- make_design_pairs(
    design_a = design_a_g,
    design_b = design_b_g
  )
  
  
  for (s in seq_len(nrow(scenario_grid))) {
    
    d_s <- scenario_grid$d[s]
    
    k1_s <- k0_g
    
    
    ## ========================================================
    ## TABLE 6 SHIFT:
    ## mu1 = mu0 (1 + d), with d <= 0
    ## ========================================================
    
    mean_1_s <- weibull_mean0_g * (1 + d_s)
    
    lambda1_s <- (
      mean_1_s / gamma(1 + 1 / k1_s)
    )
    
    var_1_s <- lambda1_s^2 * (
      gamma(1 + 2 / k1_s) -
        gamma(1 + 1 / k1_s)^2
    )
    
    
    ## ========================================================
    ## TABLE 6 from benchmark IS ALWAYS LOWER-SIDED
    ##
    ## Even d = 0 must use the LOW convention so that
    ## H0 aligns with the LCL comparison in Table 6.
    ## ========================================================
    
    mean_shift_direction_s <- "decrease"
    chart_side_used_s      <- "low"
    
    
    tmp <- design_pairs_g
    
    
    ## ========================================================
    ## Probability that an observation falls BELOW the
    ## lower critical time under H1.
    ## ========================================================
    
    p_h1_a <- pweibull(
      tmp$tcrit_low_a,
      shape = k1_s,
      scale = lambda1_s,
      lower.tail = TRUE
    )
    
    p_h1_b <- pweibull(
      tmp$tcrit_low_b,
      shape = k1_s,
      scale = lambda1_s,
      lower.tail = TRUE
    )
    
    
    tcrit_a <- tmp$tcrit_low_a
    tcrit_b <- tmp$tcrit_low_b
    
    
    ## ========================================================
    ## Signal probability:
    ##
    ## N >= Y
    ##
    ## This remains the SAME binomial-tail logic used in
    ## the original code.
    ## ========================================================
    
    alpha1_a <- pbinom(
      tmp$Y_a - 1,
      size = n_a_g,
      prob = p_h1_a,
      lower.tail = FALSE
    )
    
    alpha1_b <- pbinom(
      tmp$Y_b - 1,
      size = n_b_g,
      prob = p_h1_b,
      lower.tail = FALSE
    )
    
    
    ## Two-stage signal requires A AND B
    alpha1_global_h1 <- alpha1_a * alpha1_b
    
    
    ARL1 <- ifelse(
      alpha1_global_h1 > 0,
      1 / alpha1_global_h1,
      Inf
    )
    
    
    Avg_Sample1_per_cycle <- (
      n_a_g +
        n_b_g * alpha1_a
    )
    
    
    Expected_Alert_Time_H1_per_cycle <- (
      tcrit_a +
        alpha1_a * tcrit_b
    )
    
    
    Expected_Alert_Time_To_Signal_H1 <- (
      ARL1 *
        Expected_Alert_Time_H1_per_cycle
    )
    
    
    ## --------------------------------------------------------
    ## Store results
    ## --------------------------------------------------------
    
    tmp$d                    <- d_s
    tmp$k0                   <- k0_g
    tmp$lambda0              <- lambda0_g
    tmp$k1                   <- k1_s
    tmp$lambda1              <- lambda1_s
    
    tmp$mean_weibull_0       <- weibull_mean0_g
    tmp$var_weibull_0        <- weibull_var0_g
    tmp$mean_weibull_1       <- mean_1_s
    tmp$var_weibull_1        <- var_1_s
    
    tmp$mean_shift_direction <- mean_shift_direction_s
    tmp$chart_side_used      <- chart_side_used_s
    
    tmp$n_a                  <- n_a_g
    tmp$n_b                  <- n_b_g
    
    tmp$alpha_a_h0           <- alpha_a_g
    tmp$alpha_b_h0           <- alpha_b_g
    tmp$alpha_global_h0      <- alpha_global_h0
    tmp$ARL0_global          <- 1 / alpha_global_h0
    
    tmp$Avg_Sample0_per_cycle <- Avg_Sample0_g
    
    tmp$tcrit_a              <- tcrit_a
    tmp$tcrit_b              <- tcrit_b
    
    tmp$p_h1_a               <- p_h1_a
    tmp$p_h1_b               <- p_h1_b
    
    tmp$alpha1_a             <- alpha1_a
    tmp$alpha1_b             <- alpha1_b
    tmp$alpha1_global_h1     <- alpha1_global_h1
    
    tmp$ARL1                 <- ARL1
    
    tmp$Avg_Sample1_per_cycle <- Avg_Sample1_per_cycle
    
    tmp$Expected_Alert_Time_To_Signal_H1 <- (
      Expected_Alert_Time_To_Signal_H1
    )
    
    
    double_performance_list[[counter]] <- tmp
    
    counter <- counter + 1
  }
}


double_performance_table <- do.call(
  rbind,
  double_performance_list
)


## ============================================================
## Step 06 - Mapped scenarios
## TABLE 6 benchmark LCL comparison
## ============================================================

mapped_scenarios <- double_performance_table[, c(
  
  "n_a",
  "n_b",
  "Y_a",
  "Y_b",
  
  "alpha_a_h0",
  "alpha_b_h0",
  "Avg_Sample0_per_cycle",
  
  "p1_a",
  "p1_b",
  
  "k0",
  "lambda0",
  "lambda1",
  "d",
  
  "mean_weibull_0",
  "var_weibull_0",
  "mean_weibull_1",
  "var_weibull_1",
  
  "mean_shift_direction",
  "chart_side_used",
  
  "tcrit_low_a",
  "tcrit_low_b",
  "tcrit_high_a",
  "tcrit_high_b",
  
  "alpha1_a",
  "alpha1_b",
  "alpha1_global_h1",
  
  "ARL1"
)]


## ============================================================
## Attach Table 6 benchmark
## ============================================================

mapped_scenarios <- merge(
  mapped_scenarios,
  table6_benchmark,
  by = c("k0", "d"),
  all.x = TRUE
)


mapped_scenarios$k <- mapped_scenarios$k0

mapped_scenarios$ASS0 <- (
  mapped_scenarios$Avg_Sample0_per_cycle
)


## ============================================================
## ARL1 comparisons
## ============================================================

mapped_scenarios$ARL1_improvement_vs_table6_true_n5 <- (
  mapped_scenarios$table6_true_ARL_n5 -
    mapped_scenarios$ARL1
)

mapped_scenarios$ARL1_ratio_vs_table6_true_n5 <- (
  mapped_scenarios$ARL1 /
    mapped_scenarios$table6_true_ARL_n5
)

mapped_scenarios$better_than_table6_true_n5 <- (
  mapped_scenarios$ARL1 <=
    mapped_scenarios$table6_true_ARL_n5
)


mapped_scenarios$ARL1_improvement_vs_table6_true_n10 <- (
  mapped_scenarios$table6_true_ARL_n10 -
    mapped_scenarios$ARL1
)

mapped_scenarios$ARL1_ratio_vs_table6_true_n10 <- (
  mapped_scenarios$ARL1 /
    mapped_scenarios$table6_true_ARL_n10
)

mapped_scenarios$better_than_table6_true_n10 <- (
  mapped_scenarios$ARL1 <=
    mapped_scenarios$table6_true_ARL_n10
)


mapped_scenarios$ARL1_improvement_vs_table6_true_n30 <- (
  mapped_scenarios$table6_true_ARL_n30 -
    mapped_scenarios$ARL1
)

mapped_scenarios$ARL1_ratio_vs_table6_true_n30 <- (
  mapped_scenarios$ARL1 /
    mapped_scenarios$table6_true_ARL_n30
)

mapped_scenarios$better_than_table6_true_n30 <- (
  mapped_scenarios$ARL1 <=
    mapped_scenarios$table6_true_ARL_n30
)


## ============================================================
## Final selected columns
## ============================================================

mapped_scenarios <- mapped_scenarios[, c(
  
  "n_a",
  "n_b",
  "Y_a",
  "Y_b",
  
  "alpha_a_h0",
  "alpha_b_h0",
  
  "ASS0",
  
  "p1_a",
  "p1_b",
  
  "k",
  "lambda0",
  "lambda1",
  "d",
  
  "mean_weibull_0",
  "var_weibull_0",
  "mean_weibull_1",
  "var_weibull_1",
  
  "mean_shift_direction",
  "chart_side_used",
  
  ## np-Weibull lower limits
  "tcrit_low_a",
  "tcrit_low_b",
  
  ## Kept for reference
  #"tcrit_high_a",
  #"tcrit_high_b",
  
  "alpha1_a",
  "alpha1_b",
  "alpha1_global_h1",
  
  ## Table 6 TRUE ARL benchmark
  "table6_true_ARL_n5",
  "table6_true_ARL_n10",
  "table6_true_ARL_n30",
  
  ## Table 6 TRUE LCL benchmark
  "table6_true_LCL_n5",
  "table6_true_LCL_n10",
  "table6_true_LCL_n30",
  
  ## Table 6 asymptotic LCL
  "table6_false_LCL_n5",
  "table6_false_LCL_n10",
  "table6_false_LCL_n30",
  
  ## Proposed chart
  "ARL1"
)]


## ============================================================
## Sort
## ============================================================

mapped_scenarios <- mapped_scenarios[order(
  mapped_scenarios$ASS0,
  mapped_scenarios$k,
  mapped_scenarios$d,
  mapped_scenarios$ARL1,
  mapped_scenarios$n_a,
  mapped_scenarios$n_b,
  mapped_scenarios$Y_a,
  mapped_scenarios$Y_b
), ]


## ============================================================
## Design ID
## ============================================================

mapped_scenarios$design_id <- paste0(
  
  "ASS",
  gsub(
    "\\.",
    "p",
    formatC(
      mapped_scenarios$ASS0,
      format = "f",
      digits = 2
    )
  ),
  
  "_na",
  mapped_scenarios$n_a,
  
  "_nb",
  mapped_scenarios$n_b
)


## ============================================================
## Export
## ============================================================

write.csv(
  mapped_scenarios,
  "np_weibull_DS_mapped_scenarios_table6_n5_n10_n30.csv",
  row.names = FALSE
)


## ============================================================
## Step 07 - Final checks
## ============================================================

cat("\n============================================================\n")
cat("FINAL CHECKS\n")
cat("============================================================\n")

cat("Table 6 comparison: DECREASES in mean\n")
cat("Control side: LOWER / LCL\n")
cat("Table 6 benchmarks included: n = 5, 10 and 30\n")

cat("alpha_global_h0:", alpha_global_h0, "\n")
cat("ARL0 target:", ARL0_target, "\n")

cat(
  "ASS0 restriction:",
  ASS_min,
  "<= ASS0 <=",
  ASS_max,
  "\n"
)

cat(
  "ASS0 allowed values:",
  paste(ASS0_allowed, collapse = ", "),
  "\n"
)

cat("n_b constraint: n_b <=", n_b_max, "\n")

cat(
  "Double configurations:",
  nrow(double_param_grid),
  "\n"
)

cat(
  "Unique ASS0 values:",
  length(unique(
    double_param_grid$Avg_Sample0_per_cycle
  )),
  "\n"
)

cat(
  "ASS0 values found:",
  paste(
    sort(unique(
      double_param_grid$Avg_Sample0_per_cycle
    )),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Scenarios:",
  nrow(scenario_grid),
  "\n"
)

cat(
  "Mapped scenario rows:",
  nrow(mapped_scenarios),
  "\n"
)

cat("\nFiles written:\n")

cat(
  "- table6_benchmark_n5_n10_n30.csv\n"
)

cat(
  "- np_weibull_DS_mapped_scenarios_table6_n5_n10_n30.csv\n"
)

