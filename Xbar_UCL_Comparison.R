## ============================================================
## npˆs2_wb comparison based on UCL
## ============================================================

options(digits = 16)

## ============================================================
## Step 00 - Global parameters
## ============================================================

alpha_global_h0 <- 0.0027
ARL0_target <- 1 / alpha_global_h0

ASS_min <- 3
ASS_max <- 4

#ASS0_allowed <- seq(ASS_min, ASS_max, by = 0.25)
ASS0_allowed <- c(3.0, 3.10, 3.20)

tol_ASS0 <- 1e-10

n_a_vec <- 1:ASS_max

n_b_max <- 20
n_b_vec <- 1:n_b_max

alpha_a_vec <- c(0.1)

## Table 5 shape parameters
k0_vec <- c(0.5, 3, 5, 10, 15, 20)

## Scale-invariant reference mean for proportional shifts
mu0 <- 1

## Table 5 increases in mean
d_vec <- c(0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40)


cat("\n============================================================\n")
cat("GLOBAL SETTINGS - TABLE 5 BENCHMARK COMPARISON\n")
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
cat("Table 5 benchmark n values: 5 and 10\n")


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


get_mean_shift_direction <- function(mean_1, mean_0, tol = 1e-12) {
  
  direction <- rep("same", length(mean_1))
  direction[mean_1 < mean_0 - tol] <- "decrease"
  direction[mean_1 > mean_0 + tol] <- "increase"
  
  direction
}


make_stage_design <- function(stage, n_sample, alpha, k0, lambda0) {
  
  Y_vec <- seq_len(n_sample)
  
  p1_vec <- find_p1(
    Y        = Y_vec,
    n_sample = n_sample,
    alpha    = alpha
  )
  
  tcrit_low_vec <- qweibull(
    p1_vec,
    shape = k0,
    scale = lambda0
  )
  
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
## Increases in mean
## ============================================================

scenario_grid <- data.frame(
  d = d_vec
)

scenario_grid$scenario_id <- paste0(
  "d",
  gsub("\\.", "p", as.character(scenario_grid$d))
)


## ============================================================
## Step 03 - Benchmark values for n = 5 and n = 10
## TRUE and FALSE ARL values from Table 5
## ============================================================

table5_benchmark <- data.frame(
  k0 = c(
    rep(0.5, 7),
    rep(3, 7),
    rep(5, 7),
    rep(10, 7),
    rep(15, 7),
    rep(20, 7)
  ),
  d = rep(d_vec, times = 6),
  
  table5_true_ARL_n5 = c(
    ## delta = 0.5
    370.370, 373.134, 284.090, 240.384, 162.337, 112.107, 94.161,
    ## delta = 3
    370.370, 295.858, 124.069, 43.898, 11.030, 4.552, 2.605,
    ## delta = 5
    370.370, 207.468, 50.607, 14.667, 3.286, 1.640, 1.209,
    ## delta = 10
    370.370, 142.857, 15.133, 3.242, 1.169, 1.013, 1.000,
    ## delta = 15
    370.370, 107.526, 5.946, 1.547, 1.010, 1.000, 1.000,
    ## delta = 20
    370.370, 88.652, 3.239, 1.153, 1.000, 1.000, 1.000
  ),
  
  table5_true_ARL_n10 = c(
    ## delta = 0.5
    370.370, 322.580, 276.243, 210.084, 140.449, 99.601, 65.189,
    ## delta = 3
    370.370, 271.739, 84.745, 25.654, 5.510, 2.334, 1.484,
    ## delta = 5
    370.370, 235.849, 38.197, 8.748, 1.888, 1.159, 1.024,
    ## delta = 10
    370.370, 112.359, 8.029, 1.795, 1.014, 1.000, 1.000,
    ## delta = 15
    370.370, 97.276, 3.411, 1.128, 1.000, 1.000, 1.000,
    ## delta = 20
    370.370, 55.248, 1.817, 1.013, 1.000, 1.000, 1.000
  ),
  
  table5_true_ARL_n30 = c(
    ## delta = 0.5
    370.370, 304.878, 203.252, 178.571, 92.764, 51.493, 11.978,
    ## delta = 3
    370.370, 233.644, 37.850, 8.418, 1.760, 1.104, 1.010,
    ## delta = 5
    370.370, 225.225, 15.276, 2.765, 1.052, 1.000, 1.000,
    ## delta = 10
    370.370, 79.617, 2.679, 1.044, 1.000, 1.000, 1.000,
    ## delta = 15
    370.370, 47.619, 1.322, 1.000, 1.000, 1.000, 1.000,
    ## delta = 20
    370.370, 25.150, 1.049, 1.000, 1.000, 1.000, 1.000
  ),
  
  table5_false_ARL_n5 = c(
    ## delta = 0.5
    45.085, 42.194, 37.037, 34.722, 25.667, 19.976, 16.778,
    ## delta = 3
    295.858, 255.102, 103.092, 38.941, 10.042, 4.291, 2.496,
    ## delta = 5
    833.333, 384.615, 81.566, 21.026, 4.004, 1.797, 1.264,
    ## delta = 10
    1562.500, 574.712, 34.940, 4.926, 1.269, 1.023, 1.001,
    ## delta = 15
    2941.176, 500.000, 12.581, 1.974, 1.021, 1.000, 1.000,
    ## delta = 20
    2272.727, 438.596, 5.616, 1.276, 1.000, 1.000, 1.000
  ),
  
  table5_false_ARL_n10 = c(
    ## delta = 0.5
    51.124, 45.620, 38.819, 33.355, 24.582, 18.789, 14.216,
    ## delta = 3
    306.748, 236.966, 75.414, 23.430, 5.189, 2.255, 1.455,
    ## delta = 5
    515.463, 324.675, 46.598, 10.014, 1.995, 1.177, 1.027,
    ## delta = 10
    1136.364, 289.017, 13.401, 2.188, 1.024, 1.000, 1.000,
    ## delta = 15
    1063.830, 248.756, 4.737, 1.189, 1.000, 1.000, 1.000,
    ## delta = 20
    1515.152, 148.368, 2.316, 1.024, 1.000, 1.000, 1.000
  ),
  
  table5_true_UCL_n5 = c(
    rep(5.632, 7),
    rep(2.841, 7),
    rep(2.594, 7),
    rep(2.438, 7),
    rep(2.379, 7),
    rep(2.379, 7)
  ),
  
  table5_true_UCL_n10 = c(
    rep(4.965, 7),
    rep(2.831, 7),
    rep(2.697, 7),
    rep(2.507, 7),
    rep(2.527, 7),
    rep(2.466, 7)
  ),
  
  table5_true_UCL_n30 = c(
    rep(4.128, 7),
    rep(2.810, 7),
    rep(2.780, 7),
    rep(2.646, 7),
    rep(2.635, 7),
    rep(2.602, 7)
  ),
  
  table5_false_UCL_n30 = rep(2.780, 42),
  table5_false_UCL_n5 = rep(2.780, 42),
  table5_false_UCL_n10 = rep(2.780, 42),
  
  stringsAsFactors = FALSE
)

write.csv(
  table5_benchmark,
  "table5_benchmark_n5_n10.csv",
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
## ============================================================

double_performance_list <- list()

counter <- 1

for (g in seq_len(nrow(double_param_grid))) {
  
  k0_g      <- double_param_grid$k0[g]
  lambda0_g <- double_param_grid$lambda0[g]
  
  weibull_mean0_g <- double_param_grid$mean_weibull_0[g]
  weibull_var0_g  <- double_param_grid$var_weibull_0[g]
  
  tol_mean_g <- 1e-12 * max(1, abs(weibull_mean0_g))
  
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
    
    mean_1_s <- weibull_mean0_g * (1 + d_s)
    
    lambda1_s <- (
      mean_1_s / gamma(1 + 1 / k1_s)
    )
    
    var_1_s <- lambda1_s^2 * (
      gamma(1 + 2 / k1_s) -
        gamma(1 + 1 / k1_s)^2
    )
    
    mean_shift_direction_s <- get_mean_shift_direction(
      mean_1 = mean_1_s,
      mean_0 = weibull_mean0_g,
      tol    = tol_mean_g
    )
    
    ## Table 5 is UCL / increase-in-mean.
    ## For d = 0, force high-side convention so the row aligns with Table 5 UCL.
    if (d_s == 0) {
      mean_shift_direction_s <- "increase"
    }
    
    tmp <- design_pairs_g
    
    if (mean_shift_direction_s == "increase") {
      
      chart_side_used_s <- "high"
      
      p_h1_a <- pweibull(
        tmp$tcrit_high_a,
        shape = k1_s,
        scale = lambda1_s,
        lower.tail = FALSE
      )
      
      p_h1_b <- pweibull(
        tmp$tcrit_high_b,
        shape = k1_s,
        scale = lambda1_s,
        lower.tail = FALSE
      )
      
      tcrit_a <- tmp$tcrit_high_a
      tcrit_b <- tmp$tcrit_high_b
      
    } else {
      
      chart_side_used_s <- "low"
      
      p_h1_a <- pweibull(
        tmp$tcrit_low_a,
        shape = k1_s,
        scale = lambda1_s
      )
      
      p_h1_b <- pweibull(
        tmp$tcrit_low_b,
        shape = k1_s,
        scale = lambda1_s
      )
      
      tcrit_a <- tmp$tcrit_low_a
      tcrit_b <- tmp$tcrit_low_b
    }
    
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
    
    alpha1_global_h1 <- alpha1_a * alpha1_b
    
    ARL1 <- ifelse(
      alpha1_global_h1 > 0,
      1 / alpha1_global_h1,
      Inf
    )
    
    Avg_Sample1_per_cycle <- n_a_g + n_b_g * alpha1_a
    
    Expected_Alert_Time_H1_per_cycle <- (
      tcrit_a + alpha1_a * tcrit_b
    )
    
    Expected_Alert_Time_To_Signal_H1 <- (
      ARL1 * Expected_Alert_Time_H1_per_cycle
    )
    
    tmp$d                     <- d_s
    tmp$k0                    <- k0_g
    tmp$lambda0               <- lambda0_g
    tmp$k1                    <- k1_s
    tmp$lambda1               <- lambda1_s
    tmp$mean_weibull_0        <- weibull_mean0_g
    tmp$var_weibull_0         <- weibull_var0_g
    tmp$mean_weibull_1        <- mean_1_s
    tmp$var_weibull_1         <- var_1_s
    tmp$mean_shift_direction  <- mean_shift_direction_s
    tmp$chart_side_used       <- chart_side_used_s
    
    tmp$n_a                   <- n_a_g
    tmp$n_b                   <- n_b_g
    tmp$alpha_a_h0            <- alpha_a_g
    tmp$alpha_b_h0            <- alpha_b_g
    tmp$alpha_global_h0       <- alpha_global_h0
    tmp$ARL0_global           <- 1 / alpha_global_h0
    
    tmp$Avg_Sample0_per_cycle <- Avg_Sample0_g
    
    tmp$tcrit_a               <- tcrit_a
    tmp$tcrit_b               <- tcrit_b
    
    tmp$Expected_Alert_Time_To_Signal_H1 <- Expected_Alert_Time_To_Signal_H1
    
    tmp$p_h1_a                <- p_h1_a
    tmp$p_h1_b                <- p_h1_b
    tmp$alpha1_a              <- alpha1_a
    tmp$alpha1_b              <- alpha1_b
    tmp$alpha1_global_h1      <- alpha1_global_h1
    tmp$ARL1                  <- ARL1
    tmp$Avg_Sample1_per_cycle <- Avg_Sample1_per_cycle
    
    tmp$ASS0_in_range <- (
      tmp$Avg_Sample0_per_cycle >= ASS_min &
        tmp$Avg_Sample0_per_cycle <= ASS_max
    )
    
    tmp$ASS1_in_range <- (
      tmp$Avg_Sample1_per_cycle >= ASS_min &
        tmp$Avg_Sample1_per_cycle <= ASS_max
    )
    
    double_performance_list[[counter]] <- tmp
    counter <- counter + 1
  }
}

double_performance_table <- do.call(rbind, double_performance_list)


## ============================================================
## Step 06 - Mapped scenarios only
## Table 5 benchmark comparison, benchmark n = 5 and n = 10
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
  
  #"tcrit_a",
  #"tcrit_b",
  
  "alpha1_a",
  "alpha1_b",
  "alpha1_global_h1",
  
  #"Expected_Alert_Time_To_Signal_H1",
  
  "ARL1"
)]

mapped_scenarios <- merge(
  mapped_scenarios,
  table5_benchmark,
  by = c("k0", "d"),
  all.x = TRUE
)

mapped_scenarios$k <- mapped_scenarios$k0

mapped_scenarios$ASS0 <- mapped_scenarios$Avg_Sample0_per_cycle

mapped_scenarios$ARL1_improvement_vs_table5_true_n5 <- (
  mapped_scenarios$table5_true_ARL_n5 - mapped_scenarios$ARL1
)

mapped_scenarios$ARL1_ratio_vs_table5_true_n5 <- (
  mapped_scenarios$ARL1 / mapped_scenarios$table5_true_ARL_n5
)

mapped_scenarios$better_than_table5_true_n5 <- (
  mapped_scenarios$ARL1 <= mapped_scenarios$table5_true_ARL_n5
)

mapped_scenarios$ARL1_improvement_vs_table5_true_n10 <- (
  mapped_scenarios$table5_true_ARL_n10 - mapped_scenarios$ARL1
)

mapped_scenarios$ARL1_ratio_vs_table5_true_n10 <- (
  mapped_scenarios$ARL1 / mapped_scenarios$table5_true_ARL_n10
)

mapped_scenarios$better_than_table5_true_n10 <- (
  mapped_scenarios$ARL1 <= mapped_scenarios$table5_true_ARL_n10
)

mapped_scenarios$ARL1_improvement_vs_table5_true_n30 <- (
  mapped_scenarios$table5_true_ARL_n30 - mapped_scenarios$ARL1
)

mapped_scenarios$ARL1_ratio_vs_table5_true_n30 <- (
  mapped_scenarios$ARL1 / mapped_scenarios$table5_true_ARL_n30
)

mapped_scenarios$better_than_table5_true_n30 <- (
  mapped_scenarios$ARL1 <= mapped_scenarios$table5_true_ARL_n30
)

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
  
  "tcrit_low_a",
  "tcrit_low_b",
  "tcrit_high_a",
  "tcrit_high_b",
  
  #"tcrit_a",
  #"tcrit_b",
  
  "alpha1_a",
  "alpha1_b",
  "alpha1_global_h1",
  
  #"Expected_Alert_Time_To_Signal_H1",
  
  "table5_true_ARL_n5",
  "table5_true_ARL_n10",
  "table5_true_ARL_n30",
  #"table5_false_ARL_n5",
  #"table5_false_ARL_n10",
  "table5_true_UCL_n5",
  "table5_true_UCL_n10",
  "table5_true_UCL_n30",
  #"table5_false_UCL_n5",
  #"table5_false_UCL_n10",
  
  #"ARL1_improvement_vs_table5_true_n5",
  #"ARL1_ratio_vs_table5_true_n5",
  #"better_than_table5_true_n5",
  
  #"ARL1_improvement_vs_table5_true_n10",
  #"ARL1_ratio_vs_table5_true_n10",
  #"better_than_table5_true_n10",
  
  "ARL1"
)]

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

mapped_scenarios$design_id <- paste0(
  "ASS", gsub("\\.", "p", formatC(mapped_scenarios$ASS0, format = "f", digits = 2)),
  "_na", mapped_scenarios$n_a,
  "_nb", mapped_scenarios$n_b
)

write.csv(
  mapped_scenarios,
  "np_weibull_DS_mapped_scenarios_table5_n5_n10.csv",
  row.names = FALSE
)


## ============================================================
## Step 07 - Final checks
## ============================================================

cat("\n============================================================\n")
cat("FINAL CHECKS\n")
cat("============================================================\n")

cat("Table 5 comparison: increases in mean\n")
cat("Table 5 benchmarks included: n = 5 and n = 10\n")
cat("alpha_global_h0:", alpha_global_h0, "\n")
cat("ARL0 target:", ARL0_target, "\n")
cat("ASS0 restriction:", ASS_min, "<= ASS0 <=", ASS_max, "\n")
cat("ASS0 allowed values:", paste(ASS0_allowed, collapse = ", "), "\n")
cat("n_b constraint: n_b <=", n_b_max, "\n")
cat("Double configurations:", nrow(double_param_grid), "\n")
cat("Unique ASS0 values:", length(unique(double_param_grid$Avg_Sample0_per_cycle)), "\n")
cat("ASS0 values found:", paste(sort(unique(double_param_grid$Avg_Sample0_per_cycle)), collapse = ", "), "\n")
cat("Scenarios:", nrow(scenario_grid), "\n")
cat("Mapped scenario rows:", nrow(mapped_scenarios), "\n")

cat("\nFiles written:\n")
cat("- table5_benchmark_n5_n10.csv\n")
cat("- np_weibull_DS_mapped_scenarios_table5_n5_n10.csv\n")

