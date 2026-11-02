## ============================================================
## SPECIAL CASE - SAMPLE A NOT COLLECTED
##
## UCL Evaluation
## ASS values investigated:
## 3, 5, 7, 10, 15, 20
##
## ============================================================

options(digits = 16)

## ============================================================
## Step 00 - Global parameters
## ============================================================

alpha_global_h0 <- 0.0027

ARL0_target <- 1 / alpha_global_h0


## ------------------------------------------------------------
## Special case:
## A is not collected.
##
## Therefore:
## ASS = n_b
## ------------------------------------------------------------

ASS0_allowed <- c(
  3,
  5,
  7,
  10,
  15,
  20
)

n_b_vec <- ASS0_allowed


## ----------------------------------------------------------------
## Shape parameters used in Table 6 from Fernandes et Al benchmark
## ----------------------------------------------------------------

k0_vec <- c(
  0.5,
  3,
  5,
  10,
  15,
  20
)


## ------------------------------------------------------------
## Scale-invariant reference mean
## ------------------------------------------------------------

mu0 <- 1


## ------------------------------------------------------------
## Table 6:
## decreases in mean
## ------------------------------------------------------------

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
cat("GLOBAL SETTINGS - SPECIAL CASE: SAMPLE A NOT COLLECTED\n")
cat("============================================================\n")

cat("Chart side             = LCL / LOW\n")

cat(
  "alpha_global_h0       =",
  alpha_global_h0,
  "\n"
)

cat(
  "ARL0_target           =",
  ARL0_target,
  "\n"
)

cat(
  "Sample A collected?   = NO\n"
)

cat(
  "Sample B mandatory?   = YES\n"
)

cat(
  "alpha_b special       =",
  alpha_global_h0,
  "\n"
)

cat(
  "ASS values            =",
  paste(ASS0_allowed, collapse = ", "),
  "\n"
)

cat(
  "k0 values             =",
  paste(k0_vec, collapse = ", "),
  "\n"
)

cat(
  "d values              =",
  paste(d_vec, collapse = ", "),
  "\n"
)


## ============================================================
## Step 01 - Helper functions
## ============================================================

.validate_inputs <- function(
    Y,
    n_sample,
    alpha
) {
  
  if (
    length(n_sample) != 1 ||
    !is.finite(n_sample) ||
    n_sample <= 0 ||
    n_sample != as.integer(n_sample)
  ) {
    
    stop(
      "n_sample must be a positive integer."
    )
    
  }
  
  
  if (
    length(alpha) != 1 ||
    !is.finite(alpha) ||
    alpha <= 0 ||
    alpha >= 1
  ) {
    
    stop(
      "alpha must be in (0,1)."
    )
    
  }
  
  
  if (
    any(!is.finite(Y)) ||
    any(Y != as.integer(Y)) ||
    any(Y < 1) ||
    any(Y > n_sample)
  ) {
    
    stop(
      "Y must contain integers in [1, n_sample]."
    )
    
  }
  
  
  invisible(TRUE)
}


## ------------------------------------------------------------
## Find p1 such that
##
## P(Binomial(n,p1) >= Y) = alpha
## ------------------------------------------------------------

find_p1 <- function(
    Y,
    n_sample,
    alpha
) {
  
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


## ------------------------------------------------------------
## Build B-stage design
##
## Both low and high thresholds are calculated,
## but ONLY tcrit_low is used for Table 6.
## ------------------------------------------------------------

make_stage_design <- function(
    stage,
    n_sample,
    alpha,
    k0,
    lambda0
) {
  
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
  
  
  ## ----------------------------------------------------------
  ## H0 check
  ##
  ## For LCL:
  ## p1 = P(T <= tcrit_low | H0)
  ##
  ## Signal if N >= Y.
  ## ----------------------------------------------------------
  
  alpha_h0_check <- pbinom(
    Y_vec - 1,
    size = n_sample,
    prob = p1_vec,
    lower.tail = FALSE
  )
  
  
  ARL0_check <- ifelse(
    alpha_h0_check > 0,
    1 / alpha_h0_check,
    Inf
  )
  
  
  data.frame(
    
    stage          = stage,
    
    n_sample       = n_sample,
    
    alpha_stage    = alpha,
    
    Y              = Y_vec,
    
    p1             = p1_vec,
    
    tcrit_low      = tcrit_low_vec,
    
    tcrit_high     = tcrit_high_vec,
    
    alpha_h0_check = alpha_h0_check,
    
    ARL0_check     = ARL0_check,
    
    stringsAsFactors = FALSE
  )
}


## ============================================================
## Step 02 - Out-of-control scenarios
##
## Table 6 = decreases in mean
## ============================================================

scenario_grid <- data.frame(
  d = d_vec
)


scenario_grid$scenario_id <- paste0(
  
  "d",
  
  gsub(
    "-",
    "m",
    gsub(
      "\\.",
      "p",
      as.character(
        scenario_grid$d
      )
    )
  )
)


## ============================================================
## Step 03 - Table 6 benchmark
##
## TRUE ARL values and TRUE LCL values
##
## n = 3, 5, 10 and 30
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
  
  
  d = rep(
    d_vec,
    times = 6
  ),
  
  
  ## ==========================================================
  ## TRUE ARL - n = 3
  ## ==========================================================
  
  table6_true_ARL_n3 = c(
    
    ## k = 0.5
    179.856,
    226.244,
    267.379,
    274.725,
    373.134,
    378.787,
    370.370,
    
    ## k = 3
    7.717,
    21.394,
    56.561,
    136.239,
    214.592,
    352.112,
    370.370,
    
    ## k = 5
    1.649,
    4.941,
    19.186,
    83.472,
    166.666,
    318.471,
    370.370,
    
    ## k = 10
    1.000,
    1.037,
    2.912,
    25.720,
    93.457,
    245.098,
    370.370,
    
    ## k = 15
    1.000,
    1.000,
    1.131,
    10.307,
    55.679,
    287.356,
    370.370,
    
    ## k = 20
    1.000,
    1.000,
    1.000,
    3.757,
    29.154,
    187.969,
    370.370
  ),
  
  
  ## ==========================================================
  ## TRUE ARL - n = 5
  ## ==========================================================
  
  table6_true_ARL_n5 = c(
    
    ## k = 0.5
    126.582,
    178.571,
    210.970,
    310.559,
    326.797,
    362.318,
    370.370,
    
    ## k = 3
    2.458,
    7.074,
    25.667,
    89.766,
    153.374,
    251.256,
    370.370,
    
    ## k = 5
    1.042,
    1.941,
    8.007,
    52.083,
    144.508,
    299.401,
    370.370,
    
    ## k = 10
    1.000,
    1.000,
    1.351,
    12.431,
    72.568,
    280.898,
    370.370,
    
    ## k = 15
    1.000,
    1.000,
    1.000,
    3.606,
    28.620,
    214.592,
    370.370,
    
    ## k = 20
    1.000,
    1.000,
    1.000,
    1.555,
    13.488,
    201.612,
    370.370
  ),
  
  
  ## ==========================================================
  ## TRUE ARL - n = 10
  ## ==========================================================
  
  table6_true_ARL_n10 = c(
    
    ## k = 0.5
    66.666,
    109.409,
    179.856,
    255.102,
    349.650,
    367.647,
    370.370,
    
    ## k = 3
    1.110,
    2.174,
    8.647,
    52.137,
    129.199,
    284.090,
    370.370,
    
    ## k = 5
    1.000,
    1.034,
    2.279,
    19.305,
    82.508,
    268.817,
    370.370,
    
    ## k = 10
    1.000,
    1.000,
    1.000,
    2.900,
    21.376,
    183.823,
    370.370,
    
    ## k = 15
    1.000,
    1.000,
    1.000,
    1.260,
    9.281,
    163.934,
    370.370,
    
    ## k = 20
    1.000,
    1.000,
    1.000,
    1.009,
    4.013,
    151.057,
    370.370
  ),
  
  
  ## ==========================================================
  ## TRUE ARL - n = 30
  ## ==========================================================
  
  table6_true_ARL_n30 = c(
    
    ## k = 0.5
    19.076,
    43.554,
    91.911,
    198.412,
    236.966,
    344.827,
    370.370,
    
    ## k = 3
    1.000,
    1.006,
    1.611,
    12.764,
    64.350,
    274.725,
    370.370,
    
    ## k = 5
    1.000,
    1.000,
    1.007,
    3.395,
    25.188,
    201.612,
    370.370,
    
    ## k = 10
    1.000,
    1.000,
    1.000,
    1.033,
    4.179,
    127.226,
    370.370,
    
    ## k = 15
    1.000,
    1.000,
    1.000,
    1.000,
    1.513,
    77.881,
    370.370,
    
    ## k = 20
    1.000,
    1.000,
    1.000,
    1.000,
    1.057,
    43.706,
    370.370
  ),
  
  
  ## ==========================================================
  ## TRUE LCL - Table 6
  ## ==========================================================
  
  table6_true_LCL_n3 = c(
    
    rep(-0.770, 7),
    rep(-2.582, 7),
    rep(-2.848, 7),
    rep(-3.155, 7),
    rep(-3.317, 7),
    rep(-3.304, 7)
  ),
  
  
  table6_true_LCL_n5 = c(
    
    rep(-0.973, 7),
    rep(-2.610, 7),
    rep(-2.884, 7),
    rep(-3.179, 7),
    rep(-3.218, 7),
    rep(-3.233, 7)
  ),
  
  
  table6_true_LCL_n10 = c(
    
    rep(-1.277, 7),
    rep(-2.690, 7),
    rep(-2.866, 7),
    rep(-2.953, 7),
    rep(-3.119, 7),
    rep(-3.158, 7)
  ),
  
  
  table6_true_LCL_n30 = c(
    
    rep(-1.745, 7),
    rep(-2.770, 7),
    rep(-2.874, 7),
    rep(-3.124, 7),
    rep(-2.945, 7),
    rep(-2.926, 7)
  ),
  
  
  ## ----------------------------------------------------------
  ## Asymptotic LCL from Table 6
  ## ----------------------------------------------------------
  
  table6_false_LCL_n3  = rep(-2.780, 42),
  table6_false_LCL_n5  = rep(-2.780, 42),
  table6_false_LCL_n10 = rep(-2.780, 42),
  table6_false_LCL_n30 = rep(-2.780, 42),
  
  
  stringsAsFactors = FALSE
)


write.csv(
  
  table6_benchmark,
  
  "table6_benchmark_LCL_n3_n5_n10_n30.csv",
  
  row.names = FALSE
)


## ============================================================
## Step 04 - Special-case parameter grid
##
## A IS NOT COLLECTED
##
## n_a = 0
##
## B is mandatory:
##
## alpha_b = alpha_global_h0 = 0.0027
##
## ASS0 = n_b
## ============================================================

special_param_grid <- expand.grid(
  
  k0  = k0_vec,
  
  n_b = n_b_vec
)


special_param_grid$n_a <- 0


special_param_grid$lambda0 <- (
  
  mu0 /
    
    gamma(
      1 + 1 / special_param_grid$k0
    )
)


special_param_grid$mean_weibull_0 <- mu0


special_param_grid$var_weibull_0 <- (
  
  special_param_grid$lambda0^2 *
    
    (
      gamma(
        1 + 2 / special_param_grid$k0
      ) -
        
        gamma(
          1 + 1 / special_param_grid$k0
        )^2
    )
)


## ------------------------------------------------------------
## B is now the only decision stage.
##
## Therefore:
##
## alpha_b = alpha_global
## ------------------------------------------------------------

special_param_grid$alpha_b <- (
  alpha_global_h0
)


special_param_grid$alpha_global_h0 <- (
  alpha_global_h0
)


special_param_grid$ARL0_global <- (
  1 / alpha_global_h0
)


## ------------------------------------------------------------
## Since A was not collected and B is mandatory:
##
## ASS0 = n_b
## ------------------------------------------------------------

special_param_grid$Avg_Sample0_per_cycle <- (
  special_param_grid$n_b
)


special_param_grid$ASS0 <- (
  special_param_grid$n_b
)


## ============================================================
## Step 05 - Sweep
##
## SPECIAL CASE:
## A missing
## B mandatory
##
## LCL ONLY
## ============================================================

performance_list <- list()

counter <- 1


for (
  g in seq_len(
    nrow(special_param_grid)
  )
) {
  
  
  ## ----------------------------------------------------------
  ## Current design
  ## ----------------------------------------------------------
  
  k0_g <- (
    special_param_grid$k0[g]
  )
  
  
  lambda0_g <- (
    special_param_grid$lambda0[g]
  )
  
  
  weibull_mean0_g <- (
    special_param_grid$mean_weibull_0[g]
  )
  
  
  weibull_var0_g <- (
    special_param_grid$var_weibull_0[g]
  )
  
  
  n_b_g <- (
    special_param_grid$n_b[g]
  )
  
  
  ## ----------------------------------------------------------
  ## B is mandatory and alone preserves ARL0.
  ## ----------------------------------------------------------
  
  alpha_b_g <- (
    alpha_global_h0
  )
  
  
  ASS0_g <- n_b_g
  
  
  ## ----------------------------------------------------------
  ## Build B using alpha = 0.0027
  ## ----------------------------------------------------------
  
  design_b_g <- make_stage_design(
    
    stage = "B_mandatory_sample_A_missing",
    
    n_sample = n_b_g,
    
    alpha = alpha_b_g,
    
    k0 = k0_g,
    
    lambda0 = lambda0_g
  )
  
  
  ## ----------------------------------------------------------
  ## Each Y_b defines one possible B-only rule
  ## ----------------------------------------------------------
  
  for (
    s in seq_len(
      nrow(scenario_grid)
    )
  ) {
    
    
    d_s <- (
      scenario_grid$d[s]
    )
    
    
    k1_s <- k0_g
    
    
    ## ========================================================
    ## Table 6:
    ##
    ## mu1 = mu0 * (1 + d)
    ##
    ## d <= 0
    ## ========================================================
    
    mean_1_s <- (
      weibull_mean0_g *
        (1 + d_s)
    )
    
    
    lambda1_s <- (
      
      mean_1_s /
        
        gamma(
          1 + 1 / k1_s
        )
    )
    
    
    var_1_s <- (
      
      lambda1_s^2 *
        
        (
          gamma(
            1 + 2 / k1_s
          ) -
            
            gamma(
              1 + 1 / k1_s
            )^2
        )
    )
    
    
    ## ========================================================
    ## LCL / LOW SIDE ONLY
    ## ========================================================
    
    chart_side_used_s <- "low"
    
    
    mean_shift_direction_s <- ifelse(
      d_s < 0,
      "decrease",
      "same"
    )
    
    
    ## --------------------------------------------------------
    ## Probability that an observation falls BELOW
    ## tcrit_low under H1
    ## --------------------------------------------------------
    
    p_h1_b <- pweibull(
      
      design_b_g$tcrit_low,
      
      shape = k1_s,
      
      scale = lambda1_s,
      
      lower.tail = TRUE
    )
    
    
    ## --------------------------------------------------------
    ## Probability that B signals:
    ##
    ## P(N_b >= Y_b | H1)
    ## --------------------------------------------------------
    
    alpha1_b <- pbinom(
      
      design_b_g$Y - 1,
      
      size = n_b_g,
      
      prob = p_h1_b,
      
      lower.tail = FALSE
    )
    
    
    ## ========================================================
    ## Since there is NO sample A:
    ##
    ## global signal probability = B signal probability
    ## ========================================================
    
    alpha1_global_h1 <- alpha1_b
    
    
    ARL1 <- ifelse(
      
      alpha1_global_h1 > 0,
      
      1 / alpha1_global_h1,
      
      Inf
    )
    
    
    ## ========================================================
    ## B is mandatory:
    ##
    ## sample size is always n_b
    ## ========================================================
    
    Avg_Sample0_per_cycle <- n_b_g
    
    Avg_Sample1_per_cycle <- n_b_g
    
    
    ## --------------------------------------------------------
    ## Build output
    ## --------------------------------------------------------
    
    tmp <- data.frame(
      
      sample_A_collected = FALSE,
      
      sample_B_mandatory = TRUE,
      
      n_a = 0,
      
      n_b = n_b_g,
      
      Y_a = NA_integer_,
      
      Y_b = design_b_g$Y,
      
      alpha_a_h0 = NA_real_,
      
      alpha_b_h0 = alpha_b_g,
      
      alpha_global_h0 = alpha_global_h0,
      
      ARL0_target = ARL0_target,
      
      p1_b = design_b_g$p1,
      
      alpha_b_h0_check = design_b_g$alpha_h0_check,
      
      ARL0_check = design_b_g$ARL0_check,
      
      ASS0 = ASS0_g,
      
      Avg_Sample0_per_cycle = Avg_Sample0_per_cycle,
      
      Avg_Sample1_per_cycle = Avg_Sample1_per_cycle,
      
      k0 = k0_g,
      
      lambda0 = lambda0_g,
      
      k1 = k1_s,
      
      lambda1 = lambda1_s,
      
      d = d_s,
      
      mean_weibull_0 = weibull_mean0_g,
      
      var_weibull_0 = weibull_var0_g,
      
      mean_weibull_1 = mean_1_s,
      
      var_weibull_1 = var_1_s,
      
      mean_shift_direction = mean_shift_direction_s,
      
      chart_side_used = chart_side_used_s,
      
      tcrit_low_b = design_b_g$tcrit_low,
      
      p_h1_b = p_h1_b,
      
      alpha1_b = alpha1_b,
      
      alpha1_global_h1 = alpha1_global_h1,
      
      ARL1 = ARL1,
      
      stringsAsFactors = FALSE
    )
    
    
    performance_list[[counter]] <- tmp
    
    
    counter <- counter + 1
  }
}


special_performance_table <- do.call(
  rbind,
  performance_list
)


## ============================================================
## Step 06 - Attach Table 6 benchmarks
## ============================================================

mapped_scenarios <- merge(
  
  special_performance_table,
  
  table6_benchmark,
  
  by = c(
    "k0",
    "d"
  ),
  
  all.x = TRUE
)


mapped_scenarios$k <- (
  mapped_scenarios$k0
)


## ============================================================
## Step 07 - ARL1 comparisons
## ============================================================

mapped_scenarios$ARL1_improvement_vs_table6_true_n3 <- (
  
  mapped_scenarios$table6_true_ARL_n3 -
    mapped_scenarios$ARL1
)


mapped_scenarios$ARL1_ratio_vs_table6_true_n3 <- (
  
  mapped_scenarios$ARL1 /
    mapped_scenarios$table6_true_ARL_n3
)


mapped_scenarios$better_than_table6_true_n3 <- (
  
  mapped_scenarios$ARL1 <=
    mapped_scenarios$table6_true_ARL_n3
)


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
## Step 08 - Final selected columns
## ============================================================

mapped_scenarios <- mapped_scenarios[, c(
  
  ## Special-case description
  "sample_A_collected",
  
  "sample_B_mandatory",
  
  "n_a",
  
  "n_b",
  
  "Y_a",
  
  "Y_b",
  
  
  ## H0 calibration
  "alpha_a_h0",
  
  "alpha_b_h0",
  
  "alpha_global_h0",
  
  "ARL0_target",
  
  "alpha_b_h0_check",
  
  "ARL0_check",
  
  
  ## Sample size
  "ASS0",
  
  
  ## B design
  "p1_b",
  
  "tcrit_low_b",
  
  
  ## Weibull scenario
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
  
  
  ## H1 probabilities
  "p_h1_b",
  
  "alpha1_b",
  
  "alpha1_global_h1",
  
  
  ## Table 6 ARL benchmarks
  "table6_true_ARL_n3",
  
  "table6_true_ARL_n5",
  
  "table6_true_ARL_n10",
  
  "table6_true_ARL_n30",
  
  
  ## Table 6 LCL benchmarks
  "table6_true_LCL_n3",
  
  "table6_true_LCL_n5",
  
  "table6_true_LCL_n10",
  
  "table6_true_LCL_n30",
  
  
  ## Proposed special-case ARL1
  "ARL1"
)]


## ============================================================
## Step 09 - Sort
## ============================================================

mapped_scenarios <- mapped_scenarios[
  
  order(
    
    mapped_scenarios$ASS0,
    
    mapped_scenarios$k,
    
    mapped_scenarios$d,
    
    mapped_scenarios$ARL1,
    
    mapped_scenarios$Y_b
  ),
]


## ============================================================
## Step 10 - Design ID
## ============================================================

mapped_scenarios$design_id <- paste0(
  
  "SKIPA",
  
  "_ASS",
  mapped_scenarios$ASS0,
  
  "_nb",
  mapped_scenarios$n_b,
  
  "_Yb",
  mapped_scenarios$Y_b
)


## ============================================================
## Step 11 - Export
## ============================================================

write.csv(
  
  mapped_scenarios,
  
  "np_weibull_skipA_LCL_table6_ASS_3_5_7_10_15_20.csv",
  
  row.names = FALSE
)


## ============================================================
## Step 12 - Final checks
## ============================================================

cat("\n============================================================\n")
cat("FINAL CHECKS\n")
cat("============================================================\n")


cat(
  "Chart side: LCL / LOW\n"
)


cat(
  "Sample A collected: NO\n"
)


cat(
  "Sample B mandatory: YES\n"
)


cat(
  "alpha_global_h0:",
  alpha_global_h0,
  "\n"
)


cat(
  "ARL0 target:",
  ARL0_target,
  "\n"
)


cat(
  "alpha_b special:",
  alpha_global_h0,
  "\n"
)


cat(
  "ASS values:",
  paste(
    sort(
      unique(
        mapped_scenarios$ASS0
      )
    ),
    collapse = ", "
  ),
  "\n"
)


cat(
  "Unique n_b values:",
  paste(
    sort(
      unique(
        mapped_scenarios$n_b
      )
    ),
    collapse = ", "
  ),
  "\n"
)


cat(
  "Number of mapped rows:",
  nrow(mapped_scenarios),
  "\n"
)


## ------------------------------------------------------------
## Verify H0 calibration
## ------------------------------------------------------------

max_alpha_error <- max(
  
  abs(
    mapped_scenarios$alpha_b_h0_check -
      alpha_global_h0
  ),
  
  na.rm = TRUE
)


max_ARL0_error <- max(
  
  abs(
    mapped_scenarios$ARL0_check -
      ARL0_target
  ),
  
  na.rm = TRUE
)


cat("\nH0 CALIBRATION CHECK\n")


cat(
  "Maximum alpha error:",
  max_alpha_error,
  "\n"
)


cat(
  "Maximum ARL0 error:",
  max_ARL0_error,
  "\n"
)


cat(
  "ARL0 min:",
  min(
    mapped_scenarios$ARL0_check,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "ARL0 max:",
  max(
    mapped_scenarios$ARL0_check,
    na.rm = TRUE
  ),
  "\n"
)


cat("\nFiles written:\n")


cat(
  "- table6_benchmark_LCL_n3_n5_n10_n30.csv\n"
)


cat(
  "- np_weibull_skipA_LCL_table6_ASS_3_5_7_10_15_20.csv\n"
)

