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
## Therefore the only sample collected is B.
## ASS = n_b.
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


## ------------------------------------------------------------------
## Shape parameters used in Table 5 from Fernandes et Al. benchmark
## ------------------------------------------------------------------

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
## Table 5 from Benchmark:
## increases in mean
## ------------------------------------------------------------

d_vec <- c(
  0,
  0.01,
  0.05,
  0.10,
  0.20,
  0.30,
  0.40
)


cat("\n============================================================\n")
cat("GLOBAL SETTINGS - SPECIAL CASE: SAMPLE A NOT COLLECTED\n")
cat("============================================================\n")

cat("Chart side             = UCL / HIGH\n")

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
## Find the probability p1 such that
##
## P(Binomial(n,p1) >= Y) = alpha
##
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
## Both low and high critical values are calculated,
## but in this experiment ONLY tcrit_high is used.
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
  ## Check H0 probability for each Y
  ##
  ## By construction this should equal alpha.
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
## Table 5 = increases in mean
## ============================================================

scenario_grid <- data.frame(
  d = d_vec
)


scenario_grid$scenario_id <- paste0(
  
  "d",
  
  gsub(
    "\\.",
    "p",
    as.character(
      scenario_grid$d
    )
  )
)


## ============================================================
## Step 03 - Table 5 benchmark
##
## TRUE ARL values and TRUE UCL values
##
## n = 3, 5, 10 and 30
##
## Including n = 3 here is useful because we will now
## investigate ASS = 3 as well.
## ============================================================

table5_benchmark <- data.frame(
  
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
  
  table5_true_ARL_n3 = c(
    
    ## k = 0.5
    370.370,
    352.110,
    314.460,
    276.243,
    211.864,
    140.056,
    123.762,
    
    ## k = 3
    370.370,
    252.525,
    126.582,
    53.078,
    15.484,
    6.831,
    3.934,
    
    ## k = 5
    370.370,
    234.741,
    70.323,
    21.949,
    5.170,
    2.366,
    1.566,
    
    ## k = 10
    370.370,
    200.000,
    23.158,
    5.124,
    1.477,
    1.082,
    1.014,
    
    ## k = 15
    370.370,
    150.602,
    10.072,
    2.216,
    1.073,
    1.004,
    1.000,
    
    ## k = 20
    370.370,
    100.603,
    4.938,
    1.420,
    1.008,
    1.000,
    1.000
  ),
  
  
  ## ==========================================================
  ## TRUE ARL - n = 5
  ## ==========================================================
  
  table5_true_ARL_n5 = c(
    
    ## k = 0.5
    370.370,
    373.134,
    284.090,
    240.384,
    162.337,
    112.107,
    94.161,
    
    ## k = 3
    370.370,
    295.858,
    124.069,
    43.898,
    11.030,
    4.552,
    2.605,
    
    ## k = 5
    370.370,
    207.468,
    50.607,
    14.667,
    3.286,
    1.640,
    1.209,
    
    ## k = 10
    370.370,
    142.857,
    15.133,
    3.242,
    1.169,
    1.013,
    1.000,
    
    ## k = 15
    370.370,
    107.526,
    5.946,
    1.547,
    1.010,
    1.000,
    1.000,
    
    ## k = 20
    370.370,
    88.652,
    3.239,
    1.153,
    1.000,
    1.000,
    1.000
  ),
  
  
  ## ==========================================================
  ## TRUE ARL - n = 10
  ## ==========================================================
  
  table5_true_ARL_n10 = c(
    
    ## k = 0.5
    370.370,
    322.580,
    276.243,
    210.084,
    140.449,
    99.601,
    65.189,
    
    ## k = 3
    370.370,
    271.739,
    84.745,
    25.654,
    5.510,
    2.334,
    1.484,
    
    ## k = 5
    370.370,
    235.849,
    38.197,
    8.748,
    1.888,
    1.159,
    1.024,
    
    ## k = 10
    370.370,
    112.359,
    8.029,
    1.795,
    1.014,
    1.000,
    1.000,
    
    ## k = 15
    370.370,
    97.276,
    3.411,
    1.128,
    1.000,
    1.000,
    1.000,
    
    ## k = 20
    370.370,
    55.248,
    1.817,
    1.013,
    1.000,
    1.000,
    1.000
  ),
  
  
  ## ==========================================================
  ## TRUE ARL - n = 30
  ## ==========================================================
  
  table5_true_ARL_n30 = c(
    
    ## k = 0.5
    370.370,
    304.878,
    203.252,
    178.571,
    92.764,
    51.493,
    11.978,
    
    ## k = 3
    370.370,
    233.644,
    37.850,
    8.418,
    1.760,
    1.104,
    1.010,
    
    ## k = 5
    370.370,
    225.225,
    15.276,
    2.765,
    1.052,
    1.000,
    1.000,
    
    ## k = 10
    370.370,
    79.617,
    2.679,
    1.044,
    1.000,
    1.000,
    1.000,
    
    ## k = 15
    370.370,
    47.619,
    1.322,
    1.000,
    1.000,
    1.000,
    1.000,
    
    ## k = 20
    370.370,
    25.150,
    1.049,
    1.000,
    1.000,
    1.000,
    1.000
  ),
  
  
  ## ==========================================================
  ## TRUE UCL
  ## ==========================================================
  
  table5_true_UCL_n3 = c(
    
    rep(6.275, 7),
    rep(2.823, 7),
    rep(2.554, 7),
    rep(2.393, 7),
    rep(2.340, 7),
    rep(2.278, 7)
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
  
  
  ## Asymptotic UCL from Table 5
  table5_false_UCL_n3  = rep(2.780, 42),
  table5_false_UCL_n5  = rep(2.780, 42),
  table5_false_UCL_n10 = rep(2.780, 42),
  table5_false_UCL_n30 = rep(2.780, 42),
  
  
  stringsAsFactors = FALSE
)


write.csv(
  
  table5_benchmark,
  
  "table5_benchmark_UCL_n3_n5_n10_n30.csv",
  
  row.names = FALSE
)


## ============================================================
## Step 04 - Special-case parameter grid
##
## A IS NOT COLLECTED
##
## Therefore:
##
## n_a = 0
## alpha_a = NA
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
## Because B is now the only decision stage:
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
## A is missing
## B is mandatory
##
## UCL ONLY
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
  ## B is mandatory and alone must preserve ARL0.
  ## ----------------------------------------------------------
  
  alpha_b_g <- (
    alpha_global_h0
  )
  
  
  ASS0_g <- n_b_g
  
  
  ## ----------------------------------------------------------
  ## Build B using alpha = 0.0027
  ##
  ## This is DIFFERENT from normal double sampling,
  ## where alpha_b = alpha_global / alpha_a.
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
    ## Table 5:
    ## mu1 = mu0 * (1 + d)
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
    ## UCL / HIGH SIDE ONLY
    ## ========================================================
    
    chart_side_used_s <- "high"
    
    
    mean_shift_direction_s <- "increase"
    
    
    ## --------------------------------------------------------
    ## Probability that an observation exceeds the
    ## HIGH critical time under H1
    ## --------------------------------------------------------
    
    p_h1_b <- pweibull(
      
      design_b_g$tcrit_high,
      
      shape = k1_s,
      
      scale = lambda1_s,
      
      lower.tail = FALSE
    )
    
    
    ## --------------------------------------------------------
    ## Probability of B signaling:
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
    ## Since there is NO A:
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
    ## Since B is mandatory:
    ##
    ## sample size is always n_b under H0 and H1
    ## ========================================================
    
    Avg_Sample0_per_cycle <- n_b_g
    
    Avg_Sample1_per_cycle <- n_b_g
    
    
    ## --------------------------------------------------------
    ## Build output rows
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
      
      tcrit_high_b = design_b_g$tcrit_high,
      
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
## Step 06 - Attach Table 5 benchmarks
## ============================================================

mapped_scenarios <- merge(
  
  special_performance_table,
  
  table5_benchmark,
  
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

mapped_scenarios$ARL1_improvement_vs_table5_true_n3 <- (
  
  mapped_scenarios$table5_true_ARL_n3 -
    mapped_scenarios$ARL1
)


mapped_scenarios$ARL1_ratio_vs_table5_true_n3 <- (
  
  mapped_scenarios$ARL1 /
    mapped_scenarios$table5_true_ARL_n3
)


mapped_scenarios$better_than_table5_true_n3 <- (
  
  mapped_scenarios$ARL1 <=
    mapped_scenarios$table5_true_ARL_n3
)


mapped_scenarios$ARL1_improvement_vs_table5_true_n5 <- (
  
  mapped_scenarios$table5_true_ARL_n5 -
    mapped_scenarios$ARL1
)


mapped_scenarios$ARL1_ratio_vs_table5_true_n5 <- (
  
  mapped_scenarios$ARL1 /
    mapped_scenarios$table5_true_ARL_n5
)


mapped_scenarios$better_than_table5_true_n5 <- (
  
  mapped_scenarios$ARL1 <=
    mapped_scenarios$table5_true_ARL_n5
)


mapped_scenarios$ARL1_improvement_vs_table5_true_n10 <- (
  
  mapped_scenarios$table5_true_ARL_n10 -
    mapped_scenarios$ARL1
)


mapped_scenarios$ARL1_ratio_vs_table5_true_n10 <- (
  
  mapped_scenarios$ARL1 /
    mapped_scenarios$table5_true_ARL_n10
)


mapped_scenarios$better_than_table5_true_n10 <- (
  
  mapped_scenarios$ARL1 <=
    mapped_scenarios$table5_true_ARL_n10
)


mapped_scenarios$ARL1_improvement_vs_table5_true_n30 <- (
  
  mapped_scenarios$table5_true_ARL_n30 -
    mapped_scenarios$ARL1
)


mapped_scenarios$ARL1_ratio_vs_table5_true_n30 <- (
  
  mapped_scenarios$ARL1 /
    mapped_scenarios$table5_true_ARL_n30
)


mapped_scenarios$better_than_table5_true_n30 <- (
  
  mapped_scenarios$ARL1 <=
    mapped_scenarios$table5_true_ARL_n30
)


## ============================================================
## Step 08 - Final selected columns
## ============================================================

mapped_scenarios <- mapped_scenarios[, c(
  
  ## ----------------------------------------------------------
  ## Special-case description
  ## ----------------------------------------------------------
  
  "sample_A_collected",
  
  "sample_B_mandatory",
  
  "n_a",
  
  "n_b",
  
  "Y_a",
  
  "Y_b",
  
  
  ## ----------------------------------------------------------
  ## H0 calibration
  ## ----------------------------------------------------------
  
  "alpha_a_h0",
  
  "alpha_b_h0",
  
  "alpha_global_h0",
  
  "ARL0_target",
  
  "alpha_b_h0_check",
  
  "ARL0_check",
  
  
  ## ----------------------------------------------------------
  ## Sample size
  ## ----------------------------------------------------------
  
  "ASS0",
  
  
  ## ----------------------------------------------------------
  ## B design
  ## ----------------------------------------------------------
  
  "p1_b",
  
  "tcrit_high_b",
  
  
  ## ----------------------------------------------------------
  ## Weibull scenario
  ## ----------------------------------------------------------
  
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
  
  
  ## ----------------------------------------------------------
  ## H1 probabilities
  ## ----------------------------------------------------------
  
  "p_h1_b",
  
  "alpha1_b",
  
  "alpha1_global_h1",
  
  
  ## ----------------------------------------------------------
  ## Table 5 ARL benchmarks
  ## ----------------------------------------------------------
  
  "table5_true_ARL_n3",
  
  "table5_true_ARL_n5",
  
  "table5_true_ARL_n10",
  
  "table5_true_ARL_n30",
  
  
  ## ----------------------------------------------------------
  ## Table 5 UCL benchmarks
  ## ----------------------------------------------------------
  
  "table5_true_UCL_n3",
  
  "table5_true_UCL_n5",
  
  "table5_true_UCL_n10",
  
  "table5_true_UCL_n30",
  
  
  ## ----------------------------------------------------------
  ## Proposed special-case ARL1
  ## ----------------------------------------------------------
  
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
  
  "np_weibull_skipA_UCL_table5_ASS_3_5_7_10_15_20.csv",
  
  row.names = FALSE
)


## ============================================================
## Step 12 - Final checks
## ============================================================

cat("\n============================================================\n")
cat("FINAL CHECKS\n")
cat("============================================================\n")


cat(
  "Chart side: UCL / HIGH\n"
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
  "- table5_benchmark_UCL_n3_n5_n10_n30.csv\n"
)


cat(
  "- np_weibull_skipA_UCL_table5_ASS_3_5_7_10_15_20.csv\n"
)

