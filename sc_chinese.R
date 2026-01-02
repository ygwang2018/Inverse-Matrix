## =========================================================
## Minimal code for PAPER RESULTS ONLY
## Outputs:
##   tab_expA_markov_ou.csv, tab_expB_ma1.csv
##   fig_expA_solve_time.pdf, fig_expA_logdet_time.pdf, fig_expB_ma1_residual.pdf
## =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

## =========================
## 0) USER SETTINGS (paper)
## =========================
set.seed(1L)
out_dir <- "."

## A: runtime scaling table (display can go to 1e10; compute is capped)
n_axis_max    <- 1e10
n_vec_table   <- 10^(seq(3, 10, by = 0.5))   # 1e3 ... 1e10 (for axis + table)
max_n_compute <- 1e5                         # only compute when n <= this
R_rep_time    <- 5

rho_set <- c(0.2, 0.5, 0.8, 0.95)
h_set   <- c(1e-4, 1e-2, 1e-1, 1)

grid_irregular <- "unif"   # "unif" or "cluster"
delta_markov   <- 1

## B: MA(1) stability
n_set_ma1 <- c(50, 200, 1000)
eps_set   <- c(1e-6, 1e-8)
rho_grid_inside <- seq(-0.49, 0.49, length.out = 61)
rho_grid_boundary <- function(eps) c(0.5 - eps, 0.5 + eps, -0.5 + eps, -0.5 - eps)
k_sing_points <- c(1, 2, 3)
sing_delta <- 1e-6
J_cols_ma1 <- 30

save_tables <- TRUE
save_figs   <- TRUE


## =========================
## 1) Tridiagonal solver + logdet (symmetric)
## =========================
solve_tridiag_sym <- function(a, b, d) {
  n <- length(a)
  if (length(b) != n - 1) stop("length(b) must be n-1")
  if (length(d) != n) stop("length(d) must be n")
  
  cp <- numeric(n - 1)
  dp <- numeric(n)
  
  denom <- a[1]
  cp[1] <- b[1] / denom
  dp[1] <- d[1] / denom
  
  if (n > 2) {
    for (i in 2:(n - 1)) {
      denom <- a[i] - b[i - 1] * cp[i - 1]
      cp[i] <- b[i] / denom
      dp[i] <- (d[i] - b[i - 1] * dp[i - 1]) / denom
    }
  }
  
  denom <- a[n] - b[n - 1] * cp[n - 1]
  dp[n] <- (d[n] - b[n - 1] * dp[n - 1]) / denom
  
  x <- numeric(n)
  x[n] <- dp[n]
  for (i in (n - 1):1) x[i] <- dp[i] - cp[i] * x[i + 1]
  x
}

logdet_tridiag_sym <- function(a, b) {
  n <- length(a)
  if (length(b) != n - 1) stop("length(b) must be n-1")
  
  D0 <- 1
  D1 <- a[1]
  scale_log <- 0
  
  if (!is.finite(D1) || abs(D1) < 1e-300) stop("Nearly singular recursion at start.")
  
  for (k in 2:n) {
    Dk <- a[k] * D1 - (b[k - 1]^2) * D0
    
    s <- max(abs(Dk), abs(D1), abs(D0))
    if (s > 1e100 || s < 1e-100) {
      Dk <- Dk / s
      D1 <- D1 / s
      D0 <- D0 / s
      scale_log <- scale_log + log(s)
    }
    
    D0 <- D1
    D1 <- Dk
    if (!is.finite(D1) || abs(D1) < 1e-300) stop("Nearly singular recursion.")
  }
  
  log(abs(D1)) + scale_log
}


## =========================
## 2) Experiment A: (R_CM + H) via tearing
## =========================
gen_t_regular <- function(n) 1:n

gen_t_irregular <- function(n, type = c("unif", "cluster")) {
  type <- match.arg(type)
  if (type == "unif") {
    d <- runif(n - 1, 0.2, 1.8)
  } else {
    mix <- rbinom(n - 1, 1, 0.35)
    d <- ifelse(mix == 1, runif(n - 1, 0.02, 0.12), runif(n - 1, 0.8, 2.5))
  }
  cumsum(c(0, d))
}

gen_H <- function(n, h0, hetero = FALSE) {
  if (!hetero) return(rep(h0, n))
  i <- 1:n
  h0 * (1 + i / n)^2
}

build_Q_markov <- function(t, rho, delta = 1) {
  n <- length(t)
  td <- t^delta
  dk <- diff(td)
  g  <- rho^dk
  s  <- 1 - rho^(2 * dk)
  
  off <- -g / s
  
  diag <- numeric(n)
  diag[1] <- 1 / s[1]
  diag[n] <- 1 / s[n - 1]
  if (n > 2) {
    for (i in 2:(n - 1)) diag[i] <- 1 / s[i - 1] + 1 / s[i] - 1
  }
  
  list(diag = diag, off = off, dk = dk)
}

solve_Sigma_b_tearing <- function(t, rho, hvec, b, delta = 1) {
  Q <- build_Q_markov(t, rho, delta)
  
  Delta <- 1 / hvec
  aF    <- Q$diag + Delta
  bF    <- Q$off
  
  u <- b * Delta
  y <- solve_tridiag_sym(aF, bF, u)
  u - y * Delta
}

logdet_Sigma_tearing <- function(t, rho, hvec, delta = 1) {
  Q  <- build_Q_markov(t, rho, delta)
  dk <- Q$dk
  
  logdet_R <- sum(log(1 - rho^(2 * dk)))   # log|R_CM|
  Delta <- 1 / hvec
  aF <- Q$diag + Delta
  bF <- Q$off
  logdet_F <- logdet_tridiag_sym(aF, bF)   # log|Q+Delta|
  
  logdet_R + logdet_F + sum(log(hvec))     # + log|H|
}

one_run_A <- function(n, rho, h0, grid = c("regular", "irregular"),
                      irregular_type = "unif", hetero = FALSE, delta = 1) {
  grid <- match.arg(grid)
  b <- rnorm(n)
  t <- if (grid == "regular") gen_t_regular(n) else gen_t_irregular(n, irregular_type)
  hvec <- gen_H(n, h0, hetero)
  
  tm_solve <- system.time({
    solve_Sigma_b_tearing(t, rho, hvec, b, delta)
  })["elapsed"]
  
  tm_logdet <- system.time({
    logdet_Sigma_tearing(t, rho, hvec, delta)
  })["elapsed"]
  
  c(time_solve = as.numeric(tm_solve),
    time_logdet = as.numeric(tm_logdet))
}

run_expA <- function(n_vec, rho_set, h_set, grid, hetero,
                     irregular_type = "unif", delta = 1, R_rep = 5) {
  out <- list()
  idx <- 1
  
  for (rho in rho_set) for (h0 in h_set) for (n in n_vec) {
    
    rec <- data.frame(
      grid = grid,
      hetero = hetero,
      rho = rho,
      h = h0,
      n = n,
      mean_time_solve = NA_real_,
      mean_time_logdet = NA_real_,
      stringsAsFactors = FALSE
    )
    
    if (n <= max_n_compute) {
      runs <- replicate(
        R_rep,
        one_run_A(n, rho, h0, grid = grid, irregular_type = irregular_type,
                  hetero = hetero, delta = delta),
        simplify = TRUE
      )
      rec$mean_time_solve  <- mean(runs["time_solve", ])
      rec$mean_time_logdet <- mean(runs["time_logdet", ])
    }
    
    out[[idx]] <- rec
    idx <- idx + 1
  }
  
  do.call(rbind, out)
}


## =========================
## 3) Experiment B: MA(1) closed-form inverse stability
## =========================
apply_MA1_R <- function(x, rho) {
  n <- length(x)
  y <- x
  if (n >= 2) {
    y[1] <- x[1] + rho * x[2]
    for (i in 2:(n - 1)) y[i] <- rho * x[i - 1] + x[i] + rho * x[i + 1]
    y[n] <- rho * x[n - 1] + x[n]
  }
  y
}

inv_MA1_col_closedform <- function(n, rho, j) {
  if (rho == 0) {
    x <- numeric(n); x[j] <- 1
    return(x)
  }
  
  ## boundary limits
  if (abs(rho - 0.5) < 1e-14) {
    x <- numeric(n)
    for (i in 1:n) {
      ii <- min(i, j); jj <- max(i, j)
      x[i] <- ((-1)^(ii - jj)) * (n - jj + 1) * ii / (0.5 * (n + 1))
    }
    return(x)
  }
  if (abs(rho + 0.5) < 1e-14) {
    x <- numeric(n)
    for (i in 1:n) {
      ii <- min(i, j); jj <- max(i, j)
      x[i] <- - (n - jj + 1) * ii / ((-0.5) * (n + 1))
    }
    return(x)
  }
  
  theta <- acos(-1 / (2 * rho))
  denom <- rho * sin(theta) * sin((n + 1) * theta)
  
  sj   <- sin(j * theta)
  sNj1 <- sin((n - j + 1) * theta)
  
  x <- rep(0 + 0i, n)
  for (i in 1:n) {
    if (i <= j) x[i] <- - sNj1 * sin(i * theta) / denom
    else        x[i] <- - sin((n - i + 1) * theta) * sj / denom
  }
  
  if (max(abs(Im(x))) < 1e-10) x <- Re(x)
  x
}

ma1_residual_one <- function(n, rho, j) {
  ej <- numeric(n); ej[j] <- 1
  x  <- inv_MA1_col_closedform(n, rho, j)
  r  <- apply_MA1_R(x, rho)
  sqrt(sum((r - ej)^2))
}

rho_set_B <- function(n) {
  rb <- unlist(lapply(eps_set, rho_grid_boundary))
  rsing <- c()
  for (k in k_sing_points) {
    rk <- -1 / (2 * cos(k * pi / (n + 1)))
    rsing <- c(rsing, rk - sing_delta, rk + sing_delta)
  }
  unique(c(rho_grid_inside, rb, rsing))
}

run_expB <- function(n_set, J_cols, seed = 1L) {
  set.seed(seed)
  out <- list()
  idx <- 1
  
  for (n in n_set) {
    rhos <- rho_set_B(n)
    cols <- sort(sample(1:n, min(J_cols, n)))
    
    for (rho in rhos) {
      res <- sapply(cols, function(j) ma1_residual_one(n, rho, j))
      out[[idx]] <- data.frame(
        n = n,
        rho = rho,
        mean_residual = mean(res),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  do.call(rbind, out)
}


## =========================
## 4) RUN + SAVE TABLES
## =========================
cat("Running Experiment A...\n")
df_A <- rbind(
  run_expA(n_vec_table, rho_set, h_set, grid = "regular",   hetero = FALSE,
           delta = delta_markov, R_rep = R_rep_time),
  run_expA(n_vec_table, rho_set, h_set, grid = "regular",   hetero = TRUE,
           delta = delta_markov, R_rep = R_rep_time),
  run_expA(n_vec_table, rho_set, h_set, grid = "irregular", hetero = FALSE,
           irregular_type = grid_irregular, delta = delta_markov, R_rep = R_rep_time),
  run_expA(n_vec_table, rho_set, h_set, grid = "irregular", hetero = TRUE,
           irregular_type = grid_irregular, delta = delta_markov, R_rep = R_rep_time)
)

cat("Running Experiment B...\n")
df_B <- run_expB(n_set_ma1, J_cols_ma1, seed = 1L)

if (save_tables) {
  write.csv(df_A, file.path(out_dir, "tab_expA_markov_ou.csv"), row.names = FALSE)
  write.csv(df_B, file.path(out_dir, "tab_expB_ma1.csv"), row.names = FALSE)
}


## =========================
## 5) FIGURES (paper)
## =========================
df_A_plot <- df_A
df_A_plot$grid <- factor(df_A_plot$grid, levels = c("regular", "irregular"))

pA_time <- ggplot(
  df_A_plot,
  aes(x = n, y = mean_time_solve,
      color = factor(rho),
      linetype = factor(h),
      group = interaction(rho, h))
) +
  geom_line(na.rm = TRUE, linewidth = 0.6) +
  geom_point(na.rm = TRUE, size = 1.5) +
  scale_x_log10(limits = c(min(n_vec_table), n_axis_max)) +
  scale_y_log10() +
  facet_grid(rows = vars(grid), cols = vars(hetero), labeller = label_both) +
  labs(
    x = "Matrix dimension n",
    y = expression("Mean elapsed time for solving " * Sigma^{-1} * "b (seconds)"),
    title = expression("Experiment A: mean runtime for solving " * Sigma^{-1} * "b"),
    color = expression(rho),
    linetype = "nugget h"
  ) +
  theme_bw(base_size = 12.5) +
  theme(panel.grid.minor = element_blank())

pA_ld <- ggplot(
  df_A_plot,
  aes(x = n, y = mean_time_logdet,
      color = factor(rho),
      linetype = factor(h),
      group = interaction(rho, h))
) +
  geom_line(na.rm = TRUE, linewidth = 0.6) +
  geom_point(na.rm = TRUE, size = 1.5) +
  scale_x_log10(limits = c(min(n_vec_table), n_axis_max)) +
  scale_y_log10() +
  facet_grid(rows = vars(grid), cols = vars(hetero), labeller = label_both) +
  labs(
    x = "Matrix dimension n",
    y = expression("Mean elapsed time for computing log|" * Sigma * "| (seconds)"),
    title = expression("Experiment A: mean runtime for log|" * Sigma * "|"),
    color = expression(rho),
    linetype = "nugget h"
  ) +
  theme_bw(base_size = 12.5) +
  theme(panel.grid.minor = element_blank())

pB <- ggplot(df_B, aes(x = rho, y = mean_residual)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "gray45") +
  scale_y_log10() +
  facet_wrap(~ n, scales = "free_y") +
  labs(
    x = expression(rho),
    y = "Mean residual ||R x - e_j||_2",
    title = expression("Experiment B: MA(1) closed-form inverse stability (mean residual)")
  ) +
  theme_bw(base_size = 12.5) +
  theme(panel.grid.minor = element_blank())

if (save_figs) {
  ggsave(file.path(out_dir, "fig_expA_solve_time.pdf"),  pA_time, width = 8.6, height = 5.2)
  ggsave(file.path(out_dir, "fig_expA_logdet_time.pdf"), pA_ld,   width = 8.6, height = 5.2)
  ggsave(file.path(out_dir, "fig_expB_ma1_residual.pdf"), pB,     width = 8.0, height = 4.8)
}

print(pA_time)
print(pA_ld)
print(pB)

cat("Done.\n")
cat("Note: n > max_n_compute are kept on x-axis but not evaluated (NA).\n")

