# dea_run_all_cases2.R
# Robust, fixed pipeline for three cases: lvpico, lvnorm, fds
# - Runs Benchmarking::dea: VRS output-oriented, VRS input-oriented, CRS input-oriented
# - Extracts: phi (input-style, <=1), eff_original, per-input s_minus, per-output s_plus, peers
# - Attempts super-efficiency (VRS output-oriented) using XREF/YREF if supported
# - Writes outputs to results/<case>/ and writes results/summary_across_cases.csv
#
# Run:
# 1) Put this file in your working directory that contains:
#      inn_lvpico.csv out_lvpico.csv
#      inn_lvnorm.csv out_lvnorm.csv
#      inn_fds.csv    out_fds.csv
# 2) In R/RStudio console run:
#      source("dea_run_all_cases2.R")
#
# NOTE: this script avoids rgl/deaR and only uses Benchmarking (no OpenGL issues).
# It is defensive about 1-column inputs/outputs and about shapes of returned matrices.

# -------------------- Setup --------------------
if(!requireNamespace("Benchmarking", quietly = TRUE)) {
  install.packages("Benchmarking", repos = "https://cloud.r-project.org")
}
library(Benchmarking)

# -------------------- Config --------------------
cases <- c("lvpico", "lvnorm", "fds")
input_prefix  <- "inn_"
output_prefix <- "out_"
dmuid_col <- "DMUs"       # change if your CSV uses a different DMU id column name
results_root <- "results"
tol_lambda <- 1e-8

ensure_dir <- function(d) if(!dir.exists(d)) dir.create(d, recursive = TRUE)
safe_read <- function(path) {
  if(!file.exists(path)) stop("File not found: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# -------------------- Utility: ensure matrix with correct dims --------------------
# Try to coerce 'mat' into an nrows x ncols matrix, handling vectors and transposed shapes.
ensure_matrix_dims <- function(mat, nrows, ncols, name = "matrix") {
  if(is.null(mat)) return(NULL)
  # if atomic vector of correct length, replicate to rows if needed
  if(is.atomic(mat) && !is.matrix(mat)) {
    if(length(mat) == ncols) {
      return(matrix(rep(mat, each = nrows), nrow = nrows, ncol = ncols, byrow = FALSE,
                    dimnames = list(NULL, paste0("c", seq_len(ncols)))))
    }
    if(length(mat) == nrows) {
      return(matrix(mat, nrow = nrows, ncol = ncols, byrow = FALSE,
                    dimnames = list(NULL, paste0("c", seq_len(ncols)))))
    }
    # unexpected length -> try to coerce to matrix by column
    return(matrix(as.numeric(mat), nrow = nrows, ncol = ncols, byrow = FALSE))
  }
  # if it's a matrix, try to align dims
  if(is.matrix(mat)) {
    r <- nrow(mat); c <- ncol(mat)
    if(r == nrows && c == ncols) return(mat)
    if(r == ncols && c == nrows) return(t(mat))
    # if single row but ncols matches, replicate rows
    if(r == 1 && c == ncols) return(matrix(rep(mat, each = nrows), nrow = nrows, ncol = ncols, byrow = FALSE))
    # if single column and ncols == 1, replicate to ncols
    if(r == nrows && c == 1 && ncols == 1) return(mat)
    if(r == nrows && c == 1 && ncols > 1) return(matrix(rep(mat, times = ncols), nrow = nrows, ncol = ncols))
    # otherwise try transpose if that helps
    if(r == ncols && c == nrows) return(t(mat))
    # fallback: coerce to numeric then shape
    numeric_vec <- as.numeric(mat)
    if(length(numeric_vec) >= nrows * ncols) {
      return(matrix(numeric_vec[1:(nrows*ncols)], nrow = nrows, ncol = ncols))
    } else if(length(numeric_vec) == ncols) {
      return(matrix(rep(numeric_vec, each = nrows), nrow = nrows, ncol = ncols))
    } else if(length(numeric_vec) == nrows) {
      return(matrix(rep(numeric_vec, times = ncols), nrow = nrows, ncol = ncols))
    } else {
      stop(sprintf("Cannot coerce %s to %dx%d matrix (found %dx%d).", name, nrows, ncols, r, c))
    }
  }
  # fallback: coerce to numeric and shape
  numeric_vec <- as.numeric(mat)
  if(length(numeric_vec) < nrows * ncols) {
    # try to repeat
    numeric_vec <- rep(numeric_vec, length.out = nrows * ncols)
  }
  return(matrix(numeric_vec, nrow = nrows, ncol = ncols))
}

# -------------------- Core: write_model_results (robust) --------------------
write_model_results <- function(out_dir, model_name, X, Y, res) {
  n <- nrow(X); m <- ncol(X); s <- ncol(Y)
  dmus <- rownames(X)

  # Original efficiency from package
  eff_orig <- as.numeric(res$eff)
  # Convert to input-style phi (<=1). For output-oriented res$eff >=1 -> phi = 1/eff
  phi <- ifelse(eff_orig >= 1, 1/eff_orig, eff_orig)
  phi[abs(phi - 1) < 1e-10] <- 1

  # Try to get projected X and Y (Xopt/Yopt); coerce shapes robustly
  Xopt_mat <- NULL; Yopt_mat <- NULL
  if(!is.null(res$Xopt)) {
    Xopt_mat <- tryCatch(ensure_matrix_dims(res$Xopt, n, m, "Xopt"), error = function(e) NULL)
  }
  if(!is.null(res$Yopt)) {
    Yopt_mat <- tryCatch(ensure_matrix_dims(res$Yopt, n, s, "Yopt"), error = function(e) NULL)
  }

  # Try compute from lambda if needed
  lambda <- NULL
  if(!is.null(res$lambda)) {
    lambda <- res$lambda
    # ensure lambda is matrix with dims = n x n (rows references, cols evaluated)
    if(is.vector(lambda)) lambda <- matrix(lambda, nrow = length(lambda), ncol = 1)
    # if dims mismatch, try to coerce
    if(is.matrix(lambda)) {
      # if lambda has appropriate dims use as is; else try to pad/trim
      lr <- nrow(lambda); lc <- ncol(lambda)
      if(!(lr == n && lc == n)) {
        # if it's n x 1 or 1 x n, convert appropriately to n x n by repeating
        if(lr == n && lc == 1) lambda <- matrix(rep(lambda, n), nrow = n, ncol = n)
        else if(lr == 1 && lc == n) lambda <- matrix(rep(lambda, each = n), nrow = n, ncol = n)
        else if(lr == n && lc < n) {
          # pad columns with zeros
          new <- matrix(0, nrow = n, ncol = n)
          new[, seq_len(lc)] <- lambda
          lambda <- new
        } else if(lc == n && lr < n) {
          new <- matrix(0, nrow = n, ncol = n)
          new[seq_len(lr), ] <- lambda
          lambda <- new
        } else {
          # try to coerce to n x n by repetition
          lambda <- matrix(as.numeric(lambda), nrow = n, ncol = n, byrow = TRUE)
        }
      }
    } else {
      # not matrix -> coerce to n x n
      lambda <- matrix(as.numeric(lambda), nrow = n, ncol = n, byrow = TRUE)
    }
  }

  # If still missing projections, compute via convex combination if lambda available
  if(is.null(Xopt_mat) && !is.null(lambda)) {
    # Xopt = t( t(X) %*% lambda )  => produces n x m
    Xopt_mat <- tryCatch(ensure_matrix_dims(t(t(X) %*% lambda), n, m, "Xopt_from_lambda"), error = function(e) NULL)
  }
  if(is.null(Yopt_mat) && !is.null(lambda)) {
    Yopt_mat <- tryCatch(ensure_matrix_dims(t(t(Y) %*% lambda), n, s, "Yopt_from_lambda"), error = function(e) NULL)
  }

  # Fallback: radial projection using phi
  if(is.null(Xopt_mat)) Xopt_mat <- ensure_matrix_dims(X * phi, n, m, "Xopt_radial")
  if(is.null(Yopt_mat)) Yopt_mat <- ensure_matrix_dims(Y / pmax(phi, 1e-12), n, s, "Yopt_radial")

  # Compute slacks: s_minus = pmax(0, X - Xopt), s_plus = pmax(0, Yopt - Y)
  s_minus_mat <- pmax(0, (X - Xopt_mat))
  s_plus_mat  <- pmax(0, (Yopt_mat - Y))

  # Ensure s_minus_mat and s_plus_mat are matrices with proper dims
  s_minus_mat <- ensure_matrix_dims(s_minus_mat, n, m, "s_minus")
  s_plus_mat  <- ensure_matrix_dims(s_plus_mat, n, s, "s_plus")

  # peers from lambda: for column j, peers are rows with lambda[,j] > tol
  peers_vec <- rep(NA_character_, n)
  if(!is.null(lambda)) {
    for(j in seq_len(n)) {
      inds <- which(lambda[, j] > tol_lambda)
      if(length(inds) == 0) peers_vec[j] <- NA else peers_vec[j] <- paste(rownames(X)[inds], collapse = ",")
    }
  }

  # Assemble dataframe
  df <- data.frame(DMU = dmus, phi = phi, eff_original = eff_orig, stringsAsFactors = FALSE)
  df$s_minus_sum <- rowSums(s_minus_mat)
  df$s_plus_sum  <- rowSums(s_plus_mat)
  for(i in seq_len(m)) df[[paste0("s_minus_input", i)]] <- s_minus_mat[, i]
  for(k in seq_len(s)) df[[paste0("s_plus_output", k)]] <- s_plus_mat[, k]
  df$peers <- peers_vec

  csv_path <- file.path(out_dir, paste0("results_", model_name, ".csv"))
  write.csv(df, file = csv_path, row.names = FALSE)
  saveRDS(res, file = file.path(out_dir, paste0("model_", model_name, ".rds")))
  invisible(list(df = df, csv = csv_path))
}

# -------------------- Super-efficiency function (safe) --------------------
compute_super_eff <- function(X, Y) {
  n <- nrow(X)
  super_eff <- rep(NA_real_, n)
  names(super_eff) <- rownames(X)
  if(n <= 1) return(super_eff)

  # test whether dea() accepts XREF/YREF (use tmp variable, not underscore)
  test_ok <- tryCatch({
    tmp <- dea(X[1, , drop = FALSE], Y[1, , drop = FALSE],
               XREF = X[-1, , drop = FALSE], YREF = Y[-1, , drop = FALSE],
               RTS = "vrs", ORIENTATION = "out")
    TRUE
  }, error = function(e) FALSE)

  if(!test_ok) {
    warning("Benchmarking::dea does not accept XREF/YREF in this build. Super-efficiency will be skipped.")
    return(super_eff)
  }

  for(i in seq_len(n)) {
    Xi <- X[i, , drop = FALSE]; Yi <- Y[i, , drop = FALSE]
    Xref <- X[-i, , drop = FALSE]; Yref <- Y[-i, , drop = FALSE]
    res_i <- tryCatch({
      dea(Xi, Yi, XREF = Xref, YREF = Yref, RTS = "vrs", ORIENTATION = "out")
    }, error = function(e) {
      warning("Super-eff error for DMU ", rownames(X)[i], ": ", conditionMessage(e))
      NULL
    })
    if(!is.null(res_i)) super_eff[i] <- as.numeric(res_i$eff)
  }
  super_eff
}

# -------------------- Main loop --------------------
ensure_dir(results_root)
case_summaries <- list()

for(case in cases) {
  cat("\n===== CASE:", case, "=====\n")
  in_file <- paste0(input_prefix, case, ".csv")
  out_file <- paste0(output_prefix, case, ".csv")
  if(!file.exists(in_file)) stop("Missing file: ", in_file)
  if(!file.exists(out_file)) stop("Missing file: ", out_file)

  inn <- safe_read(in_file)
  out <- safe_read(out_file)

  # require DMU id column
  if(!(dmuid_col %in% colnames(inn)) || !(dmuid_col %in% colnames(out))) {
    stop("DMU id column '", dmuid_col, "' must be present in both CSVs for case ", case)
  }

  # set rownames and drop id column
  rownames(inn) <- as.character(inn[[dmuid_col]]); inn[[dmuid_col]] <- NULL
  rownames(out) <- as.character(out[[dmuid_col]]); out[[dmuid_col]] <- NULL

  # keep common DMUs in same order
  common <- intersect(rownames(inn), rownames(out))
  if(length(common) == 0) stop("No DMUs in common for case ", case)
  inn <- inn[common, , drop = FALSE]
  out <- out[common, , drop = FALSE]

  # Convert to numeric matrices and ensure dims
  X <- as.matrix(sapply(inn, as.numeric))
  Y <- as.matrix(sapply(out, as.numeric))
  rownames(X) <- rownames(inn); rownames(Y) <- rownames(out)

  if(any(is.na(X)) || any(is.na(Y))) stop("NA values found in inputs/outputs for case ", case)
  if(any(X < 0) || any(Y < 0)) stop("Negative values detected in data for case ", case)

  out_dir <- file.path(results_root, case)
  ensure_dir(out_dir)

  # 1) VRS output-oriented
  cat("Running VRS (vrs, output-oriented)...\n")
  res_vrs_out <- dea(X, Y, RTS = "vrs", ORIENTATION = "out")
  write_model_results(out_dir, "vrs_output", X, Y, res_vrs_out)

  # 2) VRS input-oriented
  cat("Running VRS (vrs, input-oriented)...\n")
  res_vrs_in <- dea(X, Y, RTS = "vrs", ORIENTATION = "in")
  ret_vrs_in <- write_model_results(out_dir, "vrs_input", X, Y, res_vrs_in)

  # 3) CRS input-oriented
  cat("Running CRS (crs, input-oriented)...\n")
  res_crs_in <- dea(X, Y, RTS = "crs", ORIENTATION = "in")
  ret_crs_in <- write_model_results(out_dir, "crs_input", X, Y, res_crs_in)

  # Scale-efficiency summary
  phi_vrs_in <- ret_vrs_in$df$phi
  phi_crs_in <- ret_crs_in$df$phi
  scale_eff <- phi_crs_in / phi_vrs_in
  scale_df <- data.frame(DMU = ret_vrs_in$df$DMU,
                         phi_crs = phi_crs_in,
                         phi_vrs = phi_vrs_in,
                         scale_efficiency = scale_eff,
                         stringsAsFactors = FALSE)
  write.csv(scale_df, file = file.path(out_dir, "scale_efficiency_summary.csv"), row.names = FALSE)

  # Super-efficiency (attempt)
  cat("Attempting super-efficiency (VRS output-oriented) ...\n")
  super_eff_raw <- compute_super_eff(X, Y)
  super_phi <- ifelse(!is.na(super_eff_raw) & super_eff_raw > 0, 1 / super_eff_raw, NA)
  order_non_na <- order(-super_eff_raw, na.last = TRUE)
  super_rank <- rep(NA_integer_, length(super_eff_raw))
  super_rank[order_non_na] <- seq_len(length(super_eff_raw))
  super_df <- data.frame(DMU = rownames(X),
                         super_eff_outputstyle = super_eff_raw,
                         super_phi_inputstyle = super_phi,
                         super_rank = super_rank,
                         stringsAsFactors = FALSE)
  write.csv(super_df, file = file.path(out_dir, "super_efficiency.csv"), row.names = FALSE)

  # Simple histograms
  png(file.path(out_dir, "hist_vrs_out_phi.png"), width = 800, height = 600)
  hist(1 / as.numeric(res_vrs_out$eff), main = paste0("VRS out phi (input-style) - ", case),
       xlab = "phi (<=1)", col = "lightblue")
  dev.off()

  png(file.path(out_dir, "hist_vrs_in_phi.png"), width = 800, height = 600)
  hist(as.numeric(res_vrs_in$eff), main = paste0("VRS in phi - ", case), xlab = "phi (<=1)", col = "lightgreen")
  dev.off()

  png(file.path(out_dir, "hist_crs_in_phi.png"), width = 800, height = 600)
  hist(as.numeric(res_crs_in$eff), main = paste0("CRS in phi - ", case), xlab = "phi (<=1)", col = "lightgray")
  dev.off()

  # collect summary info
  case_summary <- list(
    case = case,
    n_dmus = nrow(X),
    inputs = ncol(X),
    outputs = ncol(Y),
    vrs_out_mean_phi = mean(1 / as.numeric(res_vrs_out$eff), na.rm = TRUE),
    vrs_out_median_phi = median(1 / as.numeric(res_vrs_out$eff), na.rm = TRUE),
    vrs_out_pct_efficient = mean((1 / as.numeric(res_vrs_out$eff)) >= 1 - 1e-10) * 100,
    vrs_in_mean_phi = mean(as.numeric(res_vrs_in$eff), na.rm = TRUE),
    crs_in_mean_phi = mean(as.numeric(res_crs_in$eff), na.rm = TRUE),
    super_best_eff = if(all(is.na(super_eff_raw))) NA else max(super_eff_raw, na.rm = TRUE)
  )
  case_summaries[[case]] <- case_summary

  cat("Case", case, "completed. Results in:", out_dir, "\n")
}

# -------------------- Cross-case summary --------------------
summary_rows <- lapply(case_summaries, function(x) {
  data.frame(case = x$case,
             n_dmus = x$n_dmus,
             inputs = x$inputs,
             outputs = x$outputs,
             vrs_out_mean_phi = x$vrs_out_mean_phi,
             vrs_out_median_phi = x$vrs_out_median_phi,
             vrs_out_pct_efficient = x$vrs_out_pct_efficient,
             vrs_in_mean_phi = x$vrs_in_mean_phi,
             crs_in_mean_phi = x$crs_in_mean_phi,
             super_best_eff_outputstyle = x$super_best_eff,
             stringsAsFactors = FALSE)
})
summary_df <- do.call(rbind, summary_rows)
ensure_dir(results_root)
write.csv(summary_df, file = file.path(results_root, "summary_across_cases.csv"), row.names = FALSE)

cat("\nALL CASES COMPLETED. Check the 'results' folder for outputs.\n")