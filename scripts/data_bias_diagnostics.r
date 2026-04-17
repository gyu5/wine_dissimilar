set.seed(1234)
get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    script_dir <- dirname(normalizePath(script_path))
    return(normalizePath(file.path(script_dir, "..")))
  }
  normalizePath(file.path(getwd(), ".."))
}
read_matrix_csv <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  mat
}
extract_offdiag_finite <- function(mat) {
  if (nrow(mat) != ncol(mat)) stop("相関行列は正方行列である必要があります", call. = FALSE)
  diag(mat) <- NA_real_
  vals <- as.vector(mat)
  vals[is.finite(vals)]
}
calc_user_negative_stats <- function(corr_mat) {
  user_ids <- intersect(rownames(corr_mat), colnames(corr_mat))
  rows <- vector("list", length(user_ids))
  for (i in seq_along(user_ids)) {
    uid <- user_ids[i]
    sims <- corr_mat[uid, user_ids]
    sims <- sims[user_ids != uid]
    sims <- sims[is.finite(sims)]
    if (length(sims) == 0) next
    rows[[i]] <- data.frame(
      user_id = uid,
      n_valid_corr = length(sims),
      neg_ratio = mean(sims < 0),
      strong_neg_ratio_02 = mean(sims < -0.2),
      strong_neg_ratio_03 = mean(sims < -0.3),
      mean_corr = mean(sims)
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
calc_neighbor_negative_availability <- function(rating_mat, corr_mat, k_values, n_samples = 3000) {
  user_ids <- intersect(colnames(rating_mat), rownames(corr_mat))
  if (length(user_ids) == 0) stop("評価行列と相関行列のユーザIDが一致しません", call. = FALSE)
  rows <- list()
  for (k in k_values) {
    done <- 0L
    attempts <- 0L
    max_attempts <- n_samples * 500L
    n_cand <- c()
    n_cand_neg <- c()
    n_top <- c()
    n_top_neg <- c()
    any_top_neg <- 0L
    while (done < n_samples && attempts < max_attempts) {
      attempts <- attempts + 1L
      uid <- sample(user_ids, 1)
      u <- match(uid, colnames(rating_mat))
      rated <- which(is.finite(rating_mat[, u]))
      if (length(rated) == 0) next
      item_idx <- sample(rated, 1)
      candidates <- setdiff(which(is.finite(rating_mat[item_idx, ])), u)
      if (length(candidates) == 0) next
      cand_ids <- colnames(rating_mat)[candidates]
      ok_id <- cand_ids %in% colnames(corr_mat)
      if (!any(ok_id)) next
      cand_ids <- cand_ids[ok_id]
      sims <- corr_mat[uid, cand_ids]
      sims <- sims[is.finite(sims)]
      if (length(sims) == 0) next
      ord <- order(abs(sims), decreasing = TRUE)
      ord <- ord[seq_len(min(k, length(ord)))]
      top_sims <- sims[ord]
      n_cand <- c(n_cand, length(sims))
      n_cand_neg <- c(n_cand_neg, sum(sims < 0))
      n_top <- c(n_top, length(top_sims))
      n_top_neg <- c(n_top_neg, sum(top_sims < 0))
      if (sum(top_sims < 0) > 0) any_top_neg <- any_top_neg + 1L
      done <- done + 1L
    }
    rows[[length(rows) + 1L]] <- data.frame(
      k = k,
      samples = done,
      mean_candidates = mean(n_cand),
      mean_negative_candidates = mean(n_cand_neg),
      neg_share_candidates = sum(n_cand_neg) / sum(n_cand),
      mean_topk_size = mean(n_top),
      mean_negative_topk = mean(n_top_neg),
      neg_share_topk = sum(n_top_neg) / sum(n_top),
      any_negative_in_topk_rate = any_top_neg / done
    )
  }
  do.call(rbind, rows)
}
main <- function() {
  root <- get_project_root()
  out_dir <- file.path(root, "results", "diagnostics")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  experts <- read_matrix_csv(file.path(root, "data", "experts.csv"))
  amateurs <- read_matrix_csv(file.path(root, "data", "amateurs.csv"))
  rating_mat <- cbind(experts, amateurs)
  corr_mat <- read_matrix_csv(file.path(root, "data", "correlation.csv"))
  # 1) 相関分布（対角除外・有限値のみ）
  corr_vals <- extract_offdiag_finite(corr_mat)
  corr_summary <- data.frame(
    n_pairs = length(corr_vals),
    neg_ratio = mean(corr_vals < 0),
    pos_ratio = mean(corr_vals > 0),
    near_zero_ratio = mean(abs(corr_vals) < 0.1),
    strong_neg_ratio_02 = mean(corr_vals < -0.2),
    strong_neg_ratio_03 = mean(corr_vals < -0.3),
    strong_pos_ratio_02 = mean(corr_vals > 0.2),
    strong_pos_ratio_03 = mean(corr_vals > 0.3),
    mean_corr = mean(corr_vals),
    median_corr = median(corr_vals)
  )
  # 2) ユーザ別の負相関保有量
  user_stats <- calc_user_negative_stats(corr_mat)
  # 3) signed top-k で実際に負相関近傍を確保できるか
  k_values <- c(1, 3, 5, 10, 20, 30, 50)
  availability <- calc_neighbor_negative_availability(
    rating_mat = rating_mat,
    corr_mat = corr_mat,
    k_values = k_values,
    n_samples = 3000
  )
  write.csv(corr_summary, file.path(out_dir, "corr_global_summary.csv"), row.names = FALSE)
  write.csv(user_stats, file.path(out_dir, "corr_user_negative_stats.csv"), row.names = FALSE)
  write.csv(availability, file.path(out_dir, "signed_negative_neighbor_availability.csv"), row.names = FALSE)
  cat("=== Correlation global summary ===\n")
  print(corr_summary)
  cat("\n=== User-level negative-correlation summary ===\n")
  print(summary(user_stats[, c("neg_ratio", "strong_neg_ratio_02", "strong_neg_ratio_03", "mean_corr")]))
  cat("\n=== Negative-neighbor availability for signed top-k ===\n")
  print(availability, row.names = FALSE)
  cat("\nSaved files:\n")
  cat(file.path(out_dir, "corr_global_summary.csv"), "\n")
  cat(file.path(out_dir, "corr_user_negative_stats.csv"), "\n")
  cat(file.path(out_dir, "signed_negative_neighbor_availability.csv"), "\n")
}
main()