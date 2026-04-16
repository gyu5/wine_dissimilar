set.seed(1234)

get_project_root <- function() {
  # Rscript 実行時は --file=... からこのスクリプトの場所が取れる
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    script_dir <- dirname(normalizePath(script_path))
    return(normalizePath(file.path(script_dir, "..")))
  }
  # それ以外（対話実行など）は作業ディレクトリ基準で親を想定
  normalizePath(file.path(getwd(), ".."))
}

read_rating_matrix <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  mat
}

read_correlation_matrix <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  mat
}

predict_knn_one <- function(rating_mat, corr_mat, user_id, item_idx, k = 10, use_signed = TRUE) {
  if (!(user_id %in% colnames(rating_mat))) return(NA_real_)
  if (!(user_id %in% rownames(corr_mat))) return(NA_real_)

  u <- match(user_id, colnames(rating_mat))

  true_row <- rating_mat[item_idx, ]
  candidates <- which(is.finite(true_row))
  candidates <- setdiff(candidates, u)
  if (length(candidates) == 0) return(NA_real_)

  cand_ids <- colnames(rating_mat)[candidates]
  # 評価行列にいるが相関行列にいないユーザ（例: 134 など）を先に落とす
  okc <- cand_ids %in% colnames(corr_mat)
  if (!any(okc)) return(NA_real_)
  candidates <- candidates[okc]
  cand_ids <- cand_ids[okc]

  sims <- corr_mat[user_id, cand_ids]

  ok <- is.finite(sims)
  if (!any(ok)) return(NA_real_)
  candidates <- candidates[ok]
  cand_ids <- cand_ids[ok]
  sims <- sims[ok]

  if (use_signed) {
    ord <- order(abs(sims), decreasing = TRUE)
  } else {
    pos <- sims > 0
    if (!any(pos)) return(NA_real_)
    candidates <- candidates[pos]
    sims <- sims[pos]
    ord <- order(sims, decreasing = TRUE)
  }
  ord <- ord[seq_len(min(k, length(ord)))]

  nn <- candidates[ord]
  w <- sims[ord]

  user_mean <- mean(rating_mat[, u], na.rm = TRUE)
  nn_means <- colMeans(rating_mat[, nn, drop = FALSE], na.rm = TRUE)
  vals <- rating_mat[item_idx, nn]

  denom <- sum(abs(w))
  if (!is.finite(denom) || denom == 0) return(NA_real_)

  pred <- user_mean + sum(w * (vals - nn_means)) / denom
  pred
}

simulate_pairwise_choice <- function(rating_mat, corr_mat, k = 10, n_sim = 1000) {
  user_ids <- intersect(colnames(rating_mat), rownames(corr_mat))
  if (length(user_ids) == 0) stop("相関行列と評価行列でユーザIDが一致しません", call. = FALSE)

  correct <- 0L
  done <- 0L
  attempts <- 0L
  max_attempts <- n_sim * 200L

  while (done < n_sim) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      stop("試行回数が最大試行回数を超えました。データがsparseすぎます。", call. = FALSE)
    }

    user_id <- sample(user_ids, 1)
    u <- match(user_id, colnames(rating_mat))

    rated_items <- which(is.finite(rating_mat[, u]))
    if (length(rated_items) < 2) next

    ij <- sample(rated_items, 2, replace = FALSE)
    i <- ij[1]
    j <- ij[2]

    pred_i <- predict_knn_one(rating_mat, corr_mat, user_id, i, k = k, use_signed = TRUE)
    pred_j <- predict_knn_one(rating_mat, corr_mat, user_id, j, k = k, use_signed = TRUE)
    if (!is.finite(pred_i) || !is.finite(pred_j)) next
    if (pred_i == pred_j) next

    true_i <- rating_mat[i, u]
    true_j <- rating_mat[j, u]
    if (!is.finite(true_i) || !is.finite(true_j)) next
    if (true_i == true_j) next

    pred_pref <- pred_i > pred_j
    true_pref <- true_i > true_j

    if (identical(pred_pref, true_pref)) correct <- correct + 1L
    done <- done + 1L
  }

  list(
    n = n_sim,
    k = k,
    correct = correct,
    accuracy = correct / n_sim
  )
}

simulate_pairwise_compare <- function(rating_mat, corr_mat, k = 10, n_sim = 1000, user_ids = NULL) {
  if (is.null(user_ids)) {
    user_ids <- intersect(colnames(rating_mat), rownames(corr_mat))
  } else {
    user_ids <- intersect(user_ids, intersect(colnames(rating_mat), rownames(corr_mat)))
  }
  if (length(user_ids) == 0) stop("相関行列と評価行列でユーザIDが一致しません", call. = FALSE)

  signed_correct <- 0L
  baseline_correct <- 0L
  done <- 0L
  attempts <- 0L
  max_attempts <- n_sim * 300L

  # McNemar 用 2x2:
  # 行 = signed(incorrect/correct), 列 = baseline(incorrect/correct)
  pair_tab <- matrix(
    0L, nrow = 2, ncol = 2,
    dimnames = list(
      signed = c("incorrect", "correct"),
      baseline = c("incorrect", "correct")
    )
  )

  while (done < n_sim) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      stop("試行回数が最大試行回数を超えました。データがsparseすぎます。", call. = FALSE)
    }

    user_id <- sample(user_ids, 1)
    u <- match(user_id, colnames(rating_mat))
    rated_items <- which(is.finite(rating_mat[, u]))
    if (length(rated_items) < 2) next

    ij <- sample(rated_items, 2, replace = FALSE)
    i <- ij[1]
    j <- ij[2]

    true_i <- rating_mat[i, u]
    true_j <- rating_mat[j, u]
    if (!is.finite(true_i) || !is.finite(true_j)) next
    if (true_i == true_j) next
    true_pref <- true_i > true_j

    pred_s_i <- predict_knn_one(rating_mat, corr_mat, user_id, i, k = k, use_signed = TRUE)
    pred_s_j <- predict_knn_one(rating_mat, corr_mat, user_id, j, k = k, use_signed = TRUE)
    pred_b_i <- predict_knn_one(rating_mat, corr_mat, user_id, i, k = k, use_signed = FALSE)
    pred_b_j <- predict_knn_one(rating_mat, corr_mat, user_id, j, k = k, use_signed = FALSE)

    if (!is.finite(pred_s_i) || !is.finite(pred_s_j)) next
    if (!is.finite(pred_b_i) || !is.finite(pred_b_j)) next
    if (pred_s_i == pred_s_j || pred_b_i == pred_b_j) next

    s_ok <- (pred_s_i > pred_s_j) == true_pref
    b_ok <- (pred_b_i > pred_b_j) == true_pref

    if (s_ok) signed_correct <- signed_correct + 1L
    if (b_ok) baseline_correct <- baseline_correct + 1L

    pair_tab[ifelse(s_ok, 2, 1), ifelse(b_ok, 2, 1)] <- pair_tab[ifelse(s_ok, 2, 1), ifelse(b_ok, 2, 1)] + 1L
    done <- done + 1L
  }

  list(
    n = n_sim,
    k = k,
    signed_correct = signed_correct,
    signed_accuracy = signed_correct / n_sim,
    baseline_correct = baseline_correct,
    baseline_accuracy = baseline_correct / n_sim,
    accuracy_gain = (signed_correct - baseline_correct) / n_sim,
    pair_table = pair_tab
  )
}

compute_user_sharpness <- function(corr_mat, rating_mat) {
  user_ids <- intersect(colnames(rating_mat), rownames(corr_mat))
  rows <- vector("list", length(user_ids))
  for (i in seq_along(user_ids)) {
    uid <- user_ids[i]
    sims <- corr_mat[uid, user_ids]
    sims <- sims[user_ids != uid]
    sims <- sims[is.finite(sims)]
    if (length(sims) == 0) next

    u_col <- match(uid, colnames(rating_mat))
    user_ratings <- rating_mat[, u_col]
    user_ratings <- user_ratings[is.finite(user_ratings)]

    rows[[i]] <- data.frame(
      user_id = uid,
      neg_ratio = mean(sims < 0),
      mean_corr = mean(sims),
      rating_sd = if (length(user_ratings) >= 2) stats::sd(user_ratings) else NA_real_,
      n_rated = length(user_ratings)
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

run_stratified_eval <- function(rating_mat, corr_mat, k_values, n_sim, sharp_df, top_frac = 0.25, bottom_frac = 0.25) {
  hi_th <- as.numeric(stats::quantile(sharp_df$neg_ratio, probs = 1 - top_frac, na.rm = TRUE, type = 7))
  lo_th <- as.numeric(stats::quantile(sharp_df$neg_ratio, probs = bottom_frac, na.rm = TRUE, type = 7))

  sharp_users <- sharp_df$user_id[sharp_df$neg_ratio >= hi_th]
  mild_users <- sharp_df$user_id[sharp_df$neg_ratio <= lo_th]

  cat(sprintf("\n[層別条件] top=%.0f%%, bottom=%.0f%%\n", top_frac * 100, bottom_frac * 100))
  cat(sprintf("neg_ratio 閾値: sharp >= %.4f, mild <= %.4f\n", hi_th, lo_th))
  cat(sprintf("ユーザ数: sharp=%d, mild=%d\n", length(sharp_users), length(mild_users)))

  rows <- list()
  for (k in k_values) {
    res_sharp <- simulate_pairwise_compare(rating_mat, corr_mat, k = k, n_sim = n_sim, user_ids = sharp_users)
    res_mild <- simulate_pairwise_compare(rating_mat, corr_mat, k = k, n_sim = n_sim, user_ids = mild_users)

    mc_sharp <- mcnemar.test(res_sharp$pair_table, correct = TRUE)
    mc_mild <- mcnemar.test(res_mild$pair_table, correct = TRUE)

    rows[[length(rows) + 1L]] <- data.frame(
      group = "sharp",
      k = k,
      n = res_sharp$n,
      signed_acc = res_sharp$signed_accuracy,
      baseline_acc = res_sharp$baseline_accuracy,
      gain = res_sharp$accuracy_gain,
      mcnemar_p = mc_sharp$p.value
    )
    rows[[length(rows) + 1L]] <- data.frame(
      group = "mild",
      k = k,
      n = res_mild$n,
      signed_acc = res_mild$signed_accuracy,
      baseline_acc = res_mild$baseline_accuracy,
      gain = res_mild$accuracy_gain,
      mcnemar_p = mc_mild$p.value
    )
  }

  do.call(rbind, rows)
}

main <- function() {
  root <- get_project_root()
  amateurs <- read_rating_matrix(file.path(root, "data", "amateurs.csv"))
  experts <- read_rating_matrix(file.path(root, "data", "experts.csv"))
  rating_mat <- cbind(experts, amateurs) # experts(0-13) + amateurs(14-133) を想定

  corr_mat <- read_correlation_matrix(file.path(root, "data", "correlation.csv"))

  n_sim <- 1000
  k_values <- c(1, 3, 5, 10, 20, 30, 50, 80, 120)

  summary_rows <- list()

  cat("\n=== K sweep: Signed vs Baseline ===\n")
  for (k in k_values) {
    res_cmp <- simulate_pairwise_compare(rating_mat, corr_mat, k = k, n_sim = n_sim)
    mc <- mcnemar.test(res_cmp$pair_table, correct = TRUE)
    signed_bt <- binom.test(res_cmp$signed_correct, res_cmp$n, p = 0.5, alternative = "greater")
    baseline_bt <- binom.test(res_cmp$baseline_correct, res_cmp$n, p = 0.5, alternative = "greater")

    row <- data.frame(
      k = k,
      n = res_cmp$n,
      signed_acc = res_cmp$signed_accuracy,
      baseline_acc = res_cmp$baseline_accuracy,
      gain = res_cmp$accuracy_gain,
      signed_vs_random_p = signed_bt$p.value,
      baseline_vs_random_p = baseline_bt$p.value,
      signed_vs_baseline_mcnemar_p = mc$p.value
    )
    summary_rows[[length(summary_rows) + 1L]] <- row

    cat(sprintf(
      "k=%3d | signed=%.3f baseline=%.3f gain=%+.3f | McNemar p=%.4g\n",
      k, res_cmp$signed_accuracy, res_cmp$baseline_accuracy, res_cmp$accuracy_gain, mc$p.value
    ))
  }

  summary_df <- do.call(rbind, summary_rows)
  cat("\n=== Summary table ===\n")
  print(summary_df, row.names = FALSE)

  sharp_df <- compute_user_sharpness(corr_mat, rating_mat)
  cat("\n=== Sharpness summary (neg_ratio) ===\n")
  print(summary(sharp_df$neg_ratio))

  cat("\n=== Stratified result: top25% vs bottom25% (by neg_ratio) ===\n")
  strat_25 <- run_stratified_eval(
    rating_mat = rating_mat,
    corr_mat = corr_mat,
    k_values = k_values,
    n_sim = n_sim,
    sharp_df = sharp_df,
    top_frac = 0.25,
    bottom_frac = 0.25
  )
  print(strat_25, row.names = FALSE)

  cat("\n=== Threshold sensitivity (k=20, top/bottom %) ===\n")
  sens_fracs <- c(0.20, 0.25, 0.30, 0.40)
  sens_rows <- list()
  for (f in sens_fracs) {
    s <- run_stratified_eval(
      rating_mat = rating_mat,
      corr_mat = corr_mat,
      k_values = c(20),
      n_sim = n_sim,
      sharp_df = sharp_df,
      top_frac = f,
      bottom_frac = f
    )
    s$frac <- f
    sens_rows[[length(sens_rows) + 1L]] <- s
  }
  sens_df <- do.call(rbind, sens_rows)
  print(sens_df[, c("frac", "group", "k", "signed_acc", "baseline_acc", "gain", "mcnemar_p")], row.names = FALSE)
}

main()