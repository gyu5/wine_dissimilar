############################################
# Signed k-NN vs Baseline k-NN (Leave-One-Out)
# For sparse wine rating matrices
# rows = items (wines), cols = raters
############################################

# -----------------------------
# 1. Utility functions
# -----------------------------
load_rating_matrix <- function(path) {
  # rating matrix の読み込み:
  # - CSVは先頭列が wine id の行ラベル(ヘッダが空のことが多い)になっている
  # - 欠損は空セルで入っている想定
  # それらを頑健に処理する
  na_strings <- c("NA", "NaN", "", "·")
  df <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = na_strings
  )
  
  # よくある余計な列を削除（先頭ヘッダが空のケースにも対応）
  drop_name_cols <- intersect(names(df), c("Unnamed: 0", "X", "", "V1", "...1"))
  if (length(drop_name_cols) > 0) {
    df <- df[, !(names(df) %in% drop_name_cols), drop = FALSE]
  }
  
  # それでも先頭列が wine/item のインデックスになっている場合があるので、
  # 「先頭列が 0..(nrow-1) の完全な連番」なら落とす
  if (ncol(df) >= 2) {
    first_col_num <- suppressWarnings(as.numeric(as.character(df[[1]])))
    if (all(is.finite(first_col_num))) {
      expected <- seq_len(nrow(df)) - 1
      uniq <- sort(unique(first_col_num))
      if (length(uniq) == length(expected) && all(uniq == expected)) {
        df <- df[, -1, drop = FALSE]
      }
    }
  }
  
  # 数値化（文字列混入があっても NA 扱い）
  mat <- data.matrix(df)
  storage.mode(mat) <- "numeric"
  mat
}

# 2人のユーザーのピアソン相関
# min_overlap 未満なら NA
pairwise_pearson <- function(x, y, min_overlap = 5) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < min_overlap) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
}

# 全ユーザー間の相関行列
compute_similarity_matrix <- function(rating_mat, min_overlap = 5) {
  n_users <- ncol(rating_mat)
  sim <- matrix(NA_real_, n_users, n_users)
  
  for (i in seq_len(n_users)) {
    sim[i, i] <- NA_real_
    if (i < n_users) {
      for (j in seq.int(i + 1, n_users)) {
        s <- pairwise_pearson(rating_mat[, i], rating_mat[, j], min_overlap = min_overlap)
        sim[i, j] <- s
        sim[j, i] <- s
      }
    }
  }
  
  colnames(sim) <- colnames(rating_mat)
  rownames(sim) <- colnames(rating_mat)
  sim
}

# ユーザー平均との差分を用いた予測
# Baseline: 正の相関のみ
predict_baseline_knn <- function(user_idx, item_idx, rating_mat, sim_mat, k = 10, min_sim = 0) {
  user_mean <- mean(rating_mat[, user_idx], na.rm = TRUE)
  
  # そのアイテムを評価している他者
  candidates <- which(is.finite(rating_mat[item_idx, ]))
  candidates <- setdiff(candidates, user_idx)
  if (length(candidates) == 0) return(NA_real_)
  
  sims <- sim_mat[user_idx, candidates]
  vals <- rating_mat[item_idx, candidates]
  
  ok <- is.finite(sims) & is.finite(vals) & (sims > min_sim)
  if (!any(ok)) return(NA_real_)
  
  candidates <- candidates[ok]
  sims <- sims[ok]
  vals <- vals[ok]
  
  # 相関の大きい順に上位k
  ord <- order(sims, decreasing = TRUE)
  ord <- ord[seq_len(min(k, length(ord)))]
  
  nn <- candidates[ord]
  w <- sims[ord]
  nn_means <- colMeans(rating_mat[, nn, drop = FALSE], na.rm = TRUE)
  
  denom <- sum(abs(w))
  if (denom == 0) return(NA_real_)
  
  pred <- user_mean + sum(w * (vals[ord] - nn_means)) / denom
  pred
}

# Signed k-NN: 正負の相関をそのまま使う
# 近傍選択は |相関| 上位
predict_signed_knn <- function(user_idx, item_idx, rating_mat, sim_mat, k = 10) {
  user_mean <- mean(rating_mat[, user_idx], na.rm = TRUE)
  
  candidates <- which(is.finite(rating_mat[item_idx, ]))
  candidates <- setdiff(candidates, user_idx)
  if (length(candidates) == 0) return(NA_real_)
  
  sims <- sim_mat[user_idx, candidates]
  vals <- rating_mat[item_idx, candidates]
  
  ok <- is.finite(sims) & is.finite(vals)
  if (!any(ok)) return(NA_real_)
  
  candidates <- candidates[ok]
  sims <- sims[ok]
  vals <- vals[ok]
  
  # 絶対値相関で上位k
  ord <- order(abs(sims), decreasing = TRUE)
  ord <- ord[seq_len(min(k, length(ord)))]
  
  nn <- candidates[ord]
  w <- sims[ord]
  nn_means <- colMeans(rating_mat[, nn, drop = FALSE], na.rm = TRUE)
  
  denom <- sum(abs(w))
  if (denom == 0) return(NA_real_)
  
  pred <- user_mean + sum(w * (vals[ord] - nn_means)) / denom
  pred
}

# 予測値を評価尺度の範囲にクリップ
clip_prediction <- function(x, min_rating, max_rating) {
  if (!is.finite(x)) return(NA_real_)
  pmin(pmax(x, min_rating), max_rating)
}

# -----------------------------
# 2. Leave-one-out evaluation
# -----------------------------
evaluate_signed_vs_baseline <- function(
    rating_mat,
    k = 10,
    min_overlap = 5,
    test_per_user = 20,
    seed = 123
) {
  set.seed(seed)
  
  n_items <- nrow(rating_mat)
  n_users <- ncol(rating_mat)
  
  min_rating <- min(rating_mat, na.rm = TRUE)
  max_rating <- max(rating_mat, na.rm = TRUE)
  
  # 相関はフルデータで計算するのではなく、
  # leave-one-outごとに item_idx の評価を一時的に隠して再計算してもよいが重い。
  # まずは実用上、各テスト点で該当セルだけ隠した行列から user-row/col を再計算。
  # ここでは簡便のため、毎テスト点で sim 行列を再計算。
  
  results <- list()
  row_counter <- 1L
  
  for (u in seq_len(n_users)) {
    rated_items <- which(is.finite(rating_mat[, u]))
    if (length(rated_items) < 2) next
    
    test_items <- sample(rated_items, size = min(test_per_user, length(rated_items)), replace = FALSE)
    
    for (i in test_items) {
      true_rating <- rating_mat[i, u]
      
      # leave-one-out: 当該セルを隠す
      train_mat <- rating_mat
      train_mat[i, u] <- NA_real_
      
      sim_mat <- compute_similarity_matrix(train_mat, min_overlap = min_overlap)
      
      pred_base <- predict_baseline_knn(
        user_idx = u,
        item_idx = i,
        rating_mat = train_mat,
        sim_mat = sim_mat,
        k = k,
        min_sim = 0
      )
      
      pred_signed <- predict_signed_knn(
        user_idx = u,
        item_idx = i,
        rating_mat = train_mat,
        sim_mat = sim_mat,
        k = k
      )
      
      pred_base <- clip_prediction(pred_base, min_rating, max_rating)
      pred_signed <- clip_prediction(pred_signed, min_rating, max_rating)
      
      # signed が負相関隣人を使ったか記録
      candidates <- which(is.finite(train_mat[i, ]))
      candidates <- setdiff(candidates, u)
      sims <- sim_mat[u, candidates]
      vals <- train_mat[i, candidates]
      ok <- is.finite(sims) & is.finite(vals)
      
      neg_used <- NA
      if (any(ok)) {
        sims2 <- sims[ok]
        if (length(sims2) > 0) {
          ord <- order(abs(sims2), decreasing = TRUE)
          ord <- ord[seq_len(min(k, length(ord)))]
          neg_used <- any(sims2[ord] < 0)
        }
      }
      
      results[[row_counter]] <- data.frame(
        user = u,
        item = i,
        true_rating = true_rating,
        pred_baseline = pred_base,
        pred_signed = pred_signed,
        err_baseline = ifelse(is.finite(pred_base), abs(pred_base - true_rating), NA_real_),
        err_signed = ifelse(is.finite(pred_signed), abs(pred_signed - true_rating), NA_real_),
        sqerr_baseline = ifelse(is.finite(pred_base), (pred_base - true_rating)^2, NA_real_),
        sqerr_signed = ifelse(is.finite(pred_signed), (pred_signed - true_rating)^2, NA_real_),
        signed_better = ifelse(is.finite(pred_base) & is.finite(pred_signed),
                               abs(pred_signed - true_rating) < abs(pred_base - true_rating),
                               NA),
        negative_neighbor_used = neg_used
      )
      row_counter <- row_counter + 1L
    }
  }
  
  res_df <- do.call(rbind, results)
  
  summary_list <- list(
    n_test_cases = nrow(res_df),
    baseline_rmse = sqrt(mean(res_df$sqerr_baseline, na.rm = TRUE)),
    signed_rmse = sqrt(mean(res_df$sqerr_signed, na.rm = TRUE)),
    baseline_mae = mean(res_df$err_baseline, na.rm = TRUE),
    signed_mae = mean(res_df$err_signed, na.rm = TRUE),
    signed_better_rate = mean(res_df$signed_better, na.rm = TRUE),
    negative_neighbor_usage_rate = mean(res_df$negative_neighbor_used, na.rm = TRUE),
    detailed = res_df
  )
  
  return(summary_list)
}

# -----------------------------
# 3. Optional: user-level summary
# -----------------------------
summarize_by_user <- function(result_df, rating_mat, min_overlap = 5) {
  sim_mat <- compute_similarity_matrix(rating_mat, min_overlap = min_overlap)
  diag(sim_mat) <- NA_real_
  
  typicality <- colMeans(sim_mat, na.rm = TRUE)
  
  agg <- aggregate(
    cbind(err_baseline, err_signed, sqerr_baseline, sqerr_signed) ~ user,
    data = result_df,
    FUN = mean,
    na.rm = TRUE
  )
  
  agg$rmse_baseline <- sqrt(agg$sqerr_baseline)
  agg$rmse_signed <- sqrt(agg$sqerr_signed)
  agg$rmse_gain <- agg$rmse_baseline - agg$rmse_signed
  agg$typicality <- typicality[agg$user]
  
  agg[order(-agg$rmse_gain), ]
}

# -----------------------------
# 4. Example run
# -----------------------------
# path <- "/mnt/data/amateurs.csv"
# R <- load_rating_matrix(path)
# out <- evaluate_signed_vs_baseline(
#   rating_mat = R,
#   k = 10,
#   min_overlap = 5,
#   test_per_user = 20,
#   seed = 123
# )
#
# print(out$n_test_cases)
# print(out$baseline_rmse)
# print(out$signed_rmse)
# print(out$baseline_mae)
# print(out$signed_mae)
# print(out$signed_better_rate)
# print(out$negative_neighbor_usage_rate)
#
# by_user <- summarize_by_user(out$detailed, R, min_overlap = 5)
# head(by_user, 10)
#
# write.csv(out$detailed, "signed_knn_detailed_results.csv", row.names = FALSE)
# write.csv(by_user, "signed_knn_user_summary.csv", row.names = FALSE)