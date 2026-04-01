# opposite_knn_eval.R
# 目的: 各ユーザの好みを反転した仮想ユーザ(全員)を作り、Signed/Baseline k-NNの当たり具合を評価

source("sim.r")  # load_rating_matrix(), evaluate_signed_vs_baseline() を読み込む

invert_ratings <- function(rating_mat, strategy = c("negate", "mirror")) {
  strategy <- match.arg(strategy)
  if (strategy == "negate") {
    return(-rating_mat)
  } else {
    finite <- is.finite(rating_mat)
    mn <- min(rating_mat[finite])
    mx <- max(rating_mat[finite])
    # 大きいほど良い、という順序を反転させる
    return((mn + mx) - rating_mat)
  }
}

summarize_subset <- function(df_subset) {
  list(
    n_test_cases = nrow(df_subset),
    baseline_rmse = sqrt(mean(df_subset$sqerr_baseline, na.rm = TRUE)),
    signed_rmse   = sqrt(mean(df_subset$sqerr_signed, na.rm = TRUE)),
    baseline_mae  = mean(df_subset$err_baseline, na.rm = TRUE),
    signed_mae    = mean(df_subset$err_signed, na.rm = TRUE),
    signed_rmse_gain = (sqrt(mean(df_subset$sqerr_baseline, na.rm=TRUE)) -
                         sqrt(mean(df_subset$sqerr_signed, na.rm=TRUE))),
    signed_mae_gain = (mean(df_subset$err_baseline, na.rm=TRUE) -
                        mean(df_subset$err_signed, na.rm=TRUE)),
    signed_better_rate = mean(df_subset$signed_better, na.rm = TRUE),
    negative_neighbor_usage_rate = mean(df_subset$negative_neighbor_used, na.rm = TRUE)
  )
}

evaluate_opposite_users <- function(
  rating_mat,
  k = 10,
  min_overlap = 5,
  test_per_user = 20,
  seed = 123,
  invert_strategy = "negate"
) {
  n_orig <- ncol(rating_mat)

  rating_aug <- cbind(
    rating_mat,
    invert_ratings(rating_mat, strategy = invert_strategy)
  )

  # sim.r の評価器を、そのまま“反転ユーザ込み”で実行
  out <- evaluate_signed_vs_baseline(
    rating_mat = rating_aug,
    k = k,
    min_overlap = min_overlap,
    test_per_user = test_per_user,
    seed = seed
  )

  df <- out$detailed
  # df$user は「列番号」なので、元ユーザ数 n_orig より大きいのが反転ユーザ
  df$group <- ifelse(df$user > n_orig, "opposite", "original")

  list(
    invert_strategy = invert_strategy,
    n_original_users = n_orig,
    n_total_users = ncol(rating_aug),
    original_summary = summarize_subset(df[df$group == "original", ]),
    opposite_summary = summarize_subset(df[df$group == "opposite", ]),
    detailed = df,
    all_summary = out
  )
}

run_one <- function(path, invert_strategy = "negate") {
  rating_mat <- load_rating_matrix(path)

  res <- evaluate_opposite_users(
    rating_mat = rating_mat,
    k = 10,
    min_overlap = 5,
    test_per_user = 20,  # 重いので必要なら下げてください
    seed = 123,
    invert_strategy = invert_strategy
  )
  cat("\n====", basename(path), "====\n")
  print(res$original_summary)
  print(res$opposite_summary)
  invisible(res)
}

# まずはアマチュア
res_amateurs <- run_one("amateurs.csv", invert_strategy = "negate")

# 次にプロ（ソムリエ）
res_experts  <- run_one("experts.csv", invert_strategy = "negate")

# （任意）反転をミラー変換にしたい場合は以下に変更
# res_amateurs_mirror <- run_one("data/amateurs.csv", invert_strategy = "mirror")