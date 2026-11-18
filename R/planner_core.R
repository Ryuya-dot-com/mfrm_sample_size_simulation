# Core simulation and decision engine for MFRM sample-size & rating-design planning.
# Powered by mfrmr (JMLE and MMLE under RSM / PCM across arbitrary K >= 2 facets).

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

build_p_candidates <- function(p_range, n_points = 3L) {
  lo <- max(6L, as.integer(min(p_range)))
  hi <- max(lo, as.integer(max(p_range)))
  n_pt <- max(1L, min(6L, as.integer(n_points)))
  if (n_pt == 1L || lo == hi) return(lo)
  unique(as.integer(round(seq(lo, hi, length.out = n_pt))))
}

normalize_extra_facets <- function(extra_facets = NULL) {
  if (is.null(extra_facets) || nrow(extra_facets) == 0L) {
    return(data.frame(name = character(0), levels = integer(0), sd = numeric(0), stringsAsFactors = FALSE))
  }
  nms <- trimws(as.character(extra_facets$name))
  nms[!nzchar(nms)] <- sprintf("Facet%d", which(!nzchar(nms)))
  nms <- make.names(nms, unique = TRUE)
  reserved <- c("Person", "Rater", "Score")
  clash <- nms %in% reserved
  nms[clash] <- paste0(nms[clash], "_Facet")
  nms <- make.unique(nms, sep = "_")

  lvls <- pmax(2L, as.integer(extra_facets$levels))
  sds  <- pmax(0.05, as.numeric(extra_facets$sd))
  data.frame(name = nms, levels = lvls, sd = sds, stringsAsFactors = FALSE)
}

build_person_rater_pairs <- function(n_person, n_rater, assignment = c("rotating", "anchor", "crossed"),
                                     raters_per_person = 2L, anchor_persons = 6L) {
  assignment <- match.arg(assignment)
  persons <- sprintf("P%03d", seq_len(n_person))
  raters  <- sprintf("R%02d", seq_len(n_rater))

  if (assignment == "crossed" || raters_per_person >= n_rater) {
    return(expand.grid(Person = persons, Rater = raters, stringsAsFactors = FALSE))
  }

  if (assignment == "rotating") {
    m <- max(1L, min(as.integer(raters_per_person), n_rater))
    rows <- lapply(seq_len(n_person), function(i) {
      # Step by 1 across raters (cyclic spiral) so m >= 2 forms a connected ring even when gcd(m, n_rater) > 1
      idx <- ((i - 1L) + seq_len(m) - 1L) %% n_rater + 1L
      data.frame(Person = persons[i], Rater = raters[idx], stringsAsFactors = FALSE)
    })
    return(do.call(rbind, rows))
  }

  k <- max(0L, min(as.integer(anchor_persons), n_person))
  p1 <- raters[seq_len(max(1L, floor(n_rater / 2)))]
  p2 <- setdiff(raters, p1)
  if (length(p2) == 0L) p2 <- p1

  out <- list()
  if (k > 0L) {
    out[[length(out) + 1L]] <- expand.grid(
      Person = persons[seq_len(k)], Rater = raters, stringsAsFactors = FALSE
    )
  }
  rem <- setdiff(seq_len(n_person), seq_len(k))
  if (length(rem) > 0L) {
    half <- ceiling(length(rem) / 2)
    rem1 <- rem[seq_len(half)]
    rem2 <- setdiff(rem, rem1)
    if (length(rem1) > 0L) {
      out[[length(out) + 1L]] <- expand.grid(Person = persons[rem1], Rater = p1, stringsAsFactors = FALSE)
    }
    if (length(rem2) > 0L) {
      out[[length(out) + 1L]] <- expand.grid(Person = persons[rem2], Rater = p2, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

build_candidate_grid <- function(p_vals, r_vals, assignment = "rotating",
                                 raters_per_person = 2L, anchor_persons = 6L,
                                 extra_facets = data.frame(name = "Criterion", levels = 4L, sd = 0.30),
                                 n_cat = 5L) {
  ef <- normalize_extra_facets(extra_facets)
  grid <- expand.grid(
    n_person = sort(unique(as.integer(p_vals))),
    n_rater = sort(unique(as.integer(r_vals))),
    stringsAsFactors = FALSE
  )
  grid$assignment <- assignment
  grid$raters_per_person <- ifelse(assignment == "crossed", grid$n_rater, pmin(as.integer(raters_per_person), grid$n_rater))
  grid$anchor_persons <- ifelse(assignment == "anchor", pmin(as.integer(anchor_persons), grid$n_person), NA_integer_)
  grid$n_cat <- as.integer(n_cat)

  cell_mult <- if (nrow(ef) == 0L) 1L else as.integer(prod(ef$levels))
  pr_counts <- vapply(seq_len(nrow(grid)), function(i) {
    nrow(build_person_rater_pairs(
      n_person = grid$n_person[i],
      n_rater = grid$n_rater[i],
      assignment = grid$assignment[i],
      raters_per_person = grid$raters_per_person[i],
      anchor_persons = grid$anchor_persons[i] %||% 0L
    ))
  }, integer(1))

  grid$cell_multiplier <- cell_mult
  grid$total_ratings <- pr_counts * cell_mult
  grid$ratings_per_rater <- round(grid$total_ratings / grid$n_rater, 1)
  grid$ratings_per_person <- round(grid$total_ratings / grid$n_person, 1)
  grid$plan_id <- seq_len(nrow(grid))
  grid$plan_label <- vapply(seq_len(nrow(grid)), function(i) {
    assign_tag <- switch(
      grid$assignment[i],
      crossed  = sprintf("all %dR", grid$n_rater[i]),
      rotating = sprintf("%d/%dR", grid$raters_per_person[i], grid$n_rater[i]),
      anchor   = sprintf("%dR (k=%d)", grid$n_rater[i], grid$anchor_persons[i])
    )
    sprintf("P=%d, %s", grid$n_person[i], assign_tag)
  }, character(1))
  attr(grid, "extra_facets") <- ef
  grid
}

simulate_mfrm_dataset <- function(row_cfg, extra_facets = attr(row_cfg, "extra_facets"),
                                  model = c("RSM", "PCM"), step_facet = NULL,
                                  sd_person = 1.0, sd_rater = 0.40,
                                  step_span = 2.4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  model <- match.arg(model)
  ef <- normalize_extra_facets(extra_facets)

  n_p <- row_cfg$n_person
  n_r <- row_cfg$n_rater
  n_k <- row_cfg$n_cat
  persons <- sprintf("P%03d", seq_len(n_p))
  raters  <- sprintf("R%02d", seq_len(n_r))

  center_vec <- function(n, s) {
    v <- stats::rnorm(n, 0, s)
    v - mean(v)
  }

  theta <- stats::setNames(center_vec(n_p, sd_person), persons)
  rho   <- stats::setNames(center_vec(n_r, sd_rater), raters)

  dat <- build_person_rater_pairs(
    n_person = n_p,
    n_rater = n_r,
    assignment = row_cfg$assignment,
    raters_per_person = row_cfg$raters_per_person,
    anchor_persons = row_cfg$anchor_persons %||% 0L
  )

  truth <- list(person = theta, rater = rho, extra = list())
  for (j in seq_len(nrow(ef))) {
    fnm <- ef$name[j]
    flv <- sprintf("%s%02d", substr(fnm, 1L, 2L), seq_len(ef$levels[j]))
    f_eff <- stats::setNames(center_vec(length(flv), ef$sd[j]), flv)
    truth$extra[[fnm]] <- f_eff
    add_df <- stats::setNames(data.frame(flv, stringsAsFactors = FALSE), fnm)
    dat <- merge(dat, add_df, by = NULL)
  }

  eta <- theta[dat$Person] - rho[dat$Rater]
  for (fnm in names(truth$extra)) {
    eta <- eta - truth$extra[[fnm]][dat[[fnm]]]
  }

  all_facets <- c("Rater", ef$name)
  resolved_step <- if (model == "PCM") {
    if (!is.null(step_facet) && step_facet %in% all_facets) step_facet
    else if (nrow(ef) > 0L) ef$name[1L]
    else "Rater"
  } else {
    NULL
  }

  base_steps <- if (n_k == 2L) 0 else seq(-step_span / 2, step_span / 2, length.out = n_k - 1L)
  if (is.null(resolved_step)) {
    cum_tau <- matrix(c(0, cumsum(base_steps)), nrow = nrow(dat), ncol = n_k, byrow = TRUE)
  } else {
    sf_levels <- unique(dat[[resolved_step]])
    tau_by_lvl <- lapply(sf_levels, function(lv) {
      jit <- if (n_k == 2L) 0 else stats::rnorm(n_k - 1L, 0, 0.18)
      s <- base_steps + (jit - mean(jit))
      c(0, cumsum(s))
    })
    names(tau_by_lvl) <- sf_levels
    cum_tau <- do.call(rbind, tau_by_lvl[dat[[resolved_step]]])
  }

  k_seq <- matrix(0:(n_k - 1L), nrow = nrow(dat), ncol = n_k, byrow = TRUE)
  lp <- k_seq * eta - cum_tau
  lp <- lp - apply(lp, 1L, max)
  probs <- exp(lp)
  probs <- probs / rowSums(probs)
  u <- stats::runif(nrow(dat))
  cum_p <- t(apply(probs, 1L, cumsum))
  dat$Score <- rowSums(u > cum_p) + 1L

  list(
    data = dat,
    truth = truth,
    facets = all_facets,
    step_facet = resolved_step,
    n_cat = n_k
  )
}

evaluate_single_rep <- function(sim_obj, method = c("MML", "JML"), model = c("RSM", "PCM"), quad_points = 15L) {
  method <- match.arg(method)
  model  <- match.arg(model)
  dat    <- sim_obj$data
  n_k    <- sim_obj$n_cat

  if (length(unique(dat$Score)) < n_k) {
    return(data.frame(
      valid = FALSE, fail_reason = "unobserved_category",
      person_rel = NA_real_, person_rmse = NA_real_,
      rater_rmse = NA_real_, rater_cov95 = NA_real_, scale_ratio = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  fit <- tryCatch(
    suppressWarnings(mfrmr::fit_mfrm(
      data = dat,
      person = "Person",
      facets = sim_obj$facets,
      score = "Score",
      model = model,
      method = method,
      step_facet = sim_obj$step_facet,
      rating_min = 1L,
      rating_max = n_k,
      quad_points = as.integer(quad_points),
      maxit = 80L
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(
      valid = FALSE, fail_reason = "fit_error",
      person_rel = NA_real_, person_rmse = NA_real_,
      rater_rmse = NA_real_, rater_cov95 = NA_real_, scale_ratio = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  rd <- tryCatch(summary(fit)$readiness, error = function(e) NULL)
  fit_ok <- !is.null(rd) &&
    identical(as.character(rd$Status[rd$Domain == "Numerical"]), "pass") &&
    identical(as.character(rd$Status[rd$Domain == "Design"]), "pass_linked") &&
    identical(as.character(rd$Status[rd$Domain == "Stability"]), "pass")
  if (!fit_ok) {
    return(data.frame(
      valid = FALSE, fail_reason = "nonconverged_or_disconnected",
      person_rel = NA_real_, person_rmse = NA_real_,
      rater_rmse = NA_real_, rater_cov95 = NA_real_, scale_ratio = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  dg <- tryCatch(
    suppressWarnings(mfrmr::diagnose_mfrm(
      fit,
      interaction_pairs = list(),
      top_n_interactions = 0L,
      diagnostic_mode = "marginal_fit",
      residual_pca = "none"
    )),
    error = function(e) NULL
  )
  if (is.null(dg)) {
    return(data.frame(
      valid = FALSE, fail_reason = "diagnose_error",
      person_rel = NA_real_, person_rmse = NA_real_,
      rater_rmse = NA_real_, rater_cov95 = NA_real_, scale_ratio = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  p_tab <- fit$facets$person
  non_ext <- if ("Extreme" %in% names(p_tab)) p_tab$Extreme == "none" else rep(TRUE, nrow(p_tab))
  p_sub <- p_tab[non_ext & is.finite(p_tab$Estimate), , drop = FALSE]
  if (nrow(p_sub) >= 3L) {
    t_p <- sim_obj$truth$person[as.character(p_sub$Person)]
    e_p <- p_sub$Estimate - mean(p_sub$Estimate) + mean(t_p)
    person_rmse <- sqrt(mean((e_p - t_p)^2))
    person_rel  <- if (stats::sd(e_p) > 1e-8 && stats::sd(t_p) > 1e-8) stats::cor(e_p, t_p)^2 else NA_real_
  } else {
    person_rmse <- NA_real_
    person_rel  <- NA_real_
  }

  struct_t <- sim_obj$truth$rater[as.character(dg$measures$Level[dg$measures$Facet == "Rater"])]
  struct_e <- dg$measures$Estimate[dg$measures$Facet == "Rater"]
  struct_e <- struct_e - mean(struct_e, na.rm = TRUE)
  for (fnm in names(sim_obj$truth$extra)) {
    sub_m <- dg$measures[dg$measures$Facet == fnm, , drop = FALSE]
    if (nrow(sub_m) > 0L) {
      struct_e <- c(struct_e, sub_m$Estimate - mean(sub_m$Estimate, na.rm = TRUE))
      struct_t <- c(struct_t, sim_obj$truth$extra[[fnm]][as.character(sub_m$Level)])
    }
  }
  scale_ratio <- if (length(struct_t) >= 2L && stats::var(struct_t) > 1e-8) {
    as.numeric(stats::cov(struct_e, struct_t) / stats::var(struct_t))
  } else {
    NA_real_
  }

  r_truth <- sim_obj$truth$rater
  r_meas  <- dg$measures[dg$measures$Facet == "Rater", , drop = FALSE]
  r_est   <- r_meas$Estimate - mean(r_meas$Estimate, na.rm = TRUE)
  r_t     <- r_truth[as.character(r_meas$Level)]
  rater_rmse <- sqrt(mean((r_est - r_t)^2, na.rm = TRUE))

  rater_cov95 <- NA_real_
  if (method == "MML") {
    ci_obj <- tryCatch(suppressWarnings(mfrmr::mfrm_facet_intervals(fit, facet = "Rater")), error = function(e) NULL)
    if (!is.null(ci_obj) && !is.null(ci_obj$table) && all(ci_obj$table$CIEligible)) {
      ct <- ci_obj$table
      rt_ci <- r_truth[as.character(ct$Target)]
      rater_cov95 <- mean(ct$Lower <= rt_ci & rt_ci <= ct$Upper, na.rm = TRUE)
    }
  }
  if (is.na(rater_cov95) && all(is.finite(r_meas$SE))) {
    lo <- r_est - 1.96 * r_meas$SE
    hi <- r_est + 1.96 * r_meas$SE
    rater_cov95 <- mean(lo <= r_t & r_t <= hi, na.rm = TRUE)
  }

  data.frame(
    valid = TRUE, fail_reason = "ok",
    person_rel = person_rel, person_rmse = person_rmse,
    rater_rmse = rater_rmse, rater_cov95 = rater_cov95, scale_ratio = scale_ratio,
    stringsAsFactors = FALSE
  )
}

run_design_plan <- function(grid_df, methods = c("MML", "JML"), model = c("RSM", "PCM"),
                            step_facet = NULL, reps = 15L, sd_person = 1.0, sd_rater = 0.40,
                            target_rel = 0.80, target_rater_rmse = 0.20,
                            max_fail_rate = 0.10, quad_points = 15L,
                            max_rows_per_rep = 25000L, seed = 2026L, progress_cb = NULL) {
  model <- match.arg(model)
  methods <- intersect(c("MML", "JML"), toupper(methods))
  if (length(methods) == 0L) stop("Select at least one estimator (MML or JML).", call. = FALSE)
  if (max(grid_df$total_ratings) > max_rows_per_rep) {
    stop(sprintf(
      "Largest candidate generates %d ratings/replicate (limit: %d). Reduce facet levels or use script mode.",
      max(grid_df$total_ratings), max_rows_per_rep
    ), call. = FALSE)
  }

  ef <- normalize_extra_facets(attr(grid_df, "extra_facets"))
  reps <- max(2L, as.integer(reps))

  # shortcut: serial replicate loop inside Shiny session, switch to future.apply if reps > 200
  total_steps <- nrow(grid_df) * reps * length(methods)
  step_i <- 0L
  rep_rows <- list()

  for (g in seq_len(nrow(grid_df))) {
    cfg <- grid_df[g, , drop = FALSE]
    for (r in seq_len(reps)) {
      sim_seed <- as.integer(seed) + g * 10000L + r
      sim_obj <- simulate_mfrm_dataset(
        cfg, extra_facets = ef, model = model, step_facet = step_facet,
        sd_person = sd_person, sd_rater = sd_rater, seed = sim_seed
      )
      for (meth in methods) {
        step_i <- step_i + 1L
        if (is.function(progress_cb)) {
          progress_cb(step_i / total_steps, sprintf("%s | %s (rep %d/%d)", cfg$plan_label, meth, r, reps))
        }
        ev <- evaluate_single_rep(sim_obj, method = meth, model = model, quad_points = quad_points)
        ev$plan_id <- cfg$plan_id
        ev$plan_label <- cfg$plan_label
        ev$method <- meth
        ev$rep <- r
        rep_rows[[length(rep_rows) + 1L]] <- ev
      }
    }
  }

  raw_df <- do.call(rbind, rep_rows)
  sum_rows <- lapply(split(raw_df, list(raw_df$plan_id, raw_df$method), drop = TRUE), function(d) {
    pid <- d$plan_id[1L]
    cfg <- grid_df[grid_df$plan_id == pid, , drop = FALSE]
    valid_rate <- mean(d$valid)
    vd <- d[d$valid, , drop = FALSE]
    m_or_na <- function(x) if (nrow(vd) == 0L || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

    p_rel  <- m_or_na(vd$person_rel)
    p_rmse <- m_or_na(vd$person_rmse)
    r_rmse <- m_or_na(vd$rater_rmse)
    r_cov  <- m_or_na(vd$rater_cov95)
    s_rat  <- m_or_na(vd$scale_ratio)

    reasons <- character()
    if ((1 - valid_rate) > max_fail_rate) reasons <- c(reasons, "Category/link risk")
    if (is.na(p_rel)  || p_rel < target_rel) reasons <- c(reasons, "Low person rel.")
    if (is.na(r_rmse) || r_rmse > target_rater_rmse) reasons <- c(reasons, "High rater RMSE")
    decision <- if (length(reasons) == 0L) "Pass" else paste(reasons, collapse = "; ")

    data.frame(
      plan_id = pid,
      Plan = cfg$plan_label,
      n_person = cfg$n_person,
      n_rater = cfg$n_rater,
      total_ratings = cfg$total_ratings,
      ratings_per_rater = cfg$ratings_per_rater,
      Estimator = d$method[1L],
      ValidPct = 100 * valid_rate,
      PersonRel = p_rel,
      PersonRMSE = p_rmse,
      RaterRMSE = r_rmse,
      RaterCov95 = 100 * r_cov,
      ScaleRatio = s_rat,
      Pass = (length(reasons) == 0L),
      Decision = decision,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, sum_rows)
  summary_df <- summary_df[order(-summary_df$Pass, summary_df$total_ratings, summary_df$RaterRMSE), , drop = FALSE]
  rownames(summary_df) <- NULL

  pass_rows <- summary_df[summary_df$Pass, , drop = FALSE]
  recommended <- if (nrow(pass_rows) > 0L) pass_rows[1L, , drop = FALSE] else NULL

  list(
    summary = summary_df,
    replicates = raw_df,
    grid = grid_df,
    extra_facets = ef,
    recommended = recommended,
    targets = list(
      target_rel = target_rel,
      target_rater_rmse = target_rater_rmse,
      max_fail_rate = max_fail_rate
    )
  )
}

format_summary_kable <- function(res_obj) {
  s <- res_obj$summary
  fmt_num <- function(x, d = 3) ifelse(is.na(x), "—", sprintf(paste0("%.", d, "f"), x))
  fmt_pct <- function(x) ifelse(is.na(x), "—", sprintf("%.0f%%", x))

  badge <- ifelse(
    s$Pass,
    "<span style='background:#dff5e8;color:#0b6b3a;padding:2px 8px;border-radius:10px;font-weight:600;'>Pass</span>",
    sprintf("<span style='background:#fff0e6;color:#9a3400;padding:2px 8px;border-radius:10px;'>%s</span>", s$Decision)
  )

  disp <- data.frame(
    Plan = s$Plan,
    Estimator = s$Estimator,
    `Ratings (Per Rater)` = sprintf("%d (%.0f)", s$total_ratings, s$ratings_per_rater),
    `Valid` = fmt_pct(s$ValidPct),
    `Person r²` = fmt_num(s$PersonRel, 3),
    `Person RMSE` = fmt_num(s$PersonRMSE, 3),
    `Rater RMSE` = fmt_num(s$RaterRMSE, 3),
    `Rater 95% CI` = fmt_pct(s$RaterCov95),
    `Scale Ratio` = fmt_num(s$ScaleRatio, 2),
    Decision = badge,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  kb <- knitr::kable(
    disp,
    format = "html",
    escape = FALSE,
    align = c("l", "c", "r", "r", "r", "r", "r", "r", "r", "l")
  )
  if (requireNamespace("kableExtra", quietly = TRUE)) {
    kb <- kableExtra::kable_styling(
      kb,
      bootstrap_options = c("hover", "condensed"),
      full_width = TRUE,
      font_size = 13.5
    )
    kb <- kableExtra::add_header_above(
      kb,
      c("Design & Workload" = 3, "Examinee Precision" = 3, "Rater Calibration" = 3, "Status" = 1)
    )
    if (!is.null(res_obj$recommended)) {
      kb <- kableExtra::row_spec(kb, 1L, bold = TRUE, background = "#f4fbf7")
    }
  }
  as.character(kb)
}

generate_reproducible_code <- function(p_vals, n_rater, assignment, raters_per_person, anchor_persons,
                                       extra_facets, n_cat, methods, model, step_facet, reps,
                                       sd_person, sd_rater, target_rel, target_rater_rmse, seed) {
  ef <- normalize_extra_facets(extra_facets)
  ef_code <- if (nrow(ef) == 0L) {
    "data.frame(name = character(0), levels = integer(0), sd = numeric(0))"
  } else {
    sprintf(
      "data.frame(\n    name = c(%s),\n    levels = c(%s),\n    sd = c(%s)\n  )",
      paste(sprintf('"%s"', ef$name), collapse = ", "),
      paste(paste0(ef$levels, "L"), collapse = ", "),
      paste(sprintf("%.2f", ef$sd), collapse = ", ")
    )
  }
  step_arg <- if (identical(model, "PCM") && !is.null(step_facet)) sprintf('"%s"', step_facet) else "NULL"

  sprintf(
'library(mfrmr)
source("R/planner_core.R")

grid <- build_candidate_grid(
  p_vals = c(%s),
  r_vals = %dL,
  assignment = "%s",
  raters_per_person = %dL,
  anchor_persons = %dL,
  extra_facets = %s,
  n_cat = %dL
)

plan <- run_design_plan(
  grid_df = grid,
  methods = c(%s),
  model = "%s",
  step_facet = %s,
  reps = %dL,
  sd_person = %.2f,
  sd_rater = %.2f,
  target_rel = %.2f,
  target_rater_rmse = %.2f,
  seed = %dL
)

print(plan$summary)
',
    paste(paste0(p_vals, "L"), collapse = ", "),
    as.integer(n_rater), assignment, as.integer(raters_per_person), as.integer(anchor_persons),
    ef_code, as.integer(n_cat),
    paste(sprintf('"%s"', methods), collapse = ", "),
    model, step_arg, as.integer(reps), sd_person, sd_rater,
    target_rel, target_rater_rmse, as.integer(seed)
  )
}
