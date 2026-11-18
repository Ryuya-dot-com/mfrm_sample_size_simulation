#!/usr/bin/env Rscript
# Automated self-check for MFRM Sample-Size & Rating-Design Planner

library(mfrmr)
source("R/planner_core.R")

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop("ASSERTION FAILED: ", msg, call. = FALSE)
}

cat("1. Checking slider candidate generation & arbitrary K-facet workload (2, 3, 4, and 5 facets)...\n")
assert_true(identical(build_p_candidates(c(30, 90), 3), c(30L, 60L, 90L)), "Range slider candidate P values")

# 2-facet (0 extra)
g2 <- build_candidate_grid(30, 4, assignment = "crossed", extra_facets = NULL, n_cat = 4)
assert_true(g2$total_ratings == 120L, "2-facet crossed workload")

# 3-facet (1 extra: Criterion)
g3 <- build_candidate_grid(c(30, 60), 6, assignment = "rotating", raters_per_person = 2,
                           extra_facets = data.frame(name = "Criterion", levels = 4, sd = 0.3), n_cat = 5)
assert_true(identical(g3$total_ratings, c(240L, 480L)), "3-facet rotating workload")

# 5-facet (3 extra: Criterion x Task x Occasion)
ef5 <- data.frame(name = c("Criterion", "Task", "Occasion"), levels = c(3L, 2L, 2L), sd = c(0.3, 0.3, 0.25))
g5  <- build_candidate_grid(24, 4, assignment = "rotating", raters_per_person = 2, extra_facets = ef5, n_cat = 4)
# 24 * 2 = 48 PR pairs * (3 * 2 * 2 = 12) = 576 ratings
assert_true(g5$total_ratings == 576L, "5-facet workload calculation")

cat("2. Checking MML & JML execution across 2-, 3-, and 5-facet RSM & PCM...\n")
for (cfg in list(g2, g3[1, , drop = FALSE], g5)) {
  sim <- simulate_mfrm_dataset(cfg, extra_facets = attr(cfg, "extra_facets"), model = "RSM", seed = 101)
  ev_m <- evaluate_single_rep(sim, method = "MML", model = "RSM")
  ev_j <- evaluate_single_rep(sim, method = "JML", model = "RSM")
  assert_true(ev_m$valid && is.finite(ev_m$rater_rmse), "MML valid across arbitrary K facets")
  assert_true(ev_j$valid && is.finite(ev_j$rater_rmse), "JML valid across arbitrary K facets")
}
sim_pcm5 <- simulate_mfrm_dataset(g5, extra_facets = ef5, model = "PCM", step_facet = "Occasion", seed = 202)
ev_pcm5  <- evaluate_single_rep(sim_pcm5, method = "MML", model = "PCM")
assert_true(ev_pcm5$valid && is.finite(ev_pcm5$person_rel), "5-facet PCM with custom step_facet valid")

cat("3. Checking disconnected (k = 0) safeguard & row-explosion guardrail...\n")
g_disc <- build_candidate_grid(30, 4, assignment = "anchor", anchor_persons = 0,
                               extra_facets = data.frame(name = "Criterion", levels = 3, sd = 0.3), n_cat = 4)
sim_disc <- simulate_mfrm_dataset(g_disc, model = "RSM", seed = 303)
assert_true(!evaluate_single_rep(sim_disc, "JML", "RSM")$valid && !evaluate_single_rep(sim_disc, "MML", "RSM")$valid,
            "Disconnected k=0 design flagged")

g_huge <- build_candidate_grid(200, 20, assignment = "crossed",
                               extra_facets = data.frame(name = c("C", "T"), levels = c(5, 5), sd = c(0.3, 0.3)), n_cat = 5)
err_caught <- tryCatch({ run_design_plan(g_huge, reps = 2); FALSE }, error = function(e) grepl("limit", e$message))
assert_true(err_caught, "Row-explosion guardrail stops runaway grid")

cat("4. Checking kable HTML formatting, plotly widget, and MML vs. JML behavior...\n")
g_cmp <- build_candidate_grid(c(36, 72), 6, assignment = "rotating", raters_per_person = 2,
                              extra_facets = data.frame(name = "Criterion", levels = 3, sd = 0.3), n_cat = 4)
plan <- run_design_plan(grid_df = g_cmp, methods = c("MML", "JML"), model = "RSM",
                        reps = 6L, target_rel = 0.70, target_rater_rmse = 0.25, seed = 2026L)
kb_html <- format_summary_kable(plan)
assert_true(grepl("<table", kb_html) && grepl("Examinee Precision", kb_html), "kable grouped HTML table generated")

app_env <- new.env(parent = globalenv())
app_obj <- source("app.R", local = app_env)$value
assert_true(inherits(app_obj, "shiny.appobj"), "app.R returns shiny.appobj")
ply_obj <- app_env$build_tradeoff_plotly(plan, secondary_metric = "RaterRMSE")
assert_true(inherits(ply_obj, "plotly"), "build_tradeoff_plotly returns interactive plotly widget")

cat("5. Checking Shiny server reactivity (sliders + dynamic +Add/-Remove facets)...\n")
shiny::testServer(app_obj, {
  session$setInputs(
    p_range = c(24, 48), p_points = 2, n_cat = 4, n_rater = 4,
    assignment = "rotating", raters_per_person = 2, anchor_persons = 4,
    methods = c("MML", "JML"), model = "RSM",
    target_rel = 0.65, target_rater_rmse = 0.35,
    reps = 3, seed = 42, sd_person = 1.0, sd_rater = 0.4,
    sec_metric = "RaterRMSE"
  )
  # Add a 4th facet (Task), then verify 4-facet grid
  session$setInputs(add_facet = 1L)
  assert_true(nrow(extra_facets_df()) == 2L, "+ Add button adds a 4th facet")
  session$setInputs(run_btn = 1L)
  res <- sim_res()
  assert_true(nrow(res$summary) == 4L && nzchar(code_text()), "4-facet reactive simulation succeeded")
})

cat("ALL SELF-CHECKS PASSED CLEANLY!\n")
