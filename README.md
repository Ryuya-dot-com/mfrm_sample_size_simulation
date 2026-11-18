# MFRM Sample-Size & Rating-Design Planner

An interactive R Shiny application and scriptable simulation engine for planning **Many-Facet Rasch Measurement (MFRM)** performance assessments (e.g., L2 writing and speaking tests).

Unlike heuristic rules of thumb or cumulative-link (`clmm`) approximations, this planner generates and fits **adjacent-category Rasch models (`RSM` and `PCM`)** directly via the **`mfrmr`** engine under both **Marginal Maximum Likelihood (`MML`)** and **Joint Maximum Likelihood (`JML` / `FACETS`-type)** across **2, 3, 4, or more facets**.

---

## 1. Installation (Including Pre-CRAN `mfrmr` Setup)

### Step 1: Install required CRAN packages

```r
install.packages(c("shiny", "ggplot2", "plotly", "knitr", "kableExtra", "remotes"))
```

### Step 2: Install `mfrmr` (GitHub latest release recommended)

While `mfrmr 0.3.0` is preparing for CRAN release, install the latest version directly from the [`mfrmr` GitHub repository](https://github.com/Ryuya-dot-com/mfrmr):

```r
remotes::install_github("Ryuya-dot-com/mfrmr")
```

*(Note: The planner is also backward-compatible with `mfrmr 0.2.4` on CRAN via `install.packages("mfrmr")`, though `>= 0.3.0` from GitHub is recommended for the latest readiness and covariance-interval checks.)*

### Step 3: Launch the Shiny Planner

**Option A — Launch in one line from GitHub (no manual download needed):**

```r
shiny::runGitHub("mfrm_sample_size_similation", "Ryuya-dot-com")
```

**Option B — Clone and run locally (or from an offline peer-review archive):**

```r
shiny::runApp("path/to/app.R")
```

---

## 2. Step-by-Step User Guide

### Panel 1: Sample Size & Rater Linking (Left Sidebar)
1. **`Examinees P (candidate range)` & `Candidate P steps`**:
   Drag the two-ended slider (e.g., `30` to `90`) and choose how many candidate points to evaluate (`1` to `5`). Setting `30–90` with `3` steps evaluates `P = 30, 60, 90` side by side.
2. **`Categories K`**:
   Set the number of ordered rating scale categories (`2` to `9`, e.g., `5` for a `1–5` rubric).
3. **`Total raters R` & `Rater assignment`**:
   Choose how raters are assigned to examinees:
   - **Balanced spiral rotation (`m` raters / examinee)**: Each examinee is scored by `m` raters in a connected cyclic ring (`R1-R2`, `R2-R3`, ..., `R6-R1`), keeping workload evenly distributed even when `gcd(m, R) > 1`. The slider for `m` automatically caps at `R`.
   - **Two panels linked by `k` shared anchor examinees**: Raters split into two panels, linked by `k` common anchor examinees scored by all `R` raters (with `k` automatically capped at `min(P)`). Setting `k = 0` lets you verify how `mfrmr` flags disconnected cohorts.
   - **Fully crossed**: Every rater scores every examinee across all active facets.

### Panel 2: Additional Facets (`+ Add` / `- Remove`)
- **`Person`** and **`Rater`** are always included (`2-facet` holistic design when `0` extra facets are active).
- Click **`+ Add`** or **`- Remove`** to configure **3-, 4-, 5-, or higher-facet designs** (defaulting to `Criterion`, `Task`, `Occasion`, `Domain`, or any custom name you type) with per-facet level counts (`2–12`) and latent standard deviations (`SD` in logits).
- **Live Workload Preview**: As you move any slider, the blue summary box updates immediately with the active facet structure, candidate `P` values, **total rating count**, and **ratings per rater**—without triggering slow recalculations until you click **`Run Simulation`**. *(Interactive runs cap at `25,000` ratings/replicate to keep the browser responsive; larger designs can be run via the generated R script.)*

### Panel 3: Estimator, Scale Model, & Decision Thresholds
- **`Estimator (mfrmr)`**: Check **`MML`**, **`JML`**, or **both** to compare marginal and joint likelihood side by side on the exact same simulated datasets.
- **`Scale model`**: Choose **`RSM`** (shared step thresholds across observations) or **`PCM`** (thresholds varying across a selected facet such as `Criterion` or `Task`).
- **`Target person reliability (r²)` & `Max. rater RMSE (logits)`**: Set the minimum acceptable measurement precision for your assessment purpose (e.g., `r² >= 0.80` and `Rater RMSE <= 0.20` logits).
- Expand **`Advanced (replications & Person/Rater SDs)`** to adjust Monte Carlo replications (`5–100`, default `15` for fast interactive exploration; increase to `50–100` to finalize a design), random seed, and generating `Person` / `Rater` SDs.

---

## 3. Interpreting the Outputs & Making Decisions

### Tab 1: `Decision Summary` (`knitr::kable` Grouped Table)
Unlike `AIC` rankings (which favor the smallest sample size because log-likelihood scales with `N`), candidates are ranked by **1) meeting all user thresholds (`Pass`) and 2) requiring the fewest total human ratings**. The green banner at the top highlights the **minimum-workload passing plan**.

| Column | Meaning & Decision Guidance |
|---|---|
| **`Plan` / `Estimator`** | Candidate sample size, rater-linking topology, and estimator (`MML` or `JML`). |
| **`Ratings (Per Rater)`** | Total ratings across the study and average ratings per rater (direct human scoring cost). |
| **`Valid`** | Percentage of replicates that observed all `K` score categories, satisfied graph connectivity (`pass_linked`), and converged cleanly. Aim for **`>= 90%–100%`**. |
| **`Person r²`** | Squared correlation $r^2(\hat{\theta}, \theta)$ between estimated and true examinee ability (immune to `JML` scale inflation). Driven primarily by **ratings per examinee** ($m \times C \times T$). |
| **`Person RMSE`** | Root mean square error (logits) of non-extreme examinee ability estimates. |
| **`Rater RMSE`** | Root mean square error (logits) of centered rater severity estimates. Driven by **examinees per rater** and **anchor overlap `k`**. |
| **`Rater 95% CI`** | Empirical 95% confidence interval coverage across raters (uses full model covariance `vcov` via `mfrm_facet_intervals()` for `MML`, and conditional `SE` for `JML`). Nominal target is **`~95%`**. |
| **`Scale Ratio`** | OLS slope of estimated on true non-person facet parameters (`1.00` = unbiased logit scale; values such as `1.20` reveal `+20%` small-sample outward expansion under `JML`). |
| **`Decision`** | **`Pass`** badge if `Valid >= 90%`, `Person r² >= target`, and `Rater RMSE <= target`; otherwise lists the exact shortfall (`Low person rel.`, `High rater RMSE`, or `Category/link risk`). |

### Tab 2: `Interactive Trade-off Curves` (`plotly`)
- Hover over any point to inspect exact workload and accuracy metrics without cluttering the figure canvas.
- Switch the lower panel between **`Rater RMSE`**, **`Rater 95% CI coverage`**, and **`Structural scale ratio`**, or click **`Download Publication PNG`** for a 300-dpi white-background `ggplot2` figure.

### Tab 3: `Reproducible R Code`
- Copy or download the standalone `.R` script reproducing your exact UI settings, ideal for archiving in preregistrations/supplements or increasing `reps` to `200–500` in batch mode.

---

## 4. Scripted Batch Usage (Without GUI)

```r
library(mfrmr)
source("R/planner_core.R")

grid <- build_candidate_grid(
  p_vals = c(30L, 60L, 90L),
  r_vals = 6L,
  assignment = "rotating",
  raters_per_person = 2L,
  extra_facets = data.frame(
    name = c("Criterion", "Task"),
    levels = c(4L, 2L),
    sd = c(0.30, 0.30)
  ),
  n_cat = 5L
)

plan <- run_design_plan(
  grid_df = grid,
  methods = c("MML", "JML"),
  model = "RSM",
  reps = 50L,
  target_rel = 0.80,
  target_rater_rmse = 0.20,
  seed = 2026L
)

print(plan$summary)
```

---

## 5. Note for Double-Anonymized Peer Review

When citing or sharing this tool in a double-anonymized manuscript submission:
1. Bundle `app.R` and `R/planner_core.R` inside your anonymized OSF / Zenodo supplement (e.g., `sample_size_planner/`) rather than linking directly to a personal GitHub username in the blinded text, or use masked link text pointing to an anonymized mirror / view-only OSF link (`[MFRM Sample-Size & Rating-Design Planner](<anonymous-view-only-url>)`).
2. Reviewers can run the bundled copy offline with `shiny::runApp("sample_size_planner")`. Upon acceptance, replace the anonymized placeholder with the permanent repository link and Zenodo DOI.

---

## 6. 日本語クイックガイド（使い方要約）

1. **事前準備（GitHub 版 `mfrmr` の導入）**：
   ```r
   install.packages(c("shiny", "ggplot2", "plotly", "knitr", "kableExtra", "remotes"))
   remotes::install_github("Ryuya-dot-com/mfrmr")
   ```
2. **アプリの起動**：
   `shiny::runGitHub("mfrm_sample_size_similation", "Ryuya-dot-com")`（ローカル実行時は `shiny::runApp("app.R")`）で起動します。
3. **操作の流れ**：
   - **左サイドバー「1. Sample Size & Rater Linking」** で受験者数レンジ `P`、段階数 `K`、評価者数 `R`、割当方式（巡回ローテーション `m` 名／共通アンカー `k` 名／完全交差）をスライダーで設定します。
   - **「2. Additional Facets」** の `+ Add` / `- Remove` で `Criterion`（観点）や `Task`（課題）などの相を自由に追加・削除します（青いボックスに延べ採点数と評価者1人あたり負担が即時表示されます）。
   - **「3. Estimator & Thresholds」** で `MML` / `JML`、`RSM` / `PCM`、目標とする受験者信頼性（`r²`）と評価者 `RMSE` 上限を設定し、**`Run Simulation`** を押します。
4. **結果の読み方**：
   `Decision Summary` タブの最上段に **「目標精度を満たす最小採点コストの推奨プラン」** が表示されます。`Interactive Trade-off Curves` タブでは `plotly` 上でカーソルを合わせて詳細値を確認でき、`Reproducible R Code` タブから再現用 R スクリプトを保存できます。

---

## Automated Verification

```bash
Rscript tests/test_planner_core.R
```

## License

MIT License. See [LICENSE](LICENSE).
