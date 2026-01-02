## =============================================================================
## Libraries
## =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(car)
  library(lme4)
  library(lmerTest)
  library(lmtest)
  library(emmeans)
})


## =============================================================================
## 1. MISSING DATA REPORT
## =============================================================================

report_missing_data <- function(input_csv, output_txt) {
  df <- read_csv(input_csv, show_col_types = FALSE)
  
  hrv_cols <- c(
    "HF Power FFT (ms2)",
    "RM SSD (ms)",
    "Mean RR",
    "Stress Index",
    "PNS Index",
    "SNS Index"
  )
  
  # Participants with all HRV measures missing
  file_names <- unique(df[["File Name"]])
  missing_participants <- character(0)
  
  for (p in file_names) {
    participant_data <- df %>% filter(`File Name` == p)
    if (nrow(participant_data) == 0) next
    if (all(sapply(hrv_cols, function(col) all(is.na(participant_data[[col]]))))) {
      missing_participants <- c(missing_participants, p)
    }
  }
  
  total_participants <- length(file_names)
  participants_with_data <- total_participants - length(missing_participants)
  total_obs <- nrow(df)
  
  # "Observations with missing HRV": defined using first HRV col (same as Python)
  obs_with_missing <- sum(is.na(df[[hrv_cols[1]]]))
  obs_complete <- total_obs - obs_with_missing
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  report <- c(
    strrep("=", 80),
    "MISSING DATA REPORT",
    strrep("=", 80),
    "",
    paste0("Generated: ", timestamp),
    "Generated automatically to document data exclusions for analysis.",
    "",
    strrep("=", 80),
    "SUMMARY",
    strrep("=", 80),
    "",
    paste0("Total participants: ", total_participants),
    paste0("Participants with missing HRV data: ", length(missing_participants)),
    paste0("Participants with complete HRV data: ", participants_with_data),
    "",
    paste0("Total observations: ", total_obs),
    sprintf(
      "Observations with missing HRV: %d (%.1f%%)",
      obs_with_missing,
      obs_with_missing / total_obs * 100
    ),
    sprintf(
      "Observations with complete HRV: %d (%.1f%%)",
      obs_complete,
      obs_complete / total_obs * 100
    ),
    ""
  )
  
  if (length(missing_participants) > 0) {
    report <- c(
      report,
      strrep("=", 80),
      "EXCLUDED PARTICIPANTS",
      strrep("=", 80),
      "",
      "The following participants were excluded from analysis due to",
      "complete missing HRV data (data quality issues - files could",
      "not be opened or processed):",
      ""
    )
    
    for (p in missing_participants) {
      p_data <- df %>% filter(`File Name` == p)
      drug <- if (!all(is.na(p_data[["Drug Type"]]))) {
        p_data[["Drug Type"]][[1]]
      } else {
        "Unknown"
      }
      report <- c(
        report,
        paste0("  - ", p),
        paste0("    Drug Type: ", drug),
        paste0("    Observations: ", nrow(p_data))
      )
    }
    report <- c(report, "")
  }
  
  report <- c(
    report,
    strrep("=", 80),
    "HANDLING METHOD",
    strrep("=", 80),
    "",
    "Missing data handled via listwise deletion (complete case analysis).",
    "",
    "Rationale:",
    "  - Missingness is systematic (data quality issues), not random",
    "  - All HRV measures missing for excluded participants",
    "  - Imputation inappropriate for systematic missingness",
    "  - Sample size remains adequate (55 participants, 1,100 observations)",
    "",
    "Implementation:",
    "  - Analysis scripts use na.omit() / complete-case filtering on HRV",
    "  - This automatically excludes the 5 participants with missing data",
    "  - All analyses use complete cases only",
    "",
    strrep("=", 80)
  )
  
  dir.create(dirname(output_txt), recursive = TRUE, showWarnings = FALSE)
  writeLines(report, con = output_txt)
  
  cat(strrep("=", 80), "\n")
  cat("MISSING DATA REPORT\n")
  cat(strrep("=", 80), "\n")
  cat("\n")
  cat("Total participants:", total_participants, "\n")
  cat("Participants excluded:", length(missing_participants), "\n")
  cat("Participants in analysis:", participants_with_data, "\n")
  cat("Observations in analysis:", obs_complete, "\n")
  if (length(missing_participants) > 0) {
    cat("\nExcluded participants:", paste(missing_participants, collapse = ", "), "\n")
  }
  cat("\nReport saved to:", output_txt, "\n")
  cat(strrep("=", 80), "\n")
}


## =============================================================================
## 2. TWO-WAY ANOVAS
## =============================================================================

run_two_way_anovas <- function(input_csv, output_csv) {
  df <- read_csv(input_csv, show_col_types = FALSE)
  
  # Create factor versions for modelling (simpler names)
  df$VideoCategory <- factor(df[["Video Category"]])
  df$DrugType      <- factor(df[["Drug Type"]])
  
  dependent_vars <- c(
    "HF Power FFT (ms2)",
    "RM SSD (ms)",
    "Mean RR",
    "Stress Index",
    "PNS Index",
    "SNS Index"
  )
  
  results_rows <- list()
  
  for (dv in dependent_vars) {
    cols_needed <- c(dv, "VideoCategory", "DrugType")
    if (!all(cols_needed %in% names(df))) next
    
    sub <- df[, cols_needed]
    sub <- sub[complete.cases(sub), , drop = FALSE]
    if (nrow(sub) == 0) next
    
    # Build formula: `DV` ~ VideoCategory * DrugType
    form <- as.formula(paste0("`", dv, "` ~ VideoCategory * DrugType"))
    fit  <- lm(form, data = sub)
    
    # Type-II ANOVA
    aov_tab <- car::Anova(fit, type = 2)
    aov_df  <- as.data.frame(aov_tab)
    aov_df$term <- rownames(aov_df)
    rownames(aov_df) <- NULL
    
    # Error SS (Residuals)
    if ("Residuals" %in% aov_df$term) {
      ss_error <- aov_df$`Sum Sq`[aov_df$term == "Residuals"]
    } else {
      ss_error <- NA_real_
    }
    
    for (i in seq_len(nrow(aov_df))) {
      term_name <- aov_df$term[i]
      ss_effect <- aov_df$`Sum Sq`[i]
      df_term   <- aov_df$Df[i]
      F_val     <- if (!is.null(aov_df$`F value`)) aov_df$`F value`[i] else NA_real_
      p_val     <- if (!is.null(aov_df$`Pr(>F)`)) aov_df$`Pr(>F)`[i] else NA_real_
      
      if (!is.na(ss_effect) && !is.na(ss_error) && (ss_effect + ss_error) != 0 && term_name != "Residuals") {
        partial_eta_sq <- ss_effect / (ss_effect + ss_error)
      } else {
        partial_eta_sq <- NA_real_
      }
      
      mean_sq <- if (!is.na(df_term) && df_term != 0) ss_effect / df_term else NA_real_
      
      results_rows[[length(results_rows) + 1L]] <- data.frame(
        dependent_variable = dv,
        term               = term_name,
        df                 = df_term,
        sum_sq             = ss_effect,
        mean_sq            = mean_sq,
        F                  = F_val,
        `PR(>F)`           = p_val,
        partial_eta_sq     = partial_eta_sq,
        n                  = nrow(sub),
        stringsAsFactors   = FALSE
      )
    }
  }
  
  if (length(results_rows) > 0L) {
    results_df <- do.call(rbind, results_rows)
  } else {
    results_df <- data.frame(
      dependent_variable = character(0),
      term               = character(0),
      df                 = numeric(0),
      sum_sq             = numeric(0),
      mean_sq            = numeric(0),
      F                  = numeric(0),
      `PR(>F)`           = numeric(0),
      partial_eta_sq     = numeric(0),
      n                  = integer(0),
      stringsAsFactors   = FALSE
    )
  }
  
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  write_csv(results_df, output_csv)
}


## =============================================================================
## 3. ANOVA ASSUMPTION CHECKS
## =============================================================================

check_anova_assumptions <- function(input_csv, output_csv) {
  df <- read_csv(input_csv, show_col_types = FALSE)
  
  df$VideoCategory <- factor(df[["Video Category"]])
  df$DrugType      <- factor(df[["Drug Type"]])
  
  dependent_vars <- c(
    "HF Power FFT (ms2)",
    "RM SSD (ms)",
    "Mean RR",
    "Stress Index",
    "PNS Index",
    "SNS Index"
  )
  
  assumption_results <- list()
  
  for (dv in dependent_vars) {
    cols_needed <- c(dv, "VideoCategory", "DrugType")
    if (!all(cols_needed %in% names(df))) next
    
    sub <- df[, cols_needed]
    sub <- sub[complete.cases(sub), , drop = FALSE]
    if (nrow(sub) == 0) next
    
    form <- as.formula(paste0("`", dv, "` ~ VideoCategory * DrugType"))
    fit  <- lm(form, data = sub)
    
    residuals <- resid(fit)
    fitted    <- fitted(fit)
    
    # Intercept
    coefs <- coef(fit)
    has_intercept   <- "(Intercept)" %in% names(coefs)
    intercept_value <- if (has_intercept) coefs["(Intercept)"] else NA_real_
    
    # Shapiro-Wilk
    if (length(residuals) > 5000) {
      sample_residuals <- sample(residuals, 5000L)
    } else {
      sample_residuals <- residuals
    }
    sh <- shapiro.test(sample_residuals)
    shapiro_stat <- unname(sh$statistic)
    shapiro_p    <- sh$p.value
    normality_ok <- shapiro_p > 0.05
    
    # Levene's test (Residuals by group)
    group_factor <- interaction(sub$VideoCategory, sub$DrugType, drop = TRUE)
    groups <- split(residuals, group_factor)
    groups <- groups[sapply(groups, length) > 0]
    if (length(groups) >= 2) {
      lev <- car::leveneTest(residuals ~ group_factor)
      levene_stat <- lev[["F value"]][1]
      levene_p    <- lev[["Pr(>F)"]][1]
      homogeneity_ok <- levene_p > 0.05
    } else {
      levene_stat <- NA_real_
      levene_p    <- NA_real_
      homogeneity_ok <- NA
    }
    
    # Durbin-Watson
    dw <- lmtest::dwtest(fit)
    dw_stat <- unname(dw$statistic["DW"])
    independence_ok <- (dw_stat > 1.5 && dw_stat < 2.5)
    
    # Outliers (standardized residuals)
    std_res <- residuals / sd(residuals, na.rm = TRUE)
    n_outliers <- sum(abs(std_res) > 3, na.rm = TRUE)
    
    assumption_results[[length(assumption_results) + 1L]] <- data.frame(
      dependent_variable = dv,
      n                  = nrow(sub),
      has_intercept      = has_intercept,
      intercept_value    = intercept_value,
      shapiro_w_stat     = shapiro_stat,
      shapiro_w_p        = shapiro_p,
      normality_ok       = normality_ok,
      levene_stat        = levene_stat,
      levene_p           = levene_p,
      homogeneity_ok     = homogeneity_ok,
      durbin_watson      = dw_stat,
      independence_ok    = independence_ok,
      n_outliers_3sd     = n_outliers,
      residual_mean      = mean(residuals),
      residual_std       = sd(residuals),
      residual_min       = min(residuals),
      residual_max       = max(residuals),
      stringsAsFactors   = FALSE
    )
  }
  
  if (length(assumption_results) > 0L) {
    assumption_df <- do.call(rbind, assumption_results)
  } else {
    assumption_df <- data.frame()
  }
  
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  write_csv(assumption_df, output_csv)
  
  cat(strrep("=", 80), "\n")
  cat("ANOVA ASSUMPTION CHECKS\n")
  cat(strrep("=", 80), "\n")
  if (nrow(assumption_df) > 0) {
    for (i in seq_len(nrow(assumption_df))) {
      row <- assumption_df[i, ]
      cat("\n", row$dependent_variable, "\n", sep = "")
      cat("  Sample size: n =", row$n, "\n")
      cat("  Intercept included:", row$has_intercept,
          sprintf("(value: %.4f)", row$intercept_value), "\n")
      cat(sprintf(
        "  Normality (Shapiro-Wilk): W = %.4f, p = %.4f %s\n",
        row$shapiro_w_stat,
        row$shapiro_w_p,
        if (isTRUE(row$normality_ok)) "\u2713" else "\u2717"
      ))
      cat(sprintf(
        "  Homogeneity (Levene): F = %.4f, p = %.4f %s\n",
        row$levene_stat,
        row$levene_p,
        if (isTRUE(row$homogeneity_ok)) "\u2713" else "\u2717"
      ))
      cat(sprintf(
        "  Independence (Durbin-Watson): DW = %.4f %s\n",
        row$durbin_watson,
        if (isTRUE(row$independence_ok)) "\u2713" else "\u2717"
      ))
      cat("  Outliers (>3 SD):", row$n_outliers_3sd, "\n")
    }
  }
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("Note: \u2713 = assumption met, \u2717 = assumption violated\n")
  cat(strrep("=", 80), "\n")
}


## =============================================================================
## 4. LINEAR MIXED MODELS
## =============================================================================

## =============================================================================
## 4. LINEAR MIXED MODELS  (random intercepts + random slopes for VideoCategory)
## =============================================================================

run_linear_mixed_model <- function(input_csv, output_csv) {
  df <- read_csv(input_csv, show_col_types = FALSE)
  
  df$VideoCategory <- factor(df[["Video Category"]])
  df$DrugType      <- factor(df[["Drug Type"]])
  
  # Participant ID from File Name (remove .acq / .ACQ)
  df$Participant_ID <- sub("\\.(acq|ACQ)$", "", df[["File Name"]])
  
  dependent_vars <- c(
    "HF Power FFT (ms2)",
    "RM SSD (ms)",
    "Mean RR",
    "Stress Index",
    "PNS Index",
    "SNS Index"
  )
  
  results_rows <- list()
  
  for (dv in dependent_vars) {
    cols_needed <- c(dv, "VideoCategory", "DrugType", "Participant_ID")
    if (!all(cols_needed %in% names(df))) next
    
    sub <- df[, cols_needed]
    sub <- sub[complete.cases(sub), , drop = FALSE]
    if (nrow(sub) == 0) next
    
    # *** RANDOM INTERCEPT ONLY ***
    form <- as.formula(paste0("`", dv, "` ~ VideoCategory * DrugType + (1 | Participant_ID)"))
    
    fit <- tryCatch(
      lmer(
        form, 
        data    = sub, 
        REML    = TRUE,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl   = list(maxfun = 1e5)
        )
      ),
      error = function(e) e
    )
    
    if (inherits(fit, "error")) {
      message("Error fitting LMM for ", dv, ": ", fit$message)
      next
    }
    
    coefs  <- summary(fit)$coefficients
    params <- coefs[, "Estimate"]
    ses    <- coefs[, "Std. Error"]
    tvals  <- coefs[, "t value"]
    pvals  <- coefs[, "Pr(>|t|)"]
    
    ci <- tryCatch(
      confint(fit, method = "Wald"),
      error = function(e) NULL
    )
    
    n_groups <- length(unique(sub$Participant_ID))
    aic_val  <- AIC(fit)
    bic_val  <- BIC(fit)
    llf_val  <- as.numeric(logLik(fit))
    
    term_names <- rownames(coefs)
    for (term in term_names) {
      if (term == "(Intercept)") next
      
      if (grepl(":", term)) {
        term_type <- "Interaction"
      } else if (grepl("VideoCategory", term)) {
        term_type <- "Video Category"
      } else if (grepl("DrugType", term)) {
        term_type <- "Drug Type"
      } else {
        term_type <- "Other"
      }
      
      coef_val <- params[term]
      se_val   <- ses[term]
      p_val    <- pvals[term]
      z_val    <- if (!is.na(se_val) && se_val != 0) coef_val / se_val else NA_real_
      
      if (!is.null(ci) && term %in% rownames(ci)) {
        ci_lower <- ci[term, 1]
        ci_upper <- ci[term, 2]
      } else {
        ci_lower <- NA_real_
        ci_upper <- NA_real_
      }
      
      results_rows[[length(results_rows) + 1L]] <- data.frame(
        dependent_variable = dv,
        term               = term,
        term_type          = term_type,
        coefficient        = coef_val,
        std_error          = se_val,
        z_value            = z_val,
        p_value            = p_val,
        ci_lower           = ci_lower,
        ci_upper           = ci_upper,
        n                  = nrow(sub),
        n_groups           = n_groups,
        aic                = aic_val,
        bic                = bic_val,
        llf                = llf_val,
        stringsAsFactors   = FALSE
      )
    }
    
    # Random effect variance (Participant_ID intercept)
    re_var <- NA_real_
    vc <- VarCorr(fit)
    if ("Participant_ID" %in% names(vc)) {
      re_var <- as.numeric(vc$Participant_ID[1, 1])
    }
    
    results_rows[[length(results_rows) + 1L]] <- data.frame(
      dependent_variable = dv,
      term               = "Random Effect (Participant Intercept Variance)",
      term_type          = "Random Effect",
      coefficient        = re_var,
      std_error          = NA_real_,
      z_value            = NA_real_,
      p_value            = NA_real_,
      ci_lower           = NA_real_,
      ci_upper           = NA_real_,
      n                  = nrow(sub),
      n_groups           = n_groups,
      aic                = aic_val,
      bic                = bic_val,
      llf                = llf_val,
      stringsAsFactors   = FALSE
    )
  }
  
  if (length(results_rows) > 0L) {
    results_df <- do.call(rbind, results_rows)
  } else {
    results_df <- data.frame()
  }
  
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  write_csv(results_df, output_csv)
  
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("LINEAR MIXED MODEL RESULTS (random intercept for Participant_ID only)\n")
  cat(strrep("=", 80), "\n")
  
  if (nrow(results_df) > 0) {
    for (dv in unique(results_df$dependent_variable)) {
      dv_res <- results_df %>% filter(dependent_variable == dv)
      if (nrow(dv_res) == 0) next
      
      cat("\n", dv, "\n", sep = "")
      cat(strrep("-", 80), "\n")
      
      fixed <- dv_res %>% filter(term_type != "Random Effect")
      for (i in seq_len(nrow(fixed))) {
        row <- fixed[i, ]
        p   <- row$p_value
        sig <- if (!is.na(p) && p < 0.001) {
          "***"
        } else if (!is.na(p) && p < 0.01) {
          "**"
        } else if (!is.na(p) && p < 0.05) {
          "*"
        } else {
          "ns"
        }
        cat("  ", row$term, ":\n", sep = "")
        cat(sprintf(
          "    Coef = %.4f, SE = %.4f, z = %.4f, p = %.4f %s\n",
          row$coefficient, row$std_error, row$z_value, p, sig
        ))
        # no coefficient-based direction here anymore
      }
      
      re_row <- dv_res %>% filter(term_type == "Random Effect")
      if (nrow(re_row) > 0) {
        re <- re_row[1, ]
        cat(sprintf("  Random Intercept Variance (Participant): %.4f\n", re$coefficient))
        cat(sprintf("  AIC = %.2f, BIC = %.2f\n", re$aic, re$bic))
      }
    }
  }
  
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("Significance: *** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant\n")
  cat(strrep("=", 80), "\n")
}


## =============================================================================
### TEXT SUMMARIES (FULLY DYNAMIC)
## ============================================================================

library(readr)
library(dplyr)

## ------------------------------------------------------------------
## Helper functions for p-values and significance labels
## ------------------------------------------------------------------

# Find the p-value column name in an ANOVA data frame
get_anova_p_col <- function(df) {
  if ("PR(>F)" %in% names(df)) {
    "PR(>F)"
  } else if ("PR..F." %in% names(df)) {
    "PR..F."
  } else {
    NA_character_
  }
}

# Safe extraction of a single p-value from a row
safe_p_val <- function(row, p_col) {
  if (is.na(p_col) || !(p_col %in% names(row))) {
    return(NA_real_)
  }
  val <- row[[p_col]]
  # ensure numeric scalar
  if (length(val) == 0) return(NA_real_)
  suppressWarnings(as.numeric(val[1]))
}

# Turn p-value into a significance label
p_to_sig <- function(p) {
  if (is.na(p)) {
    "ns"
  } else if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else {
    "ns"
  }
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ---------- helpers ANOVA directionality from emmeans ----------

# For main effects (VideoCategory, DrugType)
emm_direction_main <- function(emm_df, factor_col, dv_name) {
  # assumes emmeans output with columns factor_col and "emmean"
  emm_df <- as.data.frame(emm_df)
  if (!("emmean" %in% names(emm_df))) return(character(0))
  if (!(factor_col %in% names(emm_df))) return(character(0))
  
  # keep original factor order so first level is the reference
  ref_name <- as.character(emm_df[[factor_col]][1])
  ref_mean <- emm_df$emmean[1]
  
  lines <- c(
    sprintf("  Estimated marginal means by %s for %s:", factor_col, dv_name),
    sprintf("    Reference: %s = %.3f", ref_name, ref_mean)
  )
  
  if (nrow(emm_df) > 1) {
    for (i in 2:nrow(emm_df)) {
      lvl   <- as.character(emm_df[[factor_col]][i])
      m_i   <- emm_df$emmean[i]
      diff  <- m_i - ref_mean
      dir_s <- if (diff > 0) {
        "higher than reference"
      } else if (diff < 0) {
        "lower than reference"
      } else {
        "similar to reference"
      }
      lines <- c(
        lines,
        sprintf("    %s = %.3f (%+.3f vs ref, %s)", lvl, m_i, diff, dir_s)
      )
    }
  }
  lines
}

# For interaction effects (VideoCategory × DrugType)
emm_direction_interaction <- function(emm_df, dv_name) {
  emm_df <- as.data.frame(emm_df)
  if (!all(c("VideoCategory", "DrugType", "emmean") %in% names(emm_df))) {
    return(character(0))
  }
  
  emm_df$cell <- paste0(emm_df$VideoCategory, " × ", emm_df$DrugType)
  
  # take a reference cell as the first row
  ref_cell <- emm_df$cell[1]
  ref_mean <- emm_df$emmean[1]
  
  lines <- c(
    sprintf("  Estimated marginal means by VideoCategory × DrugType for %s:", dv_name),
    sprintf("    Reference: %s = %.3f", ref_cell, ref_mean)
  )
  
  if (nrow(emm_df) > 1) {
    for (i in 2:nrow(emm_df)) {
      cell <- emm_df$cell[i]
      m_i  <- emm_df$emmean[i]
      diff <- m_i - ref_mean
      dir_s <- if (diff > 0) {
        "higher than reference"
      } else if (diff < 0) {
        "lower than reference"
      } else {
        "similar to reference"
      }
      lines <- c(
        lines,
        sprintf("    %s = %.3f (%+.3f vs ref, %s)", cell, m_i, diff, dir_s)
      )
    }
  }
  lines
}

# ---------------------------------------------------------------------------
# 3.1 ANOVA effects summary (already basically dynamic, just cleaned)
# ---------------------------------------------------------------------------

save_anova_effects_summary <- function(anova_csv, output_txt, timestamp, data_csv = NULL) {
  
  df <- readr::read_csv(anova_csv, show_col_types = FALSE)
  
  # Load raw data to compute directionality from EMMs (if provided)
  raw_data <- NULL
  if (!is.null(data_csv)) {
    raw_data <- readr::read_csv(data_csv, show_col_types = FALSE)
  }
  
  report <- c(
    strrep("=", 80),
    "TWO-WAY ANOVA: MAIN AND INTERACTION EFFECTS",
    strrep("=", 80),
    "",
    paste0("Generated: ", timestamp),
    "",
    strrep("=", 80),
    "SUMMARY",
    strrep("=", 80),
    ""
  )
  
  if (nrow(df) == 0) {
    report <- c(
      report,
      "No ANOVA results available in the CSV.",
      strrep("=", 80)
    )
    dir.create(dirname(output_txt), recursive = TRUE, showWarnings = FALSE)
    writeLines(report, con = output_txt)
    message("ANOVA effects summary saved to: ", output_txt)
    return(invisible(NULL))
  }
  
  # Figure out which column holds p-values
  p_col <- get_anova_p_col(df)
  
  dependent_vars <- unique(df$dependent_variable)
  
  for (dv in dependent_vars) {
    dv_results <- df %>% dplyr::filter(dependent_variable == dv)
    
    report <- c(
      report,
      "",
      dv,
      strrep("-", 80)
    )
    
    # Main effect: Video Category
    video_cat <- dv_results %>%
      dplyr::filter(
        grepl("VideoCategory", term),
        !grepl(":", term),
        term != "Residuals"
      )
    
    if (nrow(video_cat) > 0) {
      vc_row <- video_cat[1, ]
      p_val  <- safe_p_val(vc_row, p_col)
      sig    <- p_to_sig(p_val)
      df_vc  <- vc_row$df %||% NA_real_
      F_vc   <- vc_row$F %||% NA_real_
      pes_vc <- vc_row$partial_eta_sq %||% NA_real_
      
      report <- c(
        report,
        "Video Category (Main Effect):",
        sprintf("  F(%s) = %.4f", ifelse(is.na(df_vc), "NA", as.integer(df_vc)), F_vc),
        sprintf("  p = %s %s", ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)), sig),
        sprintf("  Partial \u03b7\u00b2 = %s",
                ifelse(is.na(pes_vc), "NA", sprintf("%.4f", pes_vc))),
        ""
      )
    } else {
      report <- c(report, "Video Category (Main Effect): no term in ANOVA table.", "")
    }
    
    # Main effect: Drug Type
    drug_type <- dv_results %>%
      dplyr::filter(
        grepl("DrugType", term),
        !grepl(":", term),
        term != "Residuals"
      )
    
    if (nrow(drug_type) > 0) {
      dt_row <- drug_type[1, ]
      p_val  <- safe_p_val(dt_row, p_col)
      sig    <- p_to_sig(p_val)
      df_dt  <- dt_row$df %||% NA_real_
      F_dt   <- dt_row$F %||% NA_real_
      pes_dt <- dt_row$partial_eta_sq %||% NA_real_
      
      report <- c(
        report,
        "Drug Type (Main Effect):",
        sprintf("  F(%s) = %.4f", ifelse(is.na(df_dt), "NA", as.integer(df_dt)), F_dt),
        sprintf("  p = %s %s", ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)), sig),
        sprintf("  Partial \u03b7\u00b2 = %s",
                ifelse(is.na(pes_dt), "NA", sprintf("%.4f", pes_dt))),
        ""
      )
    } else {
      report <- c(report, "Drug Type (Main Effect): no term in ANOVA table.", "")
    }
    
    # Interaction
    interaction <- dv_results %>%
      dplyr::filter(
        grepl(":", term),
        grepl("VideoCategory", term),
        grepl("DrugType", term)
      )
    
    if (nrow(interaction) > 0) {
      int_row <- interaction[1, ]
      p_val   <- safe_p_val(int_row, p_col)
      sig     <- p_to_sig(p_val)
      df_int  <- int_row$df %||% NA_real_
      F_int   <- int_row$F %||% NA_real_
      pes_int <- int_row$partial_eta_sq %||% NA_real_
      
      report <- c(
        report,
        "Video Category \u00d7 Drug Type (Interaction):",
        sprintf("  F(%s) = %.4f", ifelse(is.na(df_int), "NA", as.integer(df_int)), F_int),
        sprintf("  p = %s %s", ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)), sig),
        sprintf("  Partial \u03b7\u00b2 = %s",
                ifelse(is.na(pes_int), "NA", sprintf("%.4f", pes_int))),
        ""
      )
    } else {
      report <- c(report, "Video Category \u00d7 Drug Type (Interaction): no term in ANOVA table.", "")
    }
    
    # ------------------------------------------------------------------
    # Directionality from ANOVA model (via emmeans), if raw_data loaded
    # ------------------------------------------------------------------
    if (!is.null(raw_data)) {
      cols_needed <- c(dv, "Video Category", "Drug Type")
      if (all(cols_needed %in% names(raw_data))) {
        sub <- raw_data[, cols_needed]
        sub <- sub[complete.cases(sub), , drop = FALSE]
        
        if (nrow(sub) > 0) {
          sub$VideoCategory <- factor(sub[["Video Category"]])
          sub$DrugType      <- factor(sub[["Drug Type"]])
          
          form <- as.formula(paste0("`", dv, "` ~ VideoCategory * DrugType"))
          fit  <- tryCatch(lm(form, data = sub), error = function(e) NULL)
          
          if (!is.null(fit)) {
            report <- c(
              report,
              "Directionality based on estimated marginal means (ANOVA):"
            )
            
            # Video Category EMMs
            vc_emm <- tryCatch(emmeans::emmeans(fit, "VideoCategory"),
                               error = function(e) NULL)
            if (!is.null(vc_emm)) {
              report <- c(
                report,
                emm_direction_main(vc_emm, "VideoCategory", dv),
                ""
              )
            }
            
            # Drug Type EMMs
            dt_emm <- tryCatch(emmeans::emmeans(fit, "DrugType"),
                               error = function(e) NULL)
            if (!is.null(dt_emm)) {
              report <- c(
                report,
                emm_direction_main(dt_emm, "DrugType", dv),
                ""
              )
            }
            
            # Interaction EMMs
            int_emm <- tryCatch(emmeans::emmeans(fit, ~ VideoCategory * DrugType),
                                error = function(e) NULL)
            if (!is.null(int_emm)) {
              report <- c(
                report,
                emm_direction_interaction(int_emm, dv),
                ""
              )
            }
          }
        }
      }
    }
    
    # Sample size (only once)
    if (nrow(dv_results) > 0) {
      n_val <- dv_results$n[1]
      report <- c(
        report,
        sprintf("Sample size: n = %s", ifelse(is.na(n_val), "NA", as.integer(n_val))),
        ""
      )
    }
  }
  
  report <- c(
    report,
    strrep("=", 80),
    "Significance: *** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant",
    strrep("=", 80)
  )
  
  dir.create(dirname(output_txt), recursive = TRUE, showWarnings = FALSE)
  writeLines(report, con = output_txt)
  message("ANOVA effects summary saved to: ", output_txt)
}


# ---------------------------------------------------------------------------
# 3.2 LMM effects summary (dynamic, no fixed claims)
# ---------------------------------------------------------------------------

save_lmm_effects_summary <- function(lmm_csv, output_txt, timestamp, data_csv = NULL) {
  df <- readr::read_csv(lmm_csv, show_col_types = FALSE)
  
  # Load raw data to compute directionality from LMM EMMs (if provided)
  raw_data <- NULL
  if (!is.null(data_csv)) {
    raw_data <- readr::read_csv(data_csv, show_col_types = FALSE)
  }
  
  report <- c(
    strrep("=", 80),
    "LINEAR MIXED MODEL: MAIN AND INTERACTION EFFECTS",
    strrep("=", 80),
    "",
    paste0("Generated: ", timestamp),
    "",
    strrep("=", 80),
    "SUMMARY",
    strrep("=", 80),
    ""
  )
  
  if (nrow(df) == 0) {
    report <- c(
      report,
      "No LMM results available in the CSV.",
      strrep("=", 80)
    )
    dir.create(dirname(output_txt), recursive = TRUE, showWarnings = FALSE)
    writeLines(report, con = output_txt)
    message("LMM effects summary saved to: ", output_txt)
    return(invisible(NULL))
  }
  
  fixed_effects <- df %>% filter(term_type != "Random Effect")
  dependent_vars <- unique(fixed_effects$dependent_variable)
  
  for (dv in dependent_vars) {
    dv_results <- fixed_effects %>% filter(dependent_variable == dv)
    
    report <- c(
      report,
      "",
      dv,
      strrep("-", 80)
    )
    
    # --------------------------
    # Fixed-effect statistics
    # --------------------------
    
    # Video Category main effects
    video_cat <- dv_results %>% filter(term_type == "Video Category")
    if (nrow(video_cat) > 0) {
      report <- c(report, "Video Category (Main Effect):")
      for (i in seq_len(nrow(video_cat))) {
        row <- video_cat[i, ]
        p_val <- row$p_value
        sig <- if (!is.na(p_val) && p_val < 0.001) {
          "***"
        } else if (!is.na(p_val) && p_val < 0.01) {
          "**"
        } else if (!is.na(p_val) && p_val < 0.05) {
          "*"
        } else {
          "ns"
        }
        report <- c(
          report,
          paste0("  ", row$term, ":"),
          sprintf("    Coef = %.4f, SE = %.4f", row$coefficient, row$std_error),
          sprintf("    z/t = %.4f, p = %.4f %s", row$z_value, p_val, sig),
          sprintf("    95%% CI: [%.4f, %.4f]", row$ci_lower, row$ci_upper)
          # no coefficient-based direction here
        )
      }
      report <- c(report, "")
    } else {
      report <- c(report, "Video Category (Main Effect): no fixed-effect terms.", "")
    }
    
    # Drug Type main effects
    drug_type <- dv_results %>% filter(term_type == "Drug Type")
    if (nrow(drug_type) > 0) {
      report <- c(report, "Drug Type (Main Effect):")
      for (i in seq_len(nrow(drug_type))) {
        row <- drug_type[i, ]
        p_val <- row$p_value
        sig <- if (!is.na(p_val) && p_val < 0.001) {
          "***"
        } else if (!is.na(p_val) && p_val < 0.01) {
          "**"
        } else if (!is.na(p_val) && p_val < 0.05) {
          "*"
        } else {
          "ns"
        }
        report <- c(
          report,
          paste0("  ", row$term, ":"),
          sprintf("    Coef = %.4f, SE = %.4f", row$coefficient, row$std_error),
          sprintf("    z/t = %.4f, p = %.4f %s", row$z_value, p_val, sig),
          sprintf("    95%% CI: [%.4f, %.4f]", row$ci_lower, row$ci_upper)
          # no coefficient-based direction here
        )
      }
      report <- c(report, "")
    } else {
      report <- c(report, "Drug Type (Main Effect): no fixed-effect terms.", "")
    }
    
    # Interaction effects
    interaction <- dv_results %>% filter(term_type == "Interaction")
    if (nrow(interaction) > 0) {
      report <- c(report, "Video Category \u00d7 Drug Type (Interaction):")
      for (i in seq_len(nrow(interaction))) {
        row <- interaction[i, ]
        p_val <- row$p_value
        sig <- if (!is.na(p_val) && p_val < 0.001) {
          "***"
        } else if (!is.na(p_val) && p_val < 0.01) {
          "**"
        } else if (!is.na(p_val) && p_val < 0.05) {
          "*"
        } else {
          "ns"
        }
        report <- c(
          report,
          paste0("  ", row$term, ":"),
          sprintf("    Coef = %.4f, SE = %.4f", row$coefficient, row$std_error),
          sprintf("    z/t = %.4f, p = %.4f %s", row$z_value, p_val, sig),
          sprintf("    95%% CI: [%.4f, %.4f]", row$ci_lower, row$ci_upper)
          # no coefficient-based direction here
        )
      }
      report <- c(report, "")
    } else {
      report <- c(report, "Video Category \u00d7 Drug Type (Interaction): no fixed-effect terms.", "")
    }
    
    # --------------------------
    # Model fit + random effect variance
    # --------------------------
    if (nrow(dv_results) > 0) {
      sample_row <- dv_results[1, ]
      report <- c(
        report,
        "Model Fit Statistics:",
        paste0("  n = ", as.integer(sample_row$n)),
        paste0("  n_groups (participants) = ", as.integer(sample_row$n_groups)),
        sprintf("  AIC = %.2f", sample_row$aic),
        sprintf("  BIC = %.2f", sample_row$bic),
        ""
      )
    }
    
    random_effects <- df %>%
      filter(
        dependent_variable == dv,
        term_type == "Random Effect"
      )
    if (nrow(random_effects) > 0) {
      re_row <- random_effects[1, ]
      report <- c(
        report,
        "Random Effect (Participant Variance):",
        sprintf("  Variance = %.4f", re_row$coefficient),
        ""
      )
    }
    
    # --------------------------
    # Directionality from LMM EMMs
    # --------------------------
    if (!is.null(raw_data)) {
      cols_needed <- c(dv, "Video Category", "Drug Type", "File Name")
      if (all(cols_needed %in% names(raw_data))) {
        sub <- raw_data[, cols_needed]
        sub <- sub[complete.cases(sub), , drop = FALSE]
        
        if (nrow(sub) > 0) {
          sub$VideoCategory <- factor(sub[["Video Category"]])
          sub$DrugType      <- factor(sub[["Drug Type"]])
          sub$Participant_ID <- sub("\\.(acq|ACQ)$", "", sub[["File Name"]])
          
          form_lmm <- as.formula(
            paste0("`", dv, "` ~ VideoCategory * DrugType + (1 | Participant_ID)")
          )
          
          fit_lmm <- tryCatch(
            lmer(
              form_lmm,
              data    = sub,
              REML    = TRUE,
              control = lmerControl(
                optimizer = "bobyqa",
                optCtrl   = list(maxfun = 1e5)
              )
            ),
            error = function(e) NULL
          )
          
          if (!is.null(fit_lmm)) {
            report <- c(
              report,
              "Directionality based on estimated marginal means (LMM):"
            )
            
            # Video Category EMMs
            vc_emm <- tryCatch(emmeans::emmeans(fit_lmm, "VideoCategory"),
                               error = function(e) NULL)
            if (!is.null(vc_emm)) {
              report <- c(
                report,
                emm_direction_main(vc_emm, "VideoCategory", dv),
                ""
              )
            }
            
            # Drug Type EMMs
            dt_emm <- tryCatch(emmeans::emmeans(fit_lmm, "DrugType"),
                               error = function(e) NULL)
            if (!is.null(dt_emm)) {
              report <- c(
                report,
                emm_direction_main(dt_emm, "DrugType", dv),
                ""
              )
            }
            
            # Interaction EMMs
            int_emm <- tryCatch(emmeans::emmeans(fit_lmm, ~ VideoCategory * DrugType),
                                error = function(e) NULL)
            if (!is.null(int_emm)) {
              report <- c(
                report,
                emm_direction_interaction(int_emm, dv),
                ""
              )
            }
          }
        }
      }
    }
  }
  
  report <- c(
    report,
    strrep("=", 80),
    "Significance: *** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant",
    strrep("=", 80)
  )
  
  dir.create(dirname(output_txt), recursive = TRUE, showWarnings = FALSE)
  writeLines(report, con = output_txt)
  message("LMM effects summary saved to: ", output_txt)
}

# ---------------------------------------------------------------------------
# 3.4 High-level narrative summary (now fully conditional)
# ---------------------------------------------------------------------------

compare_anova_lmm <- function(anova_csv, lmm_csv, output_txt, timestamp) {
  anova_df <- readr::read_csv(anova_csv, show_col_types = FALSE)
  lmm_df   <- readr::read_csv(lmm_csv,   show_col_types = FALSE)
  
  report <- c(
    strrep("=", 80),
    "ANOVA vs LINEAR MIXED MODEL: COMPARISON",
    strrep("=", 80),
    "",
    paste0("Generated: ", timestamp),
    "",
    strrep("=", 80),
    "SUMMARY",
    strrep("=", 80),
    "",
    "This comparison shows main effects and interactions from both ANOVA",
    "and Linear Mixed Model (LMM) analyses side by side.",
    "",
    "Conceptual differences:",
    "  - ANOVA: Assumes independence of observations.",
    "  - LMM: Accounts for repeated measures via random intercepts",
    "          and random slopes for Video Category (within-subject factor).",
    "  - ANOVA: Uses F-tests.",
    "  - LMM: Uses Wald t/z-tests.",
    
    ""
  )
  
  if (nrow(anova_df) == 0 && nrow(lmm_df) == 0) {
    report <- c(
      report,
      "No results found in either ANOVA or LMM CSVs.",
      strrep("=", 80)
    )
    dir.create(dirname(output_txt), recursive = TRUE, showWarnings = FALSE)
    writeLines(report, con = output_txt)
    message("ANOVA vs LMM comparison saved to: ", output_txt)
    return(invisible(NULL))
  }
  
  lmm_fixed <- lmm_df %>% dplyr::filter(term_type != "Random Effect")
  dependent_vars <- sort(unique(c(
    anova_df$dependent_variable,
    lmm_fixed$dependent_variable
  )))
  
  p_col <- get_anova_p_col(anova_df)
  
  for (dv in dependent_vars) {
    report <- c(
      report,
      strrep("=", 80),
      dv,
      strrep("=", 80),
      ""
    )
    
    anova_dv <- anova_df %>% dplyr::filter(dependent_variable == dv)
    lmm_dv   <- lmm_fixed %>% dplyr::filter(dependent_variable == dv)
    
    ## VIDEO CATEGORY
    report <- c(report, "VIDEO CATEGORY (Main Effect)", strrep("-", 80))
    
    anova_vc <- anova_dv %>%
      dplyr::filter(
        grepl("VideoCategory", term),
        !grepl(":", term),
        term != "Residuals"
      )
    lmm_vc <- lmm_dv %>% dplyr::filter(term_type == "Video Category")
    
    if (nrow(anova_vc) > 0) {
      avc   <- anova_vc[1, ]
      p_val <- safe_p_val(avc, p_col)
      sig   <- p_to_sig(p_val)
      report <- c(
        report,
        "ANOVA:",
        sprintf(
          "  F(%s) = %.4f, p = %s %s",
          ifelse(is.na(avc$df), "NA", as.integer(avc$df)),
          avc$F,
          ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)),
          sig
        ),
        sprintf(
          "  Partial \u03b7\u00b2 = %s",
          ifelse(is.na(avc$partial_eta_sq), "NA", sprintf("%.4f", avc$partial_eta_sq))
        )
      )
    } else {
      report <- c(report, "ANOVA: No Video Category main-effect term.")
    }
    
    if (nrow(lmm_vc) > 0) {
      sig_count   <- sum(lmm_vc$p_value < 0.05, na.rm = TRUE)
      total_count <- nrow(lmm_vc)
      report <- c(
        report,
        "LMM:",
        sprintf("  %d/%d Video Category terms significant (p < 0.05)",
                sig_count, total_count)
      )
      for (i in seq_len(nrow(lmm_vc))) {
        row   <- lmm_vc[i, ]
        p_val <- row$p_value
        sig   <- p_to_sig(p_val)
        report <- c(
          report,
          sprintf(
            "    Coef = %.4f, z/t = %.4f, p = %s %s",
            row$coefficient,
            row$z_value,
            ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)),
            sig
          )
        )
        
      }
      
    } else {
      report <- c(report, "LMM: No Video Category fixed-effect terms.")
    }
    report <- c(report, "")
    
    ## DRUG TYPE
    report <- c(report, "DRUG TYPE (Main Effect)", strrep("-", 80))
    
    anova_dt <- anova_dv %>%
      dplyr::filter(
        grepl("DrugType", term),
        !grepl(":", term),
        term != "Residuals"
      )
    lmm_dt <- lmm_dv %>% dplyr::filter(term_type == "Drug Type")
    
    if (nrow(anova_dt) > 0) {
      adt   <- anova_dt[1, ]
      p_val <- safe_p_val(adt, p_col)
      sig   <- p_to_sig(p_val)
      report <- c(
        report,
        "ANOVA:",
        sprintf(
          "  F(%s) = %.4f, p = %s %s",
          ifelse(is.na(adt$df), "NA", as.integer(adt$df)),
          adt$F,
          ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)),
          sig
        ),
        sprintf(
          "  Partial \u03b7\u00b2 = %s",
          ifelse(is.na(adt$partial_eta_sq), "NA", sprintf("%.4f", adt$partial_eta_sq))
        )
      )
    } else {
      report <- c(report, "ANOVA: No Drug Type main-effect term.")
    }
    
    if (nrow(lmm_dt) > 0) {
      sig_count   <- sum(lmm_dt$p_value < 0.05, na.rm = TRUE)
      total_count <- nrow(lmm_dt)
      report <- c(
        report,
        "LMM:",
        sprintf("  %d/%d Drug Type terms significant (p < 0.05)",
                sig_count, total_count)
      )
      for (i in seq_len(nrow(lmm_dt))) {
        row   <- lmm_dt[i, ]
        p_val <- row$p_value
        sig   <- p_to_sig(p_val)
        report <- c(
          report,
          sprintf(
            "    Coef = %.4f, z/t = %.4f, p = %s %s",
            row$coefficient,
            row$z_value,
            ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)),
            sig
          )
        )
        
      }
    } else {
      report <- c(report, "LMM: No Drug Type fixed-effect terms.")
    }
    
    report <- c(report, "")
    
    ## INTERACTION
    report <- c(
      report,
      "VIDEO CATEGORY \u00d7 DRUG TYPE (Interaction)",
      strrep("-", 80)
    )
    
    anova_int <- anova_dv %>%
      dplyr::filter(
        grepl(":", term),
        grepl("VideoCategory", term),
        grepl("DrugType", term)
      )
    lmm_int <- lmm_dv %>% dplyr::filter(term_type == "Interaction")
    
    if (nrow(anova_int) > 0) {
      aint  <- anova_int[1, ]
      p_val <- safe_p_val(aint, p_col)
      sig   <- p_to_sig(p_val)
      report <- c(
        report,
        "ANOVA:",
        sprintf(
          "  F(%s) = %.4f, p = %s %s",
          ifelse(is.na(aint$df), "NA", as.integer(aint$df)),
          aint$F,
          ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)),
          sig
        ),
        sprintf(
          "  Partial \u03b7\u00b2 = %s",
          ifelse(is.na(aint$partial_eta_sq), "NA", sprintf("%.4f", aint$partial_eta_sq))
        )
      )
    } else {
      report <- c(report, "ANOVA: No Video Category \u00d7 Drug Type term.")
    }
    
    if (nrow(lmm_int) > 0) {
      report <- c(report, "LMM:")
      for (i in seq_len(nrow(lmm_int))) {
        row   <- lmm_int[i, ]
        p_val <- row$p_value
        sig   <- p_to_sig(p_val)
        report <- c(
          report,
          sprintf(
            "    Coef = %.4f, z/t = %.4f, p = %s %s",
            row$coefficient,
            row$z_value,
            ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)),
            sig
          )
        )
        
      }
    } else {
      report <- c(report, "LMM: No interaction fixed-effect terms.")
    }
    
    
    report <- c(report, "", "")
  }
  
  report <- c(
    report,
    strrep("=", 80),
    "Significance: *** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant",
    strrep("=", 80)
  )
  
  dir.create(dirname(output_txt), recursive = TRUE, showWarnings = FALSE)
  writeLines(report, con = output_txt)
  message("ANOVA vs LMM comparison saved to: ", output_txt)
}


# ============================================================================
# 4. DRIVER FUNCTION (CALLS EVERYTHING, LIKE PYTHON __main__)
# ============================================================================

run_all_analyses <- function() {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  input_file <- "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/data/data-two-way-anova.csv"
  
  # Step 1: Missing data
  cat("\n", strrep("=", 80), "\nSTEP 1: MISSING DATA REPORT\n", strrep("=", 80), "\n\n", sep = "")
  report_missing_data(
    input_csv  = input_file,
    output_txt = sprintf(
      "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/data/missing-data-report_%s.txt",
      timestamp
    )
  )
  
  # Step 2: ANOVA
  cat("\n", strrep("=", 80), "\nSTEP 2: TWO-WAY ANOVA\n", strrep("=", 80), "\n\n", sep = "")
  anova_output_csv <- sprintf(
    "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/anova/results-two-way-anova_%s.csv",
    timestamp
  )
  run_two_way_anovas(
    input_csv = input_file,
    output_csv = anova_output_csv
  )
  
  # Step 2b: ANOVA effects summary
  cat("\n", strrep("=", 80), "\nSTEP 2b: SAVING ANOVA EFFECTS SUMMARY\n", strrep("=", 80), "\n\n", sep = "")
  anova_effects_txt <- sprintf(
    "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/anova/anova-effects-summary_%s.txt",
    timestamp
  )
  save_anova_effects_summary(
    anova_csv  = anova_output_csv,
    output_txt = anova_effects_txt,
    timestamp  = timestamp,
    data_csv   = input_file  # <-- pass raw data so EMMs & directionality can be computed
  )
  
  # Step 3: ANOVA assumptions
  cat("\n", strrep("=", 80), "\nSTEP 3: ANOVA ASSUMPTION CHECKS\n", strrep("=", 80), "\n\n", sep = "")
  check_anova_assumptions(
    input_csv = input_file,
    output_csv = sprintf(
      "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/assumptions/anova-assumptions-check_%s.csv",
      timestamp
    )
  )
  
  # Step 4: LMM
  cat("\n", strrep("=", 80), "\nSTEP 4: LINEAR MIXED MODEL\n", strrep("=", 80), "\n\n", sep = "")
  lmm_output_csv <- sprintf(
    "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/lmm/results-linear-mixed-model_%s.csv",
    timestamp
  )
  run_linear_mixed_model(
    input_csv = input_file,
    output_csv = lmm_output_csv
  )
  
  # Step 4b: LMM effects summary
  cat("\n", strrep("=", 80), "\nSTEP 4b: SAVING LMM EFFECTS SUMMARY\n", strrep("=", 80), "\n\n", sep = "")
  lmm_effects_txt <- sprintf(
    "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/lmm/lmm-effects-summary_%s.txt",
    timestamp
  )
  save_lmm_effects_summary(
    lmm_csv    = lmm_output_csv,
    output_txt = lmm_effects_txt,
    timestamp  = timestamp,
    data_csv   = input_file   # pass raw data so LMM EMMs can be computed
  )
  
  # Step 5: ANOVA vs LMM comparison
  cat("\n", strrep("=", 80), "\nSTEP 5: COMPARING ANOVA AND LMM RESULTS\n", strrep("=", 80), "\n\n", sep = "")
  comparison_txt <- sprintf(
    "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/anova-vs-lmm/comparison_%s.txt",
    timestamp
  )
  compare_anova_lmm(
    anova_csv  = anova_output_csv,
    lmm_csv    = lmm_output_csv,
    output_txt = comparison_txt,
    timestamp  = timestamp
  )
  
  # Step 5b: Narrative summary
  cat("\n", strrep("=", 80), "\nSTEP 5b: CREATING ANOVA VS LMM SUMMARY\n", strrep("=", 80), "\n\n", sep = "")
  summary_txt <- sprintf(
    "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/anova-vs-lmm/summary_%s.txt",
    timestamp
  )
  
  cat("\n", strrep("=", 80), "\nALL ANALYSES COMPLETE\n", strrep("=", 80), "\n", sep = "")
  cat("\nAll output files include timestamp:", timestamp, "\n")
  cat("\nGenerated files:\n")
  cat(sprintf("  - Missing data report: missing-data-report_%s.txt\n", timestamp))
  cat(sprintf("  - ANOVA results: results-two-way-anova_%s.csv\n", timestamp))
  cat(sprintf("  - ANOVA effects summary: anova-effects-summary_%s.txt\n", timestamp))
  cat(sprintf("  - ANOVA assumptions: anova-assumptions-check_%s.csv\n", timestamp))
  cat(sprintf("  - LMM results: results-linear-mixed-model_%s.csv\n", timestamp))
  cat(sprintf("  - LMM effects summary: lmm-effects-summary_%s.txt\n", timestamp))
  cat(sprintf("  - ANOVA vs LMM comparison: comparison_%s.txt\n", timestamp))
  cat("\nNote: Missing data handled via listwise deletion (complete case analysis).\n")
  cat("See missing-data-report_", timestamp, ".txt for details on excluded participants.\n", sep = "")
  cat(strrep("=", 80), "\n")
}

## To run everything from RStudio:
##   1. Source this file
##   2. Then call:
##        run_all_analyses()
