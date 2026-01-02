## =============================================================================
## Libraries
## =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(stringr)
})

## =============================================================================
## Main function
## =============================================================================

calculate_group_means <- function(input_csv, output_csv, output_txt, timestamp) {
  # 1. Calculate raw means, SD, SEM for:
  #    - Drug Type
  #    - Video Category
  #    - Video Category × Drug Type
  # 2. Fit ANOVA and LMM and compute estimated marginal means (EMMs) for:
  #    - Video Category
  #    - Drug Type
  #    - Video Category × Drug Type
  
  df <- readr::read_csv(input_csv, show_col_types = FALSE)
  
  # HRV metrics
  hrv_metrics <- c(
    "HF Power FFT (ms2)",
    "RM SSD (ms)",
    "Mean RR",
    "Stress Index",
    "PNS Index",
    "SNS Index"
  )
  
  # Ensure categorical (raw columns)
  df[["Drug Type"]]    <- as.factor(df[["Drug Type"]])
  df[["Video Category"]] <- as.factor(df[["Video Category"]])
  
  # Alias columns with no spaces for modelling (ANOVA & LMM)
  df$DrugType      <- df[["Drug Type"]]
  df$VideoCategory <- df[["Video Category"]]
  
  # Participant_ID for LMM
  if (!("Participant_ID" %in% names(df)) && ("File Name" %in% names(df))) {
    df$Participant_ID <- sub("\\.(acq|ACQ)$", "", df[["File Name"]])
  } else if (!("Participant_ID" %in% names(df))) {
    # Fallback: use row index
    df$Participant_ID <- as.character(seq_len(nrow(df)))
  }
  
  ## ============================================================================
  ## 1) RAW GROUP MEANS
  ## ============================================================================
  
  results_rows <- list()
  
  cat(strrep("=", 80), "\n")
  cat("CALCULATING RAW MEANS BY DRUG TYPE\n")
  cat(strrep("=", 80), "\n")
  
  for (metric in hrv_metrics) {
    for (drug_type in levels(df[["Drug Type"]])) {
      subset_vals <- df[df[["Drug Type"]] == drug_type, metric, drop = TRUE]
      subset_vals <- subset_vals[!is.na(subset_vals)]
      
      if (length(subset_vals) > 0) {
        mean_val <- mean(subset_vals)
        std_val  <- sd(subset_vals)
        sem_val  <- if (length(subset_vals) > 1) std_val / sqrt(length(subset_vals)) else NA_real_
        n_val    <- length(subset_vals)
        
        results_rows[[length(results_rows) + 1L]] <- data.frame(
          grouping_variable = "Drug Type",
          group_level       = drug_type,
          metric            = metric,
          mean              = mean_val,
          std               = std_val,
          sem               = sem_val,
          n                 = n_val,
          stringsAsFactors  = FALSE
        )
      }
    }
  }
  
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("CALCULATING RAW MEANS BY VIDEO CATEGORY\n")
  cat(strrep("=", 80), "\n")
  
  for (metric in hrv_metrics) {
    for (video_cat in levels(df[["Video Category"]])) {
      subset_vals <- df[df[["Video Category"]] == video_cat, metric, drop = TRUE]
      subset_vals <- subset_vals[!is.na(subset_vals)]
      
      if (length(subset_vals) > 0) {
        mean_val <- mean(subset_vals)
        std_val  <- sd(subset_vals)
        sem_val  <- if (length(subset_vals) > 1) std_val / sqrt(length(subset_vals)) else NA_real_
        n_val    <- length(subset_vals)
        
        results_rows[[length(results_rows) + 1L]] <- data.frame(
          grouping_variable = "Video Category",
          group_level       = video_cat,
          metric            = metric,
          mean              = mean_val,
          std               = std_val,
          sem               = sem_val,
          n                 = n_val,
          stringsAsFactors  = FALSE
        )
      }
    }
  }
  
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("CALCULATING RAW MEANS BY VIDEO CATEGORY × DRUG TYPE\n")
  cat(strrep("=", 80), "\n")
  
  for (metric in hrv_metrics) {
    for (video_cat in levels(df[["Video Category"]])) {
      for (drug_type in levels(df[["Drug Type"]])) {
        subset_vals <- df[df[["Video Category"]] == video_cat &
                            df[["Drug Type"]] == drug_type, metric, drop = TRUE]
        subset_vals <- subset_vals[!is.na(subset_vals)]
        
        if (length(subset_vals) > 0) {
          mean_val <- mean(subset_vals)
          std_val  <- sd(subset_vals)
          sem_val  <- if (length(subset_vals) > 1) std_val / sqrt(length(subset_vals)) else NA_real_
          n_val    <- length(subset_vals)
          
          results_rows[[length(results_rows) + 1L]] <- data.frame(
            grouping_variable = "Video Category × Drug Type",
            group_level       = paste(video_cat, "\u00d7", drug_type),
            metric            = metric,
            mean              = mean_val,
            std               = std_val,
            sem               = sem_val,
            n                 = n_val,
            stringsAsFactors  = FALSE
          )
        }
      }
    }
  }
  
  if (length(results_rows) > 0L) {
    results_df <- dplyr::bind_rows(results_rows)
  } else {
    results_df <- data.frame(
      grouping_variable = character(0),
      group_level       = character(0),
      metric            = character(0),
      mean              = numeric(0),
      std               = numeric(0),
      sem               = numeric(0),
      n                 = integer(0),
      stringsAsFactors  = FALSE
    )
  }
  
  # Save RAW means CSV
  readr::write_csv(results_df, output_csv)
  cat("\nRaw group means saved to:", output_csv, "\n")
  
  ## ============================================================================
  ## 2) ESTIMATED MARGINAL MEANS (ANOVA & LMM)
  ## ============================================================================
  
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("CALCULATING ESTIMATED MARGINAL MEANS (ANOVA & LMM)\n")
  cat(strrep("=", 80), "\n")
  
  mm_rows <- list()
  
  for (metric in hrv_metrics) {
    cols_needed <- c(metric, "VideoCategory", "DrugType", "Participant_ID")
    if (!all(cols_needed %in% names(df))) {
      message("Skipping ", metric, ": required columns missing.")
      next
    }
    
    sub <- df[, cols_needed]
    sub <- sub[stats::complete.cases(sub), , drop = FALSE]
    if (nrow(sub) == 0L) {
      message("Skipping ", metric, ": no complete cases.")
      next
    }
    
    sub$VideoCategory <- factor(sub$VideoCategory)
    sub$DrugType      <- factor(sub$DrugType)
    
    vc_levels <- levels(sub$VideoCategory)
    dt_levels <- levels(sub$DrugType)
    
    # ANOVA model
    form_aov <- as.formula(paste0("`", metric, "` ~ VideoCategory * DrugType"))
    fit_aov  <- lm(form_aov, data = sub)
    
    # LMM model
    fit_lmm <- tryCatch(
      lmer(update(form_aov, . ~ . + (1 | Participant_ID)), data = sub, REML = TRUE),
      error = function(e) e
    )
    lmm_ok <- !inherits(fit_lmm, "error")
    
    # Cell EMMs (VideoCategory × DrugType)
    emm_cell_aov <- emmeans::emmeans(fit_aov, ~ VideoCategory * DrugType)
    emm_cell_aov_df <- as.data.frame(emm_cell_aov)
    
    if (lmm_ok) {
      emm_cell_lmm    <- emmeans::emmeans(fit_lmm, ~ VideoCategory * DrugType)
      emm_cell_lmm_df <- as.data.frame(emm_cell_lmm)
    } else {
      emm_cell_lmm_df <- NULL
    }
    
    for (i in seq_len(nrow(emm_cell_aov_df))) {
      vc <- as.character(emm_cell_aov_df$VideoCategory[i])
      dt <- as.character(emm_cell_aov_df$DrugType[i])
      cell_label <- paste(vc, "\u00d7", dt)
      
      mm_rows[[length(mm_rows) + 1L]] <- data.frame(
        model  = "ANOVA",
        effect = "VideoCategory × DrugType",
        level  = cell_label,
        metric = metric,
        emm    = emm_cell_aov_df$emmean[i],
        stringsAsFactors = FALSE
      )
      
      if (!is.null(emm_cell_lmm_df)) {
        emm_lmm_val <- emm_cell_lmm_df$emmean[
          emm_cell_lmm_df$VideoCategory == vc &
            emm_cell_lmm_df$DrugType == dt
        ]
        emm_lmm_val <- if (length(emm_lmm_val) == 1L) emm_lmm_val else NA_real_
      } else {
        emm_lmm_val <- NA_real_
      }
      
      mm_rows[[length(mm_rows) + 1L]] <- data.frame(
        model  = "LMM",
        effect = "VideoCategory × DrugType",
        level  = cell_label,
        metric = metric,
        emm    = emm_lmm_val,
        stringsAsFactors = FALSE
      )
    }
    
    # Marginal means by VideoCategory
    emm_vc_aov_df <- as.data.frame(emmeans::emmeans(fit_aov, "VideoCategory"))
    if (lmm_ok) {
      emm_vc_lmm_df <- as.data.frame(emmeans::emmeans(fit_lmm, "VideoCategory"))
    } else {
      emm_vc_lmm_df <- NULL
    }
    
    for (vc in vc_levels) {
      emm_aov_val <- emm_vc_aov_df$emmean[emm_vc_aov_df$VideoCategory == vc]
      emm_aov_val <- if (length(emm_aov_val) == 1L) emm_aov_val else NA_real_
      mm_rows[[length(mm_rows) + 1L]] <- data.frame(
        model  = "ANOVA",
        effect = "VideoCategory",
        level  = vc,
        metric = metric,
        emm    = emm_aov_val,
        stringsAsFactors = FALSE
      )
      
      if (!is.null(emm_vc_lmm_df)) {
        emm_lmm_val <- emm_vc_lmm_df$emmean[emm_vc_lmm_df$VideoCategory == vc]
        emm_lmm_val <- if (length(emm_lmm_val) == 1L) emm_lmm_val else NA_real_
      } else {
        emm_lmm_val <- NA_real_
      }
      
      mm_rows[[length(mm_rows) + 1L]] <- data.frame(
        model  = "LMM",
        effect = "VideoCategory",
        level  = vc,
        metric = metric,
        emm    = emm_lmm_val,
        stringsAsFactors = FALSE
      )
    }
    
    # Marginal means by DrugType
    emm_dt_aov_df <- as.data.frame(emmeans::emmeans(fit_aov, "DrugType"))
    if (lmm_ok) {
      emm_dt_lmm_df <- as.data.frame(emmeans::emmeans(fit_lmm, "DrugType"))
    } else {
      emm_dt_lmm_df <- NULL
    }
    
    for (dt in dt_levels) {
      emm_aov_val <- emm_dt_aov_df$emmean[emm_dt_aov_df$DrugType == dt]
      emm_aov_val <- if (length(emm_aov_val) == 1L) emm_aov_val else NA_real_
      mm_rows[[length(mm_rows) + 1L]] <- data.frame(
        model  = "ANOVA",
        effect = "DrugType",
        level  = dt,
        metric = metric,
        emm    = emm_aov_val,
        stringsAsFactors = FALSE
      )
      
      if (!is.null(emm_dt_lmm_df)) {
        emm_lmm_val <- emm_dt_lmm_df$emmean[emm_dt_lmm_df$DrugType == dt]
        emm_lmm_val <- if (length(emm_lmm_val) == 1L) emm_lmm_val else NA_real_
      } else {
        emm_lmm_val <- NA_real_
      }
      
      mm_rows[[length(mm_rows) + 1L]] <- data.frame(
        model  = "LMM",
        effect = "DrugType",
        level  = dt,
        metric = metric,
        emm    = emm_lmm_val,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(mm_rows) > 0L) {
    mm_df <- dplyr::bind_rows(mm_rows)
  } else {
    mm_df <- data.frame(
      model  = character(0),
      effect = character(0),
      level  = character(0),
      metric = character(0),
      emm    = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  # Save marginal means CSV (replace _raw.csv with _marginal_means.csv)
  mm_output_csv <- sub("_raw\\.csv$", "_marginal_means.csv", output_csv)
  readr::write_csv(mm_df, mm_output_csv)
  cat("Marginal means (ANOVA & LMM) saved to:", mm_output_csv, "\n")
  
  ## ============================================================================
  ## 3) TEXT REPORT (RAW + marginal means summary)
  ## ============================================================================
  
  report <- character(0)
  report <- c(
    report,
    strrep("=", 80),
    "GROUP MEANS SUMMARY (RAW DESCRIPTIVE MEANS)",
    strrep("=", 80),
    "",
    paste0("Generated: ", timestamp),
    "",
    "This report shows descriptive statistics (mean, standard deviation,",
    "standard error of the mean, and sample size) for each HRV metric",
    "grouped by Drug Type, by Video Category, and by Video Category × Drug Type.",
    "",
    "Estimated marginal means (EMMs) from ANOVA and LMM are also summarised",
    "below for the main effects (Video Category, Drug Type).",
    "",
    "Full EMM tables (including interaction cells) are saved in:",
    paste0("  ", mm_output_csv),
    ""
  )
  
  # RAW summary: Drug Type
  report <- c(
    report,
    strrep("=", 80),
    "MEANS BY DRUG TYPE (RAW)",
    strrep("=", 80),
    ""
  )
  
  drug_results <- results_df %>% filter(grouping_variable == "Drug Type")
  drug_types   <- sort(unique(drug_results$group_level))
  
  for (metric in hrv_metrics) {
    report <- c(report, "", metric, strrep("-", 80))
    metric_drug <- drug_results %>% filter(metric == !!metric)
    
    for (drug_type in drug_types) {
      row <- metric_drug %>% filter(group_level == !!drug_type)
      if (nrow(row) > 0) {
        r <- row[1, ]
        report <- c(
          report,
          paste0("  ", drug_type, ":"),
          sprintf("    Mean = %.4f", r$mean),
          sprintf("    SD = %.4f", r$std),
          sprintf("    SEM = %.4f", r$sem),
          paste0("    n = ", as.integer(r$n))
        )
      }
    }
    
    if (nrow(metric_drug) == 2) {
      row1 <- metric_drug[1, ]
      row2 <- metric_drug[2, ]
      diff <- row2$mean - row1$mean
      pct_diff <- if (!is.na(row1$mean) && row1$mean != 0) diff / row1$mean * 100 else NA_real_
      report <- c(
        report,
        sprintf(
          "  Difference (%s - %s): %.4f (%+.2f%%)",
          row2$group_level, row1$group_level, diff, pct_diff
        )
      )
    }
  }
  
  # RAW summary: Video Category
  report <- c(
    report,
    "",
    strrep("=", 80),
    "MEANS BY VIDEO CATEGORY (RAW)",
    strrep("=", 80),
    ""
  )
  
  video_results   <- results_df %>% filter(grouping_variable == "Video Category")
  video_categories <- sort(unique(video_results$group_level))
  
  for (metric in hrv_metrics) {
    report <- c(report, "", metric, strrep("-", 80))
    metric_video <- video_results %>% filter(metric == !!metric)
    
    for (video_cat in video_categories) {
      row <- metric_video %>% filter(group_level == !!video_cat)
      if (nrow(row) > 0) {
        r <- row[1, ]
        report <- c(
          report,
          paste0("  ", video_cat, ":"),
          sprintf("    Mean = %.4f", r$mean),
          sprintf("    SD = %.4f", r$std),
          sprintf("    SEM = %.4f", r$sem),
          paste0("    n = ", as.integer(r$n))
        )
      }
    }
    
    if (nrow(metric_video) > 1) {
      min_val   <- min(metric_video$mean, na.rm = TRUE)
      max_val   <- max(metric_video$mean, na.rm = TRUE)
      range_val <- max_val - min_val
      pct_range <- if (!is.na(min_val) && min_val != 0) range_val / min_val * 100 else NA_real_
      report <- c(
        report,
        sprintf("  Range (max - min): %.4f (%+.2f%%)", range_val, pct_range)
      )
    }
  }
  
  # RAW summary: Video Category × Drug Type
  report <- c(
    report,
    "",
    strrep("=", 80),
    "MEANS BY VIDEO CATEGORY × DRUG TYPE (RAW)",
    strrep("=", 80),
    ""
  )
  
  interaction_results <- results_df %>%
    filter(grouping_variable == "Video Category × Drug Type")
  
  for (metric in hrv_metrics) {
    report <- c(report, "", metric, strrep("-", 80))
    metric_interaction <- interaction_results %>% filter(metric == !!metric)
    
    for (video_cat in video_categories) {
      video_interactions <- metric_interaction %>%
        filter(str_detect(group_level, fixed(video_cat)))
      if (nrow(video_interactions) > 0) {
        report <- c(report, paste0("  ", video_cat, ":"))
        for (i in seq_len(nrow(video_interactions))) {
          gl <- video_interactions$group_level[i]
          parts <- strsplit(gl, "\u00d7")[[1]]
          if (length(parts) == 2) {
            drug_part <- trimws(parts[2])
            r <- video_interactions[i, ]
            report <- c(
              report,
              paste0("    ", drug_part, ":"),
              sprintf("      Mean = %.4f", r$mean),
              sprintf("      SD = %.4f", r$std),
              sprintf("      SEM = %.4f", r$sem),
              paste0("      n = ", as.integer(r$n))
            )
          }
        }
        
        if (nrow(video_interactions) == 2) {
          r1 <- video_interactions[1, ]
          r2 <- video_interactions[2, ]
          diff <- r2$mean - r1$mean
          pct_diff <- if (!is.na(r1$mean) && r1$mean != 0) diff / r1$mean * 100 else NA_real_
          report <- c(
            report,
            sprintf("    Difference (%s - %s): %.4f (%+.2f%%)",
                    video_interactions$group_level[2],
                    video_interactions$group_level[1],
                    diff, pct_diff)
          )
        }
      }
    }
  }
  
  ## ============================================================================
  ## 4) SUMMARY OF ESTIMATED MARGINAL MEANS (ANOVA & LMM)
  ## ============================================================================
  
  report <- c(
    report,
    "",
    strrep("=", 80),
    "ESTIMATED MARGINAL MEANS (MODEL-BASED: ANOVA & LMM)",
    strrep("=", 80),
    "",
    "Values below are model-based means, adjusted for the other factor.",
    "They are averaged equally over levels of the other factor (balanced design assumption).",
    ""
  )
  
  mm_main <- mm_df %>% filter(effect %in% c("VideoCategory", "DrugType"))
  
  for (metric in hrv_metrics) {
    report <- c(report, "", metric, strrep("-", 80))
    
    # Video Category EMMs
    vc_mm <- mm_main %>%
      filter(metric == !!metric, effect == "VideoCategory")
    if (nrow(vc_mm) > 0) {
      report <- c(report, "Video Category (EMMs):")
      for (level in sort(unique(vc_mm$level))) {
        row_anova <- vc_mm %>% filter(model == "ANOVA", level == !!level)
        row_lmm   <- vc_mm %>% filter(model == "LMM", level == !!level)
        emm_anova <- if (nrow(row_anova) > 0) row_anova$emm[1] else NA_real_
        emm_lmm   <- if (nrow(row_lmm) > 0)   row_lmm$emm[1]   else NA_real_
        report <- c(
          report,
          sprintf("  %s: ANOVA EMM = %.4f, LMM EMM = %.4f",
                  level, emm_anova, emm_lmm)
        )
      }
    }
    
    # Drug Type EMMs
    dt_mm <- mm_main %>%
      filter(metric == !!metric, effect == "DrugType")
    if (nrow(dt_mm) > 0) {
      report <- c(report, "", "Drug Type (EMMs):")
      for (level in sort(unique(dt_mm$level))) {
        row_anova <- dt_mm %>% filter(model == "ANOVA", level == !!level)
        row_lmm   <- dt_mm %>% filter(model == "LMM",   level == !!level)
        emm_anova <- if (nrow(row_anova) > 0) row_anova$emm[1] else NA_real_
        emm_lmm   <- if (nrow(row_lmm) > 0)   row_lmm$emm[1]   else NA_real_
        report <- c(
          report,
          sprintf("  %s: ANOVA EMM = %.4f, LMM EMM = %.4f",
                  level, emm_anova, emm_lmm)
        )
      }
    }
  }
  
  report <- c(
    report,
    "",
    strrep("=", 80),
    "NOTES",
    strrep("=", 80),
    "",
    "SD = Standard Deviation",
    "SEM = Standard Error of the Mean (SD / \u221an)",
    "n = Sample size (number of observations)",
    "",
    "Raw group means are descriptive, based on the observed data.",
    "Estimated marginal means (ANOVA & LMM) are model-based and are",
    "saved separately in the *_marginal_means.csv file.",
    strrep("=", 80)
  )
  
  # Write text report
  writeLines(report, con = output_txt)
  cat("Group means summary saved to:", output_txt, "\n")
  
  # Console summary
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("SUMMARY STATISTICS\n")
  cat(strrep("=", 80), "\n")
  cat("\nTotal groups analyzed (RAW):\n")
  cat("  Drug Type:", length(drug_types), "groups\n")
  cat("  Video Category:", length(video_categories), "groups\n")
  cat("  Video Category \u00d7 Drug Type:",
      length(drug_types) * length(video_categories), "combinations\n")
  cat("\nMarginal means file:", mm_output_csv, "\n")
  cat("\n", strrep("=", 80), "\n", sep = "")
}

## =============================================================================
## Example __main__ equivalent
## =============================================================================

# To run from R:
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
input_file <- "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/data/data-two-way-anova.csv"
output_csv <- sprintf("/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/means/group-means_%s_raw.csv", timestamp)
output_txt <- sprintf("/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/means/group-means-summary_%s.txt", timestamp)

cat("\n", strrep("=", 80), "\nCALCULATING GROUP MEANS + MARGINAL MEANS (ANOVA & LMM)\n", strrep("=", 80), "\n\n", sep = "")
calculate_group_means(
   input_csv  = input_file,
   output_csv = output_csv,
   output_txt = output_txt,
   timestamp  = timestamp
 )
 cat("\n", strrep("=", 80), "\nGROUP & MARGINAL MEANS CALCULATION COMPLETE\n", strrep("=", 80), "\n", sep = "")
 cat("\nAll output files include timestamp:", timestamp, "\n")
 cat(strrep("=", 80), "\n")
