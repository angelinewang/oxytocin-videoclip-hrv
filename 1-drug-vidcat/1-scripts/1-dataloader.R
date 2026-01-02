# Two-way ANOVA data needed for new CSV: two-way-anova.csv
#
# 1. File Name from HRV file
# 2. Sequence Label from the Prefix of the Metrics from HRV file 
# 3. Transform the Sequence Label from Sxx to Category 1, 2, 3, 4, 5
# 4. Drug Type from Randomisation File, with the File Name from HSV file 
# 5. Metrics: HF Power FFT ms, RM SSD ms, Mean RR, Stress Index, PND Index, SNS Index all from the HRV file, one per Sxx
#
# Rows: File Name per Sxx
# Columns: Sequence Label Transformed from Video Category, File Name, Drug Type, Every Metric needed
#
# Num Rows: # of files x # of sequence labels
# Num Columns: 9
#
# Maps File Names to correct Drug Types
# Removes the 5 with missing HRV Metrics

build_two_way_anova_csv <- function() {
  hrv_path <- "/Users/angwang/oxytocin-videoclip/0-data/HRV-metrics/HRV_Cosme2022_final.csv"
  rand_path <- "/Users/angwang/oxytocin-videoclip/0-data/drug-type/randomization_with_replacements_and_info.csv"
  output_path <- "/Users/angwang/oxytocin-videoclip/1-drug-vidcat/2-outputs/data/data-two-way-anova.csv"
  
  # Ensure output directory exists
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Read source data
  hrv_df <- utils::read.csv(hrv_path, stringsAsFactors = FALSE, check.names = FALSE)
  rand_df <- utils::read.csv(rand_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Map base file id (e.g., BSC_001) -> Drug Type (e.g., "Drug 1")
  # Heuristic: the first column looks like base id; the drug string appears in a column containing the word "Drug"
  base_to_drug <- character(0)
  
  for (i in seq_len(nrow(rand_df))) {
    r <- rand_df[i, , drop = FALSE]
    # Get base id from first column
    if (ncol(r) == 0) next
    
    base_id <- as.character(r[[1]])
    if (is.na(base_id) || base_id == "") next
    
    drug_val <- ""
    row_vals <- as.vector(as.matrix(r))
    for (val in row_vals) {
      if (is.character(val) && !is.na(val) && grepl("drug", tolower(val))) {
        drug_val <- val
        break
      }
    }
    # fallback: keep empty string if no drug found in row
    base_to_drug[base_id] <- drug_val
  }
  
  # Sequence to Video Category mapping
  seq_to_category <- c(
    stats::setNames(rep("Neutral (Scenery)", 4), c("S1", "S2", "S11", "S12")),
    stats::setNames(rep("Erotic couples (High arousal, High valence)", 4), c("S3", "S4", "S15", "S16")),
    stats::setNames(rep("Social Negative (Low Arousal, Low Valence)", 4), c("S5", "S6", "S13", "S14")),
    stats::setNames(rep("Social + (Low Arousal, High Valence)", 4), c("S7", "S8", "S19", "S20")),
    stats::setNames(rep("Horror (High Arousal, Low Valence)", 4), c("S9", "S10", "S17", "S18"))
  )
  
  # Desired sequences and metric column name templates
  sequences <- sprintf("S%d", 1:20)
  col_name <- function(seq, suffix) {
    paste0(seq, "_", suffix)
  }
  
  # Validate presence of FileName column
  if (!("FileName" %in% names(hrv_df))) {
    stop("Expected 'FileName' column in HRV CSV")
  }
  
  base_ids_found <- character(0)
  base_ids_mapped <- character(0)
  base_ids_unmapped <- character(0)
  
  rows_list <- list()
  row_idx <- 1L
  
  for (i in seq_len(nrow(hrv_df))) {
    row <- hrv_df[i, , drop = FALSE]
    file_name <- row[["FileName"]]
    if (is.null(file_name) || is.na(file_name)) {
      file_name <- ""
    } else {
      file_name <- as.character(file_name)
    }
    
    # Skip rows with empty file names
    if (file_name == "" || trimws(file_name) == "") next
    
    # Derive base id like BSC_001 from file name like BSC_001_IA.acq or BSC_001_ISC.acq
    # Use regex to extract BSC_XXX pattern (handles cases like BSC_034ISC.acq, BSC_001_IA.acq, BSC_001_AI.acq, or BSC_021 (text))
    pattern <- "(BSC_\\d+)"
    m <- regexpr(pattern, file_name)
    if (m[1] != -1) {
      base_id <- regmatches(file_name, m)
    } else {
      parts <- strsplit(file_name, "_", fixed = TRUE)[[1]]
      if (length(parts) >= 2) {
        base_id <- paste(parts[1], parts[2], sep = "_")
      } else {
        base_id <- parts[1]
      }
    }
    
    # Skip if we couldn't extract a valid base_id
    if (is.na(base_id) || base_id == "" || !startsWith(base_id, "BSC_")) next
    
    # Track this base_id as found
    if (!(base_id %in% base_ids_found)) {
      base_ids_found <- c(base_ids_found, base_id)
    }
    
    drug_type <- base_to_drug[[base_id]]
    if (is.null(drug_type)) drug_type <- ""
    
    # Track mapping status
    if (!is.na(drug_type) && trimws(drug_type) != "") {
      if (!(base_id %in% base_ids_mapped)) {
        base_ids_mapped <- c(base_ids_mapped, base_id)
      }
    } else {
      if (!(base_id %in% base_ids_unmapped)) {
        base_ids_unmapped <- c(base_ids_unmapped, base_id)
      }
      # Skip rows where Drug Type couldn't be determined
      next
    }
    
    for (seq in sequences) {
      # Only add row if at least one of the required columns exists for this sequence
      col_names <- c(
        "HF Power FFT (ms2)" = col_name(seq, "HFpow_FFT (ms2)"),
        "RM SSD (ms)"         = col_name(seq, "RMSSD (ms)"),
        "Mean RR"             = col_name(seq, "Mean RR (ms)"),
        "Stress Index"        = col_name(seq, "Stress index"),
        "PNS Index"           = col_name(seq, "PNS index"),
        "SNS Index"           = col_name(seq, "SNS index")
      )
      
      if (!any(col_names %in% names(hrv_df))) next
      
      get_val <- function(col_label) {
        cname <- col_names[[col_label]]
        if (!is.null(cname) && cname %in% names(hrv_df)) {
          return(row[[cname]])
        } else {
          return(NA_real_)
        }
      }
      
      hf_power     <- get_val("HF Power FFT (ms2)")
      rm_ssd       <- get_val("RM SSD (ms)")
      mean_rr      <- get_val("Mean RR")
      stress_index <- get_val("Stress Index")
      pns_index    <- get_val("PNS Index")
      sns_index    <- get_val("SNS Index")
      
      metrics <- c(hf_power, rm_ssd, mean_rr, stress_index, pns_index, sns_index)
      
      # Skip rows where ALL HRV metrics are missing/NaN
      if (all(is.na(metrics))) next
      
      video_category <- if (!is.null(seq_to_category[[seq]])) seq_to_category[[seq]] else ""
      
      rows_list[[row_idx]] <- list(
        "File Name"           = file_name,
        "Sequence Label"      = seq,
        "Video Category"      = video_category,
        "Drug Type"           = drug_type,
        "HF Power FFT (ms2)"  = hf_power,
        "RM SSD (ms)"         = rm_ssd,
        "Mean RR"             = mean_rr,
        "Stress Index"        = stress_index,
        "PNS Index"           = pns_index,
        "SNS Index"           = sns_index
      )
      row_idx <- row_idx + 1L
    }
  }
  
  # Create DataFrame with required column order
  out_cols <- c(
    "File Name",
    "Sequence Label",
    "Video Category",
    "Drug Type",
    "HF Power FFT (ms2)",
    "RM SSD (ms)",
    "Mean RR",
    "Stress Index",
    "PNS Index",
    "SNS Index"
  )
  
  # Create DataFrame with required column order
  out_cols <- c(
    "File Name",
    "Sequence Label",
    "Video Category",
    "Drug Type",
    "HF Power FFT (ms2)",
    "RM SSD (ms)",
    "Mean RR",
    "Stress Index",
    "PNS Index",
    "SNS Index"
  )
  
  if (length(rows_list) > 0L) {
    # Bind all rows into a data.frame
    out_df <- do.call(
      rbind.data.frame,
      c(rows_list, list(stringsAsFactors = FALSE))
    )
    
    # Force nice column names (with spaces) in the exact order we want
    # Assumes rows_list used fields in this order when constructing each row
    names(out_df) <- out_cols
  } else {
    out_df <- data.frame(matrix(ncol = length(out_cols), nrow = 0))
    names(out_df) <- out_cols
  }
  
  # At this point, out_df already has exactly the columns in out_cols in order,
  # so we DON'T need: out_df <- out_df[, out_cols, drop = FALSE]
  
  print(names(out_df))
  str(out_df)
  
  # Print mapping summary
  cat(strrep("=", 80), "\n")
  cat("FILE NAME TO DRUG TYPE MAPPING SUMMARY\n")
  cat(strrep("=", 80), "\n\n")
  cat("Total unique base IDs found in file names:", length(unique(base_ids_found)), "\n")
  cat("Base IDs successfully mapped to Drug Type:", length(unique(base_ids_mapped)), "\n")
  cat("Base IDs NOT found in randomization file:", length(unique(base_ids_unmapped)), "\n\n")
  
  if (length(base_ids_mapped) > 0L) {
    cat("SUCCESSFULLY MAPPED BASE IDs:\n")
    cat(strrep("-", 80), "\n")
    for (bid in sort(unique(base_ids_mapped))) {
      dt <- base_to_drug[[bid]]
      if (is.null(dt)) dt <- ""
      cat(sprintf("  %-12s \u2192 %s\n", bid, dt))
    }
    cat("\n")
  }
  
  if (length(base_ids_unmapped) > 0L) {
    cat("UNMAPPED BASE IDs (not in randomization file):\n")
    cat(strrep("-", 80), "\n")
    for (bid in sort(unique(base_ids_unmapped))) {
      cat(sprintf("  %-12s \u2192 NOT FOUND\n", bid))
      cat("\n")
    }
    cat("\n")
  }
  
  # Check for base IDs in randomization file that weren't found in HRV data
  rand_base_ids <- names(base_to_drug)
  unused_base_ids <- setdiff(rand_base_ids, base_ids_found)
  if (length(unused_base_ids) > 0L) {
    cat("BASE IDs IN RANDOMIZATION FILE BUT NOT FOUND IN HRV DATA:\n")
    cat(strrep("-", 80), "\n")
    for (bid in sort(unused_base_ids)) {
      dt <- base_to_drug[[bid]]
      if (is.null(dt)) dt <- ""
      cat(sprintf("  %-12s \u2192 %s (not in HRV file)\n", bid, dt))
    }
    cat("\n")
  }
  
  cat(strrep("=", 80), "\n\n")
  
  # Check if output file already exists and compare
  file_exists <- file.exists(output_path)
  
  if (file_exists) {
    changed <- TRUE
    tryCatch({
      existing_df <- utils::read.csv(output_path, stringsAsFactors = FALSE, check.names = FALSE)
      
      # Ensure column order matches (reorder existing_df to match out_cols)
      if (!identical(names(existing_df), out_cols)) {
        if (all(out_cols %in% names(existing_df))) {
          existing_df <- existing_df[, out_cols, drop = FALSE]
        } else {
          existing_df <- NULL
        }
      }
      
      if (!is.null(existing_df)) {
        rows_different <- nrow(existing_df) != nrow(out_df)
        if (!rows_different) {
          # Try to align types and compare
          existing_df_reset <- existing_df
          out_df_reset <- out_df
          
          for (cn in out_cols) {
            if (cn %in% names(existing_df_reset) && cn %in% names(out_df_reset)) {
              # Attempt to coerce to character to safely compare
              existing_df_reset[[cn]] <- as.character(existing_df_reset[[cn]])
              out_df_reset[[cn]] <- as.character(out_df_reset[[cn]])
            }
          }
          
          values_different <- !isTRUE(all.equal(existing_df_reset, out_df_reset, check.attributes = FALSE))
        } else {
          values_different <- TRUE
        }
        
        if (rows_different || values_different) {
          utils::write.csv(out_df, output_path, row.names = FALSE)
          cat("\u2713 Data file updated:", output_path, "\n")
          if (rows_different) {
            cat("  Row count changed:", nrow(existing_df), "\u2192", nrow(out_df), "\n")
          }
          if (values_different) {
            cat("  Values changed: Data differs from existing file\n")
          }
        } else {
          changed <- FALSE
          cat("\u2713 Data file unchanged:", output_path, "\n")
          cat("  Rows:", nrow(out_df), "(no changes detected)\n")
          cat("  Skipping file write (existing file is identical)\n")
        }
      } else {
        utils::write.csv(out_df, output_path, row.names = FALSE)
        cat("\u2713 Data file updated:", output_path, "\n")
        cat("  Column structure changed\n")
      }
    }, error = function(e) {
      cat("Warning: Could not compare with existing file:", conditionMessage(e), "\n")
      cat("Writing new file:", output_path, "\n")
      utils::write.csv(out_df, output_path, row.names = FALSE)
    })
  } else {
    # File doesn't exist, create it
    utils::write.csv(out_df, output_path, row.names = FALSE)
    cat("\u2713 Data file created:", output_path, "\n")
    cat("  Rows:", nrow(out_df), "\n")
  }
}


# Three-way ANOVA data needed for new CSV: three-way-anova.csv
#
# 1. File Name from HRV file
# 2. Sequence Label from the Prefix of the Metrics from HRV file 
# 3. Transform the Sequence Label from Sxx to Category 1, 2, 3, 4, 5
# 4. Drug Type from Randomisation File, with the File Name from HSV file 
# 5. Metrics: HF Power FFT ms, RM SSD ms, Mean RR, Stress Index, PND Index, SNS Index all from the HRV file, one per Sxx
#
# Rows: File Name per Sxx
# Columns: Valence Transformed from Sequence Label, Arousal transformed from sequence label, File Name, Drug Type, Every Metric needed
#
# Num Rows: # of files x # of sequence labels 
# Num Columns: 10
