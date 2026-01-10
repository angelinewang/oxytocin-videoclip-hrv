# -------------------------------
# Sample size summary script
# -------------------------------

# Input and output paths
input_file  <- "~/oxytocin-videoclip-hrv/1-drug-vidcat/2-outputs/data/data-two-way-anova.csv"
output_dir  <- "~/oxytocin-videoclip-hrv/1-drug-vidcat/2-outputs/data"
txt_out <- file.path(output_dir, "sample-size.txt")

# --- Helper: stop with helpful message ---
stopf <- function(...) stop(sprintf(...), call. = FALSE)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- Check input file exists ---
cat("Working directory:", getwd(), "\n")
cat("Attempting to read input file:", normalizePath(input_file, winslash = "/", mustWork = FALSE), "\n")

if (!file.exists(input_file)) {
  stopf("Input file not found: '%s'\nTip: confirm the path is correct relative to getwd().", input_file)
}

# --- Read the data safely ---
df <- tryCatch(
  read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) stopf("Failed to read CSV: %s", e$message)
)

# --- Basic sanity checks on data ---
if (!is.data.frame(df)) stopf("Read object is not a data.frame (unexpected).")
if (nrow(df) == 0) stopf("CSV was read, but it contains 0 rows.")
if (ncol(df) == 0) stopf("CSV was read, but it contains 0 columns.")

cat("Read OK. Rows:", nrow(df), "Cols:", ncol(df), "\n")
cat("Column names detected:\n")
print(names(df))

# --- Verify required columns exist exactly ---
required_cols <- c("File Name", "Drug Type")
missing <- setdiff(required_cols, names(df))

if (length(missing) > 0) {
  # Provide helpful suggestions using approximate matching
  suggestions <- lapply(missing, function(m) {
    approx <- agrep(m, names(df), value = TRUE, max.distance = 0.25)
    if (length(approx) == 0) "<no close match>" else paste(approx, collapse = ", ")
  })
  names(suggestions) <- missing
  
  cat("\nERROR: Missing required columns:\n")
  print(missing)
  cat("\nClosest matches found in your CSV:\n")
  print(suggestions)
  
  stopf(
    "Cannot proceed because required column(s) are missing.\nFix: rename columns in the CSV or update required_cols in the script."
  )
}

# --- Extract columns and validate content ---
file_names <- df[["File Name"]]
drug_types <- df[["Drug Type"]]

# Check for all-missing / all-blank
if (all(is.na(file_names)) || all(trimws(as.character(file_names)) == "")) {
  stopf("'File Name' column exists but is entirely NA/blank. Cannot compute unique counts.")
}
if (all(is.na(drug_types)) || all(trimws(as.character(drug_types)) == "")) {
  stopf("'Drug Type' column exists but is entirely NA/blank. Cannot compute drug type counts.")
}
# --- Step 1: Count unique File Names ---
# Treat blanks as missing
file_names_clean <- trimws(as.character(file_names))
file_names_clean[file_names_clean == ""] <- NA

n_unique_files <- length(unique(file_names_clean[!is.na(file_names_clean)]))
if (n_unique_files == 0) {
  stopf("After cleaning blanks, found 0 unique 'File Name' values. Something is off in the input.")
}

# --- Step 2: Keep only first row per unique File Name, then count Drug Type ---
df_unique <- df[!duplicated(file_names_clean) & !is.na(file_names_clean), , drop = FALSE]

if (nrow(df_unique) == 0) {
  stopf("After deduplicating by 'File Name', got 0 rows. Check 'File Name' values.")
}

drug_types_unique <- trimws(as.character(df_unique[["Drug Type"]]))
drug_types_unique[drug_types_unique == ""] <- NA

drug_type_counts <- sort(table(drug_types_unique, useNA = "ifany"), decreasing = TRUE)

# --- Write TXT output ---
sink(txt_out)

cat("Sample Size Summary\n")
cat("===================\n\n")

cat("Input file:\n")
cat(normalizePath(input_file, winslash = "/", mustWork = TRUE), "\n\n")

cat("Rows in raw input:", nrow(df), "\n")
cat("Rows after dedupe by 'File Name':", nrow(df_unique), "\n\n")

cat("1. Number of unique 'File Name' values:\n")
cat(n_unique_files, "\n\n")

cat("2. Samples per 'Drug Type' (first occurrence per unique 'File Name'):\n")
print(drug_type_counts)

sink()

cat("\n✅ Outputs written:\n")
cat(" - TXT:", normalizePath(txt_out, winslash = "/", mustWork = FALSE), "\n")