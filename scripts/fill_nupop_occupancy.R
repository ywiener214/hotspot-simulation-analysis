#!/usr/bin/env Rscript

# Fill nucleosome occupancy for each row in an Excel sheet using NuPoP.
#
# For each gene:
# 1) Load <Gene>.fasta from --fasta-dir
# 2) Run NuPoP once on the full gene sequence
# 3) For each row with that gene, find exact matches of Sequence_9_Long in the gene
# 4) Take occupancy at the middle base of the 9-mer (position start + 4, 1-based)
#
# Output:
# - Updates/creates occupancy column (default: Average_Nucleosome_Occupancy)
# - Adds diagnostic columns:
#   - Seq9_Match_Count
#   - Seq9_Selected_Start_1Based
#   - Seq9_Selected_Center_1Based
# - Adds nucleosome-position columns:
#   - NuPoP_Window_Left_1Based
#   - NuPoP_Window_Right_1Based
#   - NuPoP_Window_Length
#   - NuPoP_Selected146_Count
#   - NuPoP_Center_In_Selected146_1Based

suppressWarnings(suppressMessages({
  library(openxlsx)
  library(NuPoP)
}))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) return(default)
  sub(key, "", hit[[1]], fixed = TRUE)
}

required <- c("input-xlsx", "output-xlsx", "fasta-dir")
for (k in required) {
  if (is.null(get_arg(k))) {
    stop(
      paste0(
        "Missing required argument --", k, "\n",
        "Example:\n",
        "Rscript fill_nupop_occupancy.R ",
        "--input-xlsx=/path/in.xlsx ",
        "--output-xlsx=/path/out.xlsx ",
        "--fasta-dir=/path/Hotspots_FASTA"
      )
    )
  }
}

input_xlsx <- get_arg("input-xlsx")
output_xlsx <- get_arg("output-xlsx")
fasta_dir <- get_arg("fasta-dir")
sheet_name <- get_arg("sheet", "Hotspots")
gene_col <- get_arg("gene-col", "Gene")
seq9_col <- get_arg("seq9-col", "Sequence_9_Long")
occup_col <- get_arg("occup-col", "Average_Nucleosome_Occupancy")
position_col <- get_arg("position-col", "Nucleosome Position")
species <- as.integer(get_arg("species", "1"))   # 1 = human
model <- as.integer(get_arg("model", "4"))
use_chem <- tolower(get_arg("use-chem", "true")) %in% c("1", "true", "yes", "y")
multi_match_strategy <- tolower(get_arg("multi-match", "first")) # first|mean|na

if (!file.exists(input_xlsx)) stop("Input XLSX not found: ", input_xlsx)
if (!dir.exists(fasta_dir)) stop("FASTA directory not found: ", fasta_dir)
if (!(multi_match_strategy %in% c("first", "mean", "na"))) {
  stop("--multi-match must be one of: first, mean, na")
}

message("Loading workbook...")
df <- read.xlsx(input_xlsx, sheet = sheet_name, colNames = TRUE)
if (!(gene_col %in% names(df))) stop("Missing column: ", gene_col)
if (!(seq9_col %in% names(df))) stop("Missing column: ", seq9_col)

if (!(occup_col %in% names(df))) df[[occup_col]] <- NA_real_
if (!(position_col %in% names(df))) df[[position_col]] <- NA_integer_
if (!("Seq9_Match_Count" %in% names(df))) df[["Seq9_Match_Count"]] <- NA_integer_
if (!("Seq9_Selected_Start_1Based" %in% names(df))) df[["Seq9_Selected_Start_1Based"]] <- NA_integer_
if (!("Seq9_Selected_Center_1Based" %in% names(df))) df[["Seq9_Selected_Center_1Based"]] <- NA_integer_
if (!("NuPoP_Window_Left_1Based" %in% names(df))) df[["NuPoP_Window_Left_1Based"]] <- NA_integer_
if (!("NuPoP_Window_Right_1Based" %in% names(df))) df[["NuPoP_Window_Right_1Based"]] <- NA_integer_
if (!("NuPoP_Window_Length" %in% names(df))) df[["NuPoP_Window_Length"]] <- NA_integer_
if (!("NuPoP_Selected146_Count" %in% names(df))) df[["NuPoP_Selected146_Count"]] <- NA_integer_
if (!("NuPoP_Center_In_Selected146_1Based" %in% names(df))) {
  df[["NuPoP_Center_In_Selected146_1Based"]] <- NA_integer_
}

clean_seq <- function(x) {
  if (is.na(x) || is.null(x)) return(NA_character_)
  s <- toupper(gsub("[^ACGTN]", "", as.character(x)))
  if (nchar(s) == 0) return(NA_character_)
  s
}

read_fasta <- function(path) {
  if (!file.exists(path)) stop("FASTA file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  seq_lines <- lines[!grepl("^>", lines)]
  toupper(paste(seq_lines, collapse = ""))
}

all_fixed_matches_1based <- function(pattern, text) {
  # gregexpr with fixed=TRUE; returns 1-based starts.
  m <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (length(m) == 1 && m[1] == -1) integer(0) else as.integer(m)
}

compute_center_position_in_selected146 <- function(occup, center_pos, threshold = 0.1) {
  n <- length(occup)
  if (center_pos < 1 || center_pos > n) {
    return(list(
      left = NA_integer_, right = NA_integer_, window_len = NA_integer_,
      selected_count = NA_integer_, center_rank = NA_integer_, center_below_threshold = NA
    ))
  }
  if (is.na(occup[center_pos]) || occup[center_pos] <= threshold) {
    return(list(
      left = NA_integer_, right = NA_integer_, window_len = 0L,
      selected_count = 0L, center_rank = 0L, center_below_threshold = TRUE
    ))
  }

  # Expand from center until first <= threshold on each side.
  left <- center_pos
  while (left > 1 && !is.na(occup[left - 1]) && occup[left - 1] > threshold) {
    left <- left - 1
  }
  right <- center_pos
  while (right < n && !is.na(occup[right + 1]) && occup[right + 1] > threshold) {
    right <- right + 1
  }

  window_idx <- seq.int(left, right)
  window_vals <- occup[window_idx]
  window_len <- length(window_idx)
  if (window_len == 0) {
    return(list(
      left = left, right = right, window_len = 0L,
      selected_count = 0L, center_rank = 0L, center_below_threshold = FALSE
    ))
  }

  # Select top 146 probabilities from window.
  k <- min(146L, window_len)
  ord <- order(window_vals, decreasing = TRUE, na.last = NA)
  top_local <- ord[seq_len(k)]
  selected_idx <- window_idx[top_local]

  # Ensure center is included in the selected set.
  if (!(center_pos %in% selected_idx)) {
    if (k < 146L && k < window_len) {
      selected_idx <- c(selected_idx, center_pos)
    } else if (length(selected_idx) > 0) {
      # Replace the weakest selected point with center.
      sel_vals <- occup[selected_idx]
      weakest <- order(sel_vals, decreasing = FALSE, na.last = TRUE)[1]
      selected_idx[weakest] <- center_pos
    } else {
      selected_idx <- c(center_pos)
    }
  }

  # Keep genomic order.
  selected_idx <- sort(unique(as.integer(selected_idx)))
  center_rank <- match(center_pos, selected_idx)

  list(
    left = left,
    right = right,
    window_len = window_len,
    selected_count = length(selected_idx),
    center_rank = as.integer(center_rank),
    center_below_threshold = FALSE
  )
}

nupop_predict_gene <- function(fasta_path, species, model, use_chem) {
  # NuPoP writes output file in current working directory; run in temp dir.
  wd <- getwd()
  td <- tempfile("nupop_")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    setwd(wd)
    unlink(td, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  setwd(td)

  if (use_chem) {
    predNuPoP_chem(file = fasta_path, species = species, model = model)
  } else {
    predNuPoP(file = fasta_path, species = species, model = model)
  }

  base <- basename(fasta_path)
  out_file <- file.path(td, paste0(base, "_Prediction", model, ".txt"))
  if (!file.exists(out_file)) {
    # fallback in case package strips extension unexpectedly
    alt <- list.files(td, pattern = paste0("_Prediction", model, "\\.txt$"), full.names = TRUE)
    if (length(alt) == 0) stop("NuPoP output not found for ", fasta_path)
    out_file <- alt[[1]]
  }

  # Parse NuPoP output text directly for compatibility across NuPoP versions.
  pred <- tryCatch(
    read.table(out_file, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(pred)) {
    pred <- read.delim(out_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }

  nm <- names(pred)
  occ_idx <- which(tolower(nm) %in% c("occup", "occupancy"))
  if (length(occ_idx) == 0) {
    # Fallback: any column containing "occup"
    occ_idx <- grep("occup", tolower(nm), fixed = TRUE)
  }
  if (length(occ_idx) == 0) {
    stop("NuPoP output missing occupancy column in file: ", out_file)
  }
  as.numeric(pred[[occ_idx[1]]])
}

genes <- unique(na.omit(as.character(df[[gene_col]])))
genes <- genes[nzchar(genes)]
message("Genes to process: ", length(genes))

for (g in genes) {
  fasta_path <- file.path(fasta_dir, paste0(g, ".fasta"))
  if (!file.exists(fasta_path)) {
    warning("Skipping gene with missing FASTA: ", g)
    next
  }

  message("Running NuPoP for gene: ", g)
  gene_seq <- read_fasta(fasta_path)
  occup <- nupop_predict_gene(fasta_path, species = species, model = model, use_chem = use_chem)

  if (length(occup) != nchar(gene_seq)) {
    warning("Occupancy length mismatch for gene ", g, ": occup=", length(occup), ", seq=", nchar(gene_seq))
  }

  idx <- which(as.character(df[[gene_col]]) == g)
  for (i in idx) {
    seq9 <- clean_seq(df[[seq9_col]][i])
    if (is.na(seq9) || nchar(seq9) != 9) next

    starts <- all_fixed_matches_1based(seq9, gene_seq)
    df[["Seq9_Match_Count"]][i] <- length(starts)

    if (length(starts) == 0) next

    centers <- starts + 4L
    valid <- centers >= 1L & centers <= length(occup)
    centers <- centers[valid]
    starts <- starts[valid]
    if (length(centers) == 0) next

    if (length(centers) == 1) {
      chosen_center <- centers[1]
      chosen_start <- starts[1]
      val <- occup[chosen_center]
    } else if (multi_match_strategy == "mean") {
      chosen_center <- centers[1]
      chosen_start <- starts[1]
      val <- mean(occup[centers], na.rm = TRUE)
    } else if (multi_match_strategy == "na") {
      chosen_center <- NA_integer_
      chosen_start <- NA_integer_
      val <- NA_real_
    } else { # first
      chosen_center <- centers[1]
      chosen_start <- starts[1]
      val <- occup[chosen_center]
    }

    df[[occup_col]][i] <- val
    df[["Seq9_Selected_Start_1Based"]][i] <- chosen_start
    df[["Seq9_Selected_Center_1Based"]][i] <- chosen_center

    pos146 <- compute_center_position_in_selected146(occup, chosen_center, threshold = 0.1)
    df[["NuPoP_Window_Left_1Based"]][i] <- pos146$left
    df[["NuPoP_Window_Right_1Based"]][i] <- pos146$right
    df[["NuPoP_Window_Length"]][i] <- pos146$window_len
    df[["NuPoP_Selected146_Count"]][i] <- pos146$selected_count
    df[["NuPoP_Center_In_Selected146_1Based"]][i] <- pos146$center_rank
    # User-facing nucleosome position:
    # 1..146 (or <= selected count when fewer available), and 0 when center <= 0.1.
    df[[position_col]][i] <- as.integer(pos146$center_rank)
  }
}

message("Writing output: ", output_xlsx)
wb <- createWorkbook()
addWorksheet(wb, sheet_name)
writeData(wb, sheet = sheet_name, x = df)
saveWorkbook(wb, output_xlsx, overwrite = TRUE)

filled <- sum(!is.na(df[[occup_col]]))
message("Done. Filled ", filled, " occupancy values out of ", nrow(df), " rows.")
