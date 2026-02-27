#!/usr/bin/env Rscript

########################################
# Plot ONE PoPoolation2 FET file with BH-FDR (PNG)
# IMPORTANT: The value after PAIR (e.g., "1:2=") is already -log10(p).
# Robust to:
# - mixed tab/space separators
# - occasional "glued lines" (missing newline)
# - missing/irregular columns
# NON-DOWNSAMPLED: plots ALL points
########################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------
# EDIT THESE
# -----------------------------
FETFILE <- "/hb/groups/kay_lab/popoolation2/populations/LP3/p3p_p3w.fet"
FAI     <- "/hb/groups/kay_lab/popoolation2/populations/LP2/refgenome.fai"
PAIR    <- "1:2="

OUTDIR  <- "/hb/groups/kay_lab/popoolation2/populations/LP3"
OUTPNG  <- file.path(OUTDIR, "p3p_32w_fet_manhattan_FDR_ALLPOINTS.png")

FDR_ALPHA <- 0.05

# Optional filters (only applied if columns exist)
MIN_SNPS      <- 1
MIN_COV       <- -Inf
MIN_MINCOV    <- -Inf

# PNG controls
PNG_WIDTH_IN    <- 16
PNG_HEIGHT_IN   <- 6
PNG_DPI         <- 450

# Optional: cap y-axis for readability (set to Inf to disable)
Y_CAP <- Inf  # e.g., 200 or 300 if a few points dominate the y-axis

# -----------------------------
# Safety checks
# -----------------------------
stopifnot(file.exists(FETFILE))
stopifnot(file.exists(FAI))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Read .FAI and compute offsets (largest -> smallest)
# -----------------------------
fai <- fread(FAI, header = FALSE)
setnames(fai, c("CHR","scaf_len","offset0","linebases","linewidth"))
fai <- fai[, .(CHR = as.character(CHR), scaf_len = as.numeric(scaf_len))]
setorder(fai, -scaf_len)
fai[, offset := cumsum(shift(scaf_len, fill = 0))]
# -----------------------------
# Robust read-one FET file + parse PAIR column
# Expects:
# 1 CHR, 2 POS, 3 NSNPS, 4 COV, 5 AVG_MIN_COV, 6+ "i:j=value"
# where value is ALREADY -log10(p)
# -----------------------------
read_one_fet_file <- function(f) {
  dt <- suppressWarnings(
    tryCatch(
      fread(f, header = FALSE, fill = TRUE, showProgress = FALSE),
      error = function(e) NULL
    )
  )

  # Repair glued lines and re-read if needed
  if (is.null(dt) || ncol(dt) < 6) {
    lines <- readLines(f, warn = FALSE)
    lines <- gsub("(\\d)(ptg\\d)", "\\1\n\\2", lines, perl = TRUE)
    dt <- suppressWarnings(
      tryCatch(
        fread(text = paste(lines, collapse = "\n"), header = FALSE, fill = TRUE),
        error = function(e) NULL
      )
    )
  }

  if (is.null(dt) || ncol(dt) < 6) return(NULL)

  setnames(dt, 1:5, c("CHR","POS","NSNPS","COV","AVG_MIN_COV"))

  dt[, CHR := as.character(CHR)]
  dt[, POS := suppressWarnings(as.numeric(POS))]
  dt[, NSNPS := suppressWarnings(as.numeric(NSNPS))]
  dt[, COV := suppressWarnings(as.numeric(COV))]
  dt[, AVG_MIN_COV := suppressWarnings(as.numeric(AVG_MIN_COV))]

  dt <- dt[is.finite(POS)]
  if (nrow(dt) == 0) return(NULL)

  pair_cols <- names(dt)[6:ncol(dt)]
  if (length(pair_cols) == 0) return(NULL)

  has_pair <- sapply(pair_cols, function(col) any(grepl(paste0("^", PAIR), dt[[col]])))
  if (!any(has_pair)) return(NULL)
  pair_col <- pair_cols[which(has_pair)[1]]

  dt[, val := suppressWarnings(as.numeric(sub(paste0("^", PAIR), "", get(pair_col))))]
  dt <- dt[is.finite(val)]
  if (nrow(dt) == 0) return(NULL)

  dt[, .(CHR, POS, NSNPS, COV, AVG_MIN_COV, val)]
}

# -----------------------------
# Read the file
# -----------------------------
cat("Reading:", FETFILE, "\n")
dt <- read_one_fet_file(FETFILE)
stopifnot(!is.null(dt), nrow(dt) > 0)
cat("Rows (finite val only):", nrow(dt), "\n")

# -----------------------------
# Optional filters
# -----------------------------
before <- nrow(dt)
if ("NSNPS" %in% names(dt)) dt <- dt[is.na(NSNPS) | NSNPS >= MIN_SNPS]
if ("COV" %in% names(dt)) dt <- dt[is.na(COV) | COV >= MIN_COV]
if ("AVG_MIN_COV" %in% names(dt)) dt <- dt[is.na(AVG_MIN_COV) | AVG_MIN_COV >= MIN_MINCOV]
cat("Rows removed by filters:", before - nrow(dt), "\n")
stopifnot(nrow(dt) > 0)

# -----------------------------
# val is already -log10(p)
# -----------------------------
dt[, logP := val]

# Convert back to p only for BH-FDR
dt[, p := 10^(-logP)]
dt[!is.finite(p) | p <= 0, p := .Machine$double.xmin]
dt[p > 1, p := 1]

cat("\nlogP range:\n"); print(range(dt$logP, na.rm = TRUE))
cat("p range (clamped):\n"); print(range(dt$p, na.rm = TRUE))

# -----------------------------
# BH-FDR + threshold line (computed on ALL points)
# -----------------------------
dt[, q := p.adjust(p, method = "BH")]

sig <- dt[q <= FDR_ALPHA]
if (nrow(sig) == 0) {
  fdr_line <- NA_real_
  cat("No SNPs pass FDR <", FDR_ALPHA, "\n")
} else {
  p_cut <- max(sig$p, na.rm = TRUE)
  fdr_line <- -log10(p_cut)
  cat("FDR", FDR_ALPHA, "threshold p_cut =", format(p_cut, scientific = TRUE),
      " => -log10(p_cut) =", fdr_line, "\n")
}

# -----------------------------
# cumulative coordinate (from .fai)
# -----------------------------
dt <- merge(dt, fai, by = "CHR", all.x = TRUE, sort = FALSE)
stopifnot(all(!is.na(dt$offset)))
dt[, BPcum := POS + offset]

# Y cap only for plotting (optional)
dt[, logP_plot := pmin(logP, Y_CAP)]

# -----------------------------
# Plot (ALL points)
# -----------------------------
pplot <- ggplot(dt, aes(BPcum, logP_plot)) +
  geom_point(size = 0.24, alpha = 0.35) +
  labs(
    x = "Genomic position (scaffolds ordered largest → smallest)",
    y = paste0("-log10(p) from Fisher's Exact Test, pair ", PAIR)
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.margin = margin(18, 22, 18, 18),
    axis.line  = element_line(linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.9),
    axis.ticks.length = unit(3.0, "mm")
  )

if (is.finite(fdr_line)) {
  pplot <- pplot +
    geom_hline(yintercept = min(fdr_line, Y_CAP), linewidth = 0.8, linetype = "dashed")
}

# Overlay significant points so they pop
pplot <- pplot +
  geom_point(
    data = dt[q <= FDR_ALPHA],
    aes(BPcum, logP_plot),
    size = 0.38, alpha = 0.95
  )

# -----------------------------
# Save PNG with a high-quality device
# -----------------------------
cat("Saving PNG:", OUTPNG, "\n")

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    filename = OUTPNG,
    plot = pplot,
    device = ragg::agg_png,
    width = PNG_WIDTH_IN,
    height = PNG_HEIGHT_IN,
    units = "in",
    res = PNG_DPI
  )
  cat("Used ragg::agg_png\n")
} else {
  ggsave(
    filename = OUTPNG,
    plot = pplot,
    device = "png",
    type = "cairo-png",
    width = PNG_WIDTH_IN,
    height = PNG_HEIGHT_IN,
    units = "in",
    dpi = PNG_DPI
  )
  cat("Used cairo-png fallback\n")
}
cat("Done.\n")
