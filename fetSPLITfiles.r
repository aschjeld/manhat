#!/usr/bin/env Rscript

########################################
# Plot split PoPoolation2 FET files with BH-FDR (PNG)
# Key behavior:
# - BH-FDR computed on ALL points
# - Plot shows ALL significant SNPs (q <= FDR_ALPHA)
# - Thins ONLY the nonsignificant background for readability
#
# Assumption: value after PAIR is already -log10(p)
# (This matches your prior log: q99(val) ~ 3.18 -> treated as -log10(p))
########################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------
# EDIT THESE
# -----------------------------
FETDIR  <- "/hb/groups/kay_lab/popoolation2/populations/LP4/sync_chunks"
FAI     <- "/hb/groups/kay_lab/popoolation2/populations/LP2/refgenome.fai"
PAIR    <- "1:2="

OUTDIR  <- file.path(FETDIR, "plots")
OUTPNG  <- file.path(OUTDIR, "fet_manhattan_FDR_keepSigThinBg.png")

# FDR level
FDR_ALPHA <- 0.05

# Optional filters (only applied if columns exist)
MIN_SNPS      <- 1
MIN_COV       <- -Inf
MIN_MINCOV    <- -Inf

# Plot controls
MAX_PLOT_POINTS <- 5e6          # total plotted points = all sig + sampled nonsig
PNG_WIDTH_IN    <- 16
PNG_HEIGHT_IN   <- 6
PNG_DPI         <- 450

# Optional y cap for readability (set Inf to disable)
Y_CAP <- Inf  # e.g. 200 if extreme hits squash everything else

pattern <- "\\.fet$"

# -----------------------------
# FIND FILES
# -----------------------------
stopifnot(dir.exists(FETDIR))
fet_files <- list.files(path = FETDIR, pattern = pattern, full.names = TRUE)
stopifnot(length(fet_files) > 0)
cat("Found", length(fet_files), ".fet files\n")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# READ .FAI for scaffold ordering and offsets
# -----------------------------
stopifnot(file.exists(FAI))
fai <- fread(FAI, header = FALSE)
setnames(fai, c("CHR","scaf_len","offset0","linebases","linewidth"))
fai <- fai[, .(CHR = as.character(CHR), scaf_len = as.numeric(scaf_len))]
setorder(fai, -scaf_len)                         # largest -> smallest
fai[, offset := cumsum(shift(scaf_len, fill = 0))]

# -----------------------------
# Robust read-one split .FET file + parse PAIR column
# -----------------------------
read_one_fet <- function(f) {
  dt <- suppressWarnings(
    tryCatch(
      fread(f, header = FALSE, fill = TRUE, showProgress = FALSE),
      error = function(e) NULL
    )
  )

  # Repair glued lines and re-read if needed
  if (is.null(dt) || ncol(dt) < 6) {
    lines <- readLines(f, warn = FALSE)
    # "...1:2=0ptg000001l..." -> "...1:2=0\nptg000001l..."
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

  # Extract numeric value after PAIR (ASSUMED already -log10(p))
  dt[, val := suppressWarnings(as.numeric(sub(paste0("^", PAIR), "", get(pair_col))))]
  dt <- dt[is.finite(val)]
  if (nrow(dt) == 0) return(NULL)

  dt[, .(CHR, POS, NSNPS, COV, AVG_MIN_COV, val)]
}

cat("Reading split .fet files...\n")
lst <- lapply(fet_files, read_one_fet)
lst <- Filter(Negate(is.null), lst)
stopifnot(length(lst) > 0)

dt <- rbindlist(lst, use.names = TRUE, fill = TRUE)
cat("Combined rows (finite val only):", nrow(dt), "\n")

# -----------------------------
# OPTIONAL FILTERS
# -----------------------------
before <- nrow(dt)
if ("NSNPS" %in% names(dt)) dt <- dt[is.na(NSNPS) | NSNPS >= MIN_SNPS]
if ("COV" %in% names(dt)) dt <- dt[is.na(COV) | COV >= MIN_COV]
if ("AVG_MIN_COV" %in% names(dt)) dt <- dt[is.na(AVG_MIN_COV) | AVG_MIN_COV >= MIN_MINCOV]
cat("Rows removed by filters:", before - nrow(dt), "\n")
stopifnot(nrow(dt) > 0)

# -----------------------------
# IMPORTANT: val is already -log10(p)
# Convert to p only for BH-FDR
# -----------------------------
dt[, logP := val]
dt[, p := 10^(-logP)]
dt[!is.finite(p) | p <= 0, p := .Machine$double.xmin]
dt[p > 1, p := 1]

cat("\nlogP range:\n"); print(range(dt$logP, na.rm = TRUE))
cat("p range (clamped):\n"); print(range(dt$p, na.rm = TRUE))

# -----------------------------
# BH-FDR across genome (ALL points)
# -----------------------------
cat("\nComputing BH-FDR across", nrow(dt), "tests...\n")
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
# cumulative coordinate
# -----------------------------
dt <- merge(dt, fai, by = "CHR", all.x = TRUE, sort = FALSE)
stopifnot(all(!is.na(dt$offset)))
dt[, BPcum := POS + offset]

# Cap y for plotting only
dt[, logP_plot := pmin(logP, Y_CAP)]
# -----------------------------
# Plot thinning: keep ALL significant, sample ONLY nonsignificant
# -----------------------------
if (is.finite(MAX_PLOT_POINTS) && nrow(dt) > MAX_PLOT_POINTS) {
  dt_sig    <- dt[q <= FDR_ALPHA]
  dt_nonsig <- dt[q >  FDR_ALPHA]

  remaining_budget <- MAX_PLOT_POINTS - nrow(dt_sig)

  if (remaining_budget <= 0) {
    set.seed(1)
    dt_plot <- dt_sig[sample.int(nrow(dt_sig), MAX_PLOT_POINTS)]
    cat("NOTE: Significant points exceed MAX_PLOT_POINTS; sampled significant only.\n")
  } else {
    set.seed(1)
    if (nrow(dt_nonsig) > remaining_budget) {
      dt_plot <- rbind(
        dt_sig,
        dt_nonsig[sample.int(nrow(dt_nonsig), remaining_budget)],
        use.names = TRUE
      )
      cat("Plotted ALL significant + sampled nonsignificant background.\n")
    } else {
      dt_plot <- dt
    }
  }

  cat("Plotting points:", nrow(dt_plot), " (sig kept:", nrow(dt_sig), ")\n")
} else {
  dt_plot <- dt
  cat("Plotting all points (no thinning needed).\n")
}

# -----------------------------
# PLOT (PNG)
# -----------------------------
pplot <- ggplot(dt_plot, aes(BPcum, logP_plot)) +
  geom_point(size = 0.22, alpha = 0.35) +
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
  pplot <- pplot + geom_hline(yintercept = min(fdr_line, Y_CAP), linewidth = 0.8, linetype = "dashed")
}

# Overlay significant points so they pop (still monochrome)
pplot <- pplot +
  geom_point(
    data = dt_plot[q <= FDR_ALPHA],
    aes(BPcum, logP_plot),
    size = 0.38, alpha = 0.95
  )

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
