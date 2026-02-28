#!/usr/bin/env Rscript

########################################
# SNP-level PoPoolation2 FET Manhattan plot (PNG)
# Shared genome axis + rigorous SNP slimming + working labels (requires ggrepel)
#
# Designed for your last run characteristics:
# - ~69.7M SNPs total across 469 .fet files
# - BH-FDR yields millions of “significant” SNPs (too many to interpret)
# - Therefore: use SNP-level Bonferroni for “genome-wide significant SNPs”
# - Then: clump Bonferroni hits by physical distance per scaffold to get independent lead SNPs
# - Plot: sampled background + highlighted Bonferroni hits + labeled lead SNPs
#
# Assumption: value after PAIR is already -log10(p)
########################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Require ggrepel up front so you don't waste a big run with no labels
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("ggrepel is not installed in this R environment. Install: conda install -c conda-forge r-ggrepel")
}

# -----------------------------
# EDIT THESE
# -----------------------------
FETDIR   <- "/hb/groups/kay_lab/popoolation2/populations/LP4/sync_chunks"
OFFSETS  <- "/hb/groups/kay_lab/popoolation2/populations/LP2/genome_offsets.tsv"
PAIR     <- "1:2="

OUTDIR   <- file.path(FETDIR, "plots")
OUTPNG   <- file.path(OUTDIR, "fet_manhattan_SNP_BONF_clumped_labeled.png")
OUTHITS  <- file.path(OUTDIR, "fet_lead_hits_BONF.tsv")

# -----------------------------
# Multiple testing and slimming
# -----------------------------
BONF_ALPHA <- 0.05          # genome-wide Bonferroni alpha

# Clumping: keep only one lead SNP per +/- CLUMP_BP region per scaffold (for labels/table)
CLUMP_BP   <- 50000L        # try 50k or 100k
LABEL_TOP_N <- 40           # label only top N lead SNPs after clumping (0 disables)

# Optional filters (only applied if columns exist)
MIN_SNPS      <- 1
MIN_COV       <- -Inf
MIN_MINCOV    <- -Inf

# Plot controls
MAX_PLOT_POINTS <- 1e7      # sampled background points for readability
PNG_WIDTH_IN    <- 16
PNG_HEIGHT_IN   <- 6
PNG_DPI         <- 450

# Optional y cap for readability (set Inf to disable)
Y_CAP <- Inf

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
# READ genome offsets (shared coordinate system)
# -----------------------------
stopifnot(file.exists(OFFSETS))

genome <- fread(
  OFFSETS,
  sep = "\t",
  colClasses = list(
    character = "CHR",
    numeric   = c("LEN", "offset", "end")
  ),
  showProgress = FALSE
)

stopifnot(all(c("CHR","LEN","offset","end") %in% names(genome)))
genome <- genome[is.finite(end) & end > 0]
stopifnot(nrow(genome) > 0)

genome_end <- max(genome$end, na.rm = TRUE)

cat("Loaded genome offsets:", nrow(genome), "scaffolds\n")
cat("Shared genome_end:", sprintf("%.0f", genome_end), "\n")
stopifnot(genome_end > 1e6)

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
  dt[, logP := suppressWarnings(as.numeric(sub(paste0("^", PAIR), "", get(pair_col))))]
  dt <- dt[is.finite(logP)]
  if (nrow(dt) == 0) return(NULL)

  # Optional filters early (cheap)
  if ("NSNPS" %in% names(dt)) dt <- dt[is.na(NSNPS) | NSNPS >= MIN_SNPS]
  if ("COV" %in% names(dt)) dt <- dt[is.na(COV) | COV >= MIN_COV]
  if ("AVG_MIN_COV" %in% names(dt)) dt <- dt[is.na(AVG_MIN_COV) | AVG_MIN_COV >= MIN_MINCOV]
  if (nrow(dt) == 0) return(NULL)

  # Convert to p (for Bonferroni thresholding)
  dt[, p := 10^(-logP)]
  dt[!is.finite(p) | p <= 0, p := .Machine$double.xmin]
  dt[p > 1, p := 1]

  # Keep only required columns to reduce memory
  dt[, .(CHR, POS, logP, p)]
}

cat("Reading split .fet files...\n")
lst <- lapply(fet_files, read_one_fet)
lst <- Filter(Negate(is.null), lst)
stopifnot(length(lst) > 0)

dt <- rbindlist(lst, use.names = TRUE, fill = TRUE)
cat("Combined rows:", nrow(dt), "\n")

cat("\nlogP range:\n"); print(range(dt$logP, na.rm = TRUE))
cat("p range:\n"); print(range(dt$p, na.rm = TRUE))

# -----------------------------
# Map onto shared genome axis
# -----------------------------
dt <- merge(dt, genome[, .(CHR, offset)], by = "CHR", all.x = TRUE, sort = FALSE)
if (any(is.na(dt$offset))) {
  missing <- unique(dt[is.na(offset), CHR])
  stop("Some CHR values were not found in genome_offsets.tsv. Examples: ",
       paste(head(missing, 10), collapse = ", "))
}
dt[, BPcum := POS + offset]
dt[, logP_plot := pmin(logP, Y_CAP)]

cat("BPcum range:\n"); print(range(dt$BPcum, na.rm = TRUE))

# -----------------------------
# Bonferroni threshold (SNP-level, genome-wide)
# -----------------------------
m <- nrow(dt)
bonf_p <- BONF_ALPHA / m
bonf_logp <- -log10(bonf_p)

cat("\nBonferroni alpha =", BONF_ALPHA,
    " => p <=", format(bonf_p, scientific = TRUE),
    " => -log10(p) >=", bonf_logp, "\n")

sig_bonf <- dt[p <= bonf_p]
cat("Bonferroni-significant SNPs:", nrow(sig_bonf), "\n")

# -----------------------------
# Clump significant SNPs for independent lead hits (for labels/table)
# -----------------------------
lead <- sig_bonf
if (nrow(lead) > 0) {
  setorder(lead, p)  # most significant first

  if (is.finite(CLUMP_BP) && CLUMP_BP > 0) {
    out <- lead[0]

    for (chr in unique(lead$CHR)) {
      x <- lead[CHR == chr]
      kept_pos <- numeric(0)

      keep_idx <- logical(nrow(x))
      for (i in seq_len(nrow(x))) {
        pos <- x$POS[i]
        if (length(kept_pos) == 0 || all(abs(pos - kept_pos) > CLUMP_BP)) {
          keep_idx[i] <- TRUE
          kept_pos <- c(kept_pos, pos)
        }
      }
      out <- rbind(out, x[keep_idx], use.names = TRUE)
    }
    lead <- out
  }

  # Create labels and keep only top N for labeling
  lead[, label := paste0(CHR, ":", POS)]
  setorder(lead, p)

  # Write lead hits table
  fwrite(lead[, .(CHR, POS, BPcum, logP, p)], OUTHITS, sep = "\t")
  cat("Wrote lead hits:", OUTHITS, "\n")
} else {
  # still write an empty file so downstream doesn't break
  fwrite(data.table(CHR=character(), POS=integer(), BPcum=numeric(), logP=numeric(), p=numeric()),
         OUTHITS, sep="\t")
  cat("Wrote empty lead hits:", OUTHITS, "\n")
}

lead_to_label <- lead[0]
if (LABEL_TOP_N > 0 && nrow(lead) > 0) {
  lead_to_label <- lead[1:min(LABEL_TOP_N, .N)]
  cat("Labeling", nrow(lead_to_label), "lead SNPs\n")
}

# -----------------------------
# Background sampling for plotting
# -----------------------------
if (is.finite(MAX_PLOT_POINTS) && nrow(dt) > MAX_PLOT_POINTS) {
  set.seed(1)
  dt_plot <- dt[sample.int(nrow(dt), MAX_PLOT_POINTS)]
  cat("Plotting sampled background points:", nrow(dt_plot), "\n")
} else {
  dt_plot <- dt
  cat("Plotting all points (no sampling needed).\n")
}

# -----------------------------
# PLOT
# -----------------------------
pplot <- ggplot(dt_plot, aes(BPcum, logP_plot)) +
  geom_point(size = 0.22, alpha = 0.22) +
  scale_x_continuous(limits = c(0, genome_end)) +
  labs(
    x = "Genomic position (shared coordinate; scaffolds ordered largest → smallest)",
    y = paste0("-log10(p) from Fisher's Exact Test, pair ", PAIR)
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.margin = margin(18, 22, 18, 18),
    axis.line  = element_line(linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.9),
    axis.ticks.length = unit(3.0, "mm")
  )

# Bonferroni line
pplot <- pplot +
  geom_hline(yintercept = min(bonf_logp, Y_CAP), linewidth = 0.8, linetype = "dotted")

# Overlay Bonferroni significant SNPs (all of them, if any)
if (nrow(sig_bonf) > 0) {
  pplot <- pplot +
    geom_point(
      data = sig_bonf,
      aes(BPcum, logP_plot),
      size = 0.6,
      alpha = 0.95
    )
}

# Label clumped lead SNPs (top N)
if (nrow(lead_to_label) > 0) {
  pplot <- pplot +
    ggrepel::geom_text_repel(
      data = lead_to_label,
      aes(BPcum, logP_plot, label = label),
      size = 2.8,
      max.overlaps = Inf,
      min.segment.length = 0
    )
}

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
