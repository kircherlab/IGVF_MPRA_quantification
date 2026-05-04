# Argument parser
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Process mpralm variant data")
parser$add_argument("--counts", type = "character", required = TRUE, help = "Path to the counts file")
parser$add_argument("--labels", type = "character", required = TRUE, help = "Path to the labels file")
parser$add_argument("--test-label", type = "character", required = TRUE, help = "Name of the test group")
parser$add_argument("--control-label", type = "character", required = TRUE, help = "Name of the control group")
parser$add_argument("--percentile",
  type = "double", default = 0.975,
  help = "Percentile of control to test on. Default is 0.975"
)
parser$add_argument("--output", type = "character", required = TRUE, help = "Path to the output file")
parser$add_argument("--output-volcano-plot", type = "character", required = FALSE, help = "Path to store the volcano plot")
parser$add_argument("--output-density-plot", type = "character", required = FALSE, help = "Path to store the density plot")
parser$add_argument("--normalize", type = "logical", default = TRUE, help = "Whether to normalize the data (TRUE or FALSE)")
parser$add_argument("--normalize-size",
  type = "double", default = 1e9,
  help = "Scaling factor for normalization (default is 1e9)"
)

args <- parser$parse_args()

suppressPackageStartupMessages(library(mpra))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))


mpra_treat <- function( # nolint: cyclocomp_linter.
  fit, percentile = 0.975, neg_label, trend = FALSE, robust = FALSE, winsor_tail_p = c(0.05, 0.1)
) {
  # 	Check fit
  if (is(fit, "MPRASet")) {
    mpra <- attr(fit, "MArrayLM")
    mpra$logFC <- rowData(fit)$logFC
    mpra$label <- getLabel(fit)
    fit <- mpra
  }
  if (!is(fit, "MArrayLM")) stop("fit must be an MArrayLM object")
  if (is.null(fit$coefficients)) stop("coefficients not found in fit object")
  if (is.null(fit$stdev.unscaled)) stop("stdev.unscaled not found in fit object")
  if (is.null(fit$label)) stop("Your mpra fit object should contain a label column.")

  fit$lods <- NULL

  neg_mean <- mean(fit$coefficients[fit$label == neg_label])

  coefficients <- as.matrix(fit$coefficients - neg_mean)
  coefficients_neg <- as.matrix(fit$coefficients[fit$label == neg_label] - neg_mean)

  stdev_unscaled <- as.matrix(fit$stdev.unscaled)
  sigma <- fit$sigma
  df_residual <- fit$df.residual
  if (is.null(coefficients) || is.null(stdev_unscaled) || is.null(sigma) || is.null(df_residual)) {
    stop("No data, or argument is not a valid lmFit object")
  }
  if (max(df_residual) == 0) {
    stop("No residual degrees of freedom in linear model fits")
  }
  if (!any(is.finite(sigma))) {
    stop("No finite residual standard deviations")
  }
  if (trend) {
    covariate <- fit$Amean
    if (is.null(covariate)) stop("Need Amean component in fit to estimate trend")
  } else {
    covariate <- NULL
  }
  sv <- squeezeVar(sigma^2, df_residual, covariate = covariate, robust = robust, winsor.tail.p = winsor_tail_p)
  fit$df.prior <- sv$df.prior
  fit$s2.prior <- sv$var.prior
  fit$s2.post <- sv$var.post
  df_total <- df_residual + sv$df.prior
  df_pooled <- sum(df_residual, na.rm = TRUE)
  df_total <- pmin(df_total, df_pooled)
  fit$df.total <- df_total

  acoef <- abs(coefficients)
  se <- stdev_unscaled * sqrt(fit$s2.post)
  lfc_right <- quantile(coefficients_neg, percentile)
  lfc_left <- quantile(coefficients_neg, 1 - percentile)
  tstat_right <- (acoef - lfc_right) / se
  tstat_left <- (acoef - lfc_left) / se
  fit$t <- array(0, dim(coefficients), dimnames = dimnames(coefficients))
  fit$p.value <- pt(tstat_right, df = df_total, lower.tail = FALSE) + pt(tstat_left, df = df_total, lower.tail = FALSE)
  tstat_right <- pmax(tstat_right, 0)
  tstat_left <- pmax(tstat_left, 0)
  fc_up <- (coefficients > lfc_right)
  fc_down <- (coefficients < lfc_left)
  fit$t[fc_up] <- tstat_right[fc_up]
  fit$t[fc_down] <- tstat_left[fc_down]
  fit$treat.lfc_right <- lfc_right
  fit$treat.lfc_left <- lfc_left
  fit
}


# read in the data
counts_df <- read.table(args$counts, header = TRUE, sep = "\t", fill = TRUE, c("", "NA", "N/A"))
colnames(counts_df)[1] <- c("ID")
counts_df <- counts_df %>% column_to_rownames(var = "ID")

dna_elem <- counts_df[, grepl("dna", colnames(counts_df))]
colnames(dna_elem) <- gsub("dna_", "", colnames(dna_elem))
rna_elem <- counts_df[, grepl("rna", colnames(counts_df))]
colnames(rna_elem) <- gsub("rna_", "", colnames(rna_elem))

labels <- read.table(args$labels, header = FALSE, sep = "\t", col.names = c("name", "label"))

labels_vec <- as.vector(labels$label)
names(labels_vec) <- labels$name
# Use only these labels of the sequences that remained after filtering
labels_vec <- labels_vec[rownames(dna_elem)]


# create the MPRASet object
mpraset <- MPRASet(
  DNA = dna_elem,
  RNA = rna_elem,
  eid = rownames(dna_elem),
  eseq = NULL,
  barcode = NULL,
)

# create the design matrix
design <- model.matrix(~1, data = data.frame(sample = seq_len(ncol(dna_elem))))


# run the mpralm analysis
fit_elem <- mpralm(
  object = mpraset,
  design = design,
  aggregate = "none",
  normalize = args$normalize,
  normalizeSize = args$normalize_size,
  model_type = "indep_groups",
  plot <- FALSE
)

toptab_element <- topTable(fit_elem, coef = 1, number = Inf)
percentile <- args$percentile

if (!is.null(args$output_density_plot)) {
  cat("Plot density elements...\n")

  toptab_element_label <- toptab_element %>%
    rownames_to_column(var = "name") %>%
    left_join(labels, by = "name") %>%
    column_to_rownames(var = "name")

  percentile_up <- quantile(toptab_element_label$logFC[toptab_element_label$label == args$control_label], percentile)
  up_label <- paste(percentile, "th percentile of negative controls", sep = "")

  percentile_down <- quantile(toptab_element_label$logFC[toptab_element_label$label == args$control_label], 1 - percentile)
  down_label <- paste(1 - percentile, "th percentile of negative controls", sep = "")


  density_plot <- ggplot(toptab_element_label, aes(x = logFC, fill = label, y = after_stat(density))) +
    geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.1) +
    geom_density(alpha = 0.2, adjust = 1) +
    labs(x = "log2 fold change", y = "Density") +
    xlim(c(min(toptab_element_label$logFC), max(toptab_element_label$logFC))) +
    geom_vline(aes(xintercept = percentile_up, color = up_label), linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = percentile_down, color = down_label), linetype = "dashed", linewidth = 1) +
    scale_color_manual(
      values = setNames(c("green", "orange"), c(up_label, down_label)),
      guide = guide_legend(override.aes = list(linetype = "dashed"))
    ) +
    theme_minimal()

  ggsave(filename = args$output_density_plot, plot = density_plot, width = 8, height = 6)
}


# Re-evaluate
# tr <- treat(fit_elem, lfc = percentile_up)
fit_elem$label <- labels_vec
tr <- mpra_treat(fit_elem, percentile, neg_label = args$control_label)
mpra_element <- topTreat(tr, coef = 1, number = Inf)

# Make volcano plot with cutoff of FDR < 0.01
if (!is.null(args$output_volcano_plot)) {
  cat("Plot volcano...\n")
  p <- ggplot(mpra_element, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    geom_point(data = subset(mpra_element, adj.P.Val < 0.01), aes(x = logFC, y = -log10(adj.P.Val)), color = "red") +
    labs(x = "log2 fold change", y = "-log10(p-value)") +
    theme_minimal()

  ggsave(filename = args$output_volcano_plot, plot = p, width = 8, height = 6)
}


names <- c("ID", colnames(mpra_element))
mpra_element$ID <- rownames(mpra_element)
mpra_element <- mpra_element[, names]

cat("Write output to file...\n")
gzfile_output <- gzfile(args$output, "w")
write.table(mpra_element, gzfile_output, row.names = FALSE, sep = "\t", quote = FALSE)
close(gzfile_output)
