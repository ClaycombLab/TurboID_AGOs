## ============================================================================
## ppi_dotplot.R
##
##  This script was written with help from Claude AI
##
## Replicates a DotPlotViz-style "PPI dot plot":
##   - x axis  = Bait / pulldown (e.g. CSR-1b, HRDE-1, WAGO-1, PPW-1, WAGO-4)
##   - y axis  = prey gene, grouped into functional categories
##   - size    = relative average spectral count
##   - fill    = average spectral count (intensity)
##   - stroke colour = BFDR (Bayesian false discovery rate)
##
## Gene inclusion follows ProHits-viz's own rule (see their Dot Plot help,
## https://prohits-viz.org/help/analysis/dotplot): a gene is kept if it
## passes the BFDR threshold against >=1 Bait, and once kept, it can still
## appear against every Bait -- but only the Bait-prey pairs that
## individually pass a BFDR bin are actually drawn; non-passing pairs are
## left blank rather than shown faded.
##
## HOW THIS FILE IS ORGANIZED:
##   Section 1  -- EDIT THIS: your file paths and settings. This is the
##                 only section you should need to touch to run your own
##                 data.
##   Section 2  -- function definitions. You shouldn't need to edit these
##                 unless you want to change how something behaves.
##   Section 3  -- runs the pipeline using Section 1's settings and saves
##                 the plot(s). Also nothing to edit here normally.
## ============================================================================

## ---- 0. Packages -----------------------------------------------------------
## install.packages(c("ggplot2", "dplyr", "tidyr", "forcats", "svglite"))
## svglite is needed for ggsave() to write .svg files. Note on fonts: SVG
## text stays as real, editable text with a font-family CSS hint (not
## embedded font data), so "Arial" below doesn't need to be installed on
## THIS machine -- it just needs to be available wherever you later open
## or view the SVG (browser, Illustrator, Inkscape, etc.).
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)


## ============================================================================
## ---- 1. YOUR SETTINGS -- edit this section, then Source the whole file ----
## ============================================================================
##
## WHERE TO PUT YOUR FILES: the simplest setup is to put this script and
## your data files in the SAME folder, then just use plain file names below
## (no path needed) -- exactly as shown. R reads/writes relative to its
## "working directory", which is NOT necessarily the folder this script is
## in unless you set it. In RStudio: Session > Set Working Directory >
## To Source File Location. From the command line, run Rscript from that
## folder, or `cd` there first.
##
## If you'd rather keep files elsewhere, give a full path instead, e.g.
##   data_file <- "C:/Users/you/Documents/data.txt"          (Windows)
##   data_file <- "/Users/you/Documents/data.txt"             (Mac)

## Your measurements: a long-format, tab-delimited file, one row per
## (gene, Bait) pair, with (at minimum) these column names:
##   PublicName - prey gene/protein name (renamed to "gene" on load below)
##   Bait       - Bait / pulldown name
##   AvgSpec    - average spectral count (numeric)
##   BFDR       - Bayesian FDR (numeric 0-1, drives stroke colour)
## Relative spectral count (dot size) is NOT expected as a column -- it's
## computed automatically in Section 3 via add_relative_spc() (Section 2d):
## for each gene, its own maximum AvgSpec (among the Bait values that make
## it into the final plot) gets the largest dot, its own minimum gets the
## smallest, scaled linearly in between.
data_file <- "[yourdatahere].txt"

## Your category list: a tab- or comma-delimited text file where each
## COLUMN HEADER is a category name and the cells below are the gene names
## in that category (see read_gene_categories() in Section 2 for the exact
## format and how multi-category genes work). This ONLY decides which
## genes are included and which category/categories each belongs to -- it
## does NOT set row order. Row order is always determined by hierarchical
## clustering (see plot_ppi_dotplot()'s order_genes/cluster_by in Section 2)
## unless you explicitly pass gene_order yourself in Section 3.
## Set category_file to NULL if your data_file already has its own
## "category" column and you don't need a separate file.
category_file <- "[yourcategories].txt"

## TRUE  = save one plot per category, each sized to just that category's
##         genes (avoids one category's row count squishing another's).
## FALSE = save a single combined plot with every category as a facet
##         panel (the original behaviour). With use_custom_BFDR_colours
##         = TRUE, combined mode now ALSO gets genuine per-category
##         colouring (see the note there) -- the old single-palette
##         fallback is gone.
plot_separately <- T

## Where to save the finished plot(s). Files are written as SVG.
##   plot_separately = TRUE  -> used as a prefix:
##                               "<output_prefix>_<category>_<colour_tag>.svg"
##   plot_separately = FALSE -> "<output_prefix>_<colour_tag>.svg"
## colour_tag is "grey" normally, or "custom" whenever
## use_custom_BFDR_colours (below) is TRUE -- so switching palettes never
## silently overwrites the other palette's files.
#If you are plotting subsets of baits from the same data file
output_prefix <- ""

## TRUE  = draw a rounded-corner border around each panel (see
##         add_rounded_panel_border() in Section 2).
## FALSE = no border at all, just the bare panel -- e.g. if you don't
##         want every category panel "boxed in" for a combined plot.
## Applies to both plot_separately modes. In combined mode with
## use_custom_BFDR_colours = TRUE, each category panel's border matches
## that category's own accent colour, same as in per-category mode.
use_rounded_borders <- T

## BFDR stroke-colour dictionary: a named vector, BFDR bin label -> colour
## (see BFDR_colours in plot_ppi_dotplot(), Section 2f, for the full
## rules -- e.g. how a "> 0.05" entry changes whether non-passing points
## are omitted or drawn faint). This is the DEFAULT dictionary, used
## whenever use_custom_BFDR_colours (below) is FALSE.
default_BFDR_colours <- c("<= 0.01" = "black",
                          "< 0.05"  = "blue",
                          "= 0.05"  = "#809fff")

## Set TRUE to use per-category colours below instead of
## default_BFDR_colours. Output filenames automatically pick up a
## "custom" tag instead of "grey" so the two versions never collide.
## Works in BOTH plot_separately modes: in combined mode, every
## category panel gets ITS OWN colour gradient (built from a combined
## "category: bin" colour key under the hood -- see build_combined_BFDR_colours()
## in Section 2g and its use in Section 3), not a single shared palette.
use_custom_BFDR_colours <- TRUE

## Each category gets its own base colour, which becomes the STRONGEST
## bin's colour (darkened slightly first for visual headroom -- see
## category_BFDR_palette() in Section 2f-ter for the exact amounts and
## why). Weaker bins fade toward WHITE from there, so the gradient
## reads as "how much colour is left" = "how confident is this hit".
## The same colour also becomes that category's rounded panel border
## (see add_rounded_panel_border() in Section 2, and its use in
## Section 3):
##   "<= 0.01" (most significant)  -> the category's colour (darkened)
##   "< 0.05"                      -> that colour, 50% lighter
##   "= 0.05"  (least significant) -> that colour, 80% lighter
## Add a "Category Name" = "#hexcolour" pair for each category you want
## a specific colour for. Only takes effect when
## use_custom_BFDR_colours = TRUE.
category_colours <- c(
  "RNA Binding"  = "#0072b2",
  "RNA helicase" = "#009e73",
  "Nuclear Pore and Transport" = "#cc79a7",
  "Integrator Complex"  = "#ba8a80",
  "Splicing"  = "#56b4e9",
  "Germ Granule"  = "#a76add",
  "sRNA"  = "#d55e00",
  "Transcription" = "#634cd3",
  "Ribosome" ="#c44e52",
  "DNA Damage" = "#f28e2b",
  "Stress Response" = "#e69f00",
  "Mitochondrial" = "#e377c2",
  "Chromosome Associated"="#55a868",
  "Cytoskeleton"="#999933",
  "Cell Cycle"="#88ccee",
  "Kinase or Phosphatase"="#882255"
  
)

## Base colour for any category NOT listed in category_colours above.
default_category_colour <- "#809fff"

## Bait order, left to right on the x axis.
Bait_order <- c("CSR-1", "WAGO-4", 
  #"WAGO-4_dis"#, 
  "HRDE-1", "PRG-1", "WAGO-1", "PPW-2"#,
  #"ZNFX-1", "EGO-1", "RRF-1", 
  #"D2005.4"#,
  #"GLH-1","DEPS-1"
)

## Plot title (only used when plot_separately = FALSE; per-category plots
## are titled with their own category name).
#plot_title <- "testing"

## Primary filter: a gene is kept if its BFDR is <= this for at least one
## Bait (see apply_primary_filter() in Section 2).
BFDR_threshold <- 0.05

## min_abundance: gene-level inclusion filter -- see apply_min_abundance()
## in Section 2. Genuinely removes genes from the data.
## max_abundance: NOT a filter and NOT applied to the data at all --
## only sets where the fill colour scale tops out (fill_limits in
## plot_ppi_dotplot(), Section 3). A gene's true AvgSpec still drives its
## dot size and clustering regardless of this value; only its COLOUR
## gets capped at whatever colour max_abundance maps to.
min_abundance <- 0
max_abundance <- 100

## Inches of vertical space per gene row -- the main lever if your plot
## comes out too cramped or too spread out (see suggested_plot_height()).
row_height_in <- 0.29

## Inches of horizontal space per Bait column -- raise this if dots are
## overlapping column-wise (a fixed image width divided across many
## Baits can leave each column narrower than the dots drawn in it; see
## suggested_plot_width()). Should generally be at least as wide as your
## largest dot -- roughly size_range[2] in plot_ppi_dotplot(), converted
## from mm to inches (mm / 25.4).
col_width_in <- 0.33


## ============================================================================
## ---- 2. Functions -- shouldn't need to edit below here --------------------
## ============================================================================

## ---- 2a. Reading a category-list text file ----------------------------------
##
## Use this when you keep your category membership in a plain text file
## where each COLUMN HEADER is a category name and the rows beneath it are
## the gene names in that category (columns can be different lengths --
## ragged/empty cells are fine). Tab- or comma-delimited.
##
## Example file (tab-delimited):
##
##   Argonaute pathway   RNA surveillance   P-granule
##   csr-1                smg-1              pgl-1
##   wago-1                smg-2              pgl-3
##   wago-4                                    deps-1
##   csr-1
##
## Note csr-1 appears under two categories above -- that's exactly the
## multi-category case, and it's handled automatically: it just shows up
## as two rows in the long-format output below, one per category.
##
## This does double duty:
##   1. It defines which genes get plotted at all (anything NOT in this
##      file is dropped later via inner_join -- see Section 3).
##   2. Row order within each column would become the default top-to-bottom
##      gene order for that category's panel (see gene_order_list()) --
##      unused by default in Section 3, since clustering decides order
##      instead, but still available if you want it.
read_gene_categories <- function(path, delim = NULL) {
  if (is.null(delim)) {
    first_line <- readLines(path, n = 1)
    delim <- if (grepl("\t", first_line)) "\t" else ","
  }
  raw <- read.delim(path, header = TRUE, sep = delim, check.names = FALSE,
                    stringsAsFactors = FALSE, fill = TRUE,
                    colClasses = "character")
  
  long <- raw %>%
    mutate(.row = row_number()) %>%
    pivot_longer(-.row, names_to = "category", values_to = "gene") %>%
    mutate(gene = trimws(gene)) %>%
    filter(!is.na(gene), gene != "") %>%
    arrange(category, .row) %>%
    select(category, gene)
  
  n_before <- nrow(long)
  long <- distinct(long, category, gene, .keep_all = TRUE)
  n_dropped <- n_before - nrow(long)
  if (n_dropped > 0) {
    message(sprintf(
      "read_gene_categories: dropped %d duplicate (category, gene) row(s) -- the same gene was listed more than once under the same category column",
      n_dropped))
  }
  
  long
}

## Turns the long mapping table above into the named list format that
## plot_ppi_dotplot()'s `gene_order` argument expects, preserving the
## file's row order within each category. Not used by default in
## Section 3 (clustering handles order instead) -- use this yourself if
## you ever want to pin row order back to the file's order.
gene_order_list <- function(gene_categories) {
  split(gene_categories$gene, gene_categories$category)
}

## ---- 2b. Hierarchical clustering of genes by Bait-hit pattern ---------------
##
## Builds a gene x Bait matrix from one category's data and returns gene
## names in a top-to-bottom order (see plot_ppi_dotplot() for how this
## maps onto the y axis) split into two blocks:
##
##   1. Genes with MORE THAN ONE Bait hit -- hierarchically clustered by
##      profile similarity, THEN the dendrogram is rotated (without
##      changing which genes cluster together) so genes with more Bait
##      interactions, higher spectral counts, and more significant BFDR
##      trend toward the top of this block.
##   2. Genes with EXACTLY ONE Bait hit -- a single value has nothing to
##      cluster against, so instead of scattering these arbitrarily among
##      the clustered genes, they're grouped by which Bait they hit
##      (following Bait_order) and sorted by abundance within that Bait,
##      then placed at the bottom as their own block.
##
## cluster_by = "presence": binary distance on which Baits each gene hits
##              at all -- purely about which proteins it's a hit against.
## cluster_by = "AvgSpec" (default): canberra distance on the actual
##              spectral count values -- ProHits-viz's own default metric,
##              well suited to abundance data since it weights differences
##              relative to magnitude rather than in absolute terms.
## method: hclust linkage. Defaults to "ward.D2" -- ProHits-viz's own
##         default, and the correct R implementation of Ward's method
##         (Murtagh & Legendre 2014 recommend ward.D2 over the older
##         ward.D, which silently assumes pre-squared distances).
cluster_gene_order <- function(sub, cluster_by = c("AvgSpec", "presence"),
                               method = "ward.D2") {
  cluster_by <- match.arg(cluster_by)
  
  ## By this point (gene, Bait) pairs are guaranteed unique -- see the
  ## check in Section 3, right after your data is loaded -- so a plain
  ## pivot_wider() is safe here. select() first is essential, not just
  ## tidy: sub carries other columns (category, rel_spc, BFDR) that vary
  ## per row, and pivot_wider() treats any column it isn't told about as
  ## part of the id -- without narrowing to just these 3 first, it won't
  ## collapse to one row per gene, and a leftover character column
  ## (category) would force the whole matrix to character type.
  wide <- sub %>%
    select(gene, Bait, AvgSpec) %>%
    tidyr::complete(gene, Bait, fill = list(AvgSpec = 0)) %>%
    pivot_wider(names_from = Bait, values_from = AvgSpec, values_fill = 0)
  
  gene_names <- wide$gene
  mat <- as.matrix(wide[, setdiff(names(wide), "gene"), drop = FALSE])
  rownames(mat) <- gene_names
  
  ## tidyr::complete()'s fill only covers combinations it ADDED -- an
  ## AvgSpec cell that was already blank/NA in your data file (an
  ## existing row with a missing value, not a missing row) passes
  ## through untouched and would otherwise turn hit_counts, and
  ## everything downstream, into NA. Treat a missing value as no hit.
  mat[is.na(mat)] <- 0
  
  hit_counts  <- rowSums(mat > 0)
  multi_names <- gene_names[hit_counts > 1]
  single_names <- gene_names[hit_counts <= 1]
  
  cluster_subset <- function(names_subset) {
    if (length(names_subset) < 3) return(names_subset)  # nothing meaningful to cluster
    sub_mat <- mat[names_subset, , drop = FALSE]
    if (cluster_by == "presence") {
      d <- dist((sub_mat > 0) * 1, method = "binary")
    } else {
      d <- dist(sub_mat, method = "canberra")
    }
    ## dist() on all-zero rows yields NA distances -- fall back to
    ## unclustered order for this subset rather than erroring.
    if (anyNA(d)) return(names_subset)
    hc <- hclust(d, method = method)
    
    ## hclust's leaf order only reflects merge order, not "strength" --
    ## two genes can cluster tightly together yet both be weak hits, and
    ## nothing about hclust() puts strong hits near either end. Reorder
    ## the dendrogram itself (rotating branches at each merge point,
    ## WITHOUT changing which genes cluster together) so genes with more
    ## Bait interactions, higher spectral counts, and lower (more
    ## significant) BFDR trend toward the top.
    gene_stats <- sub %>%
      filter(gene %in% names_subset) %>%
      group_by(gene) %>%
      summarise(
        n_hits    = sum(AvgSpec > 0, na.rm = TRUE),
        total_spc = sum(AvgSpec, na.rm = TRUE),
        best_bfdr = min(BFDR, na.rm = TRUE),
        .groups   = "drop"
      ) %>%
      mutate(score = rank(n_hits) + rank(total_spc) + rank(-best_bfdr))
    score_lookup <- setNames(gene_stats$score, gene_stats$gene)
    
    ## There is no reorder() method for class "hclust" in base R --
    ## calling reorder() directly on hc silently falls through to
    ## reorder.default(), which is for factors and errors out on a
    ## dendrogram-shaped object. reorder.dendrogram() is the real base R
    ## method for this: it rotates branches at each merge point using an
    ## external per-leaf weight, WITHOUT changing which genes cluster
    ## together.
    dend <- as.dendrogram(hc)
    dend <- reorder(dend, score_lookup[names_subset], agglo.FUN = mean)
    ord  <- names_subset[order.dendrogram(dend)]
    
    ## reorder()'s left/right convention isn't guaranteed -- confirm
    ## high-scoring genes actually landed first (top of plot) and flip
    ## the whole block if they didn't.
    half <- ceiling(length(ord) / 2)
    if (mean(score_lookup[head(ord, half)], na.rm = TRUE) <
        mean(score_lookup[tail(ord, half)], na.rm = TRUE)) {
      ord <- rev(ord)
    }
    
    ord
  }
  
  multi_order <- cluster_subset(multi_names)
  
  single_order <- if (length(single_names) == 0) {
    character(0)
  } else {
    single_mat  <- mat[single_names, , drop = FALSE]
    which_bait  <- colnames(mat)[apply(single_mat > 0, 1, function(r) {
      idx <- which(r)
      if (length(idx) == 0) NA_integer_ else idx[1]
    })]
    strength <- apply(single_mat, 1, max)
    single_names[order(factor(which_bait, levels = colnames(mat)), -strength)]
  }
  
  c(multi_order, single_order)
}

## ---- 2c. Primary filter: gene-level inclusion, not point-level -------------
##
## Matches ProHits-viz's own filtering rule: a gene is included in the plot
## if it passes the score threshold (e.g. BFDR <= 0.05) against AT LEAST
## ONE Bait. Once included, ALL of that gene's values are kept -- but only
## the Bait-prey pairs that individually pass a BFDR bin (see
## plot_ppi_dotplot()'s BFDR_colours) are actually drawn; the rest are left
## blank. Genes that don't pass for any Bait are dropped entirely.
apply_primary_filter <- function(data, threshold = 0.05, score_col = "BFDR") {
  keep_genes <- data %>%
    filter(.data[[score_col]] <= threshold) %>%
    distinct(gene) %>%
    pull(gene)
  data %>% filter(gene %in% keep_genes)
}

## ---- 2d. Minimum / maximum abundance -----------------------------------------
##
## Minimum abundance: gene-level inclusion, same pattern as the primary
## filter above -- a gene must have an abundance value above this minimum
## for AT LEAST ONE Bait to be included; once included, all its values
## are kept. min_abundance = 0 is effectively a no-op unless you have
## real zero-count rows in your data.
apply_min_abundance <- function(data, min_abundance = 0, abundance_col = "AvgSpec") {
  keep_genes <- data %>%
    filter(.data[[abundance_col]] > min_abundance) %>%
    distinct(gene) %>%
    pull(gene)
  data %>% filter(gene %in% keep_genes)
}

## Maximum abundance: NOT used by default (see Section 3) -- unlike the
## filters above, this genuinely OVERWRITES AvgSpec, clamping values
## above max_abundance down to it. That's a real, permanent change to
## the working data for the rest of the pipeline: add_relative_spc()'s
## per-gene min/max and cluster_gene_order()'s distance calculations
## would both see the clamped value, not the true one. Visual-only
## capping (which is what Section 3 actually uses, via fill_limits +
## oob = scales::squish in plot_ppi_dotplot()) doesn't touch the data at
## all -- it just tops out the colour mapping at the plotting stage.
## Kept here in case you ever want the data itself clamped, e.g. for a
## table or export where the same numbers should match what's drawn.
cap_max_abundance <- function(data, max_abundance = 100, abundance_col = "AvgSpec") {
  data %>% mutate(!!abundance_col := pmin(.data[[abundance_col]], max_abundance))
}

## ---- 2d-bis. Relative spectral count (drives dot size) -----------------------
##
## Not a column in the raw data -- computed here per gene, from whatever
## AvgSpec values are left AFTER your filters run (primary BFDR filter,
## min/max abundance). For each gene: its own maximum AvgSpec (among the
## Baits that survive filtering) scales to 1 (largest dot), its own
## minimum scales to 0 (smallest dot), linearly in between. This is
## per-gene, not global -- a weak hit for one gene and a strong hit for
## another can both show up as a "large" dot if it's the strongest value
## THAT gene has, which is what lets you compare a gene's relative
## preference across Baits regardless of its absolute abundance.
##
## A gene with only one surviving AvgSpec value has no range to scale
## against -- it's treated as its own maximum and gets rel_spc = 1
## (drawn at the full size_range[2] diameter), regardless of that one
## value's absolute magnitude.
##
## IMPORTANT: rel_spc's OWN range is always exactly [0, 1] by construction
## (every gene's own min maps to 0, own max to 1) -- but ggplot's
## scale_radius() doesn't know that unless told explicitly (see
## plot_ppi_dotplot()'s scale_radius(..., limits = c(0, 1))). Without
## it, a category where every gene happens to be single-hit (so every
## rel_spc value is exactly 1, with no 0 anywhere in that category's
## data) has its INPUT DOMAIN auto-detected as a single point (1, 1),
## and ggplot/scales silently maps a collapsed domain to the MIDPOINT of
## the output range rather than its max -- producing visibly smaller
## dots in that category than in one with real size variation, even
## though both used identical rel_spc = 1 values and identical
## size_range settings.
add_relative_spc <- function(data, abundance_col = "AvgSpec", output_col = "rel_spc") {
  data %>%
    group_by(gene) %>%
    mutate(
      .gene_max = max(.data[[abundance_col]], na.rm = TRUE),
      .gene_min = min(.data[[abundance_col]], na.rm = TRUE),
      !!output_col := if_else(
        .gene_max > .gene_min,
        (.data[[abundance_col]] - .gene_min) / (.gene_max - .gene_min),
        1
      )
    ) %>%
    select(-.gene_max, -.gene_min) %>%
    ungroup()
}

## ---- 2e. Controlling how tightly rows/columns are packed ---------------------
##
## Within ggplot, panel_spacing (a plot_ppi_dotplot() argument) is the only
## real lever for the gap BETWEEN category panels. The gap between
## individual gene rows inside one panel isn't a separate setting -- it
## falls out of the panel's physical height vs. how many genes are in it,
## which is set by ggsave()'s height, not by anything in the ggplot object.
## This helper picks a height proportional to your gene count, so the
## vertical distance between rows stays constant regardless of how many
## genes you plot. row_height_in IS that per-gene distance, in inches --
## raise it for more breathing room, lower it to pack rows tighter. Plots
## with more category panels need a bit more base_in, since each panel's
## strip label and the spacing around it eat into the fixed overhead.
suggested_plot_height <- function(data, row_height_in, base_in = 1) {
  n_rows <- if ("category" %in% names(data)) {
    nrow(distinct(data, category, gene))   # a multi-category gene counts once per panel
  } else {
    length(unique(data$gene))
  }
  base_in + n_rows * row_height_in
}

## Companion function for width -- used ONCE in Section 3 to compute a
## single plot_width_in value shared by every category's plot (see the
## note there on why width is fixed rather than recomputed per category).
suggested_plot_width <- function(data, col_width_in, base_in = 1.5, Bait_order = NULL) {
  n_baits <- if (!is.null(Bait_order)) {
    length(Bait_order)
  } else {
    length(unique(data$Bait))
  }
  base_in + n_baits * col_width_in
}

## ---- 2f-bis. Fixed-width y-axis labels ---------------------------------------
##
## The actual fix for column-width consistency: rather than letting each
## plot's y-axis label gutter auto-size to THAT plot's own longest gene
## name (which is what varies width between categories), measure the
## longest name across the WHOLE dataset once, and force every plot's
## label gutter to that exact width regardless of what it individually
## needs. Legends being off (see plot_ppi_dotplot()'s legend.position =
## "none") removes the other major source of width variability.
##
## max_text_width_in(): renders each string via a throwaway PDF device
## (opened and closed immediately, nothing written to disk -- this is
## the standard base R trick for querying font metrics without an
## active plot window) and returns the widest one, in inches. Falls
## back to fallback_in with a message if text measurement fails for any
## reason (e.g. no working graphics device on this machine) rather than
## erroring out the whole script.
##
## Note on font accuracy: R's pdf() device doesn't have native Arial
## metrics (it isn't one of the 14 standard PDF core fonts) and
## substitutes Helvetica-equivalent widths instead. Arial was designed
## to be metrically compatible with Helvetica, so this is normally a
## very close approximation, not an exact one.
max_text_width_in <- function(texts, family = "Arial", face = "bold", size_pt = 12,
                              fallback_in = 1.5) {
  tryCatch({
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    widths <- vapply(unique(texts), function(txt) {
      grid::convertWidth(
        grid::grobWidth(grid::textGrob(
          txt, gp = grid::gpar(fontfamily = family, fontface = face, fontsize = size_pt))),
        "in", valueOnly = TRUE)
    }, numeric(1))
    max(widths, na.rm = TRUE)
  }, error = function(e) {
    message("max_text_width_in: couldn't measure text width (", conditionMessage(e),
            "); falling back to ", fallback_in, "in.")
    fallback_in
  })
}

## ---- 2f-ter. Category-coloured BFDR gradient -----------------------------
##
## Builds one category's 3-bin BFDR_colours dictionary from a single base
## colour (see category_colours in Section 1). The colour ITSELF now
## represents the strongest bin (<= 0.01), darkened slightly first for
## headroom -- some dictionary colours are already fairly light, which
## would leave the two lighter tiers below them barely distinguishable
## from each other or from white. The weaker two bins then fade toward
## WHITE (not black) by increasing amounts, so the gradient reads as
## "how much colour is left" = "how confident is this hit" at a glance.
##
## darken_amount : how much to darken base_colour before using it as the
##                 <= 0.01 anchor (0 = unchanged, 1 = pure black). 0.2
##                 default is deliberately modest -- enough headroom for
##                 the two lighter tiers without making already-dark
##                 dictionary colours (like #332288) look near-black.
## mid_lighten   : how much lighter the < 0.05 bin is than the (already
##                 darkened) anchor colour, toward white. 0.5 = "50%
##                 lighter", as specified.
## weak_lighten  : how much lighter the = 0.05 bin is, same basis.
##                 Deliberately further than mid_lighten so the palette
##                 has 3 visually distinct steps, not 2.
darken_colour <- function(colour, amount) {
  grDevices::colorRampPalette(c(colour, "black"))(101)[round(amount * 100) + 1]
}
lighten_colour <- function(colour, amount) {
  grDevices::colorRampPalette(c(colour, "white"))(101)[round(amount * 100) + 1]
}

category_BFDR_palette <- function(base_colour, darken_amount = 0.2,
                                  mid_lighten = 0.5, weak_lighten = 0.8) {
  anchor <- darken_colour(base_colour, darken_amount)
  c("<= 0.01" = anchor,
    "< 0.05"  = lighten_colour(anchor, mid_lighten),
    "= 0.05"  = lighten_colour(anchor, weak_lighten))
}

## ---- 2g. Combined "category: bin" colour key (for combined-mode plots) ------
##
## A single ggplot scale_colour_manual() can only map ONE name (like
## "<= 0.01") to ONE colour at a time -- it can't show a different
## gradient per facet panel on its own. The workaround: build a colour
## key on the combined string "Category: bin" instead of just "bin",
## with one entry per (category, bin) combination that actually appears
## in the data, each looked up from THAT category's own
## category_BFDR_palette(). plot_ppi_dotplot() then colours points by
## this combined key when colour_by_category = TRUE (see Section 2f),
## instead of by bin alone. The legend would be enormous with this many
## levels, but legends are already off (legend.position = "none"), so
## that's never actually shown anywhere.
##
## categories_present : character vector of category names that will
##                       appear in the plot (e.g. unique(plot_data$category)).
## category_colours, default_category_colour : same dictionaries as
##                       Section 1 -- a category not in category_colours
##                       falls back to default_category_colour, same as
##                       everywhere else.
build_combined_BFDR_colours <- function(categories_present, category_colours,
                                        default_category_colour) {
  out <- character(0)
  for (cat in categories_present) {
    base <- if (cat %in% names(category_colours)) category_colours[[cat]] else default_category_colour
    pal  <- category_BFDR_palette(base)
    names(pal) <- paste(cat, names(pal), sep = ": ")
    out <- c(out, pal)
  }
  out
}

## Forces a ggplot's y-axis label column to an exact physical width,
## overriding ggplot's automatic "size the gutter to this plot's own
## label content" layout. Call this on a finished plot right before
## ggsave() -- NOT earlier, since converting to a gtable here loses the
## normal $data/$layers access other code (e.g. suggested_plot_height())
## relies on. Returns a gtable (a grid grob), which ggsave() accepts
## directly, same as a ggplot object.
## Height of a single line of text at the given font spec, in inches --
## used for the ROTATED y-axis title: since it's rotated 90 degrees, its
## WIDTH footprint in the plot is just this line-height, constant
## regardless of how long the category name text itself is (unlike the
## y-axis gene labels, which run horizontal and do get wider with
## longer text). Same throwaway-PDF-device technique as
## max_text_width_in(); same Arial-via-Helvetica accuracy caveat.
measure_text_height_in <- function(family = "Arial", face = "bold", size_pt = 12,
                                   fallback_in = 0.25) {
  tryCatch({
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::convertHeight(
      grid::grobHeight(grid::textGrob("Ag", gp = grid::gpar(
        fontfamily = family, fontface = face, fontsize = size_pt))),
      "in", valueOnly = TRUE)
  }, error = function(e) {
    message("measure_text_height_in: couldn't measure text height (", conditionMessage(e),
            "); falling back to ", fallback_in, "in.")
    fallback_in
  })
}

## title_width_in (optional): also forces the rotated y-axis TITLE
## column ("ylab-l" in the gtable) to an exact width -- without this,
## ggplot still reserves its own automatic amount of space for that
## column, which silently eats into whatever's left for the panel once
## the axis-l column above is force-fixed. Usually invisible with many
## Bait columns (a small tax spread across a wide panel); becomes
## severe with very few Baits (the same small amount eaten from an
## already-narrow panel), which is what caused visible dot clipping in
## single-Bait plots. Pass measure_text_height_in()'s result here (see
## Section 3) -- do NOT skip this for low-Bait-count plots.
##
## strip_width_in (optional): also forces the category-name STRIP
## column ("strip-l" in the gtable, from facet_grid(switch = "y") --
## only present at all when length(categories) > 1, i.e. combined
## mode) to an exact width. Same failure mode as title_width_in if
## skipped, just far more severe: ggplot sizes this column
## automatically to fit the LONGEST category name with no budget
## awareness at all, which can push the panel (and everything drawn in
## it) past the edge of a canvas sized without accounting for it --
## exactly what caused dots to render off-canvas entirely in combined
## multi-category plots. Pass max_text_width_in(panel_categories, ...)
## here (see Section 3) whenever plot_separately = FALSE.
fix_axis_width <- function(p, width_in, title_width_in = NULL, strip_width_in = NULL) {
  g <- ggplotGrob(p)
  left_idx <- unique(g$layout$l[grepl("^axis-l", g$layout$name)])
  g$widths[left_idx] <- unit(width_in, "in")
  if (!is.null(title_width_in)) {
    title_idx <- unique(g$layout$l[grepl("^ylab-l", g$layout$name)])
    if (length(title_idx) > 0) g$widths[title_idx] <- unit(title_width_in, "in")
  }
  if (!is.null(strip_width_in)) {
    strip_idx <- unique(g$layout$l[grepl("^strip-l", g$layout$name)])
    if (length(strip_idx) > 0) g$widths[strip_idx] <- unit(strip_width_in, "in")
  }
  g
}

## Draws a rounded-corner border around each panel, replacing the plain
## rectangular theme(panel.border = ...) (see plot_ppi_dotplot()'s theme,
## which sets panel.border = element_blank() so the two don't double up).
## ggplot2's own border elements have no corner-radius option, so this
## overlays a grid::roundrectGrob() directly onto the gtable instead --
## same "modify the gtable after the fact" technique as fix_axis_width()
## above. Call this AFTER fix_axis_width() (order doesn't matter for
## correctness, but chaining fix_axis_width() first keeps both
## gtable-editing steps together conceptually). Works with any number of
## facet panels (one border per panel-* cell in the gtable), so it's safe
## to use in both plot_separately modes.
##
## radius: corner radius, an absolute grid::unit() -- NOT relative to
##         panel size, so it looks the same regardless of how wide/short
##         a given panel is.
## colour: EITHER a single colour (recycled for every panel -- the usual
##         case in per-category mode), OR a vector with one colour PER
##         PANEL, in the same top-to-bottom order as the category
##         factor's levels (see Section 3's combined-mode border colour
##         vector) -- lets each category's panel get its own accent
##         colour even within one combined image.
## linewidth: in the same mm-based units ggplot2 theme elements use;
##         converted to grid's own lwd units (1/96 inch) internally.
add_rounded_panel_border <- function(g, radius = unit(8, "pt"), colour = "black", linewidth = 0.6) {
  panel_idx <- which(grepl("^panel", g$layout$name))
  colour <- rep_len(colour, length(panel_idx))
  for (i in seq_along(panel_idx)) {
    li <- g$layout[panel_idx[i], ]
    g <- gtable::gtable_add_grob(
      g,
      grid::roundrectGrob(
        r  = radius,
        gp = grid::gpar(col = colour[i], fill = NA, lwd = linewidth * 96 / 25.4)
      ),
      t = li$t, l = li$l, b = li$b, r = li$r,
      z = Inf, name = paste0("round-panel-border-", i)
    )
  }
  g
}

## ---- 2f. Plotting function ---------------------------------------------------
##
## data          : long-format data frame (gene, Bait, category, AvgSpec,
##                 rel_spc, BFDR).
## categories    : character vector of category values to include, in the
##                 order you want them to appear top-to-bottom. Defaults to
##                 all categories in the order they appear in the data.
## gene_order    : optional named list of gene orders per category
##                 (list(category_name = c("geneA","geneB",...))). Takes
##                 priority over order_genes/cluster_by for any category
##                 it names -- categories not named fall back to those.
## order_genes   : "cluster" (default) orders genes within each panel via
##                 cluster_gene_order() (Section 2b): multi-Bait genes are
##                 hierarchically clustered by hit pattern at the top,
##                 single-Bait genes are grouped by which Bait they hit
##                 and placed at the bottom. "total" instead just sorts
##                 all genes by total AvgSpec, ignoring this split.
## cluster_by    : "AvgSpec" (default) clusters on the actual spectral
##                 count values using canberra distance -- ProHits-viz's
##                 own default, suited to abundance data. "presence"
##                 clusters on which Baits each gene hits (binary
##                 distance) instead, ignoring magnitude entirely. Only
##                 used when order_genes = "cluster". Linkage defaults to
##                 ward.D2 (see cluster_gene_order() in Section 2b).
## Bait_order    : optional character vector for x-axis order. Also fixes
##                 which columns are DRAWN AT ALL: scale_x_discrete() below
##                 uses drop = FALSE so every Bait in Bait_order gets a
##                 column even if this particular category has no data for
##                 it (otherwise ggplot silently drops "unused" factor
##                 levels, e.g. for a category with data for only 1 Bait).
## size_range    : range (in mm-ish -- see scale_radius()'s own internal
##                 fontsize conversion) passed to scale_radius() for dot
##                 diameter. Diameter scales LINEARLY with rel_spc
##                 (0 -> smallest, 1 -> largest), as opposed to ggplot's
##                 default area-proportional scale_size(). Paired with an
##                 explicit limits = c(0, 1) below (see add_relative_spc()
##                 in Section 2d-bis for why that's essential -- without
##                 it, categories where every gene is single-hit render
##                 visibly smaller than categories with real variation,
##                 even at identical rel_spc = 1 and identical
##                 size_range). Fixed (not derived from the data), so the
##                 largest possible dot is the same size on every plot
##                 this function produces, regardless of category.
## panel_spacing : gap between category panels (a grid::unit()).
## fill_colours  : low/mid/high colours for the AvgSpec fill gradient
## fill_limits   : fixed (min, max) the fill gradient maps to, e.g.
##                 c(0, 100). Explicit limits, rather than ggplot's
##                 default of auto-scaling each plot to its own data's
##                 min/max, are what make the SAME AvgSpec value map to
##                 the SAME colour across every plot -- without this, a
##                 category whose highest value happens to be 40 would
##                 stretch its gradient across 0-40 instead of 0-100,
##                 making it look artificially more "intense" than an
##                 identical value in a category that actually reaches
##                 100. Values above the upper limit are squished to the
##                 top colour rather than becoming NA/blank (relevant if
##                 you ever change max_abundance in Section 1 without
##                 updating this to match, or skip cap_max_abundance()
##                 in your own pipeline).
## BFDR_colours  : named vector giving the BFDR stroke-colour bins, in
##                 order from most to least significant. Only Bait-prey
##                 pairs matching one of these bins are drawn -- pairs
##                 with BFDR outside all bins are simply not plotted
##                 (blank space in that panel position), even though the
##                 gene may still appear for other Baits it does pass.
##                 Add a "> 0.05" (or similar) entry here yourself if you
##                 ever want those drawn faintly instead of omitted. When
##                 colour_by_category = TRUE, this should instead be a
##                 combined-key dictionary from build_combined_BFDR_colours()
##                 (Section 2g) -- see Section 3 for how that's built.
## colour_by_category : FALSE (default) colours points by BFDR_bin alone,
##                 same behaviour as always. TRUE colours by a combined
##                 "category: bin" key instead (built internally from the
##                 data's own category column), which is what lets a
##                 SINGLE combined multi-category plot show each
##                 category's own BFDR gradient rather than one shared
##                 palette. BFDR_colours must be a combined-key
##                 dictionary (matching keys) when this is TRUE.
## font_family   : applied to every text element (axis labels, legend,
##                 strip labels, the y-axis title) via the theme's base
##                 "text" element, which the rest inherit from unless they
##                 already set their own family/face.
## title         : NOT a plot title (there isn't one -- a long category
##                 name wrapping across multiple lines there was eating
##                 into the fixed image height and squeezing the gene
##                 rows). Used as the y-axis title instead, rotated
##                 vertically, which takes up about the same width
##                 regardless of how long the text is.
plot_ppi_dotplot <- function(data,
                             categories    = NULL,
                             gene_order    = NULL,
                             order_genes   = c("cluster", "total"),
                             cluster_by    = c("AvgSpec", "presence"),
                             Bait_order    = NULL,
                             size_range    = c(1, 8),
                             panel_spacing = unit(1, "pt"),
                             fill_colours  = c("#ebebeb", "#5c5c5c", "#292929"),
                             fill_limits   = c(0, 100),
                             BFDR_colours  = c("<= 0.01" = "black",
                                               "< 0.05"  = "blue",
                                               "= 0.05"  = "#809fff"),
                             colour_by_category = FALSE,
                             font_family   = "Arial",
                             title         = NULL) {
  
  order_genes <- match.arg(order_genes)
  cluster_by  <- match.arg(cluster_by)
  
  if (is.null(categories)) categories <- unique(data$category)
  df <- data %>% filter(category %in% categories) %>%
    mutate(category = factor(category, levels = categories))
  
  if (!is.null(Bait_order)) {
    df <- df %>% mutate(Bait = factor(Bait, levels = Bait_order))
  } else {
    df <- df %>% mutate(Bait = factor(Bait, levels = unique(data$Bait)))
  }
  
  ## order genes within each category panel
  df <- df %>%
    group_by(category) %>%
    group_modify(function(sub, key) {
      cat_name <- as.character(key$category)
      if (!is.null(gene_order) && !is.null(gene_order[[cat_name]])) {
        ord <- gene_order[[cat_name]]
      } else if (order_genes == "cluster") {
        ord <- cluster_gene_order(sub, cluster_by = cluster_by)
      } else {
        ord <- sub %>%
          group_by(gene) %>%
          summarise(tot = sum(AvgSpec, na.rm = TRUE), .groups = "drop") %>%
          arrange(desc(tot)) %>%
          pull(gene)
      }
      ord <- as.character(unlist(ord, use.names = FALSE))
      if (!is.character(ord) || length(ord) == 0) {
        stop(sprintf(
          "plot_ppi_dotplot: gene ordering for category '%s' did not come back as a plain list of gene names (got class %s, length %d). This usually means duplicate/ambiguous rows somewhere upstream in the data for this category.",
          cat_name, paste(class(ord), collapse = "/"), length(ord)))
      }
      sub %>% mutate(gene = factor(gene, levels = rev(unique(ord))))
    }) %>%
    ungroup()
  
  ## bin BFDR into the stroke-colour categories: <= 0.01, then strictly
  ## < 0.05, then whatever's left up to and including 0.05 (== 0.05 in
  ## effect, but written as the sequential <= 0.05 fallback rather than
  ## a literal == 0.05 comparison -- floating-point equality on a
  ## computed decimal like 0.05 is unreliable; this gets the same
  ## practical result safely).
  df <- df %>%
    mutate(BFDR_bin = case_when(
      BFDR <= 0.01 ~ "<= 0.01",
      BFDR < 0.05  ~ "< 0.05",
      BFDR <= 0.05 ~ "= 0.05",
      TRUE ~ NA_character_
    ))
  
  ## When colouring by category, the actual colour key becomes
  ## "Category: bin" instead of just "bin" -- this is what lets each
  ## category panel show its own gradient from ONE shared
  ## scale_colour_manual() (see build_combined_BFDR_colours(),
  ## Section 2g). BFDR_colours is expected to already be that combined
  ## dictionary in this case.
  if (colour_by_category) {
    df <- df %>% mutate(BFDR_bin = paste(as.character(category), BFDR_bin, sep = ": "))
  } else if ("> 0.05" %in% names(BFDR_colours)) {
    df <- df %>% mutate(BFDR_bin = if_else(is.na(BFDR_bin), "> 0.05", BFDR_bin))
  }
  
  n_dropped <- sum(is.na(df$BFDR_bin) | !(df$BFDR_bin %in% names(BFDR_colours)))
  if (n_dropped > 0) {
    message(sprintf(
      "plot_ppi_dotplot: dropped %d point(s) with no defined BFDR colour bin (missing BFDR, > 0.05 with no matching entry, or -- if colour_by_category -- a category not covered by BFDR_colours)",
      n_dropped))
  }
  df <- df %>%
    filter(!is.na(BFDR_bin), BFDR_bin %in% names(BFDR_colours)) %>%
    mutate(BFDR_bin = factor(BFDR_bin, levels = names(BFDR_colours)))
  
  p <- ggplot(df, aes(x = Bait, y = gene)) +
    geom_point(aes(size = rel_spc, fill = AvgSpec, colour = BFDR_bin),
               shape = 21, stroke = 1.1) +
    scale_radius(range = size_range, limits = c(0, 1),
                 name = "Relative avg.\nspectral count") +
    scale_fill_gradientn(colours = fill_colours, name = "Avg. spectral\ncount",
                         limits = fill_limits, oob = scales::squish) +
    scale_colour_manual(values = BFDR_colours, name = "BFDR", drop = FALSE) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_discrete(position = "top", drop = FALSE) +
    labs(title = NULL, x = NULL, y = title) +
    theme_minimal(base_size = 12, base_family = font_family) +
    theme(
      text               = element_text(family = font_family, face = "bold"),
      legend.position    = "none",
      panel.border       = element_blank(),  # replaced by add_rounded_panel_border() in Section 3 --
      # ggplot2 has no corner-radius option for theme borders
      panel.spacing.y    = panel_spacing,
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.placement    = "outside",
      strip.text.y.left  = if (length(categories) == 1) element_blank()
      else element_text(angle = 0, face = "bold", hjust = 1),
      strip.background   = element_blank(),
      axis.text.x.top    = element_text(angle = 90, hjust = 0, vjust = 0.5, face = "bold"),
      axis.text.y        = element_text(face = "bold"),
      axis.title.y       = element_text(face = "bold", angle = 90, hjust = 0.5,
                                        margin = margin(r = 0)),  # ggplot's default theme
      # reserves extra space between the title and the axis
      # labels next to it; zeroing this pulls the title flush
      # against the label column (which is itself forced to
      # an exact width via fix_axis_width() in Section 3)
      plot.background    = element_blank()  # removes the full-canvas white background rect
      # ggplot2 includes by default -- it spans the
      # whole image (margins included), not just the
      # visible plot, which is why deleting it in a
      # vector editor changes nothing you can see
    )
  
  p
}


## ============================================================================
## ---- 3. Run -- loads your data, builds the plot(s), saves them -----------
## ============================================================================

raw_data <- read.table(data_file, stringsAsFactors = FALSE, sep = "\t", header = TRUE) %>%
  rename(gene = PublicName)

## AvgSpec and BFDR must be numeric -- if your file has any non-numeric
## entries in these columns (stray text, "N/A", a stray space, a
## thousands separator, etc.), R silently reads the WHOLE column as
## character/text, and every downstream numeric operation on it (dist(),
## max(), -x) either coerces it back with NAs sprinkled in or fails
## outright with a confusing error far from the actual cause. Catch it
## here instead, with the actual offending values named.
check_numeric_column <- function(data, col) {
  raw_vals <- data[[col]]
  numeric_vals <- suppressWarnings(as.numeric(raw_vals))
  bad <- unique(raw_vals[is.na(numeric_vals) & !is.na(raw_vals) & raw_vals != ""])
  if (length(bad) > 0) {
    stop(sprintf(
      "Column '%s' in %s has %d non-numeric value(s) that couldn't be converted to a number, e.g.: %s\nFix these in your data file and re-run.",
      col, data_file, length(bad), paste(head(bad, 10), collapse = ", ")))
  }
  numeric_vals
}
raw_data$AvgSpec <- check_numeric_column(raw_data, "AvgSpec")
raw_data$BFDR    <- check_numeric_column(raw_data, "BFDR")

## Each (gene, Bait) pair should appear exactly once in your data --
## duplicates (e.g. un-merged replicate rows) would silently distort
## clustering and dot sizing if left in, so this stops here with the
## offending rows named rather than guessing how to collapse them.
dup_pairs <- raw_data %>% count(gene, Bait, name = "n_rows") %>% filter(n_rows > 1)
if (nrow(dup_pairs) > 0) {
  stop(sprintf(
    "Found %d duplicate (gene, Bait) pair(s) in %s -- each pair should appear exactly once.\n%s%s",
    nrow(dup_pairs), data_file,
    paste(sprintf("  %s / %s (%d rows)",
                  head(dup_pairs$gene, 15), head(dup_pairs$Bait, 15), head(dup_pairs$n_rows, 15)),
          collapse = "\n"),
    if (nrow(dup_pairs) > 15) sprintf("\n  ...and %d more", nrow(dup_pairs) - 15) else ""
  ))
}

if (is.null(category_file)) {
  ## Your data file already has its own category column -- use it as-is.
  plot_data <- raw_data
} else {
  ## Category membership comes from the text file: it decides which genes
  ## are plotted and which category/categories each belongs to. It does
  ## NOT set row order -- gene_order is deliberately left NULL below so
  ## clustering (order_genes = "cluster") determines the y-axis order.
  gene_categories <- read_gene_categories(category_file)
  plot_data <- raw_data %>%
    select(-any_of("category")) %>%
    inner_join(gene_categories, by = "gene", relationship = "many-to-many")
}
panel_categories <- unique(plot_data$category)

if (!is.null(Bait_order)) {
  ## Restrict to baits that will actually appear in the plot, BEFORE any
  ## filtering or the relative-spec calculation below -- otherwise a bait
  ## sitting in the raw data but excluded from Bait_order could still
  ## influence a gene's min/max range.
  plot_data <- plot_data %>% filter(Bait %in% Bait_order)
}

plot_data <- plot_data %>%
  apply_primary_filter(threshold = BFDR_threshold) %>%   # gene-level inclusion -- uses the FULL row set,
  # unaffected by the row-level filter below
  apply_min_abundance(min_abundance = min_abundance) %>%
  ## Drop individual bait-gene pairs whose BFDR is above the highest
  ## visible bin (0.05, matching BFDR_colours' default bins inside
  ## plot_ppi_dotplot() -- update this if you ever change those bin
  ## edges). These rows never get drawn as a dot anyway (see the BFDR
  ## binning step inside plot_ppi_dotplot()), but if left in here they'd
  ## still count toward a gene's own min/max in add_relative_spc() below
  ## -- which is exactly what was making single-Bait-hit genes come out
  ## smaller than the max diameter: a hidden, non-significant background
  ## reading at another Bait could still numerically outrank the one dot
  ## that actually gets drawn.
  filter(BFDR <= 0.05) %>%
  add_relative_spc()   # per-gene relative spectral count, computed from the
# TRUE (uncapped) AvgSpec values that survived
# filtering above -- max_abundance no longer
# touches the actual data (see the note above
# cap_max_abundance() in Section 2), only where
# the fill colour scale visually tops out.

## y_label_width_in is measured ONCE, from the longest gene name across
## the ENTIRE dataset (every gene that could appear in ANY output plot),
## at the same 12pt bold Arial everything else uses. Every plot's y-axis
## label column gets forced to this exact width via fix_axis_width()
## right before ggsave() below -- see the note on max_text_width_in() in
## Section 2 for how the measurement itself works and its accuracy
## caveat.
##
## plot_width_in is now computed directly from known quantities (label
## width + title width + panel width + a small fixed edge margin)
## rather than the old generic suggested_plot_width() guess -- that
## guess existed mainly to account for legend width, which no longer
## applies now that legends are off (legend.position = "none" in
## plot_ppi_dotplot()).
y_label_width_in <- max_text_width_in(plot_data$gene, family = "Arial", face = "bold",
                                      size_pt = 12) + 0.15   # + a little breathing room

## The rotated y-axis title's width footprint (see measure_text_height_in()
## in Section 2) -- constant regardless of category name length, since a
## rotated single line's width is just its line-height. Computed once,
## not per-category, for exactly that reason. Essential for narrow
## (especially single-Bait) plots -- see fix_axis_width()'s note on why
## leaving this unaccounted for silently shrinks the panel.
y_title_width_in <- measure_text_height_in(family = "Arial", face = "bold", size_pt = 12) + 0.05

## The category-name STRIP column's width -- ONLY present at all when
## plot_separately = FALSE (combined mode shows it; per-category mode
## sets strip.text.y.left = element_blank() since there's only ever one
## category per plot there). Measured the same way as y_label_width_in,
## just over category names instead of gene names. Leaving this
## unbudgeted is what pushed panels off-canvas entirely in combined
## mode -- see fix_axis_width()'s note on strip_width_in.
y_strip_width_in <- if (!plot_separately && length(panel_categories) > 1) {
  max_text_width_in(panel_categories, family = "Arial", face = "bold", size_pt = 12) + 0.15
} else {
  0
}

plot_width_in <- y_label_width_in + y_title_width_in + y_strip_width_in +
  length(Bait_order) * col_width_in + 0.2  # + small right-edge margin

colour_tag <- if (use_custom_BFDR_colours) "custom" else "grey"

if (plot_separately) {
  
  ## One plot per category, each sized to just its own gene count.
  for (cat in panel_categories) {
    cat_data <- plot_data %>% filter(category == cat)
    if (nrow(cat_data) == 0) {
      message(sprintf("Skipping '%s': no data left after filtering.", cat))
      next
    }
    
    ## Resolved per-category: each category's own base colour (or
    ## default_category_colour if it's not in the dictionary) drives its
    ## own black-to-colour BFDR gradient AND the rounded panel border,
    ## so the border matches that category's accent colour. When custom
    ## colouring is off, both stay at their plain defaults.
    if (use_custom_BFDR_colours) {
      base_colour          <- if (cat %in% names(category_colours)) category_colours[[cat]] else default_category_colour
      BFDR_colours_for_cat <- category_BFDR_palette(base_colour)
      border_colour        <- base_colour
    } else {
      BFDR_colours_for_cat <- default_BFDR_colours
      border_colour        <- "black"
    }
    
    p <- plot_ppi_dotplot(
      cat_data,
      categories   = cat,
      gene_order   = NULL,          # clustering decides row order
      Bait_order   = Bait_order,
      fill_limits  = c(min_abundance, max_abundance),
      BFDR_colours = BFDR_colours_for_cat,
      title        = cat
    )
    
    h <- suggested_plot_height(p$data, row_height_in = row_height_in)  # compute BEFORE fix_axis_width()
    
    safe_name <- gsub("[^A-Za-z0-9_-]+", "_", cat)
    out_path  <- paste0(output_prefix, "_", safe_name, "_", colour_tag, ".svg")
    
    p_grob <- fix_axis_width(p, width_in = y_label_width_in, title_width_in = y_title_width_in, strip_width_in = y_strip_width_in)
    if (use_rounded_borders) p_grob <- add_rounded_panel_border(p_grob, colour = border_colour)
    
    ggsave(out_path, p_grob, width = plot_width_in, height = h, bg = "transparent")
    cat("Wrote", out_path, "\n")
  }
  
} else {
  
  ## One combined plot, every category as its own facet panel. With
  ## use_custom_BFDR_colours = TRUE, this now gets genuine per-category
  ## colouring via the combined "category: bin" key (see
  ## build_combined_BFDR_colours(), Section 2g, and colour_by_category
  ## in plot_ppi_dotplot(), Section 2f) -- each panel's gradient AND
  ## border colour matches that category's own dictionary entry, same
  ## as in per-category mode, just all in one file.
  if (use_custom_BFDR_colours) {
    BFDR_colours_for_plot <- build_combined_BFDR_colours(panel_categories, category_colours, default_category_colour)
    border_colours_for_plot <- vapply(panel_categories, function(cat) {
      if (cat %in% names(category_colours)) category_colours[[cat]] else default_category_colour
    }, character(1))
  } else {
    BFDR_colours_for_plot <- default_BFDR_colours
    border_colours_for_plot <- "black"
  }
  
  p <- plot_ppi_dotplot(
    plot_data,
    categories          = panel_categories,
    gene_order          = NULL,            # clustering decides row order within each panel
    Bait_order          = Bait_order,
    fill_limits         = c(min_abundance, max_abundance),
    BFDR_colours        = BFDR_colours_for_plot,
    colour_by_category  = use_custom_BFDR_colours,
    title               = plot_title
  )
  
  h <- suggested_plot_height(p$data, row_height_in = row_height_in)  # compute BEFORE fix_axis_width()
  
  out_path <- paste0(output_prefix, "_", colour_tag, ".svg")
  
  p_grob <- fix_axis_width(p, width_in = y_label_width_in, title_width_in = y_title_width_in, strip_width_in = y_strip_width_in)
  ## border_colours_for_plot is in panel_categories order, which matches
  ## the facet_grid panel order since categories was passed as that same
  ## vector above (factor levels set in that order).
  if (use_rounded_borders) p_grob <- add_rounded_panel_border(p_grob, colour = border_colours_for_plot)
  
  ggsave(out_path, p_grob, width = plot_width_in, height = h, bg = "transparent")
  cat("Done. Wrote", out_path, "\n")
}
