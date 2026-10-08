if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# This script aims to build PAE plots for AF3 models for Sebastian's paper
# The idea is to read the json files for each model and build an individual plots
# To produce the whole set of plots I'll need a for, perhaps in bash
# Load packages
pacman::p_load(jsonlite, ggplot2, scales, optparse)
# Optpasre: This will be relevant to be able to reuse the code
opts = list(
  make_option(c("-i", "--input"), type = "character", help = "JSON from AlphaFold", metavar = "file.json"),
  make_option(c("-m", "--main"), type = "character", help = "Plot title", metavar = "PAE Plot with Chain Boundaries"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Out PNG file", metavar = "plot.png")
)

parser = OptionParser(option_list = opts)
args   = parse_args(parser)
# I'll build a function to read the json file with jsonlite
# Then the plan is to use ggplot to build the PAE plots
plot_pae = function(jfile, titulo, outFile) {
  jdat = fromJSON(jfile, simplifyVector = FALSE)
  mat  = jdat$pae
  n    = length(mat)
  if (!all(sapply(mat, length) == n)) stop("Matriz PAE no es cuadrada")
  pdat = data.frame(
    x = rep(1:n, each = n),
    y = rep(1:n, times = n),
    pae = as.numeric(unlist(mat))
  )
  # Get the chains from the json based on chain IDs
  cids   = unlist(jdat$token_chain_ids)
  bounds = which(diff(as.numeric(as.factor(cids))) != 0)
  clabs  = unique(cids)
  mids   = c(0, bounds, n)  # para etiquetas
  ctrs   = head(mids, -1) + diff(mids) / 2
  # And prepare the plot
  # En mi primenr version el pae no se ve bien por que no hay division entre las chains
  # entonces voy a agregar linias que delimiten cada cadena
  gg = ggplot(pdat, aes(x = x, y = y, fill = pae)) +
    geom_raster() +
    scale_fill_viridis_c(option = "viridis", limits = c(0, 31.75), oob = squish) +
    coord_fixed() +
    labs(title = titulo,
         x = "Residue Index", y = "Residue Index", fill = "PAE (Å)") +
    theme_minimal() +
    geom_vline(xintercept = bounds, color = "white", linetype = "dashed") +
    geom_hline(yintercept = bounds, color = "white", linetype = "dashed") +
    scale_x_continuous(breaks = ctrs, labels = clabs) +
    scale_y_continuous(breaks = ctrs, labels = clabs) +
    theme(axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size = 20)) +
    theme(axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 20, angle = 90, hjust = .5)) +
    theme(plot.title = element_text(size = 20, hjust = 0.5))
  
  # Save the plot
  ggsave(outFile, gg, width = 8, height = 7, dpi = 300)
   
  
}
plot_pae(args$input, args$main, args$output)