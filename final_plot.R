setwd("/home/isaac/Data/TURBO_DATA/PHYLO/summary")
pacman::p_load(dplyr, ggplot2, ape, ggtree, tidytree, ggpubr, reshape)


files = list.files(pattern = "best_hmm_summary.tsv")
names(files) = sub(".bes.*$", "", files)
code = read.delim("sp_code.txt")

lst_hts = list()
for (i in names(files)) {
 tmp.tab = read.delim(files[i])
 hts = sub("_.*", "",tmp.tab$target) %>% table
 ms = setdiff(code$id, names(hts))
 vls = rep(0, length(ms))
 names(vls) = ms
 hts = c(hts, vls)
 snames = code[match(names(hts), code$id),]$species
 names(hts) = snames
 hts = hts[code$species]
 lst_hts[[i]] = hts
}

mtx = do.call(rbind, lst_hts) |> t() |> as.data.frame()

tree = read.tree("worm_tree_sorted.txt")
mtx = mtx[tree$tip.label, ]

tree2 = tree

tree2$tip.label = strsplit(tree2$tip.label, "_") %>%
  sapply(function(x) paste0(x[c(1,2)], collapse = "_"))

tree2$tip.label = str_to_title(tree2$tip.label) %>%
  sub("(.).*_", paste0("\\1", ". "), .)

tree2$tip.label = paste0("italic('", tree2$tip.label, "')")

tree2$tip.label[tree2$tip.label == "italic('C. elegans')"] =
  "bolditalic('C. elegans')~'*'"

p1 = ggtree(tree2, ladderize = F) +
  geom_tiplab(
    parse = TRUE,
    align = TRUE,
    family = "Arial",
    size = 5
  ) +
  xlim_tree(1.6)

df = melt(as.matrix(mtx))
colnames(df)= c("species", "gene", "count")

df$count= as.numeric(df$count)
df$count_cat= factor(df$count, levels = 0:max(df$count))

colors= colorRampPalette(c(
  "#FFFFFF",
  "#D6EAF8",
  "#AED6F1",
  "#5DADE2",
  "#21618C",
  "#154360"
))(length(levels(df$count_cat)))

names(colors)= levels(df$count_cat)

p2 = ggplot(df, aes(x = gene, y = species, fill = count_cat)) +
  geom_tile() +
  scale_fill_manual(
    values = colors,
    name = "HMM hits",
    drop = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      face = "italic",
      size = 14
    ),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right"
  ) +
  theme(axis.text.y = element_blank())

p3 = ggarrange(p1, p2, align = "hv", widths = c(1,1))

ggsave("unchar_consrevation.pdf", device = cairo_pdf, width = 8, height = 12,
       dpi = 600, path = "/home/isaac/Data/TURBO_DATA/PHYLO/")


write.table(mtx, "worm_unch_counts.txt", sep = "\t", quote = F, col.names = T, row.names = T)
