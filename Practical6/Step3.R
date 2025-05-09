#install packegae if necessary
install.packages("ape")
install.packages("phangorn")
install.packages("ggtree")
install.packages("dplyr")
install.packages("BiocManager")
BiocManager::install("ggtree")

################################################################################
3. compare different genomes
################################################################################

library(ape)
library(ggtree)
library(ggrepel)

#  Newick tree
tree <- read.tree("mytree.tree")
rooted_tree <- midpoint(tree)

# convert in data frame
fortified_tree <- fortify(rooted_tree)

# 2 genomes (A, B)
legend_data <- data.frame(
  x = rep(0, 3),  # Position X fictive pour la légende
  y = rep(0, 3),  # Position Y fictive pour la légende
  category = c("A", "B")
  color = c("green", "pink"),
  shape = c(16, 16)
)

# circular tree
p <- ggtree(rooted_tree, layout = "circular") +
  geom_tree(linewidth = 0.05) +  # Lignes plus fines

  # A
  geom_point(aes(x = x + 0.2, y = y),
             shape = 16, size = 0.2, color = "green",
             data = subset(fortified_tree, grepl("^A", label))) +

  # B
  geom_point(aes(x = x + 0.4, y = y),
             shape = 16, size = 0.2, color = "pink",
             data = subset(fortified_tree, grepl("^B", label))) +


  # legende
  geom_point(aes(x = x, y = y, color = category, shape = category),
             size = 3, data = legend_data, show.legend = TRUE) +

  # légende
  scale_color_manual(values = setNames(legend_data$color, legend_data$category)) +
  scale_shape_manual(values = setNames(legend_data$shape, legend_data$category)) +
  labs(color = "Category", shape = "Category") +
  theme(legend.position = "right")

# Clade  labels
p <- p + geom_cladelabel(node = 28362, label = "Copia", offset = 1, color = "gray", barsize = 1, align = TRUE)+
  geom_cladelabel(node = xx, label = "Gypsy", offset = 1, color = "black", barsize = 1, align = TRUE)+
#  geom_cladelabel(node = xx, label = "Gypsy", offset = 1, color = "black", barsize = 1, align = TRUE)+
  geom_cladelabel(node = xx, label = "TAT", offset = 0.8, color = "green", barsize = 1, align = TRUE) +
  geom_cladelabel(node = xx, label = "ATHILA", offset = 0.8, color = "darkgreen", barsize = 1, align = TRUE) +
  geom_cladelabel(node = xx, label = "DEL", offset = 0.8, color = "red", barsize = 1, align = TRUE) +
  geom_cladelabel(node = xx, label = "GALADRIEL", offset = 0.8, color = "yellow", barsize = 1, align = TRUE) +
  geom_cladelabel(node = xx, label = "CRM", offset = 0.8, color = "orange", barsize = 1, align = TRUE) +
  geom_cladelabel(node = xx, label = "REINA", offset = 0.8, color = "lightgreen", barsize = 1, align = TRUE)+
  geom_cladelabel(node = xx, label = "SIRE", offset = 0.8, color = "purple", barsize = 1, align = TRUE) +
  geom_cladelabel(node = xx, label = "ORYCO", offset = 0.8, color = "pink", barsize = 1, align = TRUE) +
  geom_cladelabel(node = xx, label = "RETROFIT", offset = 0.8, color = "violet", barsize = 1, align = TRUE) +
  geom_cladelabel(node = xx, label = "TORK", offset = 0.8, color = "blue", barsize = 1, align = TRUE) +
  geom_cladelabel(node = 32255, label = "BIANCA", offset = 0.8, color = "lightblue", barsize = 1, align = TRUE)

# PDF
ggsave("my_tree_final.pdf", plot = p, device = "pdf", width = 10, height = 10, dpi = 1200)
