#install packegae if necessary
install.packages("ape")
install.packages("phangorn")
install.packages("ggtree")
install.packages("dplyr")
install.packages("BiocManager")
BiocManager::install("ggtree")

################################################################################
#Step2 See  reference RT
################################################################################
install.packages("ggrepel")
library(ggtree)
library(ggrepel)
library(dplyr)
#  Newick tree
tree <- read.tree("Gmytree.tree")
rooted_tree <- midpoint(tree)
#
rooted_tree$node <- as.factor(rooted_tree$node)
rooted_tree$parent <- as.factor(rooted_tree$parent)

# circular tree@
p <- ggtree(rooted_tree, layout = "circular") +
  geom_tree(linewidth = 0.01) +
  #
  geom_tiplab(aes(color = case_when(
    grepl("RT_", label) ~ "red",
    TRUE ~ "black"
  ), label = label), size = 0.2, offset = 0.02, align = FALSE) +
  theme(legend.position = "none")  # Supprimer la légende

# Register picture in  PDF
ggsave("mytree.with_ref.pdf", plot = p, width = 9, height = 11, dpi = 1200)
