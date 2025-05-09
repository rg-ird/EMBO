#install packegae if necessary
install.packages("ape")
install.packages("phangorn")
install.packages("ggtree")
install.packages("dplyr")
install.packages("BiocManager")
BiocManager::install("ggtree")

################################################################################
#Step 1 Get node number
################################################################################
library(ggtree)
library(phangorn)
library(ape)
library(dplyr)
# Newick tree
tree <- read.tree("mytree.tree")
#
rooted_tree <- midpoint(tree)

# see circular treewith node number
p <- ggtree(rooted_tree, layout = "circular") +
  geom_text(aes(label = node), color = "red", size = 0.5, hjust = -0.3) +
  geom_tree(color = "gray", size = 0.05, linewidth = 0.05)


#print(p)

# Register the picture in PDF
ggsave("mytree.NODE_NUMBER.pdf", plot = p, width = 9, height = 11, dpi = 1200)
