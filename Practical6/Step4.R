
library(ape)
library(ggtree)
library(dplyr)

################################################################################
4. count domains
################################################################################

# Newick
tree <- read.tree("mytree.tree")
rooted_tree <- midpoint(tree)

#
get_tips <- function(tree, node) {
  tips <- c()
  traverse_tree <- function(tree, node) {
    children <- tree$edge[tree$edge[,1] == node, 2]
    for (child in children) {
      if (child <= Ntip(tree)) {
        tips <<- c(tips, tree$tip.label[child])
      } else {
        traverse_tree(tree, child)
      }
    }
  }
  traverse_tree(tree, node)
  return(tips)
}

# get tips under node XX
tips_in_node <- get_tips(rooted_tree, xx)

# count domains depending of prefixes
observed_counts <- c(
  A = sum(grepl("^A", tips_in_node)),
  B = sum(grepl("^B", tips_in_node)),
)

# Results
print(observed_counts)

# Extract names
a_domains <- tips_in_node[grepl("^A", tips_in_node)]
b_domains <- tips_in_node[grepl("^B", tips_in_node)]

# save
writeLines(a_domains, "A_domains.txt")
writeLines(b_domains, "B_domains.txt")