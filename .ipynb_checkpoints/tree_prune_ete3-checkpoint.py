import sys
import os
from ete3 import Tree

tree_file = sys.argv[0]
species_list = sys.argv[1]

if os.path.exists(tree_file):
    tree_name = os.path.basename(tree_file)
    tree_name = os.path.splitext(tree_name)
else:
    sys.exit("Error: " + tree_file + " not found.")

t = Tree(tree_file)


# t.write(outfile = tree_name[0] + "_pruned" + tree_name[1])