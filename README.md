
## TreeSAK (A Swiss-Army-Knife for manipulating phylogenetic trees)

[![pypi licence ](https://img.shields.io/pypi/l/TreeSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/TreeSAK.svg)](https://pypi.python.org/pypi/TreeSAK) 

Contact
---

+ **Weizhi Song**, Postdoctoral Researcher
+ School of Life Sciences, The Chinese University of Hong Kong
+ E-mail: songwz03@gmail.com

    
Installation
---

1. TreeSAK has been tested on Linux/Mac, but **NOT** yet on Windows.

1. TreeSAK is implemented in python3, you can install it with pip3:

       # for the first time installation
       pip3 install TreeSAK
      
       # for updating
       pip3 install --upgrade TreeSAK

TreeSAK modules
---
   
                 ...::: TreeSAK v1.0.0 :::...

    Multiple Sequence Alignment (MSA)
       convert_align_format   ->  Convert alignment format
       OneLineAln             ->  One-line fasta format alignments
       SubsetAlnCols          ->  Subset MSA by column    

    Tree-related
       GTDB_tree_r207         ->  Infer GTDB (r207) archaeal/bacterial tree
       subset_tree            ->  Subset tree
       compare_trees          ->  Compare trees with Mantel test
       rename_leaves          ->  Rename tree leaves
       FLN                    ->  Format leaf names (e.g. remove spaces in names)
       iTOL                   ->  Prepare iTOL-compatible files for tree visualization
