
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

                 ...::: TreeSAK v1.3.1 :::...

    Marker-related
       parse_deltall_stdout   ->  Parse stdout of deltaLL.rb
       get_arCOG_seq          ->  Retrieve arCOG sequences

    Multiple Sequence Alignment
       ConvertMSA             ->  Convert alignment format
       OneLineAln             ->  One-line fasta format alignments
       SubsetAlnCols          ->  Subset MSA by column    

    Tree-related
       subset_tree            ->  Subset tree
       compare_trees          ->  Compare trees with Mantel test
       rename_leaves          ->  Rename tree leaves
       FLN                    ->  Format leaf names (e.g. remove spaces in names)

    Dating-related
       AssessCVG              ->  Assess dating convergence
