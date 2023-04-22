
## TreeSAK (A Swiss-Army-Knife for manipulating phylogenetic tree related files)

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

                 ...::: TreeSAK v1.12.0 :::...

    Marker-related
       parse_deltall_stdout   ->  Parse stdout of deltaLL.rb
       get_arCOG_seq          ->  Retrieve arCOG sequences
       AssessMarkerPA         ->  Assess Markers by P/A among groups
    
    Multiple Sequence Alignment
       ConvertMSA             ->  Convert MSA format
       fa2phy                 ->  Convert MSA format (fasta to phylip)
       OneLineAln             ->  One-line fasta format alignments
       SliceMSA               ->  Slice MSA by column 
    
    Tree-related
       GTDB_tree_r207         ->  Infer GTDB (r207) archaeal/bacterial tree
       subset_tree            ->  Subset tree
       compare_trees          ->  Compare trees with Mantel test
       rename_leaves          ->  Rename tree leaves
       FLN                    ->  Format leaf names (e.g. remove spaces in names)
    
    Dating-related
       AssessCVG              ->  Assess dating convergence
       CompareMCMC            ->  Compare MCMCTree outputs
       PlotMcmcNode           ->  distribution of node's age estimation 
       
    Dating workflow
       Marker2Tree            ->  Marker to Tree
       AssessMarkerDeltaLL    ->  Assess Markers by DeltaLL
       Dating                 ->  Perform dating