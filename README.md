
## TreeSAK (A Swiss-Army-Knife for phylogenomic analysis)

[![pypi licence ](https://img.shields.io/pypi/l/TreeSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/TreeSAK.svg)](https://pypi.python.org/pypi/TreeSAK) 

Contact
---

+ Dr Weizhi Song
+ Department of Ocean Science, Hong Kong University of Science and Technology, Hong Kong
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

                 ...::: TreeSAK v1.28.0 :::...

    Marker-related
       ExtractMarkerSeq       ->  Extract marker by blastn  
       parse_deltall_stdout   ->  Parse stdout of deltaLL.rb
       get_arCOG_seq          ->  Retrieve arCOG sequences
       AssessMarkerPA         ->  Assess Markers by P/A among groups
       SplitScore             ->  Assess markers by split score
       AssessMarkerDeltaLL    ->  Assess Markers by DeltaLL
       OMA                    ->  Prepare input files for running OMA 
       OMA2                   ->  Filter OMA predicted OGs

    Multiple Sequence Alignment
       ConvertMSA             ->  Convert MSA format
       fa2phy                 ->  Convert MSA format (fasta to phylip)
       SingleLinePhy          ->  Put sequences in single line in phylip format
       OneLineAln             ->  One-line fasta format alignments
       SliceMSA               ->  Slice MSA by column 
       AlignmentPruner        ->  Remove heterogenous sites from MSA
       BMGE                   ->  Run BMGE
       ConcateMSA             ->  Concatenate MSA
    
    Tree-related
       MarkerSeq2Tree         ->  Marker sequence to tree
       MarkerRef2Tree         ->  Marker (reference sequence) to Tree
       GTDB_tree              ->  get GTDB tree
       PMSF                   ->  run iqtree with PMSF
       subset_tree            ->  Subset tree
       compare_trees          ->  Compare trees with Mantel test
       rename_leaves          ->  Rename tree leaves
       print_leaves           ->  print out tree leaves
       FLN                    ->  Format leaf names (e.g. remove spaces in names)
       ModifyTopo             ->  Modify tree topology
       ALE                    ->  Modules for running ALE
       GeneRax                ->  Run GeneRax (to be added)    
       RootTree               ->  Root tree with outgroup leaves
       
    Dating-related
       Dating                 ->  Perform molecular dating
       AssessCVG              ->  Assess dating convergence
       CompareMCMC            ->  Compare MCMCTree outputs
       PlotMcmcNode           ->  distribution of node's age estimation 
       VisHPD95               ->  HPD95 of estimated node age
       pRTC                   ->  Perform probabilistic RTC dating
