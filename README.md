
## TreeSAK (A Swiss-Army-Knife for phylogenomic analysis)

[![pypi licence ](https://img.shields.io/pypi/l/TreeSAK.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version ](https://img.shields.io/pypi/v/TreeSAK.svg)](https://pypi.python.org/pypi/TreeSAK) 


Contact
---

+ Dr. Weizhi Song
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

                 ...::: TreeSAK v1.58.0 :::...

                 
    Marker-related
       ExtractMarkerSeq       ->  Extract marker by blastn  
       deltall                ->  Parse stdout of deltaLL.rb
       get_arCOG_seq          ->  Retrieve arCOG sequences
       AssessMarkerPA         ->  Assess Markers by P/A among groups
       SplitScore             ->  Assess markers by split score
       AssessMarkerDeltaLL    ->  Assess Markers by DeltaLL
       OMA                    ->  Prepare input files for running OMA 
       OMA2                   ->  Filter OMA predicted OGs
       filter_rename_ar53     ->  Filter rename GTDB markers
       
    Multiple Sequence Alignment
       BMGE                   ->  Run BMGE
       pruneMSA               ->  Prune MSA with alignment_pruner.pl
       recode                 ->  Recode amino acids to Dayoff 4, Dayoff 6 or SR4 categories
       fa2phy                 ->  Convert MSA format (fasta to phylip)
       phy2fa                 ->  Convert MSA format (phylip to fasta)
       SliceMSA               ->  Slice MSA by column
       ConcateMSA             ->  Concatenate MSAs
       ConvertMSA             ->  Convert MSA format
       OneLineAln             ->  One-line fasta format alignments
       SingleLinePhy          ->  Put sequences in single line in phylip format
       CS_trim                ->  (to be added) perform chi-squared trimming to reduce compositional heterogeneity
       gap_stats              ->  The percentage of gap in each sequence of a MSA

    Tree-related
       iTOL                   ->  Prepare iTOL files
       batch_itol             ->  Batch access iTOL
       iTOL_gene_tree         ->  Genome metadata to gene metadata
       iTOL_msa_stats         ->  iTOL_msa_stats
       PB                     ->  Infer tree with PhyloBayes-MPI 
       assessPB               ->  Compare PhyloBayes chains
       supertree              ->  Infer species tree from multiple gene trees 
       MarkerSeq2Tree         ->  Marker sequence to tree
       MarkerRef2Tree         ->  Marker (reference sequence) to Tree
       GTDB_tree              ->  Get GTDB tree
       subset                 ->  Subset tree
       TaxonTree              ->  Subset GTDB tree by taxon
       compare_trees          ->  Compare trees with Mantel test
       rename_leaves          ->  Rename tree leaves
       print_leaves           ->  print out tree leaves
       FLN                    ->  Format leaf names (e.g. remove spaces in names)
       ModifyTopo             ->  Modify tree topology
       GeneRax                ->  (to be added) Run GeneRax    
       ALE                    ->  Modules for running ALE
       cogTree                ->  Infer tree for individual COG function
       koTree                 ->  Infer tree for individual KEGG function
       RootTree               ->  Root tree with outgroup leaves
       RootTreeGTDB           ->  Root tree by GTDB taxonomy
       LcaToLeaves            ->  Get two leaves that define an internal node
       replace_clade          ->  Replace tree clade
       GeneTree               ->  Infer gene tree
       guide_tree             ->  Prepare guide tree for iqtree
       get_lca_id             ->  Get LCA id
       tree_fmt               ->  Change tree format
       rm_leaf                ->  Remove leaf/leaves from tree
       
    Model-related
       PMSF                   ->  run iqtree with PMSF
       PPA                    ->  (to be added) Perform Posterior Predictive Analysis (across-site)
       
    Dating-related
       dating                 ->  Perform molecular dating
       CompareMCMC            ->  Compare MCMCTree outputs
       PlotMcmcNode           ->  Distribution of node's age estimation 
       VisHPD95               ->  HPD95 of estimated node age
       pRTC                   ->  Perform probabilistic RTC dating
       mcmcTC                 ->  Adding time constraints to mcmctree tree
       mcmc2tree              ->  Get the tree with internal node from mcmctree output
       mcmctree_out           ->  Parse MCMCTree produced .out file
       
    Phylo-related stats
       PhyloBiAssoc            ->  A wrapper for binaryPGLMM test
