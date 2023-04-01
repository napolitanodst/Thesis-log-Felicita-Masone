# After extracting the table of modeling results in long format with sig_genes_extract_all(), 
# we can then use layer_boxplot() to visualize any gene of interest for any of the model types and tests.

## Visualize the expression of FREM3 

set.seed(20200206)
layer_boxplot(
  i = which(sig_genes$gene == "FREM3")[1], 
  sig_genes = sig_genes,
  sce_layer = sce_layer,
  short_title = FALSE,
  col_low_box = "palegreen3",       # Box Background color for layers with the expected lower expression 
  col_low_point = "springgreen2",
  col_high_box = "darkorange2",     # Box Background color for layers with the expected higher expression
  col_high_point = "orange1"
)

## > which(sig_genes$gene == "FREM3")
# [1]   20682   27289   45280   67000   97054  132117  155990  178303  200651  223125  244193  264955
# [13]  285734  310119  333230  340279  358618  380178  414704  424869  446827  469049  491341  513663
# [25]  535981  558846  581117  606227  625614  647928  670116  693710  717610  741493  761770  783321
# [37]  820934  847257  870359  880495  914992  937696  960136  982506 1004846 1027190 1048987 1071378
# [49] 1090930 1095034 1116989

## > sig_genes[20682,]
# DataFrame with 1 row and 12 columns
#   top     model_type      test        gene      stat        pval         fdr    gene_index      ensembl            in_rows      in_rows_top20       results
# <integer> <character> <character> <character> <numeric>   <numeric>   <numeric>  <integer>     <character>     <IntegerList>    <IntegerList>   <CharacterList>
#    20682  enrichment          WM       FREM3  -5.13202 2.10782e-06   1.43156e-05    5641   ENSG00000183090   20682,27289,45280,...     67000         Layer3_top7
