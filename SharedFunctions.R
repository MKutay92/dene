## Functions to load ##
# You will need these function when generating and using both the Pocodyte
# Damage Score and the Hepatocyte Damage Score

# damage_signature.func()

# I. Function for generating List of Damage Associated Genes 

# Output: Table of genes by average rank, in ascending order. 
# Contains the following columns:
# gene_symbol - conventional Gene names, mean_rank - average rank  and 
# direction_foldchange - common direction of fold change, where '-1' and '1'
# represent negative and positive fold changes, respectively.

damage_signature.func <- function( DElist , thresh= 0.75 ,
                                   geneFilter )
{
  require(matrixStats)
  ### *file_DElist* argument specifies an input: a .rda file, 
  # containing a list of DE tables, 
  # each table represents an experiment, and necessarily contains:
  # a column "Gene.Symbol" with conventional gene names
  # a column "log2FoldChange" with log2 fold changes, 
  # !!! experimental samples are compared against controls
  # a column "pvalue" with raw p-values
  
  ### *thresh* argument sets a threshold on proportion of studies in which 
  # a) a gene should be detected, and then
  # b) the gene should have the same direction of expression changes upon damage,
  # to be considered as a damage marker. 
  
  ### geneFilter are genes that should be considered, 
  # e.g. only genes detected in sc studies
  
  # extract gene names from all tables 
  allGenes <- Reduce( union, lapply( DElist, function(X) X$Gene.Symbol ) )
  
  # pass through a specific gene filter
  allGenes <- allGenes[ allGenes %in% geneFilter ]
  
  # create table p-value per gene for each study
  allStud_pval <- Reduce( cbind, lapply(DElist, function(X) 
    X$pvalue[match(allGenes, X$Gene.Symbol)]))
  
  # create table lfc per gene for each study
  allStud_lfc <- Reduce(cbind, lapply(DElist, function(X)
    X$log2FoldChange[match(allGenes, X$Gene.Symbol)]))
  
  rownames(allStud_pval) <- rownames(allStud_lfc ) <- allGenes
  
  
  # first sub-setting step: gene should be present in at least 75% of the studies 
  allStud_pval <- allStud_pval[rowSums( !is.na(allStud_pval)) >= thresh*ncol(allStud_pval),]
  allStud_lfc <- allStud_lfc[rownames( allStud_lfc ) %in% rownames(allStud_pval),]
  # second subs-setting step
  allStud_lfc <- allStud_lfc[rowSums(allStud_lfc > 0, na.rm = TRUE) >= thresh*ncol(allStud_lfc) |
                               rowSums(allStud_lfc <0, na.rm = TRUE )>= thresh*ncol(allStud_lfc) ,]
  allStud_pval <- allStud_pval[rownames(allStud_pval) %in% rownames(allStud_lfc
  ),]
  
  # generate ranking 
  # be aware, lowest p-value should be on top of list (rank 1)
  
  allStud_ranks <- apply( allStud_pval,2, rank )
  
  # output table: gene symbol / averaged rank / positive or negative fold change 
  final_table <- data.frame('gene_symbol' = rownames(allStud_ranks), 
                            'mean_rank' = rowMeans(allStud_ranks, na.rm = TRUE),
                            'median_rank' = rowMedians(allStud_ranks, na.rm = TRUE),
                            'direction_foldchange' = sign(rowMeans(apply(allStud_lfc,2,sign),
                                                                   na.rm = TRUE)))
  
  final_table <- final_table[ order( final_table$mean_rank, decreasing=FALSE,
                                     na.last=TRUE), ]
  
  return(final_table)
}


# II. Function for calculating Damage Score (PDS or HDS) using AUCelll

DS_calc.func <- function(exprMatrices , DSignature , 
                         geneIDname = "gene_symbol", 
                         useOrder=c( "mean_rank", "median_rank"),
                         ntop= 42 , ceilThrsh=0.05 , progStat =F ,
                         wghtd=T )
  {
  require(AUCell)
  require(GSEABase)
  
  ## exprMatrices - expression matrix with genes in rows and samples/cells
  # in columns, row names ID type should be the same as in Signature[geneIDname,]
  ## DSignature - an output of damage_signature.func - a dataframe containing columns gene_symbol, mean_rank and direction_foldchange
  ## geneIDname - a charachter or a numeric value, specifying which 
  # DSignature column to use for gene names 
  # useOrder - whether to use mean_rank or median_rank columns of DSignature for ordering genes
  # nTop - how many genes from the top of a gene signature should be used to calculate the score
  # ceilThrsh - a threshold to calculate the AUC
  # progStat - to print progress messages or not
  # genesUP and genesDOWN are charachter vectors containing the names of genes UP and DOWN regulated upon damage 
  ## wghtd - logical indicating whether score for up and down-regulated genes 
  # should be weghted by the number of genes
  
  ## evaluate choices
  useOrder <- match.arg(useOrder)
  
  # order genes according to mean or median
  DSignature <- DSignature[ order(DSignature[[useOrder]] ),]
  
  ### make a collection of gene sets  to use with AUCell
  genesUP = DSignature[ 1:ntop, geneIDname][ 
    which( DSignature$direction_foldchange[1:ntop]==1)]
  genesDOWN =  DSignature[ 1:ntop, geneIDname][ 
    which( DSignature$direction_foldchange[1:ntop]== -1)]
  
geneSets <-  GSEABase::GeneSetCollection( list( GeneSet( genesUP , 
                                               setName= "UP" ) , 
                                      GeneSet( genesDOWN , 
                                               setName= "DOWN" )) )

### uppercase row. names
exprMatrices <- as.matrix(exprMatrices)

###  1. Build gene-expression rankings for each cell  
cellRanks <- AUCell_buildRankings(  exprMatrices, nCores=1, plotStats=F, 
                                    verbose = progStat )

###  2. Calculate enrichment for the gene signatures (AUC)
# gene sets should include UP and DOWN regulated genes seperately
cells_AUC <- AUCell_calcAUC( geneSets, cellRanks , verbose = progStat,
                             # aucMaxRank parameter is a tricky one,
                             # for single cell 0.05 is a reasonable 
                             # but for bulk a bigger, e.g. 0.2 can be used
                             aucMaxRank = ceiling( ceilThrsh * nrow(cellRanks))
                             )
cells_AUC <- getAUC( cells_AUC) #this will extract results 

### calculate the damage score
# indeces specify elements of geneSets that contain UP and DOWN regulated genes, corespondingly
index1 <- 1
index2 <- 2
# account for possible empty UP or DOWN gene sets
if( isTRUE(wghtd)) {
  coef1 <- sum( genesUP%in%rownames(exprMatrices))/ntop
  coef2 <- sum( genesDOWN%in%rownames(exprMatrices))/ntop
} else coef1 <- coef2 <- 1
  
if( nrow(cells_AUC)==2 ) { 
  Dscore = cells_AUC[index1,] * coef1 - cells_AUC[index2,] * coef2
  } else if(rownames(cells_AUC)=="UP") {
    Dscore= cells_AUC[index1,] } else { Dscore= - cells_AUC[index1,] }

return(Dscore)
}
                                      

