library(data.table)
# get all genes in miRwalk
miR_130a_3p_walk <- readxl::read_xls("miR_walk/130a_threeutr.xls") # 814 
miR_132_3p_walk <- readxl::read_xls("miR_walk/132_threeutr.xls") # 1673
miR_133a_3p_walk <- readxl::read_xls("miR_walk/133a_threeutr.xls") # 2312
miR_138_5p_walk <- readxl::read_xls("miR_walk/138_threeutr.xls") # 2817
miR_629_5p_walk <- readxl::read_xls("miR_walk/629_threeutr.xls") # 2028
length(unique(miR_629_5p_walk$genesymbol))

# get all genes in miRsystem
miR_130a_3p_sys <- read.table("miRsystem/Target Gene List [hsa-miR-130a-3p].csv", sep = ",", header = T) # 838
miR_132_3p_sys <- read.table("miRsystem/Target Gene List [hsa-miR-132-3p].csv", sep = ",", header = T) # 439
miR_133a_3p_sys <- read.table("miRsystem/Target Gene List [hsa-miR-133a-3p].csv", sep = ",", header = T) # 392
miR_138_5p_sys <- read.table("miRsystem/Target Gene List [hsa-miR-138-5p].csv", sep = ",", header = T) # 413
miR_629_5p_sys <- read.table("miRsystem/Target Gene List [hsa-miR-629-5p].csv", sep = ",", header = T) # 127
length(unique(miR_629_5p_sys$TARGET_GENE))

# get all genes in miRnet
miR_130a_3p_net <- read.table("miRnet/miR-130a-3p.csv", sep = ",", header = T) # 2003 genes
miR_132_3p_net <- read.table("miRnet/miR-132-3p.csv", sep = ",", header = T) # 1119 genes
miR_133a_3p_net <- read.table("miRnet/miR-133a-3p.csv", sep = ",", header = T) # 991 genes
miR_138_5p_net <- read.table("miRnet/miR-138-5p.csv", sep = ",", header = T) # 1004 genes
miR_629_5p_net <- read.table("miRnet/miR-629-5p.csv", sep = ",", header = T) # 369 genes
length(unique(miR_130a_3p_net$Target))


# combine all genes in three websites
miR_130a_3p_3all <- unique(c(miR_130a_3p_walk$genesymbol, miR_130a_3p_sys$TARGET_GENE, miR_130a_3p_net$Target)) # 3078
miR_132_3p_3all <- unique(c(miR_132_3p_walk$genesymbol, miR_132_3p_sys$TARGET_GENE, miR_132_3p_net$Target)) # 2866
miR_133a_3p_3all <- unique(c(miR_133a_3p_walk$genesymbol, miR_133a_3p_sys$TARGET_GENE, miR_133a_3p_net$Target)) # 3373
miR_138_5p_3all <- unique(c(miR_138_5p_walk$genesymbol, miR_138_5p_sys$TARGET_GENE, miR_629_5p_net$Target)) # 3424
miR_629_5p_3all <- unique(c(miR_629_5p_walk$genesymbol, miR_629_5p_sys$TARGET_GENE, miR_629_5p_net$Target)) # 2432 
length(unique(miR_629_5p_3all))

# kegg enrichment pathway analysis 
# transfer the genes symbol into geneID
tran_sym3all <- function(miR_3all_gene_ls){
  gene_tb <- as.data.frame(AnnotationDbi::select(org.Hs.eg.db,
                                                 keys = miR_3all_gene_ls,
                                                 columns = c("ENTREZID", "SYMBOL"),
                                                 keytype = "SYMBOL"))
  return(gene_tb$ENTREZID)
}

miR_130a_3p_3all_3ID <- tran_sym3all(miR_130a_3p_3all)
miR_132_3p_3all_3ID <- tran_sym3all(miR_132_3p_3all)
miR_133a_3p_3all_3ID <- tran_sym3all(miR_133a_3p_3all)
miR_138_5p_3all_3ID <- tran_sym3all(miR_138_5p_3all)
miR_629_5p_3all_3ID <- tran_sym3all(miR_629_5p_3all)

# kegg analysis 
do_kegg <- function(miR_3ID, filename_original, p, filename_p){
  kegg_path_res <- enrichKEGG(miR_3ID, organism = "hsa")@result
  kegg_path_res2 <- subset(kegg_path_res, kegg_path_res$pvalue < p)
  write.xlsx(kegg_path_res, filename_original)
  write.xlsx(kegg_path_res2, filename_p)
  return_ls <- list(original_path = kegg_path_res,
                    pfilter_path = kegg_path_res2)
  return(return_ls)
}

miR_130a_3p_3web_path_original <- do_kegg(miR_130a_3p_3all_3ID, "miR_130a_3p_3web_allgenes_path(original).xlsx", 0.05, "miR_130a_3p_3web_allgenes_path_p0.05.xlsx")$original_path
miR_130a_3p_3web_path_p0.05 <- do_kegg(miR_130a_3p_3all_3ID, "miR_130a_3p_3web_allgenes_path.xlsx", 0.05, "miR_130a_3p_3web_allgenes_path_p0.05.xlsx")$pfilter_path
miR_132_3p_3web_path_original <- do_kegg(miR_132_3p_3all_3ID, "miR_132_3p_3web_allgenes_path(original).xlsx", 0.05, "miR_132_3p_3web_allgenes_path_p0.05.xlsx")$original_path
miR_132_3p_3web_path_p0.05 <- do_kegg(miR_132_3p_3all_3ID, "miR_132_3p_3web_allgenes_path.xlsx", 0.05, "miR_132_3p_3web_allgenes_path_p0.05.xlsx")$pfilter_path
miR_133a_3p_3web_path_original <- do_kegg(miR_133a_3p_3all_3ID, "miR_133a_3p_3web_allgenes_path(original).xlsx", 0.05, "miR_133a_3p_3web_allgenes_path_p0.05.xlsx")$original_path
miR_133a_3p_3web_path_p0.05 <- do_kegg(miR_133a_3p_3all_3ID, "miR_133a_3p_3web_allgenes_path.xlsx", 0.05, "miR_133a_3p_3web_allgenes_path_p0.05.xlsx")$pfilter_path
miR_138_5p_3web_path_original <- do_kegg(miR_138_5p_3all_3ID, "miR_138_5p_3web_allgenes_path(original).xlsx", 0.05, "miR_138_5p_3web_allgenes_path_p0.05.xlsx")$original_path
miR_138_5p_3web_path_p0.05 <- do_kegg(miR_138_5p_3all_3ID, "miR_138_5p_3web_allgenes_path.xlsx", 0.05, "miR_138_5p_3web_allgenes_path_p0.05.xlsx")$pfilter_path
miR_629_5p_3web_path_original <- do_kegg(miR_629_5p_3all_3ID, "miR_629_5p_3web_allgenes_path(original).xlsx", 0.05, "miR_629_5p_3web_allgenes_path_p0.05.xlsx")$original_path
miR_629_5p_3web_path_p0.05 <- do_kegg(miR_629_5p_3all_3ID, "miR_629_5p_3web_allgenes_path.xlsx", 0.05, "miR_629_5p_3web_allgenes_path_p0.05.xlsx")$pfilter_path

# combine all significantly enriched pathways in kegg for 5 miRNAs target all genes from three webs
# combined_3web_allgenes_path_p0.05 is the table containing all kegg esignificantly enriched pathways for 5 miRNAs (P < 0.05)
combind_sigpath <- function(miR130_path, miR132_path, miR133a_path, miR138_path, miR629_path, filename){
  combined_all_path <- rbind(miR130_path, miR132_path, miR133a_path, miR138_path, miR629_path)
  combined_all_path$miRNA <- c(rep("miR_130a_3p", nrow(miR130_path)),
                               rep("miR_132_3p", nrow(miR132_path)),
                               rep("miR_133a_3p", nrow(miR133a_path)),
                               rep("miR_138_5p", nrow(miR138_path)),
                               rep("miR_629_5p", nrow(miR629_path)))
  write.xlsx(combined_all_path,filename )
  return(combined_all_path)
}
combined_3web_allgenes_path_p0.05 <- combind_sigpath(miR_130a_3p_3web_path_p0.05, 
                                               miR_132_3p_3web_path_p0.05,
                                               miR_133a_3p_3web_path_p0.05,
                                               miR_138_5p_3web_path_p0.05,
                                               miR_629_5p_3web_path_p0.05,
                                               "combined_3web_allgenes_path_p0.05.xlsx")
View(combined_3web_allgenes_path_p0.05 )

length(unique(combined_3web_allgenes_path_p0.05$Description)) # 177 pathways 

# reorder the combined_path table, geneID list in each row
# write a function to process the table, split the geneID into one row for later transfer into gene Symbol
# combined_3web_allgenes_path_p0.05_reorder_transed  is the table reorder the combined_3web_allgenes_path_p0.05
# and transfer the geneID into gene symbol, it contains all kegg significantly (P < 0.05) pathways results
# for 5 miRNAs for all genes from 3 webs
process_table <- function(original_table) {
  new_rows <- list()
  
  # loop through each row in the original table
  for (i in 1:nrow(original_table)) {
    # split the gene IDs
    gene_ids <- unlist(strsplit(as.character(original_table[i, "geneID"]), "/"))
    
    # create a new row for each gene ID
    for (gene_id in gene_ids) {
      new_rows[[length(new_rows) + 1]] <- c(original_table[i, "Description"], gene_id, original_table[i, "miRNA"], original_table[i, "count_number"])
      }
  }
  
  # convert the list of new rows into a data frame
  new_table <- as.data.frame(do.call(rbind, new_rows))
  colnames(new_table) <- c("pathway", "geneID", "miRNA")
  trans2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = new_table$geneID,
                                        columns = c("SYMBOL"),
                                        keytype = "ENTREZID"
  )$SYMBOL
  new_table_trans2symbol <- data.frame(pathway = new_table$pathway, 
                                       genesymbol = trans2symbol,
                                       miRNA = new_table$miRNA)
                        
  return_list <- list(new_table = new_table,
                      new_table_trans2symbol = new_table_trans2symbol)
  return(return_list)
}
combined_3web_allgenes_path_p0.05_reorder <- process_table(combined_3web_allgenes_path_p0.05)$new_table
combined_3web_allgenes_path_p0.05_reorder_transed <- process_table(combined_3web_allgenes_path_p0.05)$new_table_trans2symbol

# kegg network 1: for all target genes from three web for 5 miRNAs 
# remove all record which not including "pathway" or "Pathways" in pathway name
# for combined_3web_allgenes_path_p0.05_reorder_transed table
# combined_pathname_all_pathways_kegg_p0.05 is the table only have the pathway in pathway name 
unique_paths <- unique(combined_3web_allgenes_path_p0.05_reorder_transed$pathway)
length(unique_paths) # 177 unique pathways in this kegg enrichment results for all target genes got from three websites for 5 miRNAs 
                    # and they are all significantly enriched pathways (P < 0.05)

keep_pathname <- function(combined_path, filename){
  combined_pathname_td <- data.frame(
    pathway = character(),
    genesymbol = character(),
    miRNA = character()
  )
  for (i in unique_paths){
    pathway_subset <- as.data.frame(subset(combined_path, combined_path$pathway == as.character(i)))
    if_pathname <- any(strsplit(as.character(i), " ")[[1]] %in% c("pathway", "pathways", "Pathways", "Pathway"))
    if (if_pathname){
      combined_pathname_td <- rbind(combined_pathname_td, pathway_subset)
    }
  }
  write.xlsx(combined_pathname_td, filename)
  return(combined_pathname_td)
  
}
combined_pathname_all_pathways_kegg_p0.05 <- keep_pathname(combined_3web_allgenes_path_p0.05_reorder_transed, "combined_pathname_all_pathways_kegg_p0.05.xlsx")
View(combined_pathname_all_pathways_kegg_p0.05)






# calculate the interaction miRNA number for unique genes in network1
cal_internum <- function(combined_table, filename){
  unique_genesymb <- unique(combined_table$genesymbol)
  new_combined_table <- data.frame(
    pathway = character(),
    genesymbol = character(),
    miRNA = character(),
    interaction = character()
  )
  for (gene in unique_genesymb){
    gene_subset <- subset(combined_table, combined_table$genesymbol == gene)
    gene_subset$interaction <- length(unique(gene_subset$miRNA))
    new_combined_table <- rbind(new_combined_table, gene_subset)
  }
  write.xlsx(new_combined_table, filename)
  return(new_combined_table)
  
}
# network 1: the table for all pathway including "pathway" or "Pathways" inside the pathway name for three web all target genes kegg (p < 0.05) pathways
#            combined_pathname_all_pathways_kegg_p0.05_interaction including the interaction number of miRNA, you can select the genes been targeted by more than 1 miRNA here
#            based on the interaction column number
#            combined_pathname_all_pathways_kegg_p0.05 is the table only have the pathway in pathway name 
combined_pathname_all_pathways_kegg_p0.05_interaction <- cal_internum(combined_pathname_all_pathways_kegg_p0.05, 
                                                 "combined_pathname_all_pathways_kegg_p0.05_interaction.xlsx")


View(combined_pathname_all_pathways_kegg_p0.05_interaction)
length(unique(combined_pathname_all_pathways_kegg_p0.05_interaction$pathway)) # 51 pathways including pathway name 
length(unique(combined_pathname_all_pathways_kegg_p0.05_interaction$genesymbol)) # 1210 genes


# remove the interaction number is 1, 2,3, 4 for above table
interaction_select <- function(combined_interaction, big_interaction, filename){
  pathname_interbig_tb <- subset(combined_interaction,
                                 combined_interaction$interaction > big_interaction)
  return_list <- list(sele_tb =  pathname_interbig_tb,
                      interaction_genes = unique(pathname_interbig_tb$genesymbol),
                      interaction_pathway = unique(pathname_interbig_tb$pathway),
                      interaction_miRNA = unique(pathname_interbig_tb$miRNA))
  write.xlsx(pathname_interbig_tb, filename)
  return(return_list)
  
}
pathname_interbig1 <- interaction_select(combined_pathname_all_pathways_kegg_p0.05_interaction,
                                         1, "combined_pathname_all_pathways_kegg_p0.05_interaction1.xlsx")$sele_tb
# 448 genes for pathway including "pathway", P < 0.05, interaction >1 
# 50 pathways for pathway including "pathway", P < 0.05, interaction >1 
# 5 miRNAs for pathway including "pathway", P < 0.05, interaction >1 

pathname_interbig2 <- interaction_select(combined_pathname_all_pathways_kegg_p0.05_interaction,
                                         2, "combined_pathname_all_pathways_kegg_p0.05_interaction2.xlsx")$sele_tb
# 120 genes for pathway including "pathway", P < 0.05, interaction >2
# 49 pathways for pathway including "pathway", P < 0.05, interaction >2
# 5 miRNAs  for pathway including "pathway", P < 0.05, interaction >2


pathname_interbig3 <- interaction_select(combined_pathname_all_pathways_kegg_p0.05_interaction,
                                         3, "combined_pathname_all_pathways_kegg_p0.05_interaction3.xlsx")$sele_tb
 # 30 genes for pathway including "pathway", P < 0.05, interaction >3
 # 44 pathways for pathway including "pathway", P < 0.05, interaction >3 
 # 5 miRNAs  for pathway including "pathway", P < 0.05, interaction >3

pathname_interbig4 <- interaction_select(combined_pathname_all_pathways_kegg_p0.05_interaction,
                                         4, "combined_pathname_all_pathways_kegg_p0.05_interaction4.xlsx")$sele_tb
# 2 genes for pathway including "pathway", P < 0.05, interaction >4
# 11 pathways for pathway including "pathway", P < 0.05, interaction >4 
# 5 miRNAs  for pathway including "pathway", P < 0.05, interaction >4




# generate edge table and type table for network1 and network2
net_prepare <- function(interbig_tb, edge_filename, type_filename){
  edge_pathname1 <-  data.frame(source = interbig_tb$genesymbol,
                                target = interbig_tb$pathway)
  edge_pathname1 <- unique(edge_pathname1)
  edge_pathname2 <- data.frame(source = interbig_tb$miRNA,
                               target = interbig_tb$genesymbol)
  edge_pathname2 <- unique(edge_pathname2)
  edge_pathname <- rbind(edge_pathname1, edge_pathname2)
  
  nodename_ls <- c(unique(interbig_tb$miRNA), 
                   unique(interbig_tb$genesymbol),
                   unique(interbig_tb$pathway))
  type_ls <- c(rep("miRNA", length(unique(interbig_tb$miRNA))),
               rep("gene", length(unique(interbig_tb$genesymbol))),
               rep("pathway", length(unique(interbig_tb$pathway))))
  type_pathname <- data.frame(nodename = nodename_ls,
                              type = type_ls)
  write.xlsx(edge_pathname, edge_filename)
  write.xlsx(type_pathname, type_filename)
  return_ls <- list(edgefile = edge_pathname,
                    typefile = type_pathname)
  return(return_ls)
  
}
# network 1: pathway including in pathway name (interaction >3)
# all genes from 3 websites kegg results (p < 0.05) & pathway in the pathway name & interaction >3
edge_file_net_pathname_interbig3 <- net_prepare(pathname_interbig3,
                                      "edge_pathname_interbig3_file.xlsx",
                                      "type_pathname_interbig3_file.xlsx")$edgefile
type_file_net_pathname_interbig3 <- net_prepare(pathname_interbig3, 
                                      "edge_pathname_interbig3_file.xlsx", 
                                      "type_pathname_interbig3_file.xlsx")$typefile
# don't have long term depression, it only been target by miR-133a-3p which is been excluded by using interaction selection


# network 1.1: pathway including in pathway name (interaction >3)
# all genes from 3 websites kegg results (p < 0.05) & pathway in the pathway name & interaction >1
edge_file_net_pathname_interbig1 <- net_prepare(pathname_interbig1,
                                      "edge_pathname_interbig1_file.xlsx",
                                      "type_pathname_interbig1_file.xlsx")$edgefile
type_file_net_pathname_interbig1 <- net_prepare(pathname_interbig1, 
                                      "edge_pathname_interbig1_file.xlsx", 
                                      "type_pathname_interbig1_file.xlsx")$typefile
# don't have long term depression, it only been target by miR-133a-3p which is been excluded by using interaction selection









# GO analysis (biological process)
do_go <- function(gene_ls, go_path_filename_original, p,go_path_filename_p, top_num, go_path_filename_top){
          go_path <- enrichGO(gene = na.omit(gene_ls),
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH")@result
          go_path_p <- subset(go_path, go_path$pvalue < p)
          go_path_top <- head(go_path, top_num)
          return_list <- list(go_path = go_path, 
                              go_path_p = go_path_p, 
                              go_path_top = go_path_top)
          write.xlsx(go_path, go_path_filename_original)
          write.xlsx(go_path_p, go_path_filename_p)
          write.xlsx(go_path_top, go_path_filename_top)
          return(return_list)
}

miR_130a_3p_3web_allgenes_GO_path_original <- do_go(miR_130a_3p_3all_3ID, 
                                                     "miR_130a_3p_3web_allgenes_GO_path(original).xlsx", 
                                                     0.05, "miR_130a_3p_3web_allgenes_GO_path_p0.05.xlsx",
                                                     20,"miR_130a_3p_3web_allgenes_GO_path_top20.xlsx")$go_path



miR_132_3p_3web_allgenes_GO_path_original <- do_go(miR_132_3p_3all_3ID, 
                                                    "miR_132_3p_3web_allgenes_GO_path(original).xlsx", 
                                                    0.05,"miR_132_3p_3web_allgenes_GO_path_p0.05.xlsx",
                                                    20,"miR_132_3p_3web_allgenes_GO_path_top20.xlsx")$go_path

miR_133a_3p_3web_allgenes_GO_path_original <- do_go(miR_133a_3p_3all_3ID, 
                                                    "miR_133a_3p_3web_allgenes_GO_path(original).xlsx", 
                                                    0.05,
                                                    "miR_133a_3p_3web_allgenes_GO_path_p0.05.xlsx",
                                                    20,
                                                    "miR_133a_3p_3web_allgenes_GO_path_top20.xlsx")$go_path


miR_138_5p_3web_allgenes_GO_path_original <- do_go(miR_138_5p_3all_3ID, 
                                                    "miR_138_5p_3web_allgenes_GO_path(original).xlsx", 
                                                    0.05,
                                                    "miR_138_5p_3web_allgenes_GO_path_p0.05.xlsx",
                                                    20,
                                                    "miR_138_5p_3web_allgenes_GO_path_top20.xlsx")$go_path

miR_629_5p_3web_allgenes_GO_path_original <- do_go(miR_629_5p_3all_3ID, 
                                                   "miR_629_5p_3web_allgenes_GO_path(original).xlsx", 
                                                   0.05,
                                                   "miR_629_5p_3web_allgenes_GO_path_p0.05.xlsx",
                                                   20,
                                                   "miR_629_5p_3web_allgenes_GO_path_top20.xlsx")$go_path

miR_130a_3p_3web_allgenes_GO_path_top50 <- head(miR_130a_3p_3web_allgenes_GO_path_original, 50)
miR_132_3p_3web_allgenes_GO_path_top50 <- head(miR_132_3p_3web_allgenes_GO_path_original, 50)
miR_133a_3p_3web_allgenes_GO_path_top50 <- head(miR_133a_3p_3web_allgenes_GO_path_original, 50)
miR_138_5p_3web_allgenes_GO_path_top50 <- head(miR_138_5p_3web_allgenes_GO_path_original, 50)
miR_629_5p_3web_allgenes_GO_path_top50 <- head(miR_629_5p_3web_allgenes_GO_path_original, 50)

# get common pathway names for each possible pair or three or four or five miRNA
# miR_130a_3p_3web_allgenes_GO_path_original replace it using the GO result for each miRNA (original, P < 0.05 or top20)
combine_GO_path_reorder <- function(miR_130a_3web_path, miR_132_3web_path, miR_133a_3web_path, miR_138_3web_path, miR_629_3web_path, filename){
  combined_GO_path <- rbind(miR_130a_3web_path,
                            miR_132_3web_path,
                            miR_133a_3web_path,
                            miR_138_3web_path,
                            miR_629_3web_path)
  combined_GO_path$miRNA <- c(rep("miR_130a_3p", nrow(miR_130a_3web_path)),
                              rep("miR_132_3p", nrow(miR_132_3web_path)),
                              rep("miR_133a_3p", nrow(miR_133a_3web_path)),
                              rep("miR_138_5p", nrow(miR_138_3web_path)),
                              rep("miR_629_5p", nrow(miR_629_3web_path)))
  unique_GO_path <- unique(combined_GO_path$Description)
  
  combind_GO_pathway_5miRNA <- data.frame(
    ID = as.character(),
    Description = as.character(),
    GeneRatio = as.character(),
    BgRatio = as.character(),
    pvalue = as.character(),
    p.adjust = as.character(),
    qvalue = as.character(),
    geneID = as.character(),
    Count = as.character(), 
    miRNA = as.character(),
    allmiRNA = as.character()
    
  )
  
  for (path in unique_GO_path){
    testn <- subset(combined_GO_path, combined_GO_path$Description == path)
    testn$allmiRNA <- rep(paste(strsplit(testn$miRNA, " "), collapse = "/"), nrow(testn))
    combind_GO_pathway_5miRNA <- rbind(combind_GO_pathway_5miRNA, testn)
    
    
  }
  write.xlsx(combind_GO_pathway_5miRNA, filename )
  return(combind_GO_pathway_5miRNA)
  
}

combined_GO_pathway_5miRNA_top50 <- combine_GO_path_reorder(miR_130a_3p_3web_allgenes_GO_path_top50,
                              miR_132_3p_3web_allgenes_GO_path_top50,
                              miR_133a_3p_3web_allgenes_GO_path_top50,
                              miR_138_5p_3web_allgenes_GO_path_top50 ,
                              miR_629_5p_3web_allgenes_GO_path_top50,
                              "combind_GO_top50pathway_5miRNA_allgenes_3web.xlsx")

View(combined_GO_pathway_5miRNA_top50)
# transfer the geneID into genesymbol in GO genelist
for (i in 1:nrow(combined_GO_pathway_5miRNA_top50)) {
  
  gene_ids <- unlist(strsplit(combined_GO_pathway_5miRNA_top50[i, "geneID"], "/"))
  tes2 <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = gene_ids,
                                columns = c("SYMBOL"),
                                keytype = "ENTREZID" )$SYMBOL
 
  tes3 <- paste(strsplit(tes2, " "), collapse = "/")
  gene_symbols <- tes3
  combined_GO_pathway_5miRNA_top50[i, "geneID"] <- gene_symbols
  
}

write.xlsx(combined_GO_pathway_5miRNA_top50, "combind_GO_top50pathway_5miRNA_allgenes_3web_genesymbol.xlsx")




# get the selected neuro related pathway record for GO pathways resuts above (top50)
neuro_related_GO_pathway_fortop50combined <- read.xlsx("neuro_related_GO_pathway_fortop50combined.xlsx")
View(neuro_related_GO_pathway_fortop50combined )



# select the validated genes for above top50 pathways for GO 
# (validated genes come from three webs validated genes for 5 miRNAs)
# get all validated genes in miRwalk
miR_130a_3p_walk_val <- subset(miR_130a_3p_walk, !is.na(validated)) # 41 genes
miR_132_3p_walk_val <- subset(miR_132_3p_walk, !is.na(validated)) # 34 genes
miR_133a_3p_walk_val <- subset(miR_133a_3p_walk, !is.na(validated)) # 36 genes
miR_138_5p_walk_val <- subset(miR_138_5p_walk, !is.na(validated)) # 25 genes
miR_629_5p_walk_val <- subset(miR_629_5p_walk, !is.na(validated)) # 14 genes
# get all validated genes in miRsystem
miR_130a_3p_sys_val <- subset(miR_130a_3p_sys, miR_130a_3p_sys$VALIDATION == "V") # 6 genes
miR_132_3p_sys_val <- subset(miR_132_3p_sys, miR_132_3p_sys$VALIDATION == "V") # 3 genes
miR_133a_3p_sys_val <- subset(miR_133a_3p_sys, miR_133a_3p_sys$VALIDATION == "V") # 13 genes
miR_138_5p_sys_val <- subset(miR_138_5p_sys, miR_138_5p_sys$VALIDATION == "V") # 4 genes
miR_629_5p_sys_val <- subset(miR_629_5p_sys, miR_629_5p_sys$VALIDATION == "V") # 0 genes

# combine all validated genes for each miRNA
miR_130_val <- unique(c(miR_130a_3p_walk_val$genesymbol, miR_130a_3p_sys_val$TARGET_GENE, miR_130a_3p_net$Target))
miR_132_val<- unique(c(miR_132_3p_walk_val$genesymbol, miR_132_3p_sys_val$TARGET_GENE, miR_132_3p_net$Target))
miR_133_val <- unique(c(miR_133a_3p_walk_val$genesymbol, miR_133a_3p_sys_val$TARGET_GENE, miR_133a_3p_net$Target))
miR_138_val <- unique(c(miR_138_5p_walk_val$genesymbol, miR_138_5p_sys_val$TARGET_GENE, miR_138_5p_net$Target))
miR_629_val <- unique(c(miR_629_5p_walk_val$genesymbol, miR_629_5p_sys_val$TARGET_GENE, miR_629_5p_net$Target))

# Define a function to process miRNA subsets
process_miRNA_subset <- function(miRNA_name, gene_val) {
  sub <- subset(neuro_related_GO_pathway_fortop50combined, 
                neuro_related_GO_pathway_fortop50combined$miRNA == miRNA_name)
  n <- nrow(sub)
  for (i in 1:n) {
    gene_ls <- unlist(strsplit(sub[i, "geneID"], "/"))
    val_subgene <- intersect(gene_ls, gene_val)
    sub[i, "geneID"] <- paste(val_subgene, collapse = "/")
  }
  return(sub)
}

# Define miRNA names and corresponding gene values
miRNAs <- c("miR_130a_3p", "miR_132_3p", "miR_133a_3p", "miR_138_5p", "miR_629_5p")
gene_vals <- list(miR_130_val, miR_132_val, miR_133_val, miR_138_val, miR_629_val)

# Process each miRNA subset and combine them
combined_common_validated_gene_GO_top50_neuro_related <- data.frame()
for (j in 1:length(miRNAs)) {
  miRNA <- miRNAs[j]
  gene_val <- gene_vals[[j]]
  miRNA_subset <- process_miRNA_subset(miRNA, gene_val)
  combined_common_validated_gene_GO_top50_neuro_related <- rbind(combined_common_validated_gene_GO_top50_neuro_related, miRNA_subset)
}
write.xlsx(combined_common_validated_gene_GO_top50_neuro_related, "combined_common_validated_gene_GO_top50_neuro_related.xlsx")
View(combined_common_validated_gene_GO_top50_neuro_related_test)

# get common genes appear in same pathway
# split the genes and get the pathway name, miRNA, and form a new table
for (m in 1:nrow(combined_common_validated_gene_GO_top50_neuro_related)){
  gene_list2 <- unlist(strsplit(combined_common_validated_gene_GO_top50_neuro_related[n, "geneID"], "/"))
  repeat_num <- length(gene_list2)
  
  
}

combined_common_validated_gene_GO_top50_neuro_related_reorder <- data.frame()  
for (n in 1:nrow(combined_common_validated_gene_GO_top50_neuro_related)) {
  gene_list2 <- unlist(strsplit(combined_common_validated_gene_GO_top50_neuro_related[n, "geneID"], "/"))
  repeat_num <- length(gene_list2)
  
  # Repeat Description and miRNA columns for the number of genes
  repeated_desc <- rep(combined_common_validated_gene_GO_top50_neuro_related[n, "Description"], repeat_num)
  repeated_mirna <- rep(combined_common_validated_gene_GO_top50_neuro_related[n, "miRNA"], repeat_num)
  
  # Create a new dataframe for the current row
  new_rows <- data.frame(
    Description = repeated_desc,
    geneID = gene_list2,
    miRNA = repeated_mirna
  )
  
  combined_common_validated_gene_GO_top50_neuro_related_reorder <- rbind(combined_common_validated_gene_GO_top50_neuro_related_reorder, new_rows)
}
write.xlsx(combined_common_validated_gene_GO_top50_neuro_related_reorder, "combined_common_validated_gene_GO_top50_neuro_related_reorder.xlsx")
View(combined_common_validated_gene_GO_top50_neuro_related_reorder)













