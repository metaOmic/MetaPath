##' This is the second major function in the MetaPath2.0 package which calculates the pairwise
##' cohen's kappa statistics between significant pathway analysis results.

##' @title Calculate the kappa statistics for DE Pathways
##' @param method The meta qvalue used for qvalue cutoff. It is one of the three values:
##' 'MAPE_I','MAPE_G','MAPE_P'.
##' @param q_cutoff All the pathways with qvalue smaller than q_cutoff will be used for
##' kappa statistics calculation.
##' @return The matrix of kappa statistics between each DE pathways
##' @author Chien-wei Lin and George Tseng.
##' @export
##' @examples
##' data(pathways)
##' data(MetaDEresult)
##' data(MetaDEresult)
##' MAPE2.0_result = MAPE2.0(stat='maxP',method = "MAPE", enrichment = "KS", 
##'                         size.min=15,size.max=500, MetaDE = TRUE,
##'                         meta.p = meta.res.p$meta.analysis$pval,ind.p = ind.res$p)
##' MAPE.kappa_result = MAPE.Kappa(summary = MAPE2.0_result$summary, software = MAPE2.0_result$method,
##'                                pathway = MAPE2.0_result$pathway, max_k = 20, q_cutoff = 0.1,
##'                                output_dir = tempdir())


MAPE.Kappa <- function(summary, max_k = 20, pathway, software =c("CPI","MAPE"),
                       method = c("MAPE_I","MAPE_G","MAPE_P"),
                       q_cutoff = 0.1,output_dir = getwd())
  {
  
  #ConsensusClusterPlus
  # source("http://bioconductor.org/biocLite.R")
  #biocLite("ConsensusClusterPlus")
  #biocLite("epicalc")
  #biocLite("limma")
  #biocLite("topGO")
  #library(limma)
  #library(topGO)
  #library(ConsensusClusterPlus)
  #library(gplots)
  #library(ggplot2)
  #install.packages("cluster")
  #biocLite("org.Hs.eg.db")
  #library(cluster)
  #library(org.Hs.eg.db)
  #library(irr)
  # source("./src/epicalc/R/epicalc.R")
  # source("epicalc.R")
  #install.packages("ppam_1.0.0.tar.gz", repos = NULL, type="source") #this need gfortran compiler
  #   require(ppam)
  #require(epicalc) # ?? epicalc was removed from CRAN
  
  
  #calculate kappa statistic
  # a = load("GO_gene_set_long_list.rdata"); a
  ####  a = load("GO_gene_set_ALL_list.rdata"); a
  # a = load(GO_list_input); a
  # GO_list <- GO_BP_list
  # GO_list <- GO_list_input
  # GO_list <- pathway
  
  # ?? which(names(GO_BP_list) == "GO:0043087")
  # a = unlist(read.delim("g1_top30.txt", as.is = T, header = F))
  # b = read.csv("fisher_summary.csv", as.is = T, header = T)
  ####  b = read.csv("fisher_sum_complete_BP_list.csv", as.is = T, header = T) #output of previous module
  b = summary #output of previous module
  method = match.arg(method)
  # q.cut = .05 ## ??
  q.cut = q_cutoff
  
  if (software == "MAPE"){a = rownames(b[which(b[,method] < q.cut),])} #print(length(a))
  else if (software == "CPI"){a = rownames(b[which(b$q_value_meta < q.cut),])}
  else {stop("Please check: Wrong software.")}
  #gene list of interest
  # all(a %in% names(pathway))
  P = pathway[names(pathway) %in% a] #only keep the GO list that pass q value < 0.05
  
  g = toupper(unique(unlist(P))) #all genes
  k = t(sapply(P, function(r) as.integer(g %in% toupper(r)))); colnames(k) = g #for each GO, how many overlaps
  
  #sfInit(parallel = T, cpus = Num_of_CPU); sfLibrary(irr)
  # kappa.result = t(sfSapply(1:nrow(k), function(i){ #k for overlap
  kappa.result = t(sapply(1:nrow(k), function(i){ #k for overlap
    print(i)
    # x = factor(k[i,], levels = 1:0)
    # x = k[i,] 
    y = sapply(1:nrow(k), function(j){
      #   y = sapply(11:15, function(j){
      # z = factor(k[j,], levels = 1:0)
      # z = k[j,]
      kappa2(cbind(k[i,], k[j,]))$value
      # table.k = table(x, z)
      #     print(table.k)
      # kap(table.k)$kappa
    })
    return(y)
  }))
  #  sfStop()
  kappa.result[is.na(kappa.result)] = 0
  
  rownames(kappa.result) = colnames(kappa.result) = names(P)
  #heatmap(kappa.result, scale = "none")
  
  # i = which(names(P) == "GO:0007275"); x = factor(k[i,], levels = 1:0)
  # j = which(names(P) == "GO:0030154"); z = factor(k[j,], levels = 1:0)
  
  # table(x, z)
  
  #MDS plot
  d = as.dist(1-kappa.result)
  
  fit = cmdscale(d, k = 2) # k is the number of dim
  MDS1 <- fit[,1]
  MDS2 <- fit[,2]
  # pdf(paste(output_dir,"/Secondary_files/mds_plot.pdf",sep=""))
#  pdf(paste(output_dir,"/Secondary_files/mds_plot.pdf",sep=""))
#  plot(MDS1, MDS2, xlab="Coordinate 1", ylab="Coordinate 2", main="Pathway visualization")
#  dev.off()
  #text(MDS1, MDS2, labels = rownames(kappa.result), cex=.5, xpd = NA)
  
  title=output_dir
  results.20 = ConsensusClusterPlus(d,maxK=max_k,reps=500,pItem=0.8,
                                    #results.20 = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8, ?? 
                                    title= title, clusterAlg="pam",distance="pearson",
                                    seed=1262118388.71279,plot="png")
##############  
############## output "consensus021.png" and "consensus022.png" in title
##############
  
    #   k=Num_Clusters
  #   results <- results.20
  # 
  #   results.pam <- pam(d, k, diss=T)
  # length(which(results.pam$clustering==4))  
  # names(which(results.pam$clustering==4))  
  
  
  #plot(MDS1, MDS2, xlab="Coordinate 1", ylab="Coordinate 2", main="Pathway visualization", col = results[[k]]$consensusClass)
  #text(MDS1, MDS2, labels = rownames(kappa.result), cex=.5, xpd = NA, col = results[[k]]$consensusClass)
  # results[[4]]$consensusTree$order
  
  if(0){
    kappa.order<-kappa.result[results[[4]]$consensusTree$order,results[[4]]$consensusTree$order]
    pdf("heatmap_kappa_consensus_cluster_order.pdf")
    heatmap.2(kappa.order, density.info="none", trace="none", scale = "none", margins = c(6, 4),symm=F,symbreaks=F)
    #based on delta plot, choose number of clusters k = 4
    dev.off()
  }
  #par(mar=c(0.5,0.5,0.5,0.5))
  #par(mar=c(5.1,4.1,4.1,2.1))
  # names(results[[k]])
  # results[[k]]$consensusClass
  
  #   cluster.id = function(id){
  #     id = paste("GO:", id, sep = "")
  #     results[[k]]$consensusClass[names(results[[k]]$consensusClass) == id]
  #   }
  #   
  #   names(which(results[[k]]$consensusClass == 1))
  #   names(which(results[[k]]$consensusClass == 3))
  #   names(which(results[[k]]$consensusClass == 4))
  # length(which(results[[k]]$consensusClass == 1))
  # length(which(results[[k]]$consensusClass == 3))
  # length(which(results[[k]]$consensusClass == 4))
  # length(which(results[[k]]$consensusClass == 2)) 
  
  
  
  #  cc <- list()
  #  for (i in 3:k){
  #    cc[[i]] <- results[[i]]$consensusClass
  #    write.table(cc[[i]], paste("ConsensusCluster_k",i,".txt",sep=""), quote=F, col.names=F, row.names=F)
  #  }
  #  GO:0008366
  #  GO:0042552
  #  GO:0009987
  
#  if(0){
#    penalized_k_medoids<-ppam(d,4,1,diss = T) ##?? no change when changing lambda?
#    length(which(penalized_k_medoids$clustering==-1))
#    length(which(penalized_k_medoids$clustering==2))
#    length(which(penalized_k_medoids$clustering==3))
#    #ppam  
#    
#    cat("penalized k-medoids cluster 1\n", file = "Penalized_k_medoids_clusters.csv")
#    write.table(names(which(penalized_k_medoids$clustering==1)), "Penalized_k_medoids_clusters.csv", sep=",",quote=F, append = T, row.names=F,col.names =F)
#    for (i in 2:k){
#      cat(paste("\npenalized k-medoids cluster ", i, "\n", sep = ""), file = "Penalized_k_medoids_clusters.csv", append = T)
#      write.table(names(which(penalized_k_medoids$clustering==i)), "Penalized_k_medoids_clusters.csv", sep=",",quote=F, append = T, row.names=F,col.names=F)
#    } 
#  }
  
  #par(mfrow = c(2, 1))
  #venn diagram
  
  
  #??  index.k = which(names(P) %in% c("GO:0008366", "GO:0042552", "GO:0009987"))
  #??  a <- vennCounts(t(k[index.k,]))
  #??  vennDiagram(a, mar = c(0, 0, 0, 0), cex = 1)
  
  return(list(kappa = kappa.result,method = method))
  
}


