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
##' data(Psychiatry_disease_filtered)
##' data(pathways)
##' CPI_result = MAPE2.0(arraydata = Psychiatry_diseases$expr, clinical.data = Psychiatry_diseases$clinical,
##'                     label = "response",pmtx = NULL,pathway = c(Biocarta.genesets,GOBP.genesets,GOCC.genesets,GOMF.genesets,
##'                     KEGG.genesets,Reactome.genesets),data.type ="continuous", resp.type = "twoclass",method = "CPI",
##'                     ind.method = rep("limma",length(Psychiatry_diseases$expr)),paired = rep(FALSE,length(Psychiatry_diseases$expr)),
##'                     select.group = c("CASE","CTRL"),ref.level ="CTRL",tail="abs",
##'                     enrichment = "Fisher's exact", DEgene.number = 400,stat = "AW Fisher")
##' library(irr)
##' library(ConsensusClusterPlus)
##' set.seed(15213)
##' CPI.kappa_result = MAPE.Kappa(summary = CPI_result$summary,pathway = CPI_result$pathway,
##'                               max_k = 15, q_cutoff = 0.0005,software = CPI_result$method)

MAPE.Kappa <- function(summary, max_k = 10, pathway, software =c("CPI","MAPE"),
                       method = c("MAPE_I","MAPE_G","MAPE_P"),
                       q_cutoff = 0.1,output_dir = getwd())
  {
  b = summary #output of previous module
  method = match.arg(method)
  q.cut = q_cutoff

  if (software == "MAPE"){a = rownames(b[which(b[,paste(method,"_FDR",sep="")] < q.cut),])} #print(length(a))
  else if (software == "CPI"){a = rownames(b[which(b$q_value_meta < q.cut),])}
  else {stop("Please check: Wrong software.")}

  P = pathway[names(pathway) %in% a] #only keep the GO list that pass q value < 0.05

  g = toupper(unique(unlist(P))) #all genes
  k = t(sapply(P, function(r) as.integer(g %in% toupper(r)))); colnames(k) = g #for each GO, how many overlaps

   kappa.result = t(sapply(1:nrow(k), function(i){ #k for overlap
    print(i)
    y = sapply(1:nrow(k), function(j){
      kappa2(cbind(k[i,], k[j,]))$value
    })
    return(y)
  }))

  kappa.result[is.na(kappa.result)] = 0

  rownames(kappa.result) = colnames(kappa.result) = names(P)

  #MDS plot
  d = as.dist(1-kappa.result)

  fit = cmdscale(d, k = 2) # k is the number of dim
  MDS1 <- fit[,1]
  MDS2 <- fit[,2]

  title=output_dir
  results.20 = ConsensusClusterPlus(d,maxK=max_k,reps=500,pItem=0.8,
                                    #results.20 = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8, ??
                                    title= title, clusterAlg="pam",distance="pearson",
                                    seed=1262118388.71279,plot="png")

#    kappa.order<-kappa.result[results[[4]]$consensusTree$order,results[[4]]$consensusTree$order]
#    pdf("heatmap_kappa_consensus_cluster_order.pdf")
#    heatmap.2(kappa.order, density.info="none", trace="none", scale = "none", margins = c(6, 4),symm=F,symbreaks=F)
#    dev.off()
    return(list(kappa = kappa.result,method = method))
}


