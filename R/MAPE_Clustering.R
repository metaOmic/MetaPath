##' This is the third major function in the MetaPath2.0 package which outputs the pathway
##' cluster results with text mining anotation.

##' @title Output pathway cluster results
##' @param Num_Clusters number of clusters
##' @return The pathway cluster results in a csv file.
##' @authors Zhou Fang, Xiangrui Zeng and George Tseng.
##' @export
##' @examples
##' data(pathways)
##' data(MetaDEresult)
##' data(MetaDEresult)
##' MAPE2.0_result = MAPE2.0(stat='maxP',method = "MAPE", enrichment = "KS", 
##'                         size.min=15,size.max=500, MetaDE = TRUE,
##'                         meta.p = meta.res.p$meta.analysis$pval,ind.p = ind.res$p)
##' MAPE.kappa_result = MAPE.Kappa(summary = MAPE2.0_result$summary, software = MAPE2.0_result$method,
##'                               pathway = MAPE2.0_result$pathway, max_k = 20, q_cutoff = 0.05,
##'                               output_dir = tempdir())
##' genelist = NULL
##' if (MAPE2.0_result$enrichment == "Fisher's exact" & MAPE2.0_result$method == "CPI"){
##'  genelist = MAPE2.0_result$genelist
##' }
##'
##' MAPE.Clustering(summary = MAPE2.0_result$summary,Num_Clusters=10,
##'                Num_of_gene_lists = MAPE2.0_result$Num_of_gene_lists,
##'                genelist = genelist,kappa.result = MAPE.kappa_result$kappa, 
##'                pathway = MAPE2.0_result$pathway, enrichment = MAPE2.0_result$enrichment,
##'                method = MAPE.kappa_result$method,software = MAPE2.0_result$method,
##'                output_dir = tempdir())

MAPE.Clustering <- function(summary,Num_Clusters = 3, kappa.result = kappa.result, Num_of_gene_lists, 
                            genelist = NULL,pathway = all.pathway, enrichment,method,software,
                            output_dir = getwd())
  {
  #module2
  #biocLite("Rgraphviz")
  #biocLite("AnnotationDbi")
  #install.packages("ggplot2")
  #  install.packages("gplots")
  #install.packages("shape")
  #library(cluster)
  #library(gplots)
  #library(ggplot2)
  #library(AnnotationDbi)
  #library(Rgraphviz)
  #library(shape)
  #require(ppam)
  
  #load("Hallmark.RData")
  #load("Positional.RData")
  #load("CGP.RData")
  #load("Canonical.RData")
  #load("Biocarta.RData")
  #load("KEGG.RData")
  #load("Reactome.RData")
  #load("MIR.RData")
  #load("TFT.RData")
  #load("CGN.RData")
  #load("CM.RData")
  #load("GOBP.RData")
  #load("GOCC.RData")
  #load("GOMF.RData")
  #load("Oncogenic.RData")
  #load("Immunologic.RData")
  #load("PID.RData")
  #load("CMAP_up.RData")
  #load("CMAP_down.RData")
  #load("PPI.RData")
  #load("JASPARhuman.RData")
  #load("PITA.RData")
  #load("miRanda.RData")
  #load("TargetScan.RData")
  #load("Phenocarta.RData")
  
  # load('GO_gene_set_complete_list.RData')
  k = Num_Clusters
  d = as.dist(1-kappa.result)
  # results <- results.20
  results <- pam(d, k, diss=T)
  #??pam  
  #b = read.csv("fisher_sum_complete_BP_list.csv", as.is = T, header = T)
  b = summary
#  if(Num_of_gene_lists==1){
#    b =cbind(fisher_sum_tmp,fisher_sum_tmp[,3])
#  }
  
  dir.create(paste(output_dir,"/Clustering_files",sep=""))
  sil <- silhouette(results$clustering, (1-kappa.result), diss=T)
  summary(sil)
  pdf(paste(output_dir,"/Clustering_files/silhouette_plot.pdf",sep=""))
  plot(sil, nmax= 80, cex.names=0.6)
  dev.off()
  length(results$clustering)
  dim(d)
  pdf(paste(output_dir,"/Clustering_files/silhouette_width.pdf",sep=""))
  hist(sil[,3],breaks=30)
  dev.off()  
  
  k <- k+1
  sil_cut <- 0.1
  results2 <- results$clustering
  new.kappa.result<-kappa.result
  i<-0
  while(min(sil[,3]) < sil_cut){
    results2<-results2[-which(sil[,3] == min(sil[,3]))]
    for (i in unique(results2)){
      if(length(which(results2==i))==1){
        results2<-results2[-which(results2==i)]
      }
    }
    results2temp<-results2
    for(d in 1:length(results2)){
      results2[d]<-rank(unique(results2temp))[which(unique(results2temp)==results2temp[d])]
    }
    
    new.kappa.result<-new.kappa.result[rownames(new.kappa.result)%in%names(results2),
                                       colnames(new.kappa.result)%in%names(results2)]
    sil <- silhouette(results2, (1-new.kappa.result), diss=T)
  }
  pdf(paste(output_dir,"/Clustering_files/new_silhouette_plot.pdf",sep=""))
  plot(sil, nmax= 80, cex.names=0.6)
  dev.off()
  length(results$clustering)
  dim(d)
  pdf(paste(output_dir,"/Clustering_files/new_silhouette_width.pdf",sep=""))
  hist(sil[,3],breaks=30)
  dev.off() 
  
  singleton<-results$clustering[!names(results$clustering)%in%names(results2)]
  singleton[1:length(singleton)]<- max(results2)+1
  results2<-c(results2,singleton)
  #  k=20
  # results = results.pam
  # index.cluster = results2
  #  index.cluster = results2
  
  #-------------------------------------------
  #topGO figure for each cluster
  #-------------------------------------------
  #dir.create("./topGO")
  #GO_list = pathway[pathway.type == "G"]
  
  #for(i in 1:k){
  #  print(i)
  #  e = names(which(results2 == i))
  #  e = e[e %in% names(GO_list)] #only GO term, not including KEGG, BIOCARTA, REACTOME
  #  #    e = names(which(results2==i))
  #  #input: variable e is the list of GO terms we want to look at
  #  # GO_list <- GO_list
  #  gene.interest = unique(unlist(GO_list[e]))
  #  all_vec = rep(0.01, length(gene.interest)); names(all_vec) = gene.interest
  
  #  topGenes= function(allScore) return(allScore < 0.05)
  
  #  for(onto in c("BP", "MF", "CC")){
  #    print(onto)
  #    GOdata <- new("topGOdata",
  #                  description = "Meta GO Enrichment analysis", ontology = onto,
  #                  allGenes = all_vec, geneSel = topGenes,
  #                  # nodeSize = 10,annot = annFUN.org,
  #                  nodeSize = 1,annot = annFUN.org,
  #                  mapping = "org.Hs.eg.db",
  #                  ID = "alias")
  #    sum(!e %in% usedGO(GOdata)) 
  #    e1 = e[e %in% usedGO(GOdata)] # only remain the GO in that specific ontology
  #    if(length(e1)){
  #      gene.interest = unlist(GO_list[e1])
  
  #      #p-value
  #      score = rep(1, length(usedGO(GOdata))); names(score) = usedGO(GOdata)
  #      # score[names(score) %in% e] = 0.01
  #      score[names(score) %in% e1] = b$q_value_meta[match(names(score)[names(score) %in% e1], b$X)]
  
  #      pdf(paste("./topGO/(", onto, ")Cluster", k, "_", i, ".pdf", sep = ""))
  #      #pdf(paste("Cluster3-1.pdf", sep = ""))
  #      # showSigOfNodes(GOdata, score, firstSigNodes = 10, wantedNodes = e, useInfo = 'all', )
  #      showSigOfNodes(GOdata, score, firstSigNodes = sum(score < 1), wantedNodes = e1, useInfo = c("np"))
  #      # showSigOfNodes(GOdata, score, firstSigNodes = sum(score < 1), useInfo = c("np"))
  #      dev.off()
  #    }
  #  }
  #}
  
  #================================
  #heatmap for each cluster
  #================================
  indi.heatmap.data <- list()
  for (i in 1:k){
    # print(i)
    #length(which(results2 == 1))
    e = names(which(results2 == i))
    #    e = names(which(results2==i))
    #input: variable e is the list of GO terms we want to look at
    if(0){
      GO_list <- GO_list
      gene.interest = unlist(GO_list[e])
      all_vec = rep(0.01, length(gene.interest)); names(all_vec) = gene.interest
      
      topGenes= function (allScore) 
      {
        return(allScore < 0.05)
      }
      
      
      # GOdata <- new("topGOdata",
      #                   description = "Meta GO Enrichment analysis", ontology = "BP",
      #                   allGenes = all_vec, geneSel = topGenes,
      #                   nodeSize = 10,annot = annFUN.org,
      #                                     mapping = "org.Hs.eg.db",
      #                                     ID = "symbol")
      
      GOdata <- new("topGOdata",
                    description = "Meta GO Enrichment analysis", ontology = "BP",
                    allGenes = all_vec, geneSel = topGenes,
                    nodeSize = 10,annot = annFUN.org,
                    mapping = "org.Hs.eg.db",
                    ID = "alias")
      sum(!e %in% usedGO(GOdata))
      e = e[ e %in% usedGO(GOdata)] #this kick out a lot??
    }
    
    #     indi.heatmap.data[[i]] <- matrix(0, length(e),Num_of_gene_lists)
    #     for (j in 1:length(e)){
    #       for (n in 1:Num_of_gene_lists){
    #         indi.heatmap.data[[i]][j,n] <- b[which(b$X %in% e[j]),n+1]
    #       }
    #     }
    indi.heatmap.data[[i]] = as.matrix(b[match(e, rownames(b)), (1:Num_of_gene_lists) + 2])
    rownames(indi.heatmap.data[[i]]) <- e
    colnames(indi.heatmap.data[[i]]) <- paste("group",1:Num_of_gene_lists,sep="")
  }
  e = names(results2[order(results2)])
  if (software == "CPI"){indi.heatmap.data_all = as.matrix(b[match(e, rownames(b)), (1:Num_of_gene_lists) + 2])}
  if (software == "MAPE"){indi.heatmap.data_all = as.matrix(b[match(e, rownames(b)), (1:Num_of_gene_lists)])}
  rownames(indi.heatmap.data_all) <- e
  colnames(indi.heatmap.data_all) <- paste("group",1:Num_of_gene_lists,sep="")
  
  
  #dim(indi.heatmap.data[[2]])  
  
  #========================================================================================
  #heatmap for each cluster, based on kappa similarity and put number of overlapped genes
  #========================================================================================
  dir.create(paste(output_dir,"/Clustering_files/Heatmaps",sep=""))
  order.list = list()
  #  order.list.all
  last_max <-0
  for (i in 1:k){
    print(paste("Cluster", i, sep = ""))
    e = names(which(results2 == i))
    temp_file_name <- output_dir
    if( i!=k){
      temp_file_name <- paste(output_dir,"/Clustering_files/Heatmaps/Heatmap_each_cluster_kappa_number_overlap_genes_", i, ".pdf", sep = "")
    } else{
      temp_file_name <- paste(output_dir,"/Clustering_files/Heatmaps/Heatmap_each_cluster_kappa_number_overlap_genes_Singleton", ".pdf", sep = "")
    }
    pdf(temp_file_name ,height = if (length(e)<=40) 7 
        else if (length(e)<=100) 10 
        else if (length(e)<=150) 12.5
        else if (length(e)<=200) 15
        else 16, width = if (length(e)<=40) 7 
        else if (length(e)<=100) 10
        else if (length(e)<=150) 12.5
        else if (length(e)<=200) 15
        else 16)
    # e = e[rev(order.list[[i]])]
    kappa = kappa.result[e, e] #dissimilarity
    if(length(unique(as.vector(kappa))) != 1){
      overlap.matrix = matrix(0, nrow = length(e), ncol = length(e))
      for(j in 1:length(e)){
        for(l in 1:length(e)) overlap.matrix[j, l] = length(intersect(pathway[[e[j]]], pathway[[e[l]]]))
      }
      overlap.matrix = matrix(paste(round(kappa.result[e, e], 2), "\n(", overlap.matrix, ")", sep = ""), nrow = length(e))
      
      d = as.dist(1-kappa)
      hc = hclust(d); hc = as.dendrogram(hc)
      a = heatmap.2(kappa, col = intpalette(c("#FFFFD5FF", "#FF0000FF"), numcol = 30), 
                    notecex = if (length(e)==2) -4+1/(0.5*log10(length(e)))
                    else if (length(e)==3) -3+1/(0.5*log10(length(e)))
                    else if (length(e)==4) -2+1/(0.5*log10(length(e)))
                    else if (length(e)<=40) -1+1/(0.5*log10(length(e))) 
                    else if (length(e)<=100) -0.85+1/(0.5*log10(length(e))) 
                    else if (length(e)<=150) -0.78+1/(0.5*log10(length(e)))
                    else if (length(e)<=200) -0.75+1/(0.5*log10(length(e)))
                    else -0.73+1/(0.5*log10(length(e))), 
                    symbreaks = F, symkey = F, Rowv = hc, Colv = hc, notecol="black", 
                    cellnote = overlap.matrix, density.info="none", trace="none", scale = "none", 
                    margins = c(8, 9), srtCol = 45,
                    cexRow =  if (length(e)==2) 1.5
                    else if (length(e)==3) 1.48
                    else if (length(e)==4) 1.4
                    else if (length(e)<=40) 0.2 + 1/log10(length(e)) 
                    else if (length(e)<=100) 0.2 + 1/log10(length(e)) 
                    else if (length(e)<=150) 0.07 + 1/log10(length(e))
                    else if (length(e)<=200) 0.04 + 1/log10(length(e))
                    else 0.03 + 1/log10(length(e)),
                    cexCol =  if (length(e)==2) 1.4
                    else if (length(e)==3) 1.3
                    else if (length(e)==4) 1.2
                    else if (length(e)<=40) 0.1 + 1/log10(length(e)) 
                    else if (length(e)<=100) 0 + 1/log10(length(e)) 
                    else if (length(e)<=150) 0.01 + 1/log10(length(e))
                    else if (length(e)<=200)  1/log10(length(e))
                    else -0.1+1/log10(length(e)),
                    keysize = if (length(e)<=40) 1.5 
                    else if (length(e)<=100) 1.3 
                    else if (length(e)<=150) 1
                    else 0.9 ,
                    offsetRow = if (length(e)<=40) 0.5 
                    else if (length(e)<=100) 0 
                    else if (length(e)<=150) -0.1
                    else -0.2,
                    offsetCol = if (length(e)<=40) 0.5 
                    else if (length(e)<=100) 0 
                    else if (length(e)<=150) -0.1
                    else -0.2)
      # a = heatmap.2(kappa, notecex = 0.2, Rowv = T, Colv = T, notecol="black", cellnote = matrix("Hello\n test", 33, 33), density.info="none", trace="none", scale = "none", margins = c(8, 9), cexRow = if(length(e) == 3) -0.5 + 1/log10(length(e)) else 0.2 + 1/log10(length(e)), cexCol = if(length(e) == 3) -0.5 + 1/log10(length(e)) else 0.2 + 1/log10(length(e)))
      order.list[[i]] = a$rowInd
      if (i==1){
        order.list.all <- rev(a$rowInd)
      } else {
        order.list.all <- c(order.list.all, rev(a$rowInd)+last_max)
      }
      last_max <- max(a$rowInd)+last_max
      #        print(last_max)
      dev.off()
    } else {
      order.list[[i]] = 1:length(e)
      if (i==1){
        order.list.all <- rev(1:length(e))
      } else{
        order.list.all <- c(order.list.all, rev(1:length(e)))
      }
    }
  }   
  #===========================================
  #aggregate across (GO or other DBs) terms in each cluster
  #===========================================
  #xx <- as.list(GOTERM)
  # BP.GO = names(annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol"))
  # MF.GO = names(annFUN.org("MF", mapping = "org.Hs.eg.db", ID = "symbol"))
  # CC.GO = names(annFUN.org("CC", mapping = "org.Hs.eg.db", ID = "symbol"))
  # GO2ontology = inverseList(list(BP = BP.GO, MF = MF.GO, CC = CC.GO))
  
  
  GO.cluster.20 = lapply(1:k, function(i){
      print(i) #i=1
      e = names(which(results2 == i))
      #e.type = pathway.type[match(e, names(pathway))]
      #    e = names(which(results$clustering==i))
      #input: variable e is the list of GO terms we want to look at
      if(0){
        GO_list <- GO_list
        gene.interest = unlist(GO_list[e])
        all_vec = rep(0.01, length(gene.interest)); names(all_vec) = gene.interest
        
        topGenes= function (allScore) 
        {
          return(allScore < 0.05)
        }
        
        
        # GOdata <- new("topGOdata",
        #                   description = "Meta GO Enrichment analysis", ontology = "BP",
        #                   allGenes = all_vec, geneSel = topGenes,
        #                   nodeSize = 10,annot = annFUN.org,
        #                                     mapping = "org.Hs.eg.db",
        #                                     ID = "symbol")
        
        GOdata <- new("topGOdata",
                      description = "Meta GO Enrichment analysis", ontology = "BP",
                      allGenes = all_vec, geneSel = topGenes,
                      nodeSize = 10,annot = annFUN.org,
                      mapping = "org.Hs.eg.db",
                      ID = "alias")
        sum(!e %in% usedGO(GOdata))
        
        e = e[ e %in% usedGO(GOdata)] #??
      }
      
      if(length(e)!=0){
        #  term = e
        #  GO.term = e[which(e%in%names(GOgenesets))]
        #  #GO.term = sapply(e[e.type == "G"], function(ee) Term(xx[[ee]]))
        #  term[e.type == "G"] = GO.term
        term = rep("", length(e))  
        #  onto[e.type == "G"] = Ontology(e[e.type == "G"])
        #  # f = sapply(e, function(ee) Term(xx[[ee]]))
        #  # print(sum(! e %in% names(GO2ontology)))
        #  # print(sum(is.na(Ontology(e))))
        #  # f = cbind(e, unlist(GO2ontology[e]), f)#; rownames(f) = NULL; colnames(f) = c("ID", "Term")
        f = cbind(e,term)#; rownames(f) = NULL; colnames(f) = c("ID", "Term")
        
        #   for (j in 1: length(e)){
        #       NumGeneTotal[j] <- length(GO_list[[e[j]]])
        #       for (n in 1:Num_of_gene_lists){ #j=1;n=3
        #         genelist_tmp <- list()
        #         genelist_tmp <- genelists[[n]][which((genelists[[n]][,3]==1) & (genelists[[n]][,1] %in% GO_list[[e[j]]])),1]
        #         GeneInSet[j,n] <- paste(genelist_tmp,collapse="//")
        #         NumGeneInSet[[n]][j] <- length(GeneInSet[[n]][j])
        #         if(length(GeneInSet[[n]][j]))
        #           GeneInSet[[n]][j] = 0
        #       } #make table clear
        #       #Gene <- GO_list[[e[j]]]
        #     }
        if (software == "MAPE"){
          NumGeneTotalInSet <- vector()
          for (j in 1: length(e)){
            NumGeneTotalInSet[j] <- length(pathway[[e[j]]])
          }
          m = summary[e,method]
          h = cbind(f, NumGeneTotalInSet,m)
        }
        else if (software == "CPI"){
          if (enrichment == "KS"){
            NumGeneTotalInSet <- vector()
            for (j in 1: length(e)){
              NumGeneTotalInSet[j] <- length(pathway[[e[j]]])
            }
            m = summary[e,3:ncol(summary)]
            h = cbind(f, NumGeneTotalInSet,m)
          }
          else if (enrichment == "Fisher's exact"){
            NumGeneTotalInSet <- vector()
            DEGeneInSet <- matrix("NA",length(e),Num_of_gene_lists)
            NumDEGeneInSet <- matrix(0,length(e),Num_of_gene_lists)
            for (j in 1: length(e)){
              NumGeneTotalInSet[j] <- length(pathway[[e[j]]])
              for (n in 1:Num_of_gene_lists){
                genelist_tmp <- list()
                genelist_tmp <- genelist[[n]][which(toupper(genelist[[n]]) %in% toupper(pathway[[e[j]]]))]
                DEGeneInSet[j,n] <- paste(genelist_tmp,collapse="//")
                NumDEGeneInSet[j,n] <- length(genelist_tmp)
              }
            } 
            h = cbind(f, NumGeneTotalInSet,NumDEGeneInSet,DEGeneInSet)
            rownames(h) <-NULL
            colnames(h) <-NULL
          }
        }
        #     for (n in 2:Num_of_gene_lists){
        #       h = cbind(h, NumGeneInSet[[n]]) 
        #     } 
        #     h2 = cbind(h,paste(GeneInSet[[1]],collapse="//"))
        #     for (n in 2:Num_of_gene_lists){
        #       h2 = cbind(h2, paste(GeneInSet[[n]],collapse="//")) 
        #     }
        #    colnames(f) = c("ID", "Term","","NumGeneInSet_1", "")
        h = h[rev(order.list[[i]]),]
        
        #f
      }
    })
  #   })  
  #   #### end of module 2
  #   #save.image()
  #   ##module 3 heatmap
  #   #aggregate across GO term
  #   Cluster.p = t(sapply(1:k, function(i){
  #     go.id = GO.cluster.20[[i]][,1]
  #     #colnames(b)[(1:Num_of_gene_lists)+2]
  # #    p.matrix = b[b$X %in% go.id, (1:Num_of_gene_lists)+2]
  #     p.matrix = b[rownames(b) %in% go.id, (1:Num_of_gene_lists)+2]
  #     meta.p = apply(p.matrix, 2, function(p) {
  #       Fisher = -2*sum(log(p))
  #       p.fisher = pchisq(Fisher, df = 2*length(p), lower.tail = F); p.fisher
  #       # p.fisher = 1-pchisq(Fisher, df = 2*length(p), lower.tail = T); p.fisher
  #     })
  #     names(meta.p) = paste("Study", 1:Num_of_gene_lists, sep = "")  
  #     return(meta.p)
  #   }))
  #   
  #   Cluster.p[Cluster.p < 10^-10] = 10^-10
  
  
  #========================================================================================
  #heatmap for all cluster, based on kappa similarity and put number of overlapped genes
  #========================================================================================  
  #  order.list_all = list()
  print(paste("Cluster",sep = ""))
  e = names(results2[order(results2)])
  pdf(paste(output_dir,"/Heatmap_each_cluster_kappa_number_overlap_genes_all",".pdf", sep = ""),height = if (length(e)<=40) 7 
      else if (length(e)<=100) 10 
      else if (length(e)<=150) 12.5
      else if (length(e)<=200) 15
      else 16, width = if (length(e)<=40) 7 
      else if (length(e)<=100) 11
      else if (length(e)<=150) 13.5
      else if (length(e)<=200) 16
      else 17)
  # e = e[rev(order.list[[i]])]
  e = e[order.list.all]
  kappa = kappa.result[e, e] #dissimilarity
  colorbar <- as.character( results2[order(results2)][order.list.all]+1 )  
  
  if(length(unique(as.vector(kappa))) != 1){
    overlap.matrix = matrix(0, nrow = length(e), ncol = length(e))
    for(j in 1:length(e)){
      for(l in 1:length(e)) overlap.matrix[j, l] = length(intersect(pathway[[e[j]]], pathway[[e[l]]]))
    }
    overlap.matrix = matrix(paste(round(kappa.result[e, e], 2), "\n(", overlap.matrix, ")", sep = ""), nrow = length(e))
    
    d = as.dist(1-kappa)
    ##    hc = hclust(d); hc = as.dendrogram(hc) we do not do it here.
    a = heatmap.2(kappa, col = intpalette(c("#FFFFD5FF", "#FF0000FF"), numcol = 30), 
                  notecex = if (length(e)==2) -4+1/(0.5*log10(length(e)))
                  else if (length(e)==3) -3+1/(0.5*log10(length(e)))
                  else if (length(e)==4) -2+1/(0.5*log10(length(e)))
                  else if (length(e)<=40) -1+1/(0.5*log10(length(e))) 
                  else if (length(e)<=100) -0.85+1/(0.5*log10(length(e))) 
                  else if (length(e)<=150) -0.78+1/(0.5*log10(length(e)))
                  else if (length(e)<=200) -0.75+1/(0.5*log10(length(e)))
                  else -0.73+1/(0.5*log10(length(e))), 
                  symbreaks = F, symkey = F, Rowv = F, Colv = F, notecol="black", dendrogram="none", 
                  #                  symbreaks = F, symkey = F, Rowv = hc, Colv = hc, notecol="black",                   
                  cellnote = overlap.matrix, density.info="none", trace="none", scale = "none", 
                  margins = c(10, 15), srtCol = 45,
                  RowSideColors = colorbar,
                  ColSideColors = colorbar,
                  cexRow =  if (length(e)==2) 1.5
                  else if (length(e)==3) 1.48
                  else if (length(e)==4) 1.4
                  else if (length(e)<=40) 0.2 + 1/log10(length(e)) 
                  else if (length(e)<=100) 0.2 + 1/log10(length(e)) 
                  else if (length(e)<=150) 0.07 + 1/log10(length(e))
                  else if (length(e)<=200) 0.04 + 1/log10(length(e))
                  else 0.03 + 1/log10(length(e)),
                  cexCol =  if (length(e)==2) 1.4
                  else if (length(e)==3) 1.3
                  else if (length(e)==4) 1.2
                  else if (length(e)<=40) 0.1 + 1/log10(length(e)) 
                  else if (length(e)<=100) 0 + 1/log10(length(e)) 
                  else if (length(e)<=150) 0.01 + 1/log10(length(e))
                  else if (length(e)<=200)  1/log10(length(e))
                  else -0.1+1/log10(length(e)),
                  keysize = if (length(e)<=40) 1.5 
                  else if (length(e)<=100) 1.3 
                  else if (length(e)<=150) 1
                  else 0.9 ,
                  offsetRow = if (length(e)<=40) 0.5 
                  else if (length(e)<=100) 0 
                  else if (length(e)<=150) -0.1
                  else -0.2,
                  offsetCol = if (length(e)<=40) 0.5 
                  else if (length(e)<=100) 0 
                  else if (length(e)<=150) -0.1
                  else -0.2)
    # a = heatmap.2(kappa, notecex = 0.2, Rowv = T, Colv = T, notecol="black", cellnote = matrix("Hello\n test", 33, 33), density.info="none", trace="none", scale = "none", margins = c(8, 9), cexRow = if(length(e) == 3) -0.5 + 1/log10(length(e)) else 0.2 + 1/log10(length(e)), cexCol = if(length(e) == 3) -0.5 + 1/log10(length(e)) else 0.2 + 1/log10(length(e)))
    #    order.list_all[[1]] = a$rowInd
    dev.off()
  } # else order.list_all[[1]] = 1:length(e)
  
  
  
  #================================
  #heatmap for all clusters
  #================================
  
  print(paste("Cluster", sep = ""))
  e = names(results2[order(results2)])
  e = e[order.list.all]
  indi.heatmap.data_all <- indi.heatmap.data_all[order.list.all,]
  #  colorbar<-character()
  #  colorbar[which(results2[order(results2)]%%2==0)]<-"blue"
  #  colorbar[which(results2[order(results2)]%%2==1)]<-"gray"
  
  colorbar <- as.character( results2[order(results2)][order.list.all]+1 )  
  
  pdf(paste(output_dir,"/Heatmap_clusters_all",".pdf",sep=""),
      height = if (length(e)<=40) 8.5 
      else if (length(e)<=100) 10 
      else if (length(e)<=150) 12.5
      else if (length(e)<=200) 15
      else 16, 
      width = if (length(e)<=40) 9.5 
      else if (length(e)<=100) 11 
      else if (length(e)<=150) 13.5
      else if (length(e)<=200) 16
      else 17)
  if(nrow(indi.heatmap.data_all)==1){
    image(t(-log10(indi.heatmap.data_all)))
  }else if (nrow(indi.heatmap.data_all)!=0){
    heatmap.pathway.pvalue<-indi.heatmap.data_all
    heatmap.pathway.pvalue[heatmap.pathway.pvalue < 10^-10] = 10^-10
    a = heatmap.2(-log10(heatmap.pathway.pvalue), col = (c(intpalette(c("#FFFFD5FF", "#FF0000FF"), numcol = 30),rep("#FF0000FF",30))),
                  Rowv = F, Colv = F, density.info="none", trace="none", scale = "none", margins = c(7, 20), 
                  cexRow =  if (length(e)==2) 1.5
                  else if (length(e)==3) 1.48
                  else if (length(e)==4) 1.4
                  else if (length(e)<=40) 0.2 + 1/log10(length(e)) 
                  else if (length(e)<=100) 0.2 + 1/log10(length(e)) 
                  else if (length(e)<=150) 0.07 + 1/log10(length(e))
                  else if (length(e)<=200) 0.03 + 1/log10(length(e))
                  else 0.01 + 1/log10(length(e)),
                  cexCol =  if (length(e)<=40) 0.2 + 1/log10(Num_of_gene_lists) 
                  else if (length(e)<=100) 0.3 + 1/log10(Num_of_gene_lists) 
                  else if (length(e)<=150) 0.5 + 1/log10(Num_of_gene_lists)
                  else 0.6 + 1/log10(Num_of_gene_lists),
                  RowSideColors = colorbar,
                  keysize = if (length(e)<=40) 1.5 
                  else if (length(e)<=100) 1.3 
                  else if (length(e)<=150) 1
                  else 0.9 ,
                  offsetRow = if (length(e)<=40) 0.25 
                  else if (length(e)<=100) 0 
                  else if (length(e)<=150) -0.1
                  else -0.2)
  }
  dev.off()
  
  
  #====================================================================================================
  #heatmap for each cluster, based on p-value, but GO terms are ordered by kappa cluster (above steps)
  #====================================================================================================
  #   
  #   for (i in 1:k){
  #     print(paste("Cluster", i, sep = ""))
  #     e = names(which(results2 == i))
  #     indi.heatmap.data[[i]] = indi.heatmap.data[[i]][rev(order.list[[i]]),]
  # #    colorbar<-character()
  # #    colorbar[which(results2[order(results2)]%%2==0)]<-"blue"
  # #    colorbar[which(results2[order(results2)]%%2==1)]<-"gray"
  #     pdf(paste("./Heatmap_each_cluster/Heatmap_each_cluster_p_value_",i,".pdf",sep=""),
  #         height = if (length(e)<=40) 7 
  #         else if (length(e)<=100) 10 
  #         else if (length(e)<=150) 12.5
  #         else if (length(e)<=200) 15
  #         else 16, 
  #         width = if (length(e)<=40) 7 
  #         else if (length(e)<=100) 9 
  #         else if (length(e)<=150) 11
  #         else 13)
  #     if(nrow(indi.heatmap.data[[i]])==1){
  #       image(t(-log10(indi.heatmap.data[[i]])))
  #     }else if (nrow(indi.heatmap.data[[i]])!=0){
  #       a = heatmap.2(-log10(indi.heatmap.data[[i]]), col = intpalette(c("#FFFFD5FF", "#FF0000FF"), numcol = 30),Rowv = F, Colv = F, density.info="none", trace="none", scale = "none", margins = c(8, 9), 
  #                     cexRow =  if (length(e)==2) 1.5
  #                     else if (length(e)==3) 1.48
  #                     else if (length(e)==4) 1.4
  #                     else if (length(e)<=40) 0.2 + 1/log10(length(e)) 
  #                     else if (length(e)<=100) 0.2 + 1/log10(length(e)) 
  #                     else if (length(e)<=150) 0.07 + 1/log10(length(e))
  #                     else if (length(e)<=200) 0.03 + 1/log10(length(e))
  #                     else 0.01 + 1/log10(length(e)),
  #                     cexCol =  if (length(e)<=40) 0.2 + 1/log10(Num_of_gene_lists) 
  #                     else if (length(e)<=100) 0.3 + 1/log10(Num_of_gene_lists) 
  #                     else if (length(e)<=150) 0.5 + 1/log10(Num_of_gene_lists)
  #                     else 0.6 + 1/log10(Num_of_gene_lists),
  #  #                   RowSideColors = colorbar,
  #                     keysize = if (length(e)<=40) 1.5 
  #                     else if (length(e)<=100) 1.3 
  #                     else if (length(e)<=150) 1
  #                     else 0.9 ,
  #                     offsetRow = if (length(e)<=40) 0.5 
  #                     else if (length(e)<=100) 0 
  #                     else if (length(e)<=150) -0.1
  #                     else -0.2)
  #     }
  #     dev.off()
  #   }
  #   
  #   
  
  #=============================
  ##### Text Mining
  #=============================
  data(textmtx)
  new.fmat <- descmat
  new.mat <- namemat
  pathways<-c(Hallmark.genesets,Positional.genesets,CGP.genesets,Canonical.genesets,
              Biocarta.genesets,KEGG.genesets,Reactome.genesets,
              MIR.genesets,TFT.genesets,CGN.genesets,
              CM.genesets,GOBP.genesets,GOCC.genesets,
              GOMF.genesets,Oncogenic.genesets,Immunologic.genesets,
              PID.genesets,CMAP_up.genesets,CMAP_down.genesets,
              PPI.genesets,JASPARhuman.genesets,PITA.genesets,
              miRanda.genesets,TargetScan.genesets,Phenocarta.genesets)
  pathway<-c(Biocarta.genesets,KEGG.genesets,Reactome.genesets,
             GOBP.genesets,GOCC.genesets,
             GOMF.genesets,Phenocarta.genesets)
  new.fmat<-new.fmat[,which(names(pathways)%in%names(pathway))]
  new.mat<-new.mat[,which(names(pathways)%in%names(pathway))]
  
  alpha <- 0.01
  
  Xj <- rep(0,dim(new.fmat)[2]) #The penalization on word count of each pathway
  for (j in 1:length(Xj)){
    Xj[j] <- exp(-alpha * sum(new.fmat[,j]))
  }
  
#  dir.create(paste(output_dir,"/Clustering_files/text_mining",sep=""))
  #names(Reactomegenesets)
  #which(names(pathway)=="GOBP PROTEIN TETRAMERIZATION")
  #which(names(pathway)=="BIOCARTA TID PATHWAY")
  #which(names(pathway)=="GOBP NEGATIVE REGULATION OF PROGRAMMED CELL DEATH")
  #which(names(pathway)=="GOBP ANTI APOPTOSIS")
  #which(names(pathway)=="GOBP NEGATIVE REGULATION OF APOPTOSIS")
  #which(names(pathway)=="GOBP NEGATIVE REGULATION OF DEVELOPMENTAL PROCESS")
  #for (i in 1:(length(Reactomegenesets))){
  #  if(grepl(toupper("TCA"), toupper(names(Reactomegenesets[i])))){
  #    print(i)
  #  }
  #}
  #names(Reactomegenesets[6])
  #names(Reactomegenesets[37])
  #names(Reactomegenesets[321])
  
  #which(names(pathway)==toupper("Reactome Pyruvate metabolism and Citric Acid (TCA) cycle"))
  #row.names(all.mat)[3193]
  #toupper(cluster2_all)
  #which(toupper(names(pathway))==toupper(cluster2_all[8]))
  #pathway[2178]
  tm_all <- list()
  for (tmk in 1:k){
    w.size <- rep(0,dim(new.fmat)[2]) # |w|, size of word for the full description for each pathways
    for (i in 1:dim(new.fmat)[2]){
      w.size[i] <- sum(new.fmat[,i])
    }
    cluster1 <- rep(0,length(GO.cluster.20[[tmk]][,1]))
    length(cluster1)
    for (m in 1:(length(cluster1))){
      cluster1[m] <- which(toupper(names(pathway))==toupper(GO.cluster.20[[tmk]][m,1]))
    }
    #hist(w.size, breaks=200)
    #max(w.size)
    #sum(w.size)/dim(new.fmat)[2]
    #median(w.size)
    #hist(term.size, breaks=200)
    #  max(term.size)
    #  median(term.size)
    
    #    all.size <- data.frame(c(term.size,w.size),c(rep("term",dim(new.fmat)[2]),rep("full_description",dim(new.fmat)[2])))
    #colnames(all.size) <- c("word_counts", "type")
    #ggplot(all.size, aes(x=word_counts, fill=type)) + geom_histogram(binwidth=1, alpha=.5, position="identity")
    
    # exp(-400*0.01)
    #  sort(w.size)[round(0.90*length(w.size))]
    
    #which(w.size==1)
    #which(new.fmat[,1743] == 1)
    #which(w.size > 400)
    #which(new.fmat[,1919] == 1)
    
    #exp(-100*0.01)
    allwords_idx <- 0
    for (i in 1:length(cluster1)){
      for (j in 1:dim(new.fmat)[1]){
        if( new.fmat[j,cluster1[i]] ==1 ){
          allwords_idx <- c(allwords_idx,j)
        }
        if( new.mat[j,cluster1[i]] ==1 ){
          allwords_idx <- c(allwords_idx,j)
        }
      }
    }
    if(sum(allwords_idx)!=0){
      allwords_idx <- allwords_idx[2:length(allwords_idx)]
      allwords_idx <- unique(allwords_idx)
      allwords <- row.names(new.mat)[allwords_idx]
      allwords_count <- rep(0, length(allwords_idx))
      Tw <- rep(0, length(allwords_idx))
      ii <- 1
      for (i in allwords_idx){
        for (j in cluster1){
          if ( new.fmat[i,j] == 1 ){
            Tw[ii] <- Tw[ii] + max(new.mat[i,j], Xj[j])
            allwords_count[ii] <- allwords_count[ii] + 1
          } else {
            Tw[ii] <- Tw[ii] + new.mat[i,j]
            allwords_count[ii] <- allwords_count[ii] + new.mat[i,j]
          }
        }
        ii <- ii+1
      }
      
      
      TM_perm<- 1000
      subcluster1 <- matrix(0,TM_perm,length(cluster1))
      for (i in 1:TM_perm){
        subcluster1[i,] <- sample(dim(new.fmat)[2], length(cluster1))
      }
      
      allp <- rep(1, length(allwords_idx))
      allq <- rep(1, length(allwords_idx))
      ii <- 1
      for (i in allwords_idx){
        TEw <- rep(0,TM_perm)
        for (tmperm in 1:TM_perm){
          for (j in subcluster1[tmperm,]){
            if (new.fmat[i,j] == 1){
              TEw[tmperm] <- TEw[tmperm] + max(new.mat[i,j], Xj[j])
            } else {
              TEw[tmperm] <- TEw[tmperm] + new.mat[i,j]
            }
          }
          #        if ( b%%1000 == 0){
          #          print(b)
          #        }
        }
        allp[ii] <- 1-(rank(c(Tw[ii],TEw), ties.method = "min")[1])/(length(TEw)+2)
        ii <- ii+1
      }
      
      allp[which(allp>1)]=1
      allq <- p.adjust(allp, "BH")
      #    row.names(new.fmat)[allwords_idx[which(allq<0.05)]]
      #    row.names(new.fmat)[allwords_idx[order(allq)]][1:10]
      
      #    row.names(new.fmat)[allwords_idx[order(allp)]][1:40]
      #    sort(allp)
      #    sort(allp)[1:30]
      #View(new.mat[220:230,1:3])
      
      ##
      ##Test 2: Fisher exact
      fisherp <- rep(1,length(allwords_idx))
      t=1
      for (i in allwords_idx){
        fishertable <- matrix(0,2,2)
        for (j in cluster1){
          if( (new.fmat[i,j] ==1) | (new.mat[i,j] ==1)){
            fishertable[1,1] = fishertable[1,1] + 1
          }
        }
        for (j in (setdiff((1:dim(new.mat)[2]),cluster1))){
          if( (new.fmat[i,j] ==1) | (new.mat[i,j] ==1)){
            fishertable[1,2] = fishertable[1,2] + 1
          }
        }
        fishertable[2,1] = length(cluster1)-fishertable[1,1]
        fishertable[2,2] = dim(new.mat)[2] - fishertable[1,1]-fishertable[1,2] - fishertable[2,1]
        fisherp[t] <- stats::fisher.test(fishertable,alternative='g')$p.value
        t = t+1
      }
      #    fisherp
      fisherq <- p.adjust(fisherp,"BH")
      #    row.names(new.fmat)[allwords_idx[order(fisherp)]][1:100]
      #    fisherq[order(fisherp)]
      #    sort(fisherq)
      #    allq[order(allp)]
      
      #   rank(allp,ties.method="min")
      tm_all[[tmk]] <- cbind(allwords,allp,allq,fisherp,fisherq,rank(allq,ties.method="min"),rank(fisherq,ties.method="min"),allwords_count)
      colnames(tm_all[[tmk]]) <- c("words","perm_p","perm_q","Fisher_p","Fisher_q","perm_rank","Fisher_rank","count")
      
      
#      pdf(paste(output_dir,"/Secondary_files/text_mining/compare_two_tests",tmk,".pdf",sep=""))
#      permp <- allp
#      plot(-log(fisherp),-log(permp))
#      dev.off()
    } else{
      tm_all[[tmk]] <- cbind(NA,NA,NA,NA,NA,NA,NA,NA)
      colnames(tm_all[[tmk]]) <- c("words","perm_p","perm_q","Fisher_p","Fisher_q","perm_rank","Fisher_rank","count")
    }
  }
  # End of Text Mining
  tm_filtered <- list() #filter out count 1, q 0.05
  for (i in 1:k){ 
    tm_filtered[[i]] <- tm_all[[i]][ which( (as.numeric(tm_all[[i]][,3]) < 0.05) & (tm_all[[i]][,8]>1)), ]
  }
  #tm_all[[3]]
  
  setwd(output_dir)
  
  topGO.summary = GO.cluster.20
  
  if(is.null(dim(tm_filtered[[1]]))==TRUE){
    cat("Cluster 1\n", file = "Clustering_Summary.csv",append=T)
    cat("Key words,", file = "Clustering_Summary.csv",append=T)
    write.table(tm_filtered[[1]][1], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = "Clustering_Summary.csv",append=T)
    write.table(tm_filtered[[1]][3], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("count,", file = "Clustering_Summary.csv",append=T)
    write.table(tm_filtered[[1]][8], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    write.table(topGO.summary[[1]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)    
  } else {
    cat("Cluster 1\n", file = "Clustering_Summary.csv",append=T)
    cat("Key words,", file = "Clustering_Summary.csv",append=T)
    write.table(t(tm_filtered[[1]][(order(tm_filtered[[1]][,3])[1:15]),1]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = "Clustering_Summary.csv",append=T)
    write.table(t(tm_filtered[[1]][(order(tm_filtered[[1]][,3])[1:15]),3]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("count,", file = "Clustering_Summary.csv",append=T)
    write.table(t(tm_filtered[[1]][(order(tm_filtered[[1]][,3])[1:15]),8]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    write.table(topGO.summary[[1]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
  }
  for (i in 2:(k-1)){ 
    if(is.null(dim(tm_filtered[[i]]))==TRUE){
      cat(paste("\nCluster ", i, "\n", sep = ""), file = "Clustering_Summary.csv", append = T)
      cat("Key words,", file = "Clustering_Summary.csv",append=T)
      write.table(tm_filtered[[i]][1], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("q_value,", file = "Clustering_Summary.csv",append=T)
      write.table(tm_filtered[[i]][3], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = "Clustering_Summary.csv",append=T)
      write.table(tm_filtered[[i]][8], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(topGO.summary[[i]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)    
    } else {
      cat(paste("\nCluster ", i, "\n", sep = ""), file = "Clustering_Summary.csv", append = T)
      cat("Key words,", file = "Clustering_Summary.csv",append=T)
      write.table(t(tm_filtered[[i]][(order(tm_filtered[[i]][,3])[1:15]),1]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("q_value,", file = "Clustering_Summary.csv",append=T)
      write.table(t(tm_filtered[[i]][(order(tm_filtered[[i]][,3])[1:15]),3]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = "Clustering_Summary.csv",append=T)
      write.table(t(tm_filtered[[i]][(order(tm_filtered[[i]][,3])[1:15]),8]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(topGO.summary[[i]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
    }
  }
  i=k
  cat(paste("\nSingleton Cluster", "\n", sep = ""), file = "Clustering_Summary.csv", append = T)
  #  cat("Key words,", file = "Clustering_Summary.csv",append=T)
  #  write.table(t(tm_filtered[[i]][(order(tm_filtered[[i]][,3])[1:15]),1]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F)
  #  cat("q_value,", file = "Clustering_Summary.csv",append=T)
  #  write.table(t(tm_filtered[[i]][(order(tm_filtered[[i]][,3])[1:15]),3]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F)
  #  cat("count,", file = "Clustering_Summary.csv",append=T)
  #  write.table(t(tm_filtered[[i]][(order(tm_filtered[[i]][,3])[1:15]),8]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F)
  write.table(topGO.summary[[i]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
}
