##' This is the third major function in the MetaPath2.0 package which outputs the pathway
##' cluster results with text mining anotation.

##' @title Output pathway cluster results
##' @param Num_Clusters number of clusters
##' @return The pathway cluster results in a csv file.
##' @authors Zhou Fang, Xiangrui Zeng, Wei Zong and George Tseng.
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
##' CPI.kappa_result = MAPE.Kappa(summary = CPI_result$summary,pathway = CPI_result$pathway,
##'                               max_k = 15, q_cutoff = 0.0005,software = CPI_result$method)
##' MAPE.Clustering(MAPE.Clustering(summary=CPI_result$summary,Num_Clusters = 8,
##'                                 kappa.result = CPI.kappa_result$kappa,sil_cut=0.1,
##'                                 Num_of_gene_lists=CPI_result$Num_of_gene_lists,genelist =CPI_result$genelist,
##'                                 pathway=CPI_result$pathway, enrichment=CPI_result$enrichment,
##'                                 method=CPI.kappa_result$method,software=CPI_result$method,
##'                                 n.text.permute = 10000)

MAPE.Clustering <- function(summary,Num_Clusters = 3, kappa.result = kappa.result, Num_of_gene_lists,
                            genelist = NULL, pathway, enrichment,method,software, sil_cut = 0.1,
                            n.text.permute = 1000, output_dir = getwd(),parallel=FALSE)
  {
  wd = getwd()
  k = Num_Clusters
  d = as.dist(1-kappa.result)
  results <- pam(d, k, diss=T)#K medoids

  b = summary
  dir.create(paste(output_dir,"/Clustering_files","_",Num_Clusters,"_","clusters",sep=""))
  output_clustering <- paste(output_dir,"/Clustering_files","_",Num_Clusters,"_","clusters",sep="")
  sil <- silhouette(results$clustering, (1-kappa.result), diss=T)

  k <- k+1
  k_org = k
  sil_cut <- sil_cut
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
    }#rename cluster index, so it is integer from 1 to k

    new.kappa.result<-new.kappa.result[rownames(new.kappa.result)%in%names(results2),
                                       colnames(new.kappa.result)%in%names(results2)]
    sil <- silhouette(results2, (1-new.kappa.result), diss=T)#recalculate silhoutte
  }

  singleton<-results$clustering[!names(results$clustering)%in%names(results2)]#singleton pathways
  singleton[1:length(singleton)]<- max(results2)+1#label singleton as the number above total number of clusters
  results2<-c(results2,singleton)
  results2 <- results2[order(results2)]
#  if (!all(1:k %in% unique(results$clustering))){
#    stop('Please choose a smaller number of clusters k or choose a larger FDR cutoff to select more pathways for clustering')
#  }
  k = max(results2)#k is #of clusters(except singleton)+1
  if (k != k_org){warning(paste('Due to empty cluster(s), the cluster number is reduced from',
                               k_org - 1, 'to', k-1))}#when >=one cluster being removed
  pdf(paste(output_clustering,"/silhouette_plot.pdf",sep=""))
  plot(sil, nmax= 80, cex.names=0.6)
  dev.off()
  pdf(paste(output_clustering,"/silhouette_width.pdf",sep=""))
  hist(sil[,3],breaks=30)
  dev.off()




  #================================
  #heatmap for each cluster
  #================================
  indi.heatmap.data <- list()
  for (i in 1:k){
    e = names(which(results2 == i))
    indi.heatmap.data[[i]] = as.matrix(b[match(e, rownames(b)), (1:Num_of_gene_lists) + 2])
    rownames(indi.heatmap.data[[i]]) <- e
    if (software == "CPI"){
    colnames(indi.heatmap.data[[i]]) <- colnames(summary[,-c(1,2)])
    }
    if (software == "MAPE"){
    colnames(indi.heatmap.data[[i]]) <- colnames(summary[,c(1:(ncol(summary)-6))])
    }
  }
  e = names(results2)
  if (software == "CPI"){indi.heatmap.data_all = as.matrix(b[match(e, rownames(b)), (1:Num_of_gene_lists) + 2])}
  if (software == "MAPE"){indi.heatmap.data_all = as.matrix(b[match(e, rownames(b)), (1:Num_of_gene_lists)])}
  rownames(indi.heatmap.data_all) <- e
  if (software == "CPI"){
    colnames(indi.heatmap.data_all) <- colnames(summary[,-c(1,2)])
  }
  if (software == "MAPE"){
    colnames(indi.heatmap.data_all) <- colnames(summary[,c(1:(ncol(summary)-6))])
  }
  #========================================================================================
  #heatmap for each cluster, based on kappa similarity and put number of overlapped genes
  #========================================================================================
  cat("Plotting Heatmaps...\n")
  dir.create(paste(output_clustering,"/Heatmaps",sep=""))
  for (i in 1:k){
    e = names(which(results2 == i))
    temp_file_name <- output_dir
    if( i!=k){
      temp_file_name <- paste(output_clustering,"/Heatmaps/Heatmap_each_cluster_kappa_number_overlap_genes_", i, ".pdf", sep = "")
    } else{
      temp_file_name <- paste(output_clustering,"/Heatmaps/Heatmap_each_cluster_kappa_number_overlap_genes_Singleton", ".pdf", sep = "")
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
    kappa = kappa.result[e, e]
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
      dev.off()
    }
  }
  #===========================================
  #aggregate across (GO or other DBs) terms in each cluster
  #===========================================

  cat("Performing pathway clustering...\n")
  GO.cluster.20 = lapply(1:k, function(i){
    e = names(which(results2 == i))
    if(length(e)!=0){
      term = rep("", length(e))
      f = cbind(e,term)
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
        }
      }
          rownames(h) <-NULL
          colnames(h) <-NULL
        }
      h
  })


  #========================================================================================
  #heatmap for all cluster, based on kappa similarity and put number of overlapped genes
  #========================================================================================
  e = names(results2[order(results2)])
  pdf(paste(output_clustering,"/Heatmap_each_cluster_kappa_number_overlap_genes_all",".pdf", sep = ""),height = if (length(e)<=40) 7
      else if (length(e)<=100) 10
      else if (length(e)<=150) 12.5
      else if (length(e)<=200) 15
      else 16, width = if (length(e)<=40) 7
      else if (length(e)<=100) 11
      else if (length(e)<=150) 13.5
      else if (length(e)<=200) 16
      else 17)
  kappa = kappa.result[e, e] #dissimilarity
  colorIndex = results2[order(results2)]

  set.seed(1)
  nCol = max(results2)-1
  colOption = c(sample(rainbow(nCol),nCol),"black")
  colorbar = mapvalues(colorIndex,from = 1:(nCol+1),to = colOption)

  if(length(unique(as.vector(kappa))) != 1){
    overlap.matrix = matrix(0, nrow = length(e), ncol = length(e))
    for(j in 1:length(e)){
      for(l in 1:length(e)) overlap.matrix[j, l] = length(intersect(pathway[[e[j]]], pathway[[e[l]]]))
    }
    overlap.matrix = matrix(paste(round(kappa.result[e, e], 2), "\n(", overlap.matrix, ")", sep = ""), nrow = length(e))

    d = as.dist(1-kappa)
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
    dev.off()
  }



  #================================
  #heatmap for all clusters
  #================================
  e = names(results2)
  colorIndex = results2[order(results2)]

  set.seed(1)
  nCol = max(results2)-1
  colOption = c(sample(rainbow(nCol),nCol),"black")
  colorbar = mapvalues(colorIndex,from = 1:(nCol+1),to = colOption)

  png(paste(output_clustering,"/Heatmap_clusters_all",".png",sep=""),
       height = if (length(e)<=40)  600
       else if (length(e)<=100) 700
       else if (length(e)<=150) 800
       else if (length(e)<=200) 1050
       else 1200,
       width = if (length(e)<=40) 650
       else if (length(e)<=100) 800
       else if (length(e)<=150) 1000
       else if (length(e)<=200) 1200
       else 1300)
  if(nrow(indi.heatmap.data_all)==1){
    image(t(-log10(indi.heatmap.data_all)))
  }else if (nrow(indi.heatmap.data_all)!=0){
    heatmap.pathway.pvalue<-indi.heatmap.data_all
    heatmap.pathway.pvalue[heatmap.pathway.pvalue < 10^-10] = 10^-10
    a = heatmap.2(-log10(heatmap.pathway.pvalue), col = (c(intpalette(c("#FFFFD5FF", "#FF0000FF"), numcol = 30),rep("#FF0000FF",30))),
                  Rowv = F, Colv = F, density.info="none", trace="none", scale = "none", margins = c(7, 20),
                  dendrogram = "none",
                  rowsep = as.numeric(cumsum(table(as.character(results2)))),
                  sepwidth = c(0.05,0.4),
                  cexRow =  if (length(e)==2) 1.5
                  else if (length(e)==3) 1.48
                  else if (length(e)==4) 1.4
                  else if (length(e)<=40) 0.2 + 1/log10(length(e))
                  else if (length(e)<=100) 0.2 + 1/log10(length(e))
                  else if (length(e)<=150) 0.07 + 1/log10(length(e))
                  else if (length(e)<=200) 0.03 + 1/log10(length(e))
                  else 0.01 + 1/log10(length(e)),
                  cexCol =  if (length(e)<=40) 0.2 + 0.6/log10(Num_of_gene_lists)
                  else if (length(e)<=100) 0.3 + 0.6/log10(Num_of_gene_lists)
                  else if (length(e)<=150) 0.5 + 0.6/log10(Num_of_gene_lists)
                  else 0.6 + 0.6/log10(Num_of_gene_lists),
                  RowSideColors = colorbar,
                  keysize = if (length(e)<=40) 1.5
                  else if (length(e)<=100) 1.3
                  else if (length(e)<=150) 1
                  else 0.9 ,
                  offsetRow = if (length(e)<=40) 0.25
                  else if (length(e)<=100) 0
                  else if (length(e)<=150) -0.1
                  else -0.2,
                  )
  }
  dev.off()

  pdf(paste(output_clustering,"/Heatmap_clusters_all",".pdf",sep=""),
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
                  dendrogram = "none",
                  rowsep = as.numeric(cumsum(table(as.character(results2)))),
                  sepwidth = c(0.05,0.4),
                  cexRow =  if (length(e)==2) 1.5
                  else if (length(e)==3) 1.48
                  else if (length(e)==4) 1.4
                  else if (length(e)<=40) 0.2 + 1/log10(length(e))
                  else if (length(e)<=100) 0.2 + 1/log10(length(e))
                  else if (length(e)<=150) 0.07 + 1/log10(length(e))
                  else if (length(e)<=200) 0.03 + 1/log10(length(e))
                  else 0.01 + 1/log10(length(e)),
                  cexCol =  if (length(e)<=40) 0.2 + 0.6/log10(Num_of_gene_lists)
                  else if (length(e)<=100) 0.3 + 0.6/log10(Num_of_gene_lists)
                  else if (length(e)<=150) 0.5 + 0.6/log10(Num_of_gene_lists)
                  else 0.6 + 0.6/log10(Num_of_gene_lists),
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
  #========================================================================================
  #MDS for all cluster, based on kappa similarity
  #========================================================================================
  cluster.assign = results2
  fit <- cmdscale(as.dist(1-kappa.result),k=2)
  x <- fit[,1]
  y <- fit[,2]
  xlimit <- ifelse(abs(min(x))>abs(max(x)),abs(min(x)),abs(max(x)))
  ylimit <- ifelse(abs(min(y))>abs(max(y)),abs(min(y)),abs(max(y)))
  xcenter <- tapply(x,as.factor(cluster.assign),mean)
  ycenter <- tapply(y,as.factor(cluster.assign),mean)


  C <- length(unique(cluster.assign))
  unique.color <- rainbow(C)
  unique.shape <- 1:C
  sizes <- shapes <- colors <- cluster.assign
  for(i in 1:(C)){
    colors[cluster.assign==i] <- unique.color[i]
    shapes[cluster.assign==i] <- unique.shape[i]
    sizes[cluster.assign==i] <- 2
  }
  if(!is.null(singleton)){
    colors[cluster.assign== max(cluster.assign)] <- "gray50"
    shapes[cluster.assign== max(cluster.assign)] <- 20
    sizes[cluster.assign== max(cluster.assign)] <- 2
  }

  pdf(paste(output_clustering,"/MDS_clusters_kappa_numbers_all",".pdf",sep=""))
      p <- ggplot() +
        ggtitle("") +
        xlab("Coordinate 1") + ylab("Coordinate 2") +
        xlim(c(-xlimit,xlimit)) + ylim(c(-ylimit,ylimit)) +
        geom_point(aes(x, y), shape=shapes,
                   color = colors ,size=sizes) +
        geom_point(aes(xcenter,ycenter),
                   shape=unique.shape, color = unique.color,
                   size =5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 15, hjust=0.5,face="bold"),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12))
      print(p)
      dev.off()



  #=============================
  ##### Text Mining
  #=============================
      cat("Performing Text Mining Analysis...\n")
      hashtbAll = hashtb
      hashtb = hashtb[hashtb[,3]%in%which(pathways %in% names(pathway)),]
      tmk = list()
      nperm = n.text.permute
      if (nrow(hashtb) == 0){
        for (i in 1:(k-1)){
          tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)
        }
      }else{
        paraFunc = function(i){
          e = results2[results2 == i]#names are names of pathways in ith cluster
          e = which(pathways %in% names(e))#pathways index of pathways in ith cluster
          hashcl = hashtb[hashtb [,3]%in%e,]#noun phrases appeared in ith cluster
          hashcl = hashcl[duplicated(hashcl[,2]) | duplicated(hashcl[,2], fromLast=TRUE),]
          if (nrow(hashcl) != 0){
            hashf = hashcl
            hashf[,2] = 1
            hashf = aggregate(hashf[,c("row","value")],by = hashf["phrase"],FUN=sum)
            #aggregate nouns in cluster i: noun;#appearance in cluster;sum of penalized value
            rownames(hashf) = hashf[,"phrase"]
            hashf = hashf[,-1]
            colnames(hashf) = c("count","sum") #count: count noun occurance in each cluster
            hashap = hashcl#to record noun phrases name in this cluster
            hashap[,c(2,3,4)] = 0 #dim as hashcl
            mperm = matrix(nrow = nrow(hashf),ncol = nperm)
            for (j in 1:nperm){
              subtb = hashtbAll[hashtbAll[,3]%in%sample(1:length(pathways),length(e)),]#consider permutating from all possible pathways
              subtb = rbind(subtb,hashap)#make sure rownames of (subtb) includes NPs appeared in original pathways
              subtb = subtb[subtb$phrase %in% hashap$phrase,]#test NPs included in original pathways
              subtb[,2] = 1
              subtb = aggregate(subtb[,c("row","value")],by = subtb["phrase"],FUN=sum)
              rownames(subtb) = subtb[,"phrase"]
              subtb = subtb[,-1]
              colnames(subtb) = c("count","sum")
              mperm[,j] = subtb[,2]
            }
            hashf[,"p-value"] = apply(cbind(hashf[,2],mperm),1,
                                      function(x)((nperm + 2)-rank(x)[1])/(nperm + 1))
            hashf[,"q-vlaue"] = p.adjust(hashf[,"p-value"],method = "BH")
            tmk = hashf[order(hashf[,3],-hashf[,2]),]#order first by pvalue, then statistic T
          }
          else {tmk = matrix(NA,nrow = 1,ncol = 4)}
          return(tmk)
        }
        if (parallel==TRUE){
          tmk = mclapply(1:(k-1),paraFunc,mc.cores = 10)
        }else{
          tmk = lapply(1:(k-1),paraFunc)
        }
      }
      # End of Text Mining
      tm_filtered <- list() #filter out count 1, q 0.05
      for (i in 1:(k-1)){
        tm_filtered[[i]] <- tmk[[i]][ which((as.numeric(tmk[[i]][,4]) < 0.05)), ]
      }

  setwd(output_clustering)

  topGO.summary = GO.cluster.20
  if(file.exists("Clustering_Summary.csv")){file.remove("Clustering_Summary.csv")}
  if (software == "MAPE"){
    cat("ID,Term,NumGeneTotalInSet,MAPE_p-value",file="Clustering_Summary.csv",append=T)
    cat("\n",file="Clustering_Summary.csv",append=T)
  }
  else if (software == "CPI"){
    if (enrichment == "KS"){
      cat("ID,Term,NumGeneTotalInSet",file="Clustering_Summary.csv",append=T)
      for (n in 1:Num_of_gene_lists)
        cat(",Study_p-value",colnames(summary)[n+2],file="Clustering_Summary.csv",append=T)
      cat("\n",file="Clustering_Summary.csv",append=T)
    }
    else if (enrichment == "Fisher's exact"){
      cat("ID,Term,NumGeneTotalInSet",file="Clustering_Summary.csv",append=T)
      for (n in 1:Num_of_gene_lists)
        cat(",NumDEGeneInSet_",colnames(summary)[n+2],file="Clustering_Summary.csv",append=T)
      for (n in 1:Num_of_gene_lists)
        cat(",DEGeneInSet_",colnames(summary)[n+2],file="Clustering_Summary.csv",append=T)
      cat("\n",file="Clustering_Summary.csv",append=T)
    }
  }




  if(is.null(dim(tm_filtered[[1]]))==TRUE){
    cat("Cluster 1\n", file = "Clustering_Summary.csv",append=T)
    cat("Key words,", file = "Clustering_Summary.csv",append=T)
    write.table(tm_filtered[[1]][1,], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = "Clustering_Summary.csv",append=T)
    write.table(tm_filtered[[1]][4,], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("count,", file = "Clustering_Summary.csv",append=T)
    write.table(tm_filtered[[1]][1,], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    write.table(topGO.summary[[1]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
  } else {
    cat("Cluster 1\n", file = "Clustering_Summary.csv",append=T)
    cat("Key words,", file = "Clustering_Summary.csv",append=T)
    write.table(t(as.character(rownames(tm_filtered[[1]])[1:15])), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = "Clustering_Summary.csv",append=T)
     if (dim(tm_filtered[[1]])[1] == 0){
      write.table(t(tm_filtered[[1]][1:15]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = "Clustering_Summary.csv",append=T)
      write.table(t(tm_filtered[[1]][1:15]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(topGO.summary[[1]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)}

      else {
      write.table(t(tm_filtered[[1]][1:15,4]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = "Clustering_Summary.csv",append=T)
      write.table(t(tm_filtered[[1]][1:15,1]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(topGO.summary[[1]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
      }
  }
  for (i in 2:(k-1)){
    if(is.null(dim(tm_filtered[[i]]))==TRUE){
      cat(paste("\nCluster ", i, "\n", sep = ""), file = "Clustering_Summary.csv", append = T)
      cat("Key words,", file = "Clustering_Summary.csv",append=T)
      write.table(tm_filtered[[i]][1,], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("q_value,", file = "Clustering_Summary.csv",append=T)
      write.table(tm_filtered[[i]][4,], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = "Clustering_Summary.csv",append=T)
      write.table(tm_filtered[[i]][1,], "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(topGO.summary[[i]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
    } else {
      cat(paste("\nCluster ", i, "\n", sep = ""), file = "Clustering_Summary.csv", append = T)
      cat("Key words,", file = "Clustering_Summary.csv",append=T)
      write.table(t(as.character(rownames(tm_filtered[[i]])[1:15])), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("q_value,", file = "Clustering_Summary.csv",append=T)

      if (dim(tm_filtered[[i]])[1] == 0){
      write.table(t(tm_filtered[[i]][1:15]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = "Clustering_Summary.csv",append=T)
      write.table(t(tm_filtered[[i]][1:15]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(topGO.summary[[i]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)}

      else {
      write.table(t(tm_filtered[[i]][1:15,4]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = "Clustering_Summary.csv",append=T)
      write.table(t(tm_filtered[[i]][1:15,1]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(topGO.summary[[i]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
      }
    }
  }
  i=k
  cat(paste("\nSingleton Term", "\n", sep = ""), file = "Clustering_Summary.csv", append = T)
  write.table(topGO.summary[[i]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)

  #=============================
  ##### Hierarchical Clustering
  #=============================
  pdf("CPI_hclust_average.pdf",width = 7.5,height = 7.5)
  cluster.num = 1:length(unique(results2))
  par(mfrow=c(3,3))
  for (c in 1:length(unique(results2))){
    pathways = names(results2)[which(results2==c)]
    pmat = summary[pathways,-which(colnames(summary) %in% c("q_value_meta","p_value_meta"))]
    d = dist(t(pmat))
    hc1 = hclust(d, method = "average")

    title = paste("Cluster ",cluster.num[c],sep="")
    plot(hc1, cex = 1.2, hang = -1,main= title,xlab = NULL,cex.main=2)

  }
  dev.off()

  setwd(wd)
}
