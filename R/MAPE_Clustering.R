##' This is the third major function in the MetaPath2.0 package which outputs the pathway
##' cluster results with text mining anotation.

##' @title Output pathway cluster results
##' @param Num_Clusters number of clusters
##' @return The pathway cluster results in a csv file.
##' @authors Zhou Fang, Xiangrui Zeng and George Tseng.
##' @export
##' @examples
##' data(clinical)
##' data(Leukemia_v2)
##' data(pathways)
##' select.group <- c('inv(16)','t(15;17)')
##' ref.level <- "inv(16)"
##' data.type <- "continuous"
##' ind.method <- c('limma','limma','limma')
##' resp.type <- "twoclass"
##' MAPE2.0_result = MAPE2.0(arraydata = Leukemia,clinical.data = clinical,label = "label",
##'                         resp.type=resp.type,stat='maxP',method = "CPI", enrichment = "Fisher's exact", 
##'                         DEgene.number = 400,size.min=15,size.max=500,data.type=data.type,
##'                         ind.method=ind.method,ref.level=ref.level,select.group=select.group, 
##'                         paired = paired)
##' MAPE.kappa_result = MAPE.Kappa(summary = MAPE2.0_result$summary, software = MAPE2.0_result$method,
##'                                pathway = MAPE2.0_result$pathway, max_k = 10, q_cutoff = 0.0005,
##'                                output_dir = tempdir())
##' MAPE.Clustering(summary = MAPE2.0_result$summary,Num_Clusters=5,
##'                Num_of_gene_lists = MAPE2.0_result$Num_of_gene_lists,
##'                genelist = MAPE2.0_result$genelist,kappa.result = MAPE.kappa_result$kappa, 
##'                pathway = MAPE2.0_result$pathway, enrichment = MAPE2.0_result$enrichment,
##'                method = MAPE.kappa_result$method,software = MAPE2.0_result$method,
##'                output_dir = getwd())

MAPE.Clustering <- function(summary,Num_Clusters = 3, kappa.result = kappa.result, Num_of_gene_lists, 
                            genelist = NULL, pathway, enrichment,method,software,
                            output_dir = getwd())
  {
  wd = getwd()
  k = Num_Clusters
  d = as.dist(1-kappa.result)
  results <- pam(d, k, diss=T)
  
  b = summary
  dir.create(paste(output_dir,"/Clustering_files","_",Num_Clusters,"_","clusters",sep=""))
  output_clustering <- paste(output_dir,"/Clustering_files","_",Num_Clusters,"_","clusters",sep="")
  sil <- silhouette(results$clustering, (1-kappa.result), diss=T)
  
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
  
  singleton<-results$clustering[!names(results$clustering)%in%names(results2)]
  singleton[1:length(singleton)]<- max(results2)+1
  results2<-c(results2,singleton)
  results2 <- results2[order(results2)]
#  if (!all(1:k %in% unique(results$clustering))){
#    stop('Please choose a smaller number of clusters k or choose a larger FDR cutoff to select more pathways for clustering')
#  }
  k = max(results2)
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
  colorbar <- as.character( results2[order(results2)]+1 ) 
  
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
  colorbar <- as.character(results2+1 )   
  
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
  
  
  
  #=============================
  ##### Text Mining
  #=============================
  cat("Performing Text Mining Analysis...\n")
  data(hashtb)
  hashtb = hashtb[hashtb [,2]%in%which(pathways %in% names(pathway)),]
  tmk = list()
  nperm = 1000
  if (nrow(hashtb) == 0){
    for (i in 1:(k-1)){
      tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)
    }
  }
  else{
  for (i in 1:(k-1)){
    e = results2[results2 == i]
    e = which(pathways %in% names(e))
    hashcl = hashtb[hashtb [,2]%in%e,]
    hashcl = hashcl[duplicated(hashcl[,1]) | duplicated(hashcl[,1], fromLast=TRUE),]
    if (nrow(hashcl) != 0){
      hashf = hashcl
      hashf[,1] = 1
      hashf = aggregate(hashf[,-2] ~ rownames(hashf),data=hashf, FUN=sum)
      rownames(hashf) = hashf[,1]
      hashf = hashf[,-1]
      colnames(hashf) = c("count","sum")
      hashap = hashcl
      hashap[] = 0 
      mperm = matrix(nrow = nrow(hashf),ncol = nperm)
      for (j in 1:nperm){
        subtb = hashtb[hashtb [,2]%in%sample(1:length(pathways),length(e)),]
        subtb = rbind(subtb,hashap)
        subtb = subtb[rownames(subtb) %in% rownames(hashap),]
        subtb[,1] = 1
        subtb = aggregate(subtb[,-2] ~ rownames(subtb),data=subtb, FUN=sum)
        rownames(subtb) = subtb[,1]
        subtb = subtb[,-1]
        colnames(subtb) = c("count","sum")
        subtb = subtb[rownames(subtb) %in% rownames(hashf),]
        mperm[,j] = subtb[,2]
      }
      hashf[,"p-value"] = apply(cbind(hashf[,2],mperm),1,
                                function(x)((nperm + 2)-rank(x)[1])/(nperm + 1))
      hashf[,"q-vlaue"] = p.adjust(hashf[,"p-value"],method = "BH")
      tmk[[i]] = hashf[order(hashf[,3],-hashf[,2]),]
    }
    else {tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)}
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
        cat(",Study_p-value",colnames(summary[,n+2]),file="Clustering_Summary.csv",append=T)
      cat("\n",file="Clustering_Summary.csv",append=T)
    }
    else if (enrichment == "Fisher's exact"){
      cat("ID,Term,NumGeneTotalInSet",file="Clustering_Summary.csv",append=T)
      for (n in 1:Num_of_gene_lists)
        cat(",NumDEGeneInSet_",colnames(summary[,n+2]),file="Clustering_Summary.csv",append=T)
      for (n in 1:Num_of_gene_lists)
        cat(",DEGeneInSet_",colnames(summary[,n+2]),file="Clustering_Summary.csv",append=T)
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
    write.table(t(rownames(tm_filtered[[1]])[1:15]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = "Clustering_Summary.csv",append=T)
    write.table(t(tm_filtered[[1]][1:15,4]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("count,", file = "Clustering_Summary.csv",append=T)
    write.table(t(tm_filtered[[1]][1:15,1]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    write.table(topGO.summary[[1]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
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
      write.table(t(rownames(tm_filtered[[i]])[1:15]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("q_value,", file = "Clustering_Summary.csv",append=T)
      write.table(t(tm_filtered[[i]][1:15,4]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = "Clustering_Summary.csv",append=T)
      write.table(t(tm_filtered[[i]][1:15,1]), "Clustering_Summary.csv", sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(topGO.summary[[i]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
    }
  }
  i=k
  cat(paste("\nSingleton Term", "\n", sep = ""), file = "Clustering_Summary.csv", append = T)
  write.table(topGO.summary[[i]], "Clustering_Summary.csv", sep=",",quote=T, append = T, row.names=F,col.names=F)
  setwd(wd)
  
}
