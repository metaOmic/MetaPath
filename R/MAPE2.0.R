##' This is the first major function in the MetaPath2.0 package which combines the Meta-analysis
##' for Pathway Enrichment (MAPE) methods introduced by Shen and Tseng (2010) and the
##' Comparative Pathway Integrator (CPI) method introduced by Fang and Tseng (2016).
##' The default function is CPI which performs MAPE_P (integrating multiple studies at pathway level) 
##' with Adaptively Weighted Fisher's method as Meta-analysis statistics. The other function, MAPE,
##' automatically performs MAPE_G (integrating multiple studies at gene level), MAPE_P and
##' MAPE_I (a hybrid method integrating MAEP_G and MAPE_P methods).
##'
##' For MAPE, in the simulation and real data analyses in the paper, MAPE_G and MAPE_P have complementary
##' advantages and detection power depending on the data structure. In general, the integrative
##' form of MAPE_I is recommended to use. In the case that MAPE_G (or MAPE_P) detects almost
##' none pathway, the integrative MAPE_I does not improve performance and MAPE_P (or MAPE_G)
##' should be used.
##' @title Perform the Meta-Analysis for Pathway Enrichment
##' @param arraydata The arraydata is a list of microarray data sets. Each microarray data set can be
##' either an ExpressionSet or a list. If the microarray data set is a list, then it includes
##' five elements as follows: 1)x-exprs data 2)y- the phenotype of interests
##' 3)z- censoring.status if applicable. 1 stands for the event occurred and 0 stands
##' for censored. 4)geneid 5)samplename If the microarray data set is in an ExpressionSet
##' format, the users need to 1) store the phenotype of interests in the slot
##' 'label'. 2) store the censor data is the slot 'censoring.status' if applicable
##' @param pathway The pathway databases
##' @param resp.type The phenotype of interest.It is one of the four values: 'twoclass','multiclass',
##' 'continuous','survival'.
##' @param method The method for overall analysis. It is one of the two values: 'CPI','MAPE'.
##' @param enrichment The method for pathway enrichment analysis. It is one of the two values: 'KS',
##' 'Fisher's exact'.
##' @param DEgene.number For Fisher's exact method, the number of differentially expressed genes
##' identified from each studies.
##' @param stat The meta-analysis statistic to be used to combine multiple studies. For MAPE, it is
##' one of the five values: 'minP','maxP','rth','Fisher','AW Fisher'.For CPI, AW Fisher's method is the only option.
##' @param rth.value The value of the rth statistic if the meta-anlaysis statistic is 'rth'.
##' For example,rth.value = 2.
##' @param permute Whether to use permutation to calculate p-values. By default, it is false.
##' @param nperm Number of permutations to be performed.
##' @param size.min The minimum size of pathways to be considered. The default value is 15.
##' @param size max The maximum size of pathways to be considered. The default value is 500.
##' @param qvalue.cal The method to calculate the q-values if permute is true. The default method is to 
##' calcuate the q-values based on the permutation method. If qvalue.cal='estimate', the q-values
##' were estimated based on the Storey's method.
##' @param clinical.data Clinical data files
##' @param label Label selected from clinical.data
##' @param data.type It is one of the two values: 'continuous','discrete'.
##' @param pmtx Option for uploading p-value matrix.
##' @return The qvalue and pvalue of each pathway including Meta qvalue and Meta pvalue.
##' @author Kui Shen, Xiangrui Zeng and George Tseng.
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
##' paired <- rep(FALSE,length(Leukemia))
##' MAPE2.0_result_CPI = MAPE2.0(arraydata = Leukemia,clinical.data = clinical,label = "label",
##'                         resp.type=resp.type,stat='maxP',method = "CPI", enrichment = "Fisher's exact", 
##'                         DEgene.number = 400,size.min=15,size.max=500,data.type=data.type,
##'                         ind.method=ind.method,ref.level=ref.level,paired=paired,select.group=select.group)
##' 
##' MAPE2.0_result_MAPE = MAPE2.0(arraydata = Leukemia,clinical.data = clinical,label = "label",
##'                         resp.type=resp.type,stat='minP',method = "MAPE", enrichment = "KS", 
##'                         DEgene.number = 400,size.min=15,size.max=500,data.type=data.type,
##'                         ind.method=ind.method,ref.level=ref.level,paired=paired,select.group=select.group)




MAPE2.0<-function (arraydata = NULL, clinical.data = NULL, label = NULL,pmtx = NULL,
                   pathway = c(Biocarta.genesets,GOBP.genesets,GOCC.genesets,GOMF.genesets,
                               KEGG.genesets,Reactome.genesets), 
                    data.type = c("continuous", "discrete"), 
                    mixed = NULL, mix.type = NULL,
                    covariate = NULL,
                    ref.level, paired, ind.method, select.group,tail="abs",
                    resp.type = c("twoclass", "multiclass", "continuous", "survival"),
                    method = c("CPI","MAPE"), enrichment = c("KS","Fisher's exact"), DEgene.number = 200,
                    stat = c("Fisher","maxP", "minP", "rth","AW Fisher"), rth.value = NULL, permute = F, 
                    nperm = 500, size.min = 15, size.max = 500, 
                    qvalue.cal = c("estimate","permutation"))
{ 
  if(is.null(pmtx)){
  
   if(!is.null(mixed)){
    if(mixed==T) {
    	   data.type <- mix.type
    }
   } 
  	
  ind.res <- Indi.DE.Analysis (data = arraydata,clin.data= clinical.data, 
                           data.type=data.type,resp.type = resp.type,
                           response = label ,covariate = covariate,tail = tail,
                           ind.method=ind.method,select.group = select.group,
                           ref.level=ref.level,paired=paired)
  ind.p = ind.res$p
  rm(arraydata)
  }
  else {
    ind.p = pmtx
  }
  
  
  if (method == "CPI") {
    enrichment = match.arg(enrichment)
    stat = "AW Fisher"
    enrichment = match.arg(enrichment)
    pathway = pathway[which(sapply(pathway,length) >= size.min & sapply(pathway,length) <= size.max)]
    gene.in.DB = unique(unlist(pathway))
    set.name = names(pathway)
    gene.in.array = rownames(ind.p)
    gene.common = intersect(gene.in.array, gene.in.DB)
    DB.matrix = matrix(0, length(set.name), length(gene.common))
    rownames(DB.matrix) = set.name
    colnames(DB.matrix) = gene.common
    colnames(DB.matrix) = toupper(colnames(DB.matrix))
    for (t1 in 1:length(set.name)) {
      gene = toupper(intersect(pathway[[t1]], gene.common))
      DB.matrix[set.name[t1], gene] = 1
    }
    
      if (enrichment == "Fisher's exact") {
        cat("Performing MAPE_P analysis...\n")
        MAP_SET.obj = MAPE_P_Exact_DE (ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                       size.max = size.max, stat = stat, rth.value = rth.value,
                                       DEgene.number = DEgene.number)
      }
      else if (enrichment == "KS"){
          cat("Performing MAPE_P analysis...\n")
          MAP_SET.obj = MAPE_P_KS_DE (ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                      size.max = size.max, stat = stat, rth.value = rth.value)
      }
      else {
        stop("Please check: Wrong permutation methods.")
      }

    summary = data.frame(list(q_value_meta = MAP_SET.obj$qvalue.meta,p_value_meta = MAP_SET.obj$pvalue.meta))
    for (i in 1:ncol(ind.p)){
      summary[colnames(ind.p)[i]] = MAP_SET.obj$pvalue.set.study[,i]
    }
    
    if(enrichment == "KS"){return(list(summary = summary,pathway = pathway,Num_of_gene_lists = ncol(ind.p),
                                       method = method,enrichment = enrichment))}
    else if (enrichment == "Fisher's exact")
    {return(list(summary = summary,pathway = pathway,Num_of_gene_lists = ncol(ind.p),
                 method = method,enrichment = enrichment,genelist = MAP_SET.obj$genelist))}
  }
  
  
  else if (method == "MAPE"){
    enrichment = match.arg(enrichment)
    qvalue.cal = match.arg(qvalue.cal)
    stat = match.arg(stat)
    enrichment = match.arg(enrichment)
    pathway = pathway[which(sapply(pathway,length) >= size.min & sapply(pathway,length) <= size.max)]
    gene.in.DB = unique(unlist(pathway))
    set.name = names(pathway)
    gene.in.array = rownames(ind.p)
    gene.common = intersect(gene.in.array, gene.in.DB)
    DB.matrix = matrix(0, length(set.name), length(gene.common))
    rownames(DB.matrix) = set.name
    colnames(DB.matrix) = gene.common
    colnames(DB.matrix) = toupper(colnames(DB.matrix))
    for (t1 in 1:length(set.name)) {
      gene = toupper(intersect(pathway[[t1]], gene.common))
      DB.matrix[set.name[t1], gene] = 1
    }
      if (permute == T){
      if (enrichment == "KS"){
          cat("Performing MAPE_P analysis...\n")
          MAP_SET.obj = MAPE_P_KS_gene (ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                      size.max = size.max, stat = stat, rth.value = rth.value,nperm = nperm)
          cat("Performing MAPE_G analysis...\n")
          MAP_GENE.obj = MAPE_G_KS_gene (ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                         size.max = size.max, stat = stat, rth.value = rth.value,nperm = nperm)
          cat("Performing MAPE_I analysis...\n")
          MAP_I.obj = MAPE_I_KS(MAP_GENE.obj = MAP_GENE.obj, 
                                MAP_SET.obj = MAP_SET.obj)
      }
      }
      else {
        if (enrichment == "Fisher's exact") {
          cat("Performing MAPE_P analysis...\n")
          MAP_SET.obj = MAPE_P_Exact_DE (ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                         size.max = size.max, stat = stat, rth.value = rth.value,
                                         DEgene.number = DEgene.number)
          cat("Performing MAPE_G analysis...\n")
          MAP_GENE.obj = MAPE_G_Exact_DE (ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                          size.max = size.max, stat = stat, rth.value = rth.value,
                                          DEgene.number = DEgene.number)
          cat("Performing MAPE_I analysis...\n")
          MAP_I.obj = MAPE_I (MAP_GENE.obj = MAP_GENE.obj, 
                                MAP_SET.obj = MAP_SET.obj)
        }
        
        if (enrichment == "KS"){
        cat("Performing MAPE_P analysis...\n")
        MAP_SET.obj = MAPE_P_KS_DE (ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                                  size.max = size.max, stat = stat, rth.value = rth.value)
        cat("Performing MAPE_G analysis...\n")
        MAP_GENE.obj = MAPE_G_KS_DE (ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                    size.max = size.max, stat = stat, rth.value = rth.value)
        cat("Performing MAPE_I analysis...\n")
        MAP_I.obj = MAPE_I (MAP_GENE.obj = MAP_GENE.obj, 
                              MAP_SET.obj = MAP_SET.obj)
        }
      }
      study.no = ncol(ind.p)
      study.name = colnames(ind.p)
      common.set.name = rownames(MAP_I.obj$qvalue.meta)
      qvalue.all = matrix(NA, length(common.set.name), (study.no + 
                                                          3))
      rownames(qvalue.all) = common.set.name
      colnames(qvalue.all) = c(colnames(ind.p), "MAPE_P", "MAPE_G", 
                               "MAPE_I")
      qvalue.all[, 1:study.no] = MAP_SET.obj$qvalue.set.study[common.set.name,]
      qvalue.all[, (study.no + 1)] = MAP_SET.obj$qvalue.meta[common.set.name,]
      qvalue.all[, (study.no + 2)] = MAP_GENE.obj$qvalue.meta[common.set.name,]
      qvalue.all[, (study.no + 3)] = MAP_I.obj$qvalue.meta[common.set.name,]
      qvalue.all = as.data.frame(qvalue.all)
      pvalue.all = matrix(NA, length(common.set.name), (study.no + 
                                                          3))
      rownames(pvalue.all) = common.set.name
      colnames(pvalue.all) = c(colnames(ind.p), "MAPE_P", "MAPE_G", 
                               "MAPE_I")
      pvalue.all[, 1:study.no] = MAP_SET.obj$pvalue.set.study[common.set.name,]
      pvalue.all[, (study.no + 1)] = MAP_SET.obj$pvalue.meta[common.set.name,]
      pvalue.all[, (study.no + 2)] = MAP_GENE.obj$pvalue.meta[common.set.name,]
      pvalue.all[, (study.no + 3)] = MAP_I.obj$pvalue.meta[common.set.name,]
      pvalue.all = as.data.frame(pvalue.all)
      pvalue.all = pvalue.all[!apply(pvalue.all[,1:ncol(ind.p)],1,sum) == ncol(ind.p),]
      qvalue.all = qvalue.all[rownames(pvalue.all),]
      if (qvalue.cal == "estimate") {
        for (t1 in 1:ncol(pvalue.all)) {
          qvalue.all[, t1] = p.adjust(pvalue.all[, t1], 
                                      "BH")
        }
      }
#      output_dir<-paste(output.dir,Sys.time(),sep=" ")
#      dir.create(output_dir)
#      dir.create(paste(output_dir,"/Secondary_files",sep="")) 
      summary = pvalue.all
      summary["MAPE_P_FDR"] = qvalue.all[,(ncol(ind.p)+1)]
      summary["MAPE_G_FDR"] = qvalue.all[,(ncol(ind.p)+2)]
      summary["MAPE_I_FDR"] = qvalue.all[,(ncol(ind.p)+3)]
      return(list(summary = summary, qvalue = qvalue.all, pvalue = pvalue.all,pathway = pathway,
                  Num_of_gene_lists = ncol(ind.p),method = method,enrichment = enrichment
                  ))
    }
  else {stop("Please check: Wrong method.") }
}
