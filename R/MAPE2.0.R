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
##' @param permutation The options for using sample permutation or gene permutation when performing
##' enrichment analysis if permute is true. it is one of the two values: 'gene' and 'sample'. The default
##' option is sample permutation.
##' @param nperm Number of permutations to be performed.
##' @param size.min The minimum size of pathways to be considered. The default value is 15.
##' @param size max The maximum size of pathways to be considered. The default value is 500.
##' @param qvalue.cal The method to calculate the q-values if permute is true. The default method is to 
##' calcuate the q-values based on the permutation method. If qvalue.cal='estimate', the q-values
##' were estimated based on the Storey's method.
##' @param MetaDE Whether to use MetaDE's result to do pathway analysis. By default, it is false.
##' @return The qvalue and pvalue of each pathway including Meta qvalue and Meta pvalue.
##' @author Zhou Fang, Kui Shen, Xiangrui Zeng and George Tseng.
##' @export
##' @examples
##' data(clinical)
##' data(Leukemia_v2)
##' data(pathways)
##' MAPE2.0_result = MAPE2.0(arraydata = Leukemia,clinical.data = clinical,label = "label",
##'                         resp.type="multiclass",stat='maxP',method = "CPI", enrichment = "Fisher's exact", 
##'                         DEgene.number = 400,size.min=15,size.max=500)
##' 
##' data(MetaDEresult)
##' MAPE2.0_result = MAPE2.0(stat='maxP',method = "MAPE", enrichment = "KS", 
##'                         size.min=15,size.max=500, MetaDE = TRUE,
##'                         meta.p = meta.res.p$meta.analysis$pval,ind.p = ind.res$p)



MAPE2.0<-function (arraydata = NULL, clinical.data = NULL, label = NULL,censoring = NULL,
                   pathway = c(Biocarta.genesets,GOBP.genesets,GOCC.genesets,GOMF.genesets,
                               KEGG.genesets,Reactome.genesets), 
                    resp.type = c("twoclass", "multiclass", "continuous", "survival"),
                    method = c("CPI","MAPE"), enrichment = c("KS","Fisher's exact"), DEgene.number = 200,
                    stat = c("Fisher","maxP", "minP", "rth","AW Fisher"), rth.value = NULL, permute = F, 
                    permutation = c("sample", "gene"), nperm = 500, size.min = 15, size.max = 500, 
                    qvalue.cal = c("estimate","permutation"),MetaDE = F,meta.p = NULL,ind.p = NULL) 
{
  if (MetaDE == T) {
    if (is.null(meta.p) || is.null(ind.p)){
      stop("Please run MetaDE first to input MetaDE result.")
    }
  }
  
  if (MetaDE == T) {method = "MAPE"}
  
  if (method == "CPI") {
    method = match.arg(method)
    enrichment = match.arg(enrichment)
    stat = match.arg(stat)
    resp.type = match.arg(resp.type)
    permutation = match.arg(permutation)
    qvalue.cal = match.arg(qvalue.cal)
    data.type = class(arraydata[[1]])
    study = list()
    for (t1 in 1:length(arraydata)) {
      exprs = arraydata[[t1]]
      if (resp.type != "survival") {
        pheno = clinical.data[[t1]][label]
        rownames(pheno) = colnames(exprs)
        colnames(pheno) = "label"
      }
      else {
        yy = clinical.data[[t1]][label]
        zz = clinical.data[[t1]][censoring]
        pheno = as.data.frame(cbind(yy, zz))
        rownames(pheno) = arraydata[[t1]]$samplename
        colnames(pheno) = c("label", "censoring.status")
      }
      phenoData = new("AnnotatedDataFrame", data = pheno)
      study[[t1]] <- new("ExpressionSet", exprs = exprs, 
                         phenoData = phenoData)
    }
    if (is.null(names(arraydata))) {
      names(study) = paste("study", 1:length(arraydata), 
                           sep = "")
    }
    else {
      names(study) = names(arraydata)
    }
    label = "label"
    censoring.status = "censoring.status"
    rm(arraydata)
#    for (study.no in 1:length(study)) {
#      if (sum(is.na(exprs(study[[study.no]]))) > 0) {
#        exprs(study[[study.no]]) = impute::impute.knn(exprs(study[[study.no]]), 
#                                                      k = knn.neighbors)$data
#      }
#    }
    gene.in.array = featureNames(study[[1]])
    if (length(study) > 1) {
      for (t1 in 2:length(study)) {
        gene.in.array = unique(gene.in.array, featureNames(study[[t1]]))
      }
    }
    #    load("Hallmark.RData")
    #    load("Positional.RData")
    #    load("CGP.RData")
    #    load("Canonical.RData")
    #load("Biocarta.RData")
    #load("KEGG.RData")
    #load("Reactome.RData")
    #    load("MIR.RData")
    #    load("TFT.RData")
    #    load("CGN.RData")
    #    load("CM.RData")
    #load("GOBP.RData")
    #load("GOCC.RData")
    #load("GOMF.RData")
    #    load("Oncogenic.RData")
    #    load("Immunologic.RData")
    #    load("PID.RData")
    #    load("CMAP_up.RData")
    #    load("CMAP_down.RData")
    #    load("PPI.RData")
    #    load("JASPARhuman.RData")
    #    load("PITA.RData")
    #    load("miRanda.RData")
    #    load("TargetScan.RData")
    #load("Phenocarta.RData")
    #pathway<-c(Biocarta.genesets,KEGG.genesets,Reactome.genesets,GOBP.genesets,GOCC.genesets,
    #           GOMF.genesets,Phenocarta.genesets)
    pathway = pathway[which(sapply(pathway,length) >= size.min & sapply(pathway,length) <= size.max)]
    gene.in.DB = unique(unlist(pathway))
    set.name = names(pathway)
    gene.common = intersect(gene.in.array, gene.in.DB)
    DB.matrix = matrix(0, length(set.name), length(gene.common))
    rownames(DB.matrix) = set.name
    colnames(DB.matrix) = gene.common
    for (t1 in 1:length(set.name)) {
      gene = intersect(pathway[[t1]], gene.common)
      DB.matrix[set.name[t1], gene] = 1
    }
    colnames(DB.matrix) = toupper(colnames(DB.matrix))
    ##  keep.idx = (apply(DB.matrix, 1, sum) >= size.min & apply(DB.matrix, 
    ##                                                           1, sum) <= size.max)
    ##  DB.matrix = DB.matrix[keep.idx, ]
    if (length(study) == 1) {
      madata = study[[1]]
      testlabel = madata[[label]]
      if (resp.type == "survival") 
        censoring = madata[[censoring.status]]
      if (enrichment == "Fisher's exact") {
        enrich = EnrichmentC_Exact(madata = madata, label = testlabel,DEgene.number = DEgene.number,
                                   censoring = censoring, DB.matrix = DB.matrix, 
                                   size.min = size.min, size.max = size.max, nperm = nperm, 
                                   resp.type = resp.type)
      }
      else if (enrichment == "KS"){
        if (permute == F){
          enrich = Enrichment_KS (madata = madata, label = testlabel,
                                      censoring = censoring, DB.matrix = DB.matrix, 
                                      size.min = size.min, size.max = size.max, 
                                      resp.type = resp.type)
        }
        else if (permute == T) {
        if (permutation == "gene") {
          enrich = Enrichment_KS_gene(madata = madata, label = testlabel,
                                      censoring = censoring, DB.matrix = DB.matrix, 
                                      size.min = size.min, size.max = size.max, nperm = nperm, 
                                      resp.type = resp.type)
        }
        else if (permutation == "sample") {
          enrich = Enrichment_KS_sample(madata = madata, label = testlabel, 
                                        censoring = censoring, DB.matrix = DB.matrix, 
                                        size.min = size.min, size.max = size.max, nperm = nperm, 
                                        resp.type = resp.type)
        }
        }
      }
      else {
        stop("Please check: Wrong permutation methods.")
      }
      pvalue = as.data.frame(enrich$pvalue.set.0)
      colnames(pvalue) = "pvalue"
      qvalue = as.data.frame(enrich$qvalue.set.0)
      colnames(qvalue) = "qvalue"
      return(list(qvalue = qvalue, pvalue = pvalue))
    }
    else {
      if (enrichment == "Fisher's exact") {
        cat("Performing MAPE_P analysis...\n")
        MAP_SET.obj = CPI_Exact(study = study, label = label,DEgene.number = DEgene.number, 
                                censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                size.min = size.min, size.max = size.max, nperm = nperm, 
                                stat = stat, rth.value = rth.value, resp.type)
      }
      else if (enrichment == "KS"){
        if (permute == F){
          cat("Performing MAPE_P analysis...\n")
          MAP_SET.obj = CPI_KS(study = study, label = label, 
                                    censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                    size.min = size.min, size.max = size.max,  
                                    stat = stat, rth.value = rth.value, resp.type)
        }
        if (permute == T){
        if (permutation == "gene") {
          cat("Performing MAPE_P analysis...\n")
          MAP_SET.obj = CPI_gene_KS(study = study, label = label, 
                                    censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                    size.min = size.min, size.max = size.max, nperm = nperm, 
                                    stat = stat, rth.value = rth.value, resp.type)
        }
        else if (permutation == "sample") {
          cat("Performing MAPE_P analysis...\n")
          MAP_SET.obj = CPI_sample_KS(study = study, label = label, 
                                      censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                      size.min = size.min, size.max = size.max, nperm = nperm, 
                                      stat = stat, rth.value = rth.value, resp.type)
        }
       }
      }
      else {
        stop("Please check: Wrong permutation methods.")
      }
    }  
#    output_dir<-paste("./Output",Sys.time(),sep=" ")
#    dir.create(output_dir)
#    dir.create(paste(output_dir,"/Secondary_files",sep=""))
    
    if(enrichment == "KS"){return(list(summary = MAP_SET.obj$summary,pathway = pathway,Num_of_gene_lists = length(study),
                                       method = method,enrichment = enrichment))}
    else if (enrichment == "Fisher's exact")
    {return(list(summary = MAP_SET.obj$summary,pathway = pathway,Num_of_gene_lists = length(study),
                 method = method,enrichment = enrichment,genelist = MAP_SET.obj$genelist))}
  }
  
  
  else if (method == "MAPE"){
    if (MetaDE == T){
      method = match.arg(method)
      enrichment = match.arg(enrichment)
      stat = match.arg(stat)
      gene.in.array = rownames(ind.p)
      #    load("Hallmark.RData")
      #    load("Positional.RData")
      #    load("CGP.RData")
      #    load("Canonical.RData")
      #load("Biocarta.RData")
      #load("KEGG.RData")
      #load("Reactome.RData")
      #    load("MIR.RData")
      #    load("TFT.RData")
      #    load("CGN.RData")
      #    load("CM.RData")
      #load("GOBP.RData")
      #load("GOCC.RData")
      #load("GOMF.RData")
      #    load("Oncogenic.RData")
      #    load("Immunologic.RData")
      #    load("PID.RData")
      #    load("CMAP_up.RData")
      #    load("CMAP_down.RData")
      #    load("PPI.RData")
      #    load("JASPARhuman.RData")
      #    load("PITA.RData")
      #    load("miRanda.RData")
      #    load("TargetScan.RData")
      #load("Phenocarta.RData")
      #pathway<-c(Biocarta.genesets,KEGG.genesets,Reactome.genesets,GOBP.genesets,GOCC.genesets,
      #           GOMF.genesets,Phenocarta.genesets)
      pathway = pathway[which(sapply(pathway,length) >= size.min & sapply(pathway,length) <= size.max)]
      gene.in.DB = unique(unlist(pathway))
      set.name = names(pathway)
      gene.common = intersect(gene.in.array, gene.in.DB)
      DB.matrix = matrix(0, length(set.name), length(gene.common))
      rownames(DB.matrix) = set.name
      colnames(DB.matrix) = gene.common
      for (t1 in 1:length(set.name)) {
        gene = intersect(pathway[[t1]], gene.common)
        DB.matrix[set.name[t1], gene] = 1
      }
      colnames(DB.matrix) = toupper(colnames(DB.matrix))
      if (enrichment == "KS"){
      cat("Performing MAPE_P analysis...\n")
      MAP_SET.obj = MAPE_P_KS_DE(ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                 size.max = size.max, stat = stat, rth.value = rth.value)
      cat("Performing MAPE_G analysis...\n")
      MAP_GENE.obj = MAPE_G_KS_DE(meta.p = meta.p,DB.matrix,
                                  gene.common, size.min = size.min, size.max = size.max)
      cat("Performing MAPE_I analysis...\n")
      MAP_I.obj = MAPE_I (MAP_GENE.obj = MAP_GENE.obj, 
                          MAP_SET.obj = MAP_SET.obj, study = study)
      }
      else if (enrichment == "Fisher's exact" ){
      cat("Performing MAPE_P analysis...\n")
      MAP_SET.obj = MAPE_P_Exact_DE(ind.p = ind.p, DB.matrix, size.min = size.min, gene.common,
                                  size.max = size.max, stat = stat, rth.value = rth.value,
                                  DEgene.number = DEgene.number)
      cat("Performing MAPE_G analysis...\n")
      MAP_GENE.obj = MAPE_G_Exact_DE(meta.p = meta.p,DB.matrix,
                                  gene.common, size.min = size.min, size.max = size.max,
                                  DEgene.number = DEgene.number)
      cat("Performing MAPE_I analysis...\n")
      MAP_I.obj = MAPE_I (MAP_GENE.obj = MAP_GENE.obj, 
                          MAP_SET.obj = MAP_SET.obj, study = study)  
      }
      study.no = ncol(ind.p)
      study.name = colnames(ind.p)
      common.set.name = rownames(MAP_I.obj$qvalue.meta)
      qvalue.all = matrix(NA, length(common.set.name), (study.no + 
                                                          3))
      rownames(qvalue.all) = common.set.name
      colnames(qvalue.all) = c(study.name, "MAPE_P", "MAPE_G", 
                               "MAPE_I")
      qvalue.all[, 1:study.no] = MAP_SET.obj$qvalue.set.study[common.set.name,]
      qvalue.all[, (study.no + 1)] = MAP_SET.obj$qvalue.meta[common.set.name,]
      qvalue.all[, (study.no + 2)] = MAP_GENE.obj$qvalue.meta[common.set.name,]
      qvalue.all[, (study.no + 3)] = MAP_I.obj$qvalue.meta[common.set.name,]
      qvalue.all = as.data.frame(qvalue.all)
      pvalue.all = matrix(NA, length(common.set.name), (study.no + 
                                                          3))
      rownames(pvalue.all) = common.set.name
      colnames(pvalue.all) = c(study.name, "MAPE_P", "MAPE_G", 
                               "MAPE_I")
      pvalue.all[, 1:study.no] = MAP_SET.obj$pvalue.set.study[common.set.name,]
      pvalue.all[, (study.no + 1)] = MAP_SET.obj$pvalue.meta[common.set.name,]
      pvalue.all[, (study.no + 2)] = MAP_GENE.obj$pvalue.meta[common.set.name,]
      pvalue.all[, (study.no + 3)] = MAP_I.obj$pvalue.meta[common.set.name,]
      pvalue.all = as.data.frame(pvalue.all)
      pvalue.all = pvalue.all[!apply(pvalue.all[,1:study.no],1,sum) == study.no,]
      qvalue.all = qvalue.all[!apply(pvalue.all[,1:study.no],1,sum) == study.no,]
#      output_dir<-paste("./Output",Sys.time(),sep=" ")
#      dir.create(output_dir)
#      dir.create(paste(output_dir,"/Secondary_files",sep="")) 
      return(list(summary = qvalue.all, qvalue = qvalue.all, pvalue = pvalue.all,pathway = pathway,
                  Num_of_gene_lists = study.no ,method = method,enrichment = enrichment))
    
    }
    else {
    method = match.arg(method)
    enrichment = match.arg(enrichment)
    stat = match.arg(stat)
    resp.type = match.arg(resp.type)
    permutation = match.arg(permutation)
    qvalue.cal = match.arg(qvalue.cal)
    data.type = class(arraydata[[1]])
    study = list()
    for (t1 in 1:length(arraydata)) {
      exprs = arraydata[[t1]]
      if (resp.type != "survival") {
        pheno = clinical.data[[t1]][label] 
        rownames(pheno) = colnames(exprs)
        colnames(pheno) = "label"
      }
      else {
        yy = clinical.data[[t1]][label]
        zz = clinical.data[[t1]][censoring]
        pheno = as.data.frame(cbind(yy, zz))
        rownames(pheno) = arraydata[[t1]]$samplename
        colnames(pheno) = c("label", "censoring.status")
      }
      phenoData = new("AnnotatedDataFrame", data = pheno)
      study[[t1]] <- new("ExpressionSet", exprs = exprs, 
                         phenoData = phenoData)
    }
    if (is.null(names(arraydata))) {
      names(study) = paste("study", 1:length(arraydata), 
                           sep = "")
    }
    else {
      names(study) = names(arraydata)
    }
    label = "label"
    censoring.status = "censoring.status"
    rm(arraydata)
#    for (study.no in 1:length(study)) {
#      if (sum(is.na(exprs(study[[study.no]]))) > 0) {
#        exprs(study[[study.no]]) = impute::impute.knn(exprs(study[[study.no]]), 
#                                                      k = knn.neighbors)$data
#      }
#    }
    gene.in.array = featureNames(study[[1]])
    if (length(study) > 1) {
      for (t1 in 2:length(study)) {
        gene.in.array = unique(gene.in.array, featureNames(study[[t1]]))
      }
    }
    #    load("Hallmark.RData")
    #    load("Positional.RData")
    #    load("CGP.RData")
    #    load("Canonical.RData")
    #load("Biocarta.RData")
    #load("KEGG.RData")
    #load("Reactome.RData")
    #    load("MIR.RData")
    #    load("TFT.RData")
    #    load("CGN.RData")
    #    load("CM.RData")
    #load("GOBP.RData")
    #load("GOCC.RData")
    #load("GOMF.RData")
    #    load("Oncogenic.RData")
    #    load("Immunologic.RData")
    #    load("PID.RData")
    #    load("CMAP_up.RData")
    #    load("CMAP_down.RData")
    #    load("PPI.RData")
    #    load("JASPARhuman.RData")
    #    load("PITA.RData")
    #    load("miRanda.RData")
    #    load("TargetScan.RData")
    #load("Phenocarta.RData")
    #pathway<-c(Biocarta.genesets,KEGG.genesets,Reactome.genesets,GOBP.genesets,GOCC.genesets,
    #           GOMF.genesets,Phenocarta.genesets)
    pathway = pathway[which(sapply(pathway,length) >= size.min & sapply(pathway,length) <= size.max)]
    gene.in.DB = unique(unlist(pathway))
    set.name = names(pathway)
    gene.common = intersect(gene.in.array, gene.in.DB)
    DB.matrix = matrix(0, length(set.name), length(gene.common))
    rownames(DB.matrix) = set.name
    colnames(DB.matrix) = gene.common
    for (t1 in 1:length(set.name)) {
      gene = intersect(pathway[[t1]], gene.common)
      DB.matrix[set.name[t1], gene] = 1
    }
    colnames(DB.matrix) = toupper(colnames(DB.matrix))
    #  keep.idx = (apply(DB.matrix, 1, sum) >= size.min & apply(DB.matrix, 
    #                                                           1, sum) <= size.max)
    #  DB.matrix = DB.matrix[keep.idx, ]
    if (length(study) == 1) {
      madata = study[[1]]
      testlabel = madata[[label]]
      if (resp.type == "survival") 
        censoring = madata[[censoring.status]]
      if (permute == T){
      if (enrichment == "Fisher's exact") {
        enrich = EnrichmentM_Exact(madata = madata, label = testlabel,DEgene.number = DEgene.number,
                                   censoring = censoring, DB.matrix = DB.matrix, 
                                   size.min = size.min, size.max = size.max, nperm = nperm, 
                                   resp.type = resp.type)
      }
      else if (enrichment == "KS"){
        if (permutation == "gene") {
          enrich = Enrichment_KS_gene(madata = madata, label = testlabel, 
                                      censoring = censoring, DB.matrix = DB.matrix, 
                                      size.min = size.min, size.max = size.max, nperm = nperm, 
                                      resp.type = resp.type)
        }
        else if (permutation == "sample") {
          enrich = Enrichment_KS_sample(madata = madata, label = testlabel, 
                                        censoring = censoring, DB.matrix = DB.matrix, 
                                        size.min = size.min, size.max = size.max, nperm = nperm, 
                                        resp.type = resp.type)
        }
        else {
          stop("Please check: Wrong permutation methods.")
        }
      }
      }
      else {
        enrich = Enrichment_KS (madata = madata, label = testlabel, 
                                    censoring = censoring, DB.matrix = DB.matrix, 
                                    size.min = size.min, size.max = size.max, 
                                    resp.type = resp.type)
      }
      pvalue = as.data.frame(enrich$pvalue.set.0)
      colnames(pvalue) = "pvalue"
      qvalue = as.data.frame(enrich$qvalue.set.0)
      colnames(qvalue) = "qvalue"
      return(list(qvalue = qvalue, pvalue = pvalue))
    }
    else {
      if (permute == T){
      if (enrichment == "KS"){
        if (permutation == "gene") {
          cat("Performing MAPE_P analysis...\n")
          MAP_SET.obj = MAPE_P_gene_KS(study = study, label = label, 
                                       censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                       size.min = size.min, size.max = size.max, nperm = nperm, 
                                       stat = stat, rth.value = rth.value, resp.type)
          cat("Performing MAPE_G analysis...\n")
          MAP_GENE.obj = MAPE_G_gene_KS(study = study, label = label, 
                                        censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                        size.min = size.min, size.max = size.max, nperm = nperm, 
                                        stat = stat, rth.value = rth.value, resp.type)
          cat("Performing MAPE_I analysis...\n")
          MAP_I.obj = MAPE_I_KS(MAP_GENE.obj = MAP_GENE.obj, 
                                MAP_SET.obj = MAP_SET.obj, study = study)
        }
        else if (permutation == "sample") {
          cat("Performing MAPE_P analysis...\n")
          MAP_SET.obj = MAPE_P_sample_KS(study = study, label = label, 
                                         censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                         size.min = size.min, size.max = size.max, nperm = nperm, 
                                         stat = stat, rth.value = rth.value, resp.type)
          cat("Performing MAPE_G analysis...\n")
          MAP_GENE.obj = MAPE_G_sample_KS(study = study, label = label, 
                                          censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                          size.min = size.min, size.max = size.max, nperm = nperm, 
                                          stat = stat, rth.value = rth.value, resp.type)
          cat("Performing MAPE_I analysis...\n")
          MAP_I.obj = MAPE_I_KS(MAP_GENE.obj = MAP_GENE.obj, 
                                MAP_SET.obj = MAP_SET.obj, study = study)
        }
        else {
          stop("Please check: Wrong permutation methods.")
        }
      }
      }
      else {
        if (enrichment == "Fisher's exact") {
          cat("Performing MAPE_P analysis...\n")
          MAP_SET.obj = MAPE_P_Exact(study = study, label = label,DEgene.number = DEgene.number, 
                                     censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                     size.min = size.min, size.max = size.max, nperm = nperm, 
                                     stat = stat, rth.value = rth.value, resp.type)
          cat("Performing MAPE_G analysis...\n")
          MAP_GENE.obj = MAPE_G_Exact(study = study, label = label,DEgene.number = DEgene.number, 
                                      censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                      size.min = size.min, size.max = size.max, nperm = nperm, 
                                      stat = stat, rth.value = rth.value, resp.type)
          cat("Performing MAPE_I analysis...\n")
          MAP_I.obj = MAPE_I (MAP_GENE.obj = MAP_GENE.obj, 
                                MAP_SET.obj = MAP_SET.obj, study = study)
        }
        
        if (enrichment == "KS"){
        cat("Performing MAPE_P analysis...\n")
        MAP_SET.obj = MAPE_P_KS(study = study, label = label, 
                                     censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                     size.min = size.min, size.max = size.max,  
                                     stat = stat, rth.value = rth.value, resp.type)
        cat("Performing MAPE_G analysis...\n")
        MAP_GENE.obj = MAPE_G_KS(study = study, label = label, 
                                      censoring.status = censoring.status, DB.matrix = DB.matrix, 
                                      size.min = size.min, size.max = size.max,  
                                      stat = stat, rth.value = rth.value, resp.type)
        cat("Performing MAPE_I analysis...\n")
        MAP_I.obj = MAPE_I (MAP_GENE.obj = MAP_GENE.obj, 
                              MAP_SET.obj = MAP_SET.obj, study = study)
        }
      }
      study.no = length(study)
      study.name = names(study)
      common.set.name = rownames(MAP_I.obj$qvalue.meta)
      qvalue.all = matrix(NA, length(common.set.name), (study.no + 
                                                          3))
      rownames(qvalue.all) = common.set.name
      colnames(qvalue.all) = c(names(study), "MAPE_P", "MAPE_G", 
                               "MAPE_I")
      qvalue.all[, 1:study.no] = MAP_SET.obj$qvalue.set.study[common.set.name,]
      qvalue.all[, (study.no + 1)] = MAP_SET.obj$qvalue.meta[common.set.name,]
      qvalue.all[, (study.no + 2)] = MAP_GENE.obj$qvalue.meta[common.set.name,]
      qvalue.all[, (study.no + 3)] = MAP_I.obj$qvalue.meta[common.set.name,]
      qvalue.all = as.data.frame(qvalue.all)
      pvalue.all = matrix(NA, length(common.set.name), (study.no + 
                                                          3))
      rownames(pvalue.all) = common.set.name
      colnames(pvalue.all) = c(names(study), "MAPE_P", "MAPE_G", 
                               "MAPE_I")
      pvalue.all[, 1:study.no] = MAP_SET.obj$pvalue.set.study[common.set.name,]
      pvalue.all[, (study.no + 1)] = MAP_SET.obj$pvalue.meta[common.set.name,]
      pvalue.all[, (study.no + 2)] = MAP_GENE.obj$pvalue.meta[common.set.name,]
      pvalue.all[, (study.no + 3)] = MAP_I.obj$pvalue.meta[common.set.name,]
      pvalue.all = as.data.frame(pvalue.all)
      pvalue.all = pvalue.all[!apply(pvalue.all[,1:length(study)],1,sum) == length(study),]
      qvalue.all = qvalue.all[!apply(pvalue.all[,1:length(study)],1,sum) == length(study),]
      if (qvalue.cal == "estimate") {
        for (t1 in 1:ncol(pvalue.all)) {
          qvalue.all[, t1] = p.adjust(pvalue.all[, t1], 
                                      "BH")
        }
      }
#      output_dir<-paste(output.dir,Sys.time(),sep=" ")
#      dir.create(output_dir)
#      dir.create(paste(output_dir,"/Secondary_files",sep="")) 
      return(list(summary = qvalue.all, qvalue = qvalue.all, pvalue = pvalue.all,pathway = pathway,
                  Num_of_gene_lists = length(study),method = method,enrichment = enrichment
                  ))
    }
    }
  }
  else {stop("Please check: Wrong method.") }
}
