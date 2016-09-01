## Internal functions of MAPE
##
## @title Internal functions
## @author Kui Shen, Xiangrui Zeng and George Tseng.

###########################
cor.func <- function (x, y) 
{
  n <- length(y)
  xbar <- x %*% rep(1/n, n)
  sxx <- ((x - as.vector(xbar))^2) %*% rep(1, n)
  sxy <- (x - as.vector(xbar)) %*% (y - mean(y))
  syy <- sum((y - mean(y))^2)
  numer <- sxy/sxx
  sd <- sqrt((syy/sxx - numer^2)/(n - 2))
  tt <- numer/sd
  return(list(tt = tt, numer = numer, sd = sd))
}


###########################
coxfunc <- function (x, y, censoring.status) 
{
  junk <- coxscor(x, y, censoring.status)
  scor <- junk$scor
  sd <- sqrt(coxvar(x, y, censoring.status, coxstuff.obj = junk$coxstuff.obj))
  tt <- scor/sd
  return(list(tt = tt, numer = scor, sd = sd))
}


###########################
cox.perm.sample <- function (expr, testgroup, censoring, nperm) 
{
  obs = coxfunc(expr, testgroup, censoring)$tt
  perms = matrix(NA, nrow(expr), nperm)
  perms.mtx = matrix(NA, nrow = nperm, ncol = length(testgroup))
  for (i in 1:nperm) {
    perms.mtx[i, ] = sample(1:length(testgroup), size = length(testgroup))
  }
  for (t1 in 1:nperm) {
    perms[, t1] = coxfunc(expr, testgroup[perms.mtx[t1, ]], 
                          censoring)$tt
  }
  rownames(perms) = rownames(expr)
  colnames(perms) = paste("B", 1:nperm, sep = "")
  obs = abs(obs)
  perms = abs(perms)
  return(list(obs = obs, perms = perms, perms.mtx = perms.mtx))
}


###########################
coxscor <- function (x, y, ic, offset = rep(0, length(y))) 
{
  n <- length(y)
  nx <- nrow(x)
  yy <- y + (ic == 0) * (1e-05)
  otag <- order(yy)
  y <- y[otag]
  ic <- ic[otag]
  x <- x[, otag, drop = F]
  offset <- offset[otag]
  a <- coxstuff(x, y, ic, offset = offset)
  nf <- a$nf
  fail.times <- a$fail.times
  s <- a$s
  d <- a$d
  dd <- a$dd
  nn <- a$nn
  nno <- a$nno
  w <- rep(0, nx)
  for (i in (1:nf)) {
    w <- w + s[, i]
    oo <- (1:n)[y >= fail.times[i]]
    r <- rowSums(x[, oo, drop = F] * exp(offset[oo]))
    w <- w - (d[i]/nno[i]) * r
  }
  return(list(scor = w, coxstuff.obj = a))
}


#########################
coxstuff <- function (x, y, ic, offset = rep(0, length(y))) 
{
  fail.times <- unique(y[ic == 1])
  nf <- length(fail.times)
  n <- length(y)
  nn <- rep(0, nf)
  nno <- rep(0, nf)
  for (i in 1:nf) {
    nn[i] <- sum(y >= fail.times[i])
    nno[i] <- sum(exp(offset)[y >= fail.times[i]])
  }
  s <- matrix(0, ncol = nf, nrow = nrow(x))
  d <- rep(0, nf)
  for (i in 1:nf) {
    o <- (1:n)[(y == fail.times[i]) & (ic == 1)]
    d[i] <- length(o)
  }
  oo <- match(y, fail.times)
  oo[ic == 0] <- NA
  oo[is.na(oo)] <- max(oo[!is.na(oo)]) + 1
  s <- t(rowsum(t(x), oo))
  if (ncol(s) > nf) {
    s <- s[, -ncol(s)]
  }
  dd <- rep(0, n)
  for (j in 1:nf) {
    dd[(y == fail.times[j]) & (ic == 1)] <- d[j]
  }
  return(list(fail.times = fail.times, s = s, d = d, dd = dd, 
              nf = nf, nn = nn, nno = nno))
}


##########################
coxvar <- function (x, y, ic, offset = rep(0, length(y)), coxstuff.obj = NULL) 
{
  nx <- nrow(x)
  n <- length(y)
  yy <- y + (ic == 0) * (1e-06)
  otag <- order(yy)
  y <- y[otag]
  ic <- ic[otag]
  x <- x[, otag, drop = F]
  offset <- offset[otag]
  if (is.null(coxstuff.obj)) {
    coxstuff.obj <- coxstuff(x, y, ic, offset = offset)
  }
  nf <- coxstuff.obj$nf
  fail.times <- coxstuff.obj$fail.times
  s <- coxstuff.obj$s
  d <- coxstuff.obj$d
  dd <- coxstuff.obj$dd
  nn <- coxstuff.obj$nn
  nno <- coxstuff.obj$nno
  x2 <- x^2
  oo <- (1:n)[y >= fail.times[1]]
  sx <- (1/nno[1]) * rowSums(x[, oo] * exp(offset[oo]))
  s <- (1/nno[1]) * rowSums(x2[, oo] * exp(offset[oo]))
  w <- d[1] * (s - sx * sx)
  for (i in 2:nf) {
    oo <- (1:n)[y >= fail.times[i - 1] & y < fail.times[i]]
    sx <- (1/nno[i]) * (nno[i - 1] * sx - rowSums(x[, oo, 
                                                    drop = F] * exp(offset[oo])))
    s <- (1/nno[i]) * (nno[i - 1] * s - rowSums(x2[, oo, 
                                                   drop = F] * exp(offset[oo])))
    w <- w + d[i] * (s - sx * sx)
  }
  return(w)
}


#######################
reg.perm.sample <- function (expr, testgroup, nperm) 
{
  obs = cor.func(expr, testgroup)$tt[, 1]
  perms = matrix(NA, nrow(expr), nperm)
  perms.mtx = matrix(NA, nrow = nperm, ncol = length(testgroup))
  for (i in 1:nperm) {
    perms.mtx[i, ] = sample(1:length(testgroup), size = length(testgroup))
  }
  for (t1 in 1:nperm) {
    perms[, t1] = cor.func(expr, testgroup[perms.mtx[t1, 
                                                     ]])$tt[, 1]
  }
  rownames(perms) = rownames(expr)
  colnames(perms) = paste("B", 1:nperm, sep = "")
  obs = abs(obs)
  perms = abs(perms)
  return(list(obs = obs, perms = perms, perms.mtx = perms.mtx))
}


######################
Tperm.sample <- function (x, fac, nperm) 
{
  obs = genefilter::rowttests(x, fac, tstatOnly = T)$statistic
  names(obs) = rownames(x)
  perms = matrix(NA, nrow(x), nperm)
  for (t1 in 1:nperm) {
    perms[, t1] = genefilter::rowttests(x, sample(fac), tstatOnly = T)$statistic
  }
  rownames(perms) = rownames(x)
  colnames(perms) = paste("B", 1:nperm, sep = "")
  obs = abs(obs)
  perms = abs(perms)
  return(list(obs = obs, perms = perms))
}


#######################
F.perm.sample <- function(x, fac, nperm) 
{
  obs = genefilter::rowFtests(x, fac, var.equal = TRUE)$statistic
  names(obs)=rownames(x)
  perms= matrix(NA,  nrow(x),  nperm)
  
  for(t1 in 1:nperm){
    perms[,t1]=genefilter::rowFtests(x, sample(fac), var.equal = TRUE)$statistic
  }
  rownames(perms)=rownames(x)
  colnames(perms)=paste('B',1:nperm,sep="")
  
  obs=abs(obs); perms=abs(perms)
  
  return(list(obs=obs, perms=perms))
}


######################
pqvalues.compute <- function (Stat.0, Stat.B, Stat.type, PI0 = 1) 
{
  Stat.0 = as.matrix(abs(Stat.0))
  Stat.B = as.matrix(abs(Stat.B))
  if (nrow(Stat.0) != nrow(Stat.B)) 
    stop("# of rows of Stat.0 and Stat.B should be same")
  B = ncol(Stat.B)
  G = nrow(Stat.B)
  if (Stat.type == "Tstat") {
    Stat.all = cbind(Stat.0, Stat.B)
    count.0 = apply(Stat.0, 1, function(x) sum(x <= Stat.0, 
                                               na.rm = T))
    Stat.rank.all = rank(-as.vector(Stat.all), ties.method = "max")
    pvalue.0 = as.matrix(Stat.rank.all[1:G] - count.0)/(B * 
                                                          G)
    if (is.null(PI0)) {
      PI0 = sum(pvalue.0 >= 0.5)/(0.5 * G)
    }
    else {
      PI0 = 1
    }
    qvalue.0 = PI0 * pvalue.0 * G/(apply(Stat.0, 1, function(x) sum(x <= 
                                                                      Stat.0, na.rm = T)))
    qvalue.0 = ifelse(qvalue.0 <= 1, qvalue.0, 1)
    Stat.rank = rank(-as.vector(Stat.B), ties.method = "max")
    pvalue.B = matrix(Stat.rank/(B * G), G, B)
    rownames(pvalue.B) = rownames(Stat.B)
    colnames(pvalue.B) = colnames(Stat.B)
  }
  else if (Stat.type == "Pvalue") {
    Stat.all = cbind(Stat.0, Stat.B)
    count.0 = apply(Stat.0, 1, function(x) sum(x >= Stat.0))
    Stat.rank.all = rank(as.vector(Stat.all), ties.method = "max")
    pvalue.0 = as.matrix(Stat.rank.all[1:G] - count.0)/(B * 
                                                          G)
    if (is.null(PI0)) {
      PI0 = sum(pvalue.0 >= 0.5)/(0.5 * G)
    }
    else {
      PI0 = 1
    }
    qvalue.0 = PI0 * pvalue.0 * G/(apply(Stat.0, 1, function(x) sum(x >= 
                                                                      Stat.0, na.rm = T)))
    qvalue.0 = ifelse(qvalue.0 <= 1, qvalue.0, 1)
    Stat.rank = rank(as.vector(Stat.B), ties.method = "max")
    pvalue.B = matrix(Stat.rank/(B * G), G, B)
    rownames(pvalue.B) = rownames(Stat.B)
    colnames(pvalue.B) = colnames(Stat.B)
  }
  else {
    stop("wrong stat.type")
  }
  return(list(pvalue.0 = pvalue.0, pvalue.B = pvalue.B, qvalue.0 = qvalue.0, 
              PI0 = PI0))
}


##########################
EnrichmentC_Exact<- function (madata, label, censoring = NULL, DB.matrix, DEgene.number = DEgene.number,
                              size.min = 15,size.max = 500, nperm = 500, gene.pvalues = NULL, 
                              resp.type = NULL) 
{
  if (!is.null(madata)) {
    genes.in.study = featureNames(madata)
    set2allgenes.mtx = DB.matrix
    gene.common = intersect(featureNames(madata), colnames(set2allgenes.mtx))
    if (nrow(set2allgenes.mtx) > 1) {
      set2allgenes.mtx = as.matrix(set2allgenes.mtx[, gene.common])
    }
    else {
      set2allgenes.mtx = t(as.matrix(set2allgenes.mtx[, 
                                                      gene.common]))
    }
    madata = madata[gene.common, ]
    if (resp.type == "twoclass") {
      tstat = genefilter::rowttests(exprs(madata), as.factor(label), 
                                    tstatOnly = F)
      p.values = (tstat$p.value)
      names(p.values) = rownames(tstat)
      gene.name.sort = names(sort(p.values, decreasing = F))
    }
    else if (resp.type == "multiclass") {
      tstat = genefilter::rowFtests(exprs(madata), as.factor(label), 
                                    var.equal = TRUE)
      p.values = (tstat$p.value)
      names(p.values) = rownames(tstat)
      gene.name.sort = names(sort(p.values, decreasing = F))
    }
    else if (resp.type == "survival") {
      if (is.null(censoring)) {
        stop("Error: censoring status is null")
      }
      cox.out <- coxfunc(exprs(madata), label, censoring)
      scores <- abs(cox.out$tt)
      gene.name.sort = names(sort(scores, decreasing = T))
    }
    else if (resp.type == "continuous") {
      cor.out <- cor.func(exprs(madata), label)
      scores <- abs(cor.out$tt[, 1])
      gene.name.sort = names(sort(scores, decreasing = T))
    }
    else {
      stop("Error: Wrong input augument for resp.type")
    }
  }
  DEgene = gene.name.sort[1:DEgene.number]
  pvalue.set.0 = matrix(NA,nrow = nrow(DB.matrix),ncol = 1)
  rownames(pvalue.set.0) = rownames(DB.matrix)
  for (i in 1:nrow(DB.matrix)){
    count_table<-matrix(0,2,2)
    p_value <-NA
    ####in the gene list and in the pathway
    count_table[1,1]<-sum(DEgene %in% colnames(DB.matrix[,DB.matrix[i,]==1]))
    ####in the gene list but not in the pathway
    count_table[1,2]<-length(DEgene)-count_table[1,1]
    ####not in the gene list but in the pathway
    count_table[2,1]<-sum(genes.in.study%in% colnames(DB.matrix[,DB.matrix[i,]==1]))
    ####not in the gene list and not in the pathway
    count_table[2,2]<-length(genes.in.study)-count_table[2,1]       
    if(length(count_table)==4){
      pvalue.set.0[i,1] <- fisher.test(count_table, alternative="greater")$p}
  }
  pvalue.set.0 = pvalue.set.0[pvalue.set.0[,1]!=1,,drop=FALSE]
  qvalue.set.0 = matrix(p.adjust(pvalue.set.0[,1], "BH"),ncol = 1)
  return(list(pvalue.set.0 = pvalue.set.0,qvalue.set.0 = qvalue.set.0,DEgene = DEgene))
}


########################
CPI_Exact <- function (study, label, censoring.status, DB.matrix, size.min = 15,
                     DEgene.number = DEgene.number, size.max = 500, nperm = 500, 
                     stat = NULL, rth.value = NULL,resp.type) 
{
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  out = list()
  for (t1 in 1:length(study)) {
    madata = study[[t1]]
    testlabel = madata[[label]]
    out[[t1]] = list()
    if (resp.type == "survival") {
      censoring = madata[[censoring.status]]
    }
    out[[t1]] = EnrichmentC_Exact(madata = madata, label = testlabel,DEgene.number = DEgene.number, 
                                  censoring = censoring, DB.matrix = DB.matrix, size.min = size.min, 
                                  size.max = size.max, nperm = nperm, resp.type = resp.type)
  }
  set.common = rownames(out[[1]]$pvalue.set.0)
  for (t1 in 2:length(study)) {
    set.common = intersect(set.common, rownames(out[[t1]]$pvalue.set.0))
  }
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  pvalue.0.mtx = matrix(NA, length(set.common), length(study))
  for (t1 in 1:length(study)) {
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.set.0[set.common, 
                                                ]
  }
  genelist = list()
  for (t1 in 1:length(study)) {
    genelist[[t1]] = out[[t1]]$DEgene
  }
  rownames(pvalue.0.mtx) = set.common
  rm(out)
  
  p_value_meta = aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues
  q_value_meta = p.adjust(p_value_meta, "BH")
  summary<-data.frame(q_value_meta = q_value_meta,p_value_meta = p_value_meta,
                      p_data = pvalue.0.mtx)
  return(list(summary = summary,genelist = genelist))
}


####################
CPI_gene_KS <- function (study, label, censoring.status, DB.matrix, size.min = 15, 
                       size.max = 500, nperm = 500, stat = NULL, rth.value = NULL, 
                       resp.type) 
{if (is.null(names(study))) 
  names(study) = paste("study.", 1:length(study), sep = "")
out = list()
for (t1 in 1:length(study)) {
  madata = study[[t1]]
  testlabel = madata[[label]]
  out[[t1]] = list()
  if (resp.type == "survival") {
    censoring = madata[[censoring.status]]
  }
  out[[t1]] = Enrichment_KS_gene(madata = madata, label = testlabel, 
                                 censoring = censoring, DB.matrix = DB.matrix, size.min = size.min, 
                                 size.max = size.max, nperm = nperm, resp.type = resp.type)
}
set.common = rownames(out[[1]]$pvalue.set.0)
for (t1 in 2:length(study)) {
  set.common = intersect(set.common, rownames(out[[t1]]$pvalue.set.0))
}
if (is.null(names(study))) 
  names(study) = paste("study.", 1:length(study), sep = "")
pvalue.B.array = array(data = NA, dim = c(length(set.common), 
                                          nperm, length(study)))
dimnames(pvalue.B.array) = list(set.common, paste("perm", 
                                                  1:nperm, sep = ""), names(study))
pvalue.0.mtx = matrix(NA, length(set.common), length(study))
qvalue.0.mtx = matrix(NA, length(set.common), length(study))
for (t1 in 1:length(study)) {
  pvalue.B.array[, , t1] = out[[t1]]$pvalue.set.B[set.common,]
  pvalue.0.mtx[, t1] = out[[t1]]$pvalue.set.0[set.common,]
  qvalue.0.mtx[, t1] = out[[t1]]$qvalue.set.0[set.common,]
}
rownames(qvalue.0.mtx) = set.common
rownames(pvalue.0.mtx) = set.common
rm(out)

p_value_meta = aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues
q_value_meta = p.adjust(p_value_meta, "BH")
summary<-data.frame(q_value_meta = q_value_meta,p_value_meta = p_value_meta,
                    p_data = pvalue.0.mtx)
return(list(summary = summary))
}


#########################
CPI_sample_KS <- function (study, label, censoring.status = NULL, DB.matrix, size.min = 15, 
                           size.max = 500, nperm = 500, stat, rth.value = NULL, resp.type) 
{
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  out = list()
  for (t1 in 1:length(study)) {
    madata = study[[t1]]
    testlabel = madata[[label]]
    out[[t1]] = list()
    if (resp.type == "survival") {
      censoring = madata[[censoring.status]]
    }
    out[[t1]] = Enrichment_KS_sample(madata = madata, label = testlabel, 
                                     censoring = censoring, DB.matrix = DB.matrix, size.min = size.min, 
                                     size.max = size.max, nperm = nperm, resp.type = resp.type)
  }
  common.pathway = rownames(out[[1]]$pvalue.set.0)
  for (t1 in 1:length(study)) {
    common.pathway = intersect(common.pathway, rownames(out[[t1]]$pvalue.set.0))
  }
  pvalue.B.array = array(data = NA, dim = c(length(common.pathway), 
                                            nperm, length(study)))
  pvalue.0.mtx = matrix(NA, length(common.pathway), length(study))
  qvalue.0.mtx = matrix(NA, length(common.pathway), length(study))
  rownames(pvalue.0.mtx) = common.pathway
  colnames(pvalue.0.mtx) = names(study)
  rownames(qvalue.0.mtx) = common.pathway
  colnames(qvalue.0.mtx) = names(study)
  dimnames(pvalue.B.array) = list(common.pathway, paste("perm", 
                                                        1:nperm, sep = ""), names(study))
  for (t1 in 1:length(study)) {
    pvalue.B.array[, , t1] = out[[t1]]$pvalue.set.B[common.pathway, 
                                                    ]
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.set.0[common.pathway, 
                                                ]
    qvalue.0.mtx[, t1] = out[[t1]]$qvalue.set.0[common.pathway, 
                                                ]
  }
  
  p_value_meta = aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues
  q_value_meta = p.adjust(p_value_meta, "BH") 
  
  summary<-data.frame(q_value_meta = q_value_meta,p_value_meta = p_value_meta,
                      p_data = pvalue.0.mtx)
  
  return(list(summary = summary))
}


############################
EnrichmentM_Exact<- function (madata, label, censoring = NULL, DB.matrix = DB.matrix, 
                              DEgene.number = DEgene.number,
                              size.min = 15,size.max = 500, nperm = 500, gene.pvalues = NULL, 
                              resp.type = NULL,gene.name.sort,genes.in.study) 
{
  if (!is.null(madata)) {
    genes.in.study = featureNames(madata)
    set2allgenes.mtx = DB.matrix
    gene.common = intersect(featureNames(madata), colnames(set2allgenes.mtx))
    if (nrow(set2allgenes.mtx) > 1) {
      set2allgenes.mtx = as.matrix(set2allgenes.mtx[, gene.common])
    }
    else {
      set2allgenes.mtx = t(as.matrix(set2allgenes.mtx[, 
                                                      gene.common]))
    }
    madata = madata[gene.common, ]
    if (resp.type == "twoclass") {
      tstat = genefilter::rowttests(exprs(madata), as.factor(label), 
                                    tstatOnly = F)
      p.values = (tstat$p.value)
      names(p.values) = rownames(tstat)
      gene.name.sort = names(sort(p.values, decreasing = F))
    }
    else if (resp.type == "multiclass") {
      tstat = genefilter::rowFtests(exprs(madata), as.factor(label), 
                                    var.equal = TRUE)
      p.values = (tstat$p.value)
      names(p.values) = rownames(tstat)
      gene.name.sort = names(sort(p.values, decreasing = F))
    }
    else if (resp.type == "survival") {
      if (is.null(censoring)) {
        stop("Error: censoring status is null")
      }
      cox.out <- coxfunc(exprs(madata), label, censoring)
      scores <- abs(cox.out$tt)
      gene.name.sort = names(sort(scores, decreasing = T))
    }
    else if (resp.type == "continuous") {
      cor.out <- cor.func(exprs(madata), label)
      scores <- abs(cor.out$tt[, 1])
      gene.name.sort = names(sort(scores, decreasing = T))
    }
    else {
      stop("Error: Wrong input augument for resp.type")
    }
  }
  DEgene = gene.name.sort[1:DEgene.number]
  pvalue.set.0 = matrix(NA,nrow = nrow(DB.matrix),ncol = 1)
  rownames(pvalue.set.0) = rownames(DB.matrix)
  for (i in 1:nrow(DB.matrix)){
    count_table<-matrix(0,2,2)
    p_value <-NA
    ####in the gene list and in the pathway
    count_table[1,1]<-sum(DEgene %in% colnames(DB.matrix[,DB.matrix[i,]==1]))
    ####in the gene list but not in the pathway
    count_table[1,2]<-length(DEgene)-count_table[1,1]
    ####not in the gene list but in the pathway
    count_table[2,1]<-sum(genes.in.study%in% colnames(DB.matrix[,DB.matrix[i,]==1]))
    ####not in the gene list and not in the pathway
    count_table[2,2]<-length(genes.in.study)-count_table[2,1]       
    if(length(count_table)==4){
      pvalue.set.0[i,1] <- fisher.test(count_table, alternative="greater")$p}
  }
  pvalue.set.0 = pvalue.set.0[pvalue.set.0[,1]!=1,,drop=FALSE]
  qvalue.set.0 = matrix(p.adjust(pvalue.set.0[,1], "BH"),ncol = 1)
  rownames(qvalue.set.0) = rownames(pvalue.set.0)
  return(list(pvalue.set.0 = pvalue.set.0,qvalue.set.0 = qvalue.set.0,DEgene = DEgene))
  #      qvalue.set.0 = matrix(p.adjust(pvalue.set.0[,1], "BH"),ncol = 1)
  #      rownames(qvalue.set.0) = rownames(DB.matrix)
  #      return(list(pvalue.set.0 = pvalue.set.0,qvalue.set.0 = qvalue.set.0,
  #                  pvalue.set.B = pvalue.set.B))
}


######################
MAPE_P_Exact<-function (study, label, censoring.status, DB.matrix, size.min = 15,
                        DEgene.number = DEgene.number, size.max = 500, nperm = 500, 
                        stat = NULL, rth.value = NULL,resp.type) 
{
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  out = list()
  for (t1 in 1:length(study)) {
    madata = study[[t1]]
    testlabel = madata[[label]]
    out[[t1]] = list()
    if (resp.type == "survival") {
      censoring = madata[[censoring.status]]
    }
    out[[t1]] = EnrichmentM_Exact(madata = madata, label = testlabel, DEgene.number = DEgene.number,
                                  censoring = censoring, DB.matrix = DB.matrix, size.min = size.min, 
                                  size.max = size.max, nperm = nperm, resp.type = resp.type)
  }
  set.common = rownames(out[[1]]$pvalue.set.0)
  for (t1 in 2:length(study)) {
    set.common = intersect(set.common, rownames(out[[t1]]$pvalue.set.0))
  }
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  pvalue.0.mtx = matrix(NA, length(set.common), length(study))
  qvalue.0.mtx = matrix(NA, length(set.common), length(study))
  for (t1 in 1:length(study)) {
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.set.0[set.common, ]
    qvalue.0.mtx[, t1] = out[[t1]]$qvalue.set.0[set.common,]
  }
  rownames(qvalue.0.mtx) = set.common
  rownames(pvalue.0.mtx) = set.common
  rm(out)
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,n,1)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,1,n)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,rth.value,(n - rth.value + 1))
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "Fisher") {
    DF = 2 * length(study)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) (-2 * sum(log(x)))))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pchisq(P.0,DF,lower.tail = F)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "AW Fisher") {
    pvalue.meta = matrix(aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues,ncol = 1)
    rownames(pvalue.meta) = set.common
    qvalue.meta = p.adjust(pvalue.meta, "BH")
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  qvalue.meta = matrix(qvalue.meta,ncol = 1)
  rownames(qvalue.meta) = rownames(pvalue.meta)
  return(list(pvalue.meta = pvalue.meta, qvalue.meta = qvalue.meta, 
              qvalue.set.study = qvalue.0.mtx, pvalue.set.study = pvalue.0.mtx))
}


######################
MAPE_G_Exact = function (study, label, censoring.status, DB.matrix, size.min = 15, 
                         DEgene.number = DEgene.number,size.max = 500, nperm = 500, 
                         stat = NULL,rth.value = NULL,resp.type) 
{
  gene.common = featureNames(study[[1]])
  for (t1 in 2:length(study)) {
    gene.common = intersect(gene.common, featureNames(study[[t1]]))
  }
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  madata = list()
  for (t1 in 1:length(study)) {
    madata[[t1]] = study[[t1]][gene.common, ]
  }
  Tstat.p = list()
  out = list()
  for (t1 in 1:length(study)) {
    Tstat.p[[t1]] = list()
    x = exprs(madata[[t1]])
    testlabel = madata[[t1]][[label]]
    if (resp.type == "twoclass") {
      tstat = genefilter::rowttests(x, as.factor(testlabel), 
                                    tstatOnly = F)
      Tstat.p[[t1]] = (tstat$p.value)
    }
    else if (resp.type == "multiclass") {
      tstat = genefilter::rowFtests(x, as.factor(testlabel), 
                                    var.equal = TRUE)
      Tstat.p[[t1]] = (tstat$p.value)
    }
    else if (resp.type == "survival") {
      if (is.null(censoring.status)) {
        stop("Error: censoring.status is null")
      }
      censoring = madata[[t1]][[censoring.status]]
      Tstat.p[[t1]] = cox.perm.sample(expr = x, testgroup = testlabel, 
                                      censoring = censoring, nperm = 500)
    }
    else if (resp.type == "continuous") {
      Tstat.p[[t1]] = reg.perm.sample(expr = x, testgroup = testlabel, 
                                      nperm = 500)
    }
    else {
      stop("Error: Wrong input augument for resp.type")
    }
    out[[t1]] = list()
    if(resp.type %in% c("survival" , "continous")){
      out[[t1]] = pqvalues.compute(Tstat.p[[t1]]$obs, Tstat.p[[t1]]$perms, 
                                   Stat.type = "Tstat")
    }
    else {
      out[[t1]]$pvalue.0 = Tstat.p[[t1]]
      out[[t1]]$qvalue.0 = p.adjust(Tstat.p[[t1]],"BH")
    }
  }
  pvalue.0.mtx = matrix(NA, length(gene.common), length(study))
  qvalue.0.mtx = matrix(NA, length(gene.common), length(study))
  for (t1 in 1:length(study)) {
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.0
    qvalue.0.mtx[, t1] = out[[t1]]$qvalue.0
  }
  rownames(qvalue.0.mtx) = gene.common
  rownames(pvalue.0.mtx) = gene.common
  rm(out)
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,n,1)
  }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,1,n)
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,rth.value,(n - rth.value + 1))
  }
  else if (stat == "Fisher") {
    DF = 2 * length(study)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) (-2 * sum(log(x)))))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pchisq(P.0,DF,lower.tail = F)
  }
  else if (stat == "AW Fisher") {
    pvalue.meta = matrix(aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues,ncol = 1)
    rownames(pvalue.meta) = gene.common
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  gene.pvalues = as.vector(pvalue.meta)
  names(gene.pvalues) = gene.common
  gene.qvalues = p.adjust(gene.pvalues,"BH")
  gene.name.sort = names(sort(gene.pvalues, decreasing = F))
  meta.set = EnrichmentM_Exact(madata = NULL, label = NULL, DEgene.number = DEgene.number,
                               DB.matrix = DB.matrix, size.min = size.min, size.max = size.max, 
                               nperm = nperm, gene.pvalues = gene.pvalues,
                               gene.name.sort = gene.name.sort ,genes.in.study = gene.common)
  return(list(pvalue.meta = meta.set$pvalue.set.0, qvalue.meta = meta.set$qvalue.set.0, 
              gene.meta.qvalues = gene.qvalues, 
              gene.meta.pvalues = gene.pvalues, gene.study.qvalues = qvalue.0.mtx, 
              gene.study.pvalues = pvalue.0.mtx))
}


#########################
MAPE_P_gene_KS = function (study, label, censoring.status, DB.matrix, size.min = 15, 
                           size.max = 500, nperm = 500, stat = NULL, rth.value = NULL, 
                           resp.type) 
{
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  out = list()
  for (t1 in 1:length(study)) {
    madata = study[[t1]]
    testlabel = madata[[label]]
    out[[t1]] = list()
    if (resp.type == "survival") {
      censoring = madata[[censoring.status]]
    }
    out[[t1]] = Enrichment_KS_gene(madata = madata, label = testlabel, 
                                   censoring = censoring, DB.matrix = DB.matrix, size.min = size.min, 
                                   size.max = size.max, nperm = nperm, resp.type = resp.type)
  }
  set.common = rownames(out[[1]]$pvalue.set.0)
  for (t1 in 2:length(study)) {
    set.common = intersect(set.common, rownames(out[[t1]]$pvalue.set.0))
  }
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  pvalue.B.array = array(data = NA, dim = c(length(set.common), 
                                            nperm, length(study)))
  dimnames(pvalue.B.array) = list(set.common, paste("perm", 
                                                    1:nperm, sep = ""), names(study))
  pvalue.0.mtx = matrix(NA, length(set.common), length(study))
  qvalue.0.mtx = matrix(NA, length(set.common), length(study))
  for (t1 in 1:length(study)) {
    pvalue.B.array[, , t1] = out[[t1]]$pvalue.set.B[set.common,]
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.set.0[set.common, ]
    qvalue.0.mtx[, t1] = out[[t1]]$qvalue.set.0[set.common,]
  }
  rownames(qvalue.0.mtx) = set.common
  rownames(pvalue.0.mtx) = set.common
  rm(out)
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), max)
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), min)
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), function(x) sort(x)[rth.value])
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else if (stat == "Fisher") {
    DF = 2 * length(study)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) pchisq(-2 * 
                                                                sum(log(x)), DF, lower.tail = T)))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), function(x) pchisq(-2 * 
                                                              sum(log(x)), DF, lower.tail = T))
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  meta.out = pqvalues.compute(P.0, P.B, Stat.type = "Pvalue")
  return(list(pvalue.meta = meta.out$pvalue.0, qvalue.meta = meta.out$qvalue.0, 
              pvalue.meta.B = meta.out$pvalue.B, qvalue.set.study = qvalue.0.mtx, 
              pvalue.set.study = pvalue.0.mtx))
}


########################
MAPE_G_gene_KS = function (study, label, censoring.status, DB.matrix, size.min = 15, 
                           size.max = 500, nperm = 500, stat = NULL, rth.value = NULL, 
                           resp.type) 
{
  gene.common = featureNames(study[[1]])
  for (t1 in 2:length(study)) {
    gene.common = intersect(gene.common, featureNames(study[[t1]]))
  }
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  madata = list()
  for (t1 in 1:length(study)) {
    madata[[t1]] = study[[t1]][gene.common, ]
  }
  Tstat.perm = list()
  out = list()
  for (t1 in 1:length(study)) {
    Tstat.perm[[t1]] = list()
    x = exprs(madata[[t1]])
    testlabel = madata[[t1]][[label]]
    if (resp.type == "twoclass") {
      Tstat.perm[[t1]] = Tperm.sample(x = x, fac = as.factor(testlabel), 
                                      nperm = nperm)
    }
    else if (resp.type == "multiclass") {
      Tstat.perm[[t1]] = F.perm.sample(x = x, fac = as.factor(testlabel), 
                                       nperm = nperm)
    }
    else if (resp.type == "survival") {
      if (is.null(censoring.status)) {
        stop("Error: censoring.status is null")
      }
      censoring = madata[[t1]][[censoring.status]]
      Tstat.perm[[t1]] = cox.perm.sample(expr = x, testgroup = testlabel, 
                                         censoring = censoring, nperm = nperm)
    }
    else if (resp.type == "continuous") {
      Tstat.perm[[t1]] = reg.perm.sample(expr = x, testgroup = testlabel, 
                                         nperm = nperm)
    }
    else {
      stop("Error: Wrong input augument for resp.type")
    }
    out[[t1]] = list()
    out[[t1]] = pqvalues.compute(Tstat.perm[[t1]]$obs, Tstat.perm[[t1]]$perms, 
                                 Stat.type = "Tstat")
  }
  pvalue.B.array = array(data = NA, dim = c(length(gene.common), nperm, length(study)))
  dimnames(pvalue.B.array) = list(gene.common, paste("perm", 
                                                     1:nperm, sep = ""), names(study))
  pvalue.0.mtx = matrix(NA, length(gene.common), length(study))
  qvalue.0.mtx = matrix(NA, length(gene.common), length(study))
  for (t1 in 1:length(study)) {
    pvalue.B.array[, , t1] = out[[t1]]$pvalue.B
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.0
    qvalue.0.mtx[, t1] = out[[t1]]$qvalue.0
  }
  rownames(qvalue.0.mtx) = gene.common
  rownames(pvalue.0.mtx) = gene.common
  rm(out)
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), max)
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), min)
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), function(x) sort(x)[rth.value])
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else if (stat == "Fisher") {
    DF = 2 * length(study)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) pchisq(-2 * 
                                                                sum(log(x)), DF, lower.tail = T)))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), function(x) pchisq(-2 * 
                                                              sum(log(x)), DF, lower.tail = T))
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  meta.out = pqvalues.compute(P.0, P.B, Stat.type = "Pvalue")
  gene.qvalues = as.vector(meta.out$qvalue.0)
  names(gene.qvalues) = rownames(meta.out$qvalue.0)
  gene.pvalues = as.vector(meta.out$pvalue.0)
  names(gene.pvalues) = rownames(meta.out$pvalue.0)
  meta.set = Enrichment_KS_gene(madata = NULL, label = NULL, 
                                DB.matrix = DB.matrix, size.min = size.min, size.max = size.max, 
                                nperm = nperm, gene.pvalues = gene.qvalues)
  return(list(pvalue.meta = meta.set$pvalue.set.0, qvalue.meta = meta.set$qvalue.set.0, 
              pvalue.meta.B = meta.set$pvalue.set.B, gene.meta.qvalues = gene.qvalues, 
              gene.meta.pvalues = gene.pvalues, gene.study.qvalues = qvalue.0.mtx, 
              gene.study.pvalues = pvalue.0.mtx))
}


######################
MAPE_G_sample_KS = function (study, label, censoring.status = NULL, DB.matrix, size.min = 15, 
                             size.max = 500, nperm = 500, stat, rth.value = NULL, resp.type) 
{
  gene.common = featureNames(study[[1]])
  for (t1 in 2:length(study)) {
    gene.common = intersect(gene.common, featureNames(study[[t1]]))
  }
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  madata = list()
  for (t1 in 1:length(study)) {
    madata[[t1]] = study[[t1]][gene.common, ]
  }
  Tstat.perm = list()
  out = list()
  for (t1 in 1:length(study)) {
    Tstat.perm[[t1]] = list()
    x = exprs(madata[[t1]])
    testlabel = madata[[t1]][[label]]
    if (resp.type == "twoclass") {
      Tstat.perm[[t1]] = Tperm.sample(x = x, fac = as.factor(testlabel), 
                                      nperm = nperm)
    }
    else if (resp.type == "multiclass") {
      Tstat.perm[[t1]] = F.perm.sample(x = x, fac = as.factor(testlabel), 
                                       nperm = nperm)
    }
    else if (resp.type == "survival") {
      if (is.null(censoring.status)) {
        stop("Error: censoring.aus is null")
      }
      censoring = madata[[t1]][[censoring.status]]
      Tstat.perm[[t1]] = cox.perm.sample(expr = x, testgroup = testlabel, 
                                         censoring = censoring, nperm = nperm)
    }
    else if (resp.type == "continuous") {
      Tstat.perm[[t1]] = reg.perm.sample(expr = x, testgroup = testlabel, 
                                         nperm = nperm)
    }
    else {
      stop("Error: Wrong input augument for resp.type")
    }
    out[[t1]] = list()
    out[[t1]] = pqvalues.compute(Tstat.perm[[t1]]$obs, Tstat.perm[[t1]]$perms, 
                                 Stat.type = "Tstat")
  }
  pvalue.B.array = array(data = NA, dim = c(length(gene.common), 
                                            nperm, length(study)))
  dimnames(pvalue.B.array) = list(gene.common, paste("perm", 
                                                     1:nperm, sep = ""), names(study))
  pvalue.0.mtx = matrix(NA, length(gene.common), length(study))
  qvalue.0.mtx = matrix(NA, length(gene.common), length(study))
  for (t1 in 1:length(study)) {
    pvalue.B.array[, , t1] = out[[t1]]$pvalue.B
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.0
    qvalue.0.mtx[, t1] = out[[t1]]$qvalue.0
  }
  rownames(qvalue.0.mtx) = gene.common
  rownames(pvalue.0.mtx) = gene.common
  rm(out)
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), max)
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), min)
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), function(x) sort(x)[rth.value])
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else if (stat == "Fisher") {
    DF = 2 * length(study)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) pchisq(-2 * 
                                                                sum(log(x)), DF, lower.tail = T)))
    rownames(P.0) = rownames(qvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), function(x) pchisq(-2 * 
                                                              sum(log(x)), DF, lower.tail = T))
    rownames(P.B) = rownames(qvalue.0.mtx)
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  colnames(P.0) = "perm0"
  colnames(P.B) = paste("perm", 1:ncol(P.B), sep = "")
  meta.out = pqvalues.compute(P.0, P.B, Stat.type = "Pvalue")
  gene.qvalues = as.vector(meta.out$qvalue.0)
  names(gene.qvalues) = rownames(meta.out$qvalue.0)
  gene.pvalues = as.vector(meta.out$pvalue.0)
  names(gene.pvalues) = rownames(meta.out$pvalue.0)
  score.0 = P.0
  score.B = P.B
  genes.in.study = rownames(score.0)
  set2allgenes.mtx = DB.matrix
  gene.common = intersect(genes.in.study, colnames(set2allgenes.mtx))
  set2allgenes.mtx = set2allgenes.mtx[, gene.common, drop = F]
  score.0 = score.0[gene.common, , drop = F]
  score.B = score.B[gene.common, , drop = F]
  set.length = apply(set2allgenes.mtx, 1, sum)
  idx.1 = which(set.length >= size.min)
  idx.2 = which(set.length <= size.max)
  set.idx = intersect(idx.1, idx.2)
  if (length(set.idx) < 1) {
    stop("no gene sets satisfying size.min<=size<=size.max")
  }
  else {
    set2allgenes.mtx = set2allgenes.mtx[set.idx, , drop = F]
  }
  gene.name.sort = names(sort(abs(score.0[, 1]), decreasing = F))
  order.mtx.1 = set2allgenes.mtx[, gene.name.sort, drop = F]
  order.mtx.0 = (1 - order.mtx.1)
  order.mtx.1 = t(apply(order.mtx.1, 1, function(x) x/sum(x)))
  order.mtx.0 = -t(apply(order.mtx.0, 1, function(x) x/sum(x)))
  order.mtx = order.mtx.0 + order.mtx.1
  ES.0 = as.matrix(apply(t(apply(order.mtx, 1, cumsum)), 1, 
                         max))
  ES.B = matrix(NA, nrow(ES.0), ncol(score.B))
  for (t1 in 1:ncol(score.B)) {
    gene.name.sort = names(sort(abs(score.B[, t1]), decreasing = F))
    order.mtx.1 = set2allgenes.mtx[, gene.name.sort, drop = F]
    order.mtx.0 = (1 - order.mtx.1)
    order.mtx.1 = t(apply(order.mtx.1, 1, function(x) x/sum(x)))
    order.mtx.0 = -t(apply(order.mtx.0, 1, function(x) x/sum(x)))
    order.mtx = order.mtx.0 + order.mtx.1
    order.cumsum = t(apply(order.mtx, 1, cumsum))
    ES.B[, t1] = apply(order.cumsum, 1, max)
  }
  rownames(ES.B) = rownames(order.mtx)
  N.X = apply(set2allgenes.mtx, 1, sum)
  N.Y = ncol(set2allgenes.mtx) - N.X
  N = N.X * N.Y/(N.X + N.Y)
  ES.pval.0 = exp(-2 * N * ES.0^2)
  ES.pval.B = exp(-2 * N * ES.B^2)
  enrich.out = pqvalues.compute(ES.pval.0, ES.pval.B, Stat.type = "Pvalue")
  colnames(enrich.out$pvalue.0) = "MAPE_G_sample"
  colnames(enrich.out$qvalue.0) = "MAPE_G_sample"
  colnames(enrich.out$pvalue.B) = paste("perm", 1:ncol(enrich.out$pvalue.B), 
                                        sep = "")
  return(list(pvalue.meta = enrich.out$pvalue.0, qvalue.meta = enrich.out$qvalue.0, 
              pvalue.meta.B = enrich.out$pvalue.B))
}


#######################
MAPE_P_sample_KS <- function (study, label, censoring.status = NULL, DB.matrix, size.min = 15, 
                              size.max = 500, nperm = 500, stat, rth.value = NULL, resp.type) 
{
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  out = list()
  for (t1 in 1:length(study)) {
    madata = study[[t1]]
    testlabel = madata[[label]]
    out[[t1]] = list()
    if (resp.type == "survival") {
      censoring = madata[[censoring.status]]
    }
    out[[t1]] = Enrichment_KS_sample(madata = madata, label = testlabel, 
                                     censoring = censoring, DB.matrix = DB.matrix, size.min = size.min, 
                                     size.max = size.max, nperm = nperm, resp.type = resp.type)
  }
  common.pathway = rownames(out[[1]]$pvalue.set.0)
  for (t1 in 1:length(study)) {
    common.pathway = intersect(common.pathway, rownames(out[[t1]]$pvalue.set.0))
  }
  pvalue.B.array = array(data = NA, dim = c(length(common.pathway), 
                                            nperm, length(study)))
  pvalue.0.mtx = matrix(NA, length(common.pathway), length(study))
  qvalue.0.mtx = matrix(NA, length(common.pathway), length(study))
  rownames(pvalue.0.mtx) = common.pathway
  colnames(pvalue.0.mtx) = names(study)
  rownames(qvalue.0.mtx) = common.pathway
  colnames(qvalue.0.mtx) = names(study)
  dimnames(pvalue.B.array) = list(common.pathway, paste("perm", 
                                                        1:nperm, sep = ""), names(study))
  for (t1 in 1:length(study)) {
    pvalue.B.array[, , t1] = out[[t1]]$pvalue.set.B[common.pathway, 
                                                    ]
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.set.0[common.pathway, 
                                                ]
    qvalue.0.mtx[, t1] = out[[t1]]$qvalue.set.0[common.pathway, 
                                                ]
  }
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    rownames(P.0) = rownames(pvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), max)
    rownames(P.B) = rownames(pvalue.0.mtx)
  }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    rownames(P.0) = rownames(pvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), min)
    rownames(P.B) = rownames(pvalue.0.mtx)
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    rownames(P.0) = rownames(pvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), function(x) sort(x)[rth.value])
    rownames(P.B) = rownames(pvalue.0.mtx)
  }
  else if (stat == "Fisher") {
    DF = 2 * length(study)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) pchisq(-2 * 
                                                                sum(log(x)), DF, lower.tail = T)))
    rownames(P.0) = rownames(pvalue.0.mtx)
    P.B = apply(pvalue.B.array, c(1, 2), function(x) pchisq(-2 * 
                                                              sum(log(x)), DF, lower.tail = T))
    rownames(P.B) = rownames(pvalue.0.mtx)
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  colnames(P.0) = "perm0"
  colnames(P.B) = paste("perm", 1:ncol(P.B), sep = "")
  meta.out = pqvalues.compute(P.0, P.B, Stat.type = "Pvalue")
  colnames(meta.out$pvalue.0) = "MAPE_P_sample"
  colnames(meta.out$qvalue.0) = "MAPE_P_sample"
  return(list(pvalue.meta = meta.out$pvalue.0, qvalue.meta = meta.out$qvalue.0, 
              pvalue.meta.B = meta.out$pvalue.B, pvalue.set.study = pvalue.0.mtx, 
              qvalue.set.study = qvalue.0.mtx))
}


######################
MAPE_I_KS <- function (MAP_GENE.obj, MAP_SET.obj, study) 
{
  set.common = intersect(rownames(MAP_GENE.obj$qvalue.meta), 
                         rownames(MAP_SET.obj$qvalue.meta))
  nperm = ncol(MAP_SET.obj$pvalue.meta.B)
  pvalue.set.B.array = array(data = NA, dim = c(length(set.common), 
                                                nperm, 2))
  dimnames(pvalue.set.B.array) = list(set.common, paste("perm", 
                                                        1:nperm, sep = ""), c(1,2))
  pvalue.set.0.mtx = matrix(NA, length(set.common), 2)
  pvalue.set.B.array[, , 1] = MAP_SET.obj$pvalue.meta.B[set.common,]
  pvalue.set.B.array[, , 2] = MAP_GENE.obj$pvalue.meta.B[set.common,]
  pvalue.set.0.mtx[, 1] = MAP_SET.obj$pvalue.meta[set.common,]
  pvalue.set.0.mtx[, 2] = MAP_GENE.obj$pvalue.meta[set.common,]
  rownames(pvalue.set.0.mtx) = set.common
  minP.0 = as.matrix(apply(pvalue.set.0.mtx, 1, min))
  rownames(minP.0) = set.common
  minP.B = apply(pvalue.set.B.array, c(1, 2), min)
  rownames(minP.B) = set.common
  meta.out = pqvalues.compute(minP.0, minP.B, Stat.type = "Pvalue")
  return(list(pvalue.meta = meta.out$pvalue.0, qvalue.meta = meta.out$qvalue.0))
}


########################
MAPE_I <- function (MAP_GENE.obj, MAP_SET.obj, study) 
{
  set.common = intersect(rownames(MAP_GENE.obj$qvalue.meta), 
                         rownames(MAP_SET.obj$qvalue.meta))
  
  pvalue.set.0.mtx = matrix(NA, length(set.common), 2)
  pvalue.set.0.mtx[, 1] = MAP_SET.obj$pvalue.meta[set.common,]
  pvalue.set.0.mtx[, 2] = MAP_GENE.obj$pvalue.meta[set.common,]
  rownames(pvalue.set.0.mtx) = set.common
  minP.0 = as.matrix(apply(pvalue.set.0.mtx, 1, min))
  rownames(minP.0) = set.common
  pvalue.meta = pbeta(minP.0,1,2)
  rownames(pvalue.meta) = rownames(minP.0)
  qvalue.meta = p.adjust(pvalue.meta,"BH")
  qvalue.meta = matrix(qvalue.meta,ncol = 1)
  rownames(qvalue.meta) = rownames(pvalue.meta)
  return(list(pvalue.meta = pvalue.meta, qvalue.meta = qvalue.meta))
}



########################
Enrichment_KS_gene <- function (madata, label, censoring = NULL, DB.matrix, size.min = 15, 
                                size.max = 500, nperm = 500, gene.pvalues = NULL, resp.type = NULL) 
{
  if (!is.null(madata)) {
    genes.in.study = featureNames(madata)
    set2allgenes.mtx = DB.matrix
    gene.common = intersect(featureNames(madata), colnames(set2allgenes.mtx))
    if (nrow(set2allgenes.mtx) > 1) {
      set2allgenes.mtx = as.matrix(set2allgenes.mtx[, gene.common])
    }
    else {
      set2allgenes.mtx = t(as.matrix(set2allgenes.mtx[, 
                                                      gene.common]))
    }
    madata = madata[gene.common, ]
    if (resp.type == "twoclass") {
      tstat = genefilter::rowttests(exprs(madata), as.factor(label), 
                                    tstatOnly = F)
      p.values = (tstat$p.value)
      names(p.values) = rownames(tstat)
      gene.name.sort = names(sort(p.values, decreasing = F))
    }
    else if (resp.type == "multiclass") {
      tstat = genefilter::rowFtests(exprs(madata), as.factor(label), 
                                    var.equal = TRUE)
      p.values = (tstat$p.value)
      names(p.values) = rownames(tstat)
      gene.name.sort = names(sort(p.values, decreasing = F))
    }
    else if (resp.type == "survival") {
      if (is.null(censoring)) {
        stop("Error: censoring status is null")
      }
      cox.out <- coxfunc(exprs(madata), label, censoring)
      scores <- abs(cox.out$tt)
      gene.name.sort = names(sort(scores, decreasing = T))
    }
    else if (resp.type == "continuous") {
      cor.out <- cor.func(exprs(madata), label)
      scores <- abs(cor.out$tt[, 1])
      gene.name.sort = names(sort(scores, decreasing = T))
    }
    else {
      stop("Error: Wrong input augument for resp.type")
    }
  }
  if (!is.null(gene.pvalues)) {
    if ((!is.vector(gene.pvalues)) | (is.null(names(gene.pvalues)))) 
      stop("gene.pvalues should be a vector with gene names")
    genes.in.study = names(gene.pvalues)
    set2allgenes.mtx = DB.matrix
    gene.common = intersect(genes.in.study, colnames(set2allgenes.mtx))
    if (nrow(set2allgenes.mtx) > 1) {
      set2allgenes.mtx = as.matrix(set2allgenes.mtx[, gene.common])
    }
    else {
      set2allgenes.mtx = t(as.matrix(set2allgenes.mtx[, 
                                                      gene.common]))
    }
    gene.pvalues = gene.pvalues[gene.common]
    gene.name.sort = names(sort(gene.pvalues, decreasing = F))
  }
  set.length = apply(set2allgenes.mtx, 1, sum)
  idx.1 = which(set.length >= size.min)
  idx.2 = which(set.length <= size.max)
  set.idx = intersect(idx.1, idx.2)
  if (length(set.idx) < 1) 
    stop("no gene sets satisfying size.min<=size<=size.max")
  if (length(set.idx) > 1) {
    set2allgenes.mtx = set2allgenes.mtx[set.idx, ]
  }
  else {
    set2allgenes.mtx = t(as.matrix(set2allgenes.mtx[set.idx, 
                                                    ]))
  }
  
  if (nrow(set2allgenes.mtx) > 1) {
    order.mtx.1 = (set2allgenes.mtx[, gene.name.sort])
  }
  else {
    order.mtx.1 = t(as.matrix(set2allgenes.mtx[, gene.name.sort]))
  }
  order.mtx.0 = (1 - order.mtx.1)
  order.mtx.1 = t(apply(order.mtx.1, 1, function(x) x/sum(x)))
  order.mtx.0 = -t(apply(order.mtx.0, 1, function(x) x/sum(x)))
  order.mtx = order.mtx.0 + order.mtx.1
  ES.0 = as.matrix(apply(t(apply(order.mtx, 1, cumsum)), 1, 
                         max))
  ES.B = matrix(NA, nrow(ES.0), nperm)
  for (t1 in 1:nperm) {
    if (nrow(order.mtx) > 1) {
      order.mtx.perm = order.mtx[, sample(ncol(order.mtx))]
    }
    else {
      order.mtx.perm = t(as.matrix(order.mtx[, sample(ncol(order.mtx))]))
    }
    order.cumsum = t(apply(order.mtx.perm, 1, cumsum))
    ES.B[, t1] = apply(order.cumsum, 1, max)
  }
  rownames(ES.B) = rownames(order.mtx)
  N.X = apply(set2allgenes.mtx, 1, sum)
  N.Y = ncol(set2allgenes.mtx) - N.X
  N = N.X * N.Y/(N.X + N.Y)
  enrich.out = pqvalues.compute(ES.0, ES.B, Stat.type = "Tstat")
  return(list(pvalue.set.0 = enrich.out$pvalue.0, pvalue.set.B = enrich.out$pvalue.B, 
              qvalue.set.0 = enrich.out$qvalue.0))
}


######################
Enrichment_KS_sample <- function (madata, label, censoring = NULL, DB.matrix, size.min = 15, 
                                  size.max = 500, nperm = 500, resp.type = NULL) 
{
  genes.in.study = featureNames(madata)
  set2allgenes.mtx = DB.matrix
  gene.common = intersect(featureNames(madata), colnames(set2allgenes.mtx))
  set2allgenes.mtx = as.matrix(set2allgenes.mtx[, gene.common], 
                               drop = F)
  madata = madata[gene.common, ]
  set.length = apply(set2allgenes.mtx, 1, sum)
  idx.1 = which(set.length >= size.min)
  idx.2 = which(set.length <= size.max)
  set.idx = intersect(idx.1, idx.2)
  if (length(set.idx) < 1) {
    stop("no gene sets satisfying size.min<=size<=size.max")
  }
  else {
    set2allgenes.mtx = set2allgenes.mtx[set.idx, , drop = F]
  }
  x = exprs(madata)
  testlabel = label
  if (resp.type == "twoclass") {
    Tstat.perm = Tperm.sample(x = x, fac = as.factor(testlabel), 
                              nperm = nperm)
  }
  else if (resp.type == "multiclass") {
    Tstat.perm = F.perm.sample(x = x, fac = as.factor(testlabel), 
                               nperm = nperm)
  }
  else if (resp.type == "survival") {
    if (is.null(censoring)) {
      stop("Error: censoring.status is null")
    }
    Tstat.perm = cox.perm.sample(expr = x, testgroup = testlabel, 
                                 censoring = censoring, nperm = nperm)
  }
  else if (resp.type == "continuous") {
    Tstat.perm = reg.perm.sample(expr = x, testgroup = testlabel, 
                                 nperm = nperm)
  }
  else {
    stop("Error: Wrong input augument for resp.type")
  }
  score.0 = Tstat.perm$obs
  score.B = Tstat.perm$perms
  gene.name.sort = names(sort(abs(score.0), decreasing = T))
  order.mtx.1 = set2allgenes.mtx[, gene.name.sort, drop = F]
  order.mtx.0 = (1 - order.mtx.1)
  order.mtx.1 = t(apply(order.mtx.1, 1, function(x) x/sum(x)))
  order.mtx.0 = -t(apply(order.mtx.0, 1, function(x) x/sum(x)))
  order.mtx = order.mtx.0 + order.mtx.1
  ES.0 = as.matrix(apply(t(apply(order.mtx, 1, cumsum)), 1, 
                         max))
  ES.B = matrix(NA, nrow(ES.0), ncol(score.B))
  for (t1 in 1:ncol(score.B)) {
    gene.name.sort = names(sort(abs(score.B[, t1]), decreasing = T))
    order.mtx.1 = set2allgenes.mtx[, gene.name.sort, drop = F]
    order.mtx.0 = (1 - order.mtx.1)
    order.mtx.1 = t(apply(order.mtx.1, 1, function(x) x/sum(x)))
    order.mtx.0 = -t(apply(order.mtx.0, 1, function(x) x/sum(x)))
    order.mtx = order.mtx.0 + order.mtx.1
    order.cumsum = t(apply(order.mtx, 1, cumsum))
    ES.B[, t1] = apply(order.cumsum, 1, max)
  }
  rownames(ES.B) = rownames(order.mtx)
  N.X = apply(set2allgenes.mtx, 1, sum)
  N.Y = ncol(set2allgenes.mtx) - N.X
  N = N.X * N.Y/(N.X + N.Y)
  enrich.out = pqvalues.compute(ES.0, ES.B, Stat.type = "Tstat")
  return(list(pvalue.set.0 = enrich.out$pvalue.0, pvalue.set.B = enrich.out$pvalue.B, 
              qvalue.set.0 = enrich.out$qvalue.0))
}


######################
Enrichment_KS <- function (madata, label, censoring = NULL, DB.matrix, size.min = 15, 
                                size.max = 500, gene.pvalues = NULL, resp.type = NULL) 
{
  if (!is.null(madata)) {
    genes.in.study = featureNames(madata)
    set2allgenes.mtx = DB.matrix
    gene.common = intersect(featureNames(madata), colnames(set2allgenes.mtx))
    if (nrow(set2allgenes.mtx) > 1) {
      set2allgenes.mtx = as.matrix(set2allgenes.mtx[, gene.common])
    }
    else {
      set2allgenes.mtx = t(as.matrix(set2allgenes.mtx[, 
                                                      gene.common]))
    }
    madata = madata[gene.common, ]
    if (resp.type == "twoclass") {
      tstat = genefilter::rowttests(exprs(madata), as.factor(label), 
                                    tstatOnly = F)
      p.values = (tstat$p.value)
      names(p.values) = rownames(tstat)
      gene.name.sort = names(sort(p.values, decreasing = F))
    }
    else if (resp.type == "multiclass") {
      tstat = genefilter::rowFtests(exprs(madata), as.factor(label), 
                                    var.equal = TRUE)
      p.values = (tstat$p.value)
      names(p.values) = rownames(tstat)
      gene.name.sort = names(sort(p.values, decreasing = F))
    }
    else if (resp.type == "survival") {
      if (is.null(censoring)) {
        stop("Error: censoring status is null")
      }
      cox.out <- coxfunc(exprs(madata), label, censoring)
      scores <- abs(cox.out$tt)
      gene.name.sort = names(sort(scores, decreasing = T))
    }
    else if (resp.type == "continuous") {
      cor.out <- cor.func(exprs(madata), label)
      scores <- abs(cor.out$tt[, 1])
      gene.name.sort = names(sort(scores, decreasing = T))
    }
    else {
      stop("Error: Wrong input augument for resp.type")
    }
  }
  if (!is.null(gene.pvalues)) {
    if ((!is.vector(gene.pvalues)) | (is.null(names(gene.pvalues)))) 
      stop("gene.pvalues should be a vector with gene names")
    genes.in.study = names(gene.pvalues)
    set2allgenes.mtx = DB.matrix
    gene.common = intersect(genes.in.study, colnames(set2allgenes.mtx))
    if (nrow(set2allgenes.mtx) > 1) {
      set2allgenes.mtx = as.matrix(set2allgenes.mtx[, gene.common])
    }
    else {
      set2allgenes.mtx = t(as.matrix(set2allgenes.mtx[, 
                                                      gene.common]))
    }
    gene.pvalues = gene.pvalues[gene.common]
    gene.name.sort = names(sort(gene.pvalues, decreasing = F))
  }
  set.length = apply(set2allgenes.mtx, 1, sum)
  idx.1 = which(set.length >= size.min)
  idx.2 = which(set.length <= size.max)
  set.idx = intersect(idx.1, idx.2)
  if (length(set.idx) < 1) 
    stop("no gene sets satisfying size.min<=size<=size.max")
  if (length(set.idx) > 1) {
    set2allgenes.mtx = set2allgenes.mtx[set.idx, ]
  }
  else {
    set2allgenes.mtx = t(as.matrix(set2allgenes.mtx[set.idx, 
                                                    ]))
  }
  if (nrow(set2allgenes.mtx) > 1) {
    order.mtx.1 = (set2allgenes.mtx[, gene.name.sort])
  }
  else {
    order.mtx.1 = t(as.matrix(set2allgenes.mtx[, gene.name.sort]))
  }
  order.mtx.0 = (1 - order.mtx.1)
  n_hit = rowSums(order.mtx.1)
  n_miss = rowSums(order.mtx.0)
  n_genes = ncol(order.mtx.1)
  nn = (n_hit*n_miss)/n_genes
  order.mtx.1 = t(apply(order.mtx.1, 1, function(x) x/sum(x)))
  order.mtx.0 = -t(apply(order.mtx.0, 1, function(x) x/sum(x)))
  order.mtx = order.mtx.0 + order.mtx.1
  ES.0 = as.matrix(apply(t(apply(order.mtx, 1, cumsum)), 1, 
                         max))
  pvalue.0 = matrix(exp(-2*nn*(ES.0[,1]^2)), ncol = 1)
  qvalue.0 = matrix(p.adjust(pvalue.0[,1],"BH"), ncol = 1)
  rownames(pvalue.0) = rownames(ES.0)
  rownames(qvalue.0) = rownames(ES.0)
  return(list(pvalue.set.0 = pvalue.0, qvalue.set.0 = qvalue.0))
}



##########################
MAPE_P_KS = function (study, label, censoring.status, DB.matrix, size.min = 15, 
                           size.max = 500, stat = NULL, rth.value = NULL, 
                           resp.type) 
{
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  out = list()
  for (t1 in 1:length(study)) {
    madata = study[[t1]]
    testlabel = madata[[label]]
    out[[t1]] = list()
    if (resp.type == "survival") {
      censoring = madata[[censoring.status]]
    }
    out[[t1]] = Enrichment_KS (madata = madata, label = testlabel, 
                               censoring = censoring, DB.matrix = DB.matrix, size.min = size.min, 
                               size.max = size.max, resp.type = resp.type)
  }
  set.common = rownames(out[[1]]$pvalue.set.0)
  for (t1 in 2:length(study)) {
    set.common = intersect(set.common, rownames(out[[t1]]$pvalue.set.0))
  }
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  pvalue.0.mtx = matrix(NA, length(set.common), length(study))
  qvalue.0.mtx = matrix(NA, length(set.common), length(study))
  for (t1 in 1:length(study)) {
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.set.0[set.common, ]
    qvalue.0.mtx[, t1] = out[[t1]]$qvalue.set.0[set.common,]
  }
  rownames(qvalue.0.mtx) = set.common
  rownames(pvalue.0.mtx) = set.common
  rm(out)
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,n,1)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
    }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,1,n)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,rth.value,(n - rth.value + 1))
    qvalue.meta = p.adjust(pvalue.meta,"BH")
    }
  else if (stat == "Fisher") {
    DF = 2 * length(study)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) (-2 * sum(log(x)))))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pchisq(P.0,DF,lower.tail = F)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "AW Fisher") {
    pvalue.meta = matrix(aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues,ncol = 1)
    rownames(pvalue.meta) = set.common
    qvalue.meta = p.adjust(pvalue.meta, "BH")
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  qvalue.meta = matrix(qvalue.meta,ncol = 1)
  rownames(qvalue.meta) = rownames(pvalue.meta)
  return(list(pvalue.meta = pvalue.meta, qvalue.meta = qvalue.meta, 
              qvalue.set.study = qvalue.0.mtx, pvalue.set.study = pvalue.0.mtx))
}



##########################
MAPE_G_KS = function (study, label, censoring.status, DB.matrix, size.min = 15, 
                           size.max = 500, stat = NULL, rth.value = NULL, 
                           resp.type) 
{
  gene.common = featureNames(study[[1]])
  for (t1 in 2:length(study)) {
    gene.common = intersect(gene.common, featureNames(study[[t1]]))
  }
  if (is.null(names(study))) 
    names(study) = paste("study.", 1:length(study), sep = "")
  madata = list()
  for (t1 in 1:length(study)) {
    madata[[t1]] = study[[t1]][gene.common, ]
  }
  Tstat.p = list()
  out = list()
  for (t1 in 1:length(study)) {
    Tstat.p[[t1]] = list()
    x = exprs(madata[[t1]])
    testlabel = madata[[t1]][[label]]
    if (resp.type == "twoclass") {
      tstat = genefilter::rowttests(x, as.factor(testlabel), 
                                    tstatOnly = F)
      Tstat.p[[t1]] = (tstat$p.value)
    }
    else if (resp.type == "multiclass") {
      tstat = genefilter::rowFtests(x, as.factor(testlabel), 
                                    var.equal = TRUE)
      Tstat.p[[t1]] = (tstat$p.value)
    }
    else if (resp.type == "survival") {
      if (is.null(censoring.status)) {
        stop("Error: censoring.status is null")
      }
      censoring = madata[[t1]][[censoring.status]]
      Tstat.p[[t1]] = cox.perm.sample(expr = x, testgroup = testlabel, 
                                         censoring = censoring, nperm = 500)
    }
    else if (resp.type == "continuous") {
      Tstat.p[[t1]] = reg.perm.sample(expr = x, testgroup = testlabel, 
                                         nperm = 500)
    }
    else {
      stop("Error: Wrong input augument for resp.type")
    }
    out[[t1]] = list()
    if(resp.type %in% c("survival" , "continous")){
    out[[t1]] = pqvalues.compute(Tstat.p[[t1]]$obs, Tstat.p[[t1]]$perms, 
                                 Stat.type = "Tstat")
    }
    else {
    out[[t1]]$pvalue.0 = Tstat.p[[t1]]
    out[[t1]]$qvalue.0 = p.adjust(Tstat.p[[t1]],"BH")
    }
  }
  pvalue.0.mtx = matrix(NA, length(gene.common), length(study))
  qvalue.0.mtx = matrix(NA, length(gene.common), length(study))
  for (t1 in 1:length(study)) {
    pvalue.0.mtx[, t1] = out[[t1]]$pvalue.0
    qvalue.0.mtx[, t1] = out[[t1]]$qvalue.0
  }
  rownames(qvalue.0.mtx) = gene.common
  rownames(pvalue.0.mtx) = gene.common
  rm(out)
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,n,1)
  }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,1,n)
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,rth.value,(n - rth.value + 1))
  }
  else if (stat == "Fisher") {
    DF = 2 * length(study)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) (-2 * sum(log(x)))))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pchisq(P.0,DF,lower.tail = F)
  }
  else if (stat == "AW Fisher") {
    pvalue.meta = matrix(aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues,ncol = 1)
    rownames(pvalue.meta) = gene.common
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  gene.pvalues = as.vector(pvalue.meta)
  names(gene.pvalues) = gene.common
  gene.qvalues = p.adjust(gene.pvalues,"BH")
  meta.set = Enrichment_KS(madata = NULL, label = NULL, 
                                DB.matrix = DB.matrix, size.min = size.min, size.max = size.max, 
                                gene.pvalues = gene.pvalues)
  return(list(pvalue.meta = meta.set$pvalue.set.0, qvalue.meta = meta.set$qvalue.set.0, 
              gene.meta.qvalues = gene.qvalues, 
              gene.meta.pvalues = gene.pvalues, gene.study.qvalues = qvalue.0.mtx, 
              gene.study.pvalues = pvalue.0.mtx))
}



######################
MAPE_P_KS_DE = function (ind.p = ind.p, DB.matrix, size.min = 15, gene.common,
                      size.max = 500, stat = NULL, rth.value = NULL) 
{
  out = list()
  out2 = list()
  for (t1 in 1:ncol(ind.p)) {
      gene.name.sort = names(sort(ind.p[,t1],decreasing = F))
      gene.name.sort = gene.name.sort[gene.name.sort%in%gene.common]
      gene.name.sort = toupper(gene.name.sort)
      set2allgenes.mtx = DB.matrix
      order.mtx.1 = (set2allgenes.mtx[, gene.name.sort])
      order.mtx.0 = (1 - order.mtx.1)
      n_hit = rowSums(order.mtx.1)
      n_miss = rowSums(order.mtx.0)
      n_genes = ncol(order.mtx.1)
      nn = (n_hit*n_miss)/n_genes
      order.mtx.1 = t(apply(order.mtx.1, 1, function(x) x/sum(x)))
      order.mtx.0 = -t(apply(order.mtx.0, 1, function(x) x/sum(x)))
      order.mtx = order.mtx.0 + order.mtx.1
      ES.0 = as.matrix(apply(t(apply(order.mtx, 1, cumsum)), 1, 
                             max))
      pvalue.0 = matrix(exp(-2*nn*(ES.0[,1]^2)), ncol = 1)
      rownames(pvalue.0) = rownames(ES.0)
      out[[t1]] = pvalue.0 
      out2[[t1]] = as.matrix(p.adjust(pvalue.0,"BH"),ncol = 1)
      rownames(out2[[t1]]) = rownames(pvalue.0)
    }
  set.common = rownames(out[[1]])
  pvalue.0.mtx = matrix(NA, length(set.common), ncol(ind.p))
  qvalue.0.mtx = matrix(NA, length(set.common), ncol(ind.p))
  for (t1 in 1:ncol(ind.p)) {
    pvalue.0.mtx[, t1] = out[[t1]][set.common, ]
    qvalue.0.mtx[, t1] = out2[[t1]][set.common,]
  }
  rownames(qvalue.0.mtx) = set.common
  rownames(pvalue.0.mtx) = set.common
  qvalue.0.mtx = qvalue.0.mtx[complete.cases(qvalue.0.mtx),]
  pvalue.0.mtx = pvalue.0.mtx[complete.cases(pvalue.0.mtx),]
  rm(out)
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,n,1)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,1,n)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,rth.value,(n - rth.value + 1))
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "Fisher") {
    DF = 2 * ncol(ind.p)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) (-2 * sum(log(x)))))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pchisq(P.0,DF,lower.tail = F)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "AW Fisher") {
    pvalue.meta = matrix(aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues,ncol = 1)
    rownames(pvalue.meta) = rownames(pvalue.0.mtx)
    qvalue.meta = p.adjust(pvalue.meta, "BH")
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  qvalue.meta = matrix(qvalue.meta,ncol = 1)
  rownames(qvalue.meta) = rownames(pvalue.meta)
  return(list(pvalue.meta = pvalue.meta, qvalue.meta = qvalue.meta, 
              qvalue.set.study = qvalue.0.mtx, pvalue.set.study = pvalue.0.mtx))
}

########################
MAPE_G_KS_DE = function (meta.p = meta.p,DB.matrix,gene.common,
                         size.min = 15, size.max = 500) 
{
  gene.name.sort = names(sort(meta.p[,1],decreasing = F))
  gene.name.sort = gene.name.sort[gene.name.sort%in%gene.common]
  gene.name.sort = toupper(gene.name.sort)
  set2allgenes.mtx = DB.matrix
  order.mtx.1 = (set2allgenes.mtx[, gene.name.sort])
  order.mtx.0 = (1 - order.mtx.1)
  n_hit = rowSums(order.mtx.1)
  n_miss = rowSums(order.mtx.0)
  n_genes = ncol(order.mtx.1)
  nn = (n_hit*n_miss)/n_genes
  order.mtx.1 = t(apply(order.mtx.1, 1, function(x) x/sum(x)))
  order.mtx.0 = -t(apply(order.mtx.0, 1, function(x) x/sum(x)))
  order.mtx = order.mtx.0 + order.mtx.1
  ES.0 = as.matrix(apply(t(apply(order.mtx, 1, cumsum)), 1, 
                         max))
  pvalue.0 = matrix(exp(-2*nn*(ES.0[,1]^2)), ncol = 1)
  rownames(pvalue.0) = rownames(ES.0)
  pvalue.set.0 = pvalue.0 
  qvalue.set.0 = as.matrix(p.adjust(pvalue.0,"BH"),ncol = 1)
  rownames(qvalue.set.0) = rownames(pvalue.0)
  return(list(pvalue.meta = pvalue.set.0, qvalue.meta = qvalue.set.0))
}


####################
CPI_KS <- function (study, label, censoring.status, DB.matrix, size.min = 15, 
                         size.max = 500, stat = NULL, rth.value = NULL, resp.type) 
{if (is.null(names(study))) 
  names(study) = paste("study.", 1:length(study), sep = "")
out = list()
for (t1 in 1:length(study)) {
  madata = study[[t1]]
  testlabel = madata[[label]]
  out[[t1]] = list()
  if (resp.type == "survival") {
    censoring = madata[[censoring.status]]
  }
  out[[t1]] = Enrichment_KS (madata = madata, label = testlabel, 
                             censoring = censoring, DB.matrix = DB.matrix, size.min = size.min, 
                             size.max = size.max, resp.type = resp.type)
}
set.common = rownames(out[[1]]$pvalue.set.0)
for (t1 in 2:length(study)) {
  set.common = intersect(set.common, rownames(out[[t1]]$pvalue.set.0))
}
if (is.null(names(study))) 
  names(study) = paste("study.", 1:length(study), sep = "")

pvalue.0.mtx = matrix(NA, length(set.common), length(study))
qvalue.0.mtx = matrix(NA, length(set.common), length(study))
for (t1 in 1:length(study)) {
  pvalue.0.mtx[, t1] = out[[t1]]$pvalue.set.0[set.common,]
  qvalue.0.mtx[, t1] = out[[t1]]$qvalue.set.0[set.common,]
}
rownames(qvalue.0.mtx) = set.common
rownames(pvalue.0.mtx) = set.common
rm(out)

p_value_meta = aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues
q_value_meta = p.adjust(p_value_meta, "BH")
summary<-data.frame(q_value_meta = q_value_meta,p_value_meta = p_value_meta,
                    p_data = pvalue.0.mtx)
return(list(summary = summary))
}



#########################
MAPE_P_Exact_DE = function (ind.p = ind.p, DB.matrix, size.min = 15, gene.common,
                         size.max = 500, stat = NULL, rth.value = NULL, DEgene.number) 
{
  out = list()
  out2 = list()
  for (t1 in 1:ncol(ind.p)) {
    genes.in.study = names(ind.p[,t1])
    gene.name.sort = names(sort(ind.p[,t1],decreasing = F))
    gene.name.sort = gene.name.sort[gene.name.sort%in%gene.common]
    gene.name.sort = toupper(gene.name.sort)
    DEgene = gene.name.sort[1:DEgene.number]
    pvalue.0 = matrix(NA,nrow = nrow(DB.matrix),ncol = 1)
    rownames(pvalue.0) = rownames(DB.matrix)
    for (i in 1:nrow(DB.matrix)){
      count_table<-matrix(0,2,2)
      p_value <-NA
      ####in the gene list and in the pathway
      count_table[1,1]<-sum(DEgene %in% colnames(DB.matrix[,DB.matrix[i,]==1]))
      ####in the gene list but not in the pathway
      count_table[1,2]<-length(DEgene)-count_table[1,1]
      ####not in the gene list but in the pathway
      count_table[2,1]<-sum(genes.in.study%in% colnames(DB.matrix[,DB.matrix[i,]==1]))
      ####not in the gene list and not in the pathway
      count_table[2,2]<-length(genes.in.study)-count_table[2,1]       
      if(length(count_table)==4){
        pvalue.0[i,1] <- fisher.test(count_table, alternative="greater")$p}
    }
#    pvalue.0 = pvalue.0[pvalue.0[,1]!=1,,drop=FALSE]
    qvalue.0 = matrix(p.adjust(pvalue.0[,1], "BH"),ncol = 1)
    rownames(qvalue.0) = rownames(pvalue.0)
    out[[t1]] = pvalue.0 
    out2[[t1]] = as.matrix(p.adjust(pvalue.0,"BH"),ncol = 1)
    rownames(out2[[t1]]) = rownames(pvalue.0)
  }
  set.common = rownames(out[[1]])
  pvalue.0.mtx = matrix(NA, length(set.common), ncol(ind.p))
  qvalue.0.mtx = matrix(NA, length(set.common), ncol(ind.p))
  for (t1 in 1:ncol(ind.p)) {
    pvalue.0.mtx[, t1] = out[[t1]][set.common, ]
    qvalue.0.mtx[, t1] = out2[[t1]][set.common,]
    pvalue.0.mtx = pvalue.0.mtx
  }
  rownames(qvalue.0.mtx) = set.common
  rownames(pvalue.0.mtx) = set.common
  rm(out)
  if (stat == "maxP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, max))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,n,1)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "minP") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, min))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,1,n)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "rth") {
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) sort(x)[rth.value]))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pbeta(P.0,rth.value,(n - rth.value + 1))
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "Fisher") {
    DF = 2 * ncol(ind.p)
    P.0 = as.matrix(apply(pvalue.0.mtx, 1, function(x) (-2 * sum(log(x)))))
    n = ncol(pvalue.0.mtx)
    pvalue.meta = pchisq(P.0,DF,lower.tail = F)
    qvalue.meta = p.adjust(pvalue.meta,"BH")
  }
  else if (stat == "AW Fisher") {
    pvalue.meta = matrix(aw.fisher.pvalue(pvalue.0.mtx, method="original", weight.matrix=T)$pvalues,ncol = 1)
    rownames(pvalue.meta) = set.common
    qvalue.meta = p.adjust(pvalue.meta, "BH")
  }
  else {
    stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher")
  }
  qvalue.meta = matrix(qvalue.meta,ncol = 1)
  rownames(qvalue.meta) = rownames(pvalue.meta)
  return(list(pvalue.meta = pvalue.meta, qvalue.meta = qvalue.meta, 
              qvalue.set.study = qvalue.0.mtx, pvalue.set.study = pvalue.0.mtx))
}

########################
MAPE_G_Exact_DE = function (meta.p = meta.p,DB.matrix,gene.common,
                         size.min = 15, size.max = 500,DEgene.number) 
{
  gene.name.sort = names(sort(meta.p[,1],decreasing = F))
  genes.in.study = gene.name.sort
  gene.name.sort = gene.name.sort[gene.name.sort%in%gene.common]
  gene.name.sort = toupper(gene.name.sort)
  DEgene = gene.name.sort[1:DEgene.number]
  pvalue.0 = matrix(NA,nrow = nrow(DB.matrix),ncol = 1)
  rownames(pvalue.0) = rownames(DB.matrix)
  for (i in 1:nrow(DB.matrix)){
    count_table<-matrix(0,2,2)
    p_value <-NA
    ####in the gene list and in the pathway
    count_table[1,1]<-sum(DEgene %in% colnames(DB.matrix[,DB.matrix[i,]==1]))
    ####in the gene list but not in the pathway
    count_table[1,2]<-length(DEgene)-count_table[1,1]
    ####not in the gene list but in the pathway
    count_table[2,1]<-sum(genes.in.study%in% colnames(DB.matrix[,DB.matrix[i,]==1]))
    ####not in the gene list and not in the pathway
    count_table[2,2]<-length(genes.in.study)-count_table[2,1]       
    if(length(count_table)==4){
      pvalue.0[i,1] <- fisher.test(count_table, alternative="greater")$p}
  }
#  pvalue.0 = pvalue.0[pvalue.0[,1]!=1,,drop=FALSE]
  qvalue.0 = matrix(p.adjust(pvalue.0[,1], "BH"),ncol = 1)
  rownames(qvalue.0) = rownames(pvalue.0)
  pvalue.set.0 = pvalue.0 
  qvalue.set.0 = as.matrix(p.adjust(pvalue.0,"BH"),ncol = 1)
  rownames(qvalue.set.0) = rownames(pvalue.0)
  return(list(pvalue.meta = pvalue.set.0, qvalue.meta = qvalue.set.0))
}
