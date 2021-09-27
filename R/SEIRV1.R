###################SEIR#############
     # SEIR: Selector-Embedded Iterative Regression
     # A fast and extremely powerful R package for GWAS
     # Written by Mengjin Zhu
#############################################################################################################################
     # Object: Main function to perform GWAS
     # Authors: Mengjin Zhu
     # Last update: Aug 5, 2020
########################################################################################################################
#' Selector-Embedded Iterative Regression
#'
#' @param Y Phenotypic vector with length n
#' @param X n by m dataframe, matrix or big.matrix of genotype matrix
#' @param GM m by 3 dataframe or matrix for SNP name, chromosome and BP
#' @param CV n by t dataframe or matrix for t covariates
#' @param maxStep maxStep Maximum number of iteration steps
#' @param selector selector A method for Variable selection
#' @param fun.selection  fun.selection The core function for parameter calculation, lm.vec is more faster than fastLm
#' @param extraction A method for calculating p values
#' @param X.check logical
#' @param chunk.num The "pieces" number of genotype matrix when considering memory size
#' @param file.output output GWAS results, Manhattan and QQ plot
#' @param plot.style plot.style
#' @param cutOff The threshold
#'
#' @return List - includes the suggested SNPs, the results with normal 4-column format (SNP, chromosome, BP, p-values), the matrix of coefficients, SE, z-values, p-values, and the matrices of coefficients, SE, z-values, p-values from each iterative step
#' @export

#' @examples 
#'          data(geno)
#'          data(phe)
#'          data(map)
#'          SEIR(Y=phe,X=geno,GM=map,CV=NULL,maxStep=10,
#'               selector="stepwise",fun.selection="fastLm",
#'               extraction="min",X.check="FALSE",chunk.num = NULL,
#'               file.output=FALSE,plot.style="SEIR",cutOff=0.05)
SEIR <- function(Y, X, GM = NULL, CV = NULL, maxStep = 10, selector = c("L0", "L1", "MCP", "stepwise", "bess", "MARS","SCAD", "lars","vif","SES","FBED"),
                       fun.selection = c("fastLm", "lm.vec"), extraction = c("bagging-like", "1se", "min", "mean"),
					   X.check = c(TRUE, FALSE), chunk.num = NULL,file.output=TRUE,plot.style="CMplot",cutOff=0.05)
{
    # Object: Perform GWAS with SEIR function.
    # Input: y ? Phenotypic vector with length n.
    # Input: X - n by m dataframe, matrix or big.matrix of genotype matrix.
    # Input: GM - m by 3 dataframe or matrix for SNP name, chromosome and BP.
    # Input: CV - n by t dataframe or matrix for t covariates.
    # Input: X.check ? logical. Input X is pre-treated or not.
    # Input: Chunk.num is the "pieces" number of genotype matrix when considering memory size.
    # Output: List - includes the suggested SNPs, the results with normal 4-column format (SNP, chromosome, BP, p-values), the matrix of coefficients, SE, z-values, p-values, and the matrices of coefficients, SE, z-values, p-values from each iterative step.
    # Requirement: y, X and CV have same taxa order. X and GM have the same SNP order.
    # Requirement: y's missing data(NA)will be directly removed.
    # Null GM is allowed for those species without SNP annotations.
    # fun.selection is the core function for parameter calculation, lm.vec is more faster than fastLm.
    # extraction is only used for fastLm, because lm.vec performs one-round calculation.
############################################################################
    # require(Rfast)
    # require(RcppEigen)
  SEIR.0000 <- function(){
      SEIR.Version="SEIR v1.00, May 21, 2017"
      return(SEIR.Version)
}
    print("--------------------- Welcome to SEIR ----------------------------")
    echo=TRUE
    # SEIR.Version=SEIR.0000()
    print(paste('R is using', memory.size(), 'MB out of limit', memory.limit(), 'MB'))
    print("SEIR Started...")

    if(is.null(Y)){stop("Y must be required")}
    if(is.null(X)){stop("genotype matrix must be required")}
    if (!is.matrix(X)){N=dim(X)[1]; M=dim(X)[2]; X = as.matrix(X, N, M); rm(N); rm(M)}
    if (is.null(colnames(X))){colnames(X) <- c(paste('SNP', 1:ncol(X), sep=''))}
    y = as.vector(Y[,2])
    name.of.trait=colnames(Y)[2]
    missing=is.na(y) #index for missing phenotype
    if (length(which(missing==TRUE))!=0){
    y=y[!missing]
    X<- X[!missing, ]
    if(!is.null(CV)){CV=CV[!missing,]}
   }
    rm(missing)

    if (X.check == TRUE){
    varSNP= Rfast::colVars(X)
    indexSNP=which(varSNP!=0)
    X<-X[, indexSNP]
    if (!is.null(GM)){GM<-GM[indexSNP,]}
    rm(varSNP)
    rm(indexSNP)
  }

  lm.vec <- compiler::cmpfun(function(y, X, CV){
      n = length(y)
      k = ncol(CV)
      U1 <- crossprod(CV, y)
      U2 <- solve(crossprod(CV), U1) #b.cv
      yt <- y - CV %*% U2
      err.cv <- sqrt(diag(chol2inv(Rfast::cholesky(crossprod(CV)))*as.numeric(crossprod(yt)/(n-k))))
      zvalue.cv = U2/err.cv
      p.cv = 2 * pnorm(-abs(zvalue.cv))
      U3 <- crossprod(CV, X)
      U4 <- solve(crossprod(CV), U3)
      Str <- X - CV %*% U4
      Str2 <- colSums(Str ^ 2)
      b <- as.vector(crossprod(yt, Str) / Str2)
      sig <- (sum(yt ^ 2) - b ^ 2 * Str2) / (n - k - 2)
      err <- sqrt(sig * (1 / Str2))
      zvalue <- b / err
      p <- 2 * pnorm(-abs(zvalue))
      return(list(beta.cv =U2, se.cv = err.cv, zvalue.cv = zvalue.cv, p.cv = p.cv, beta = b, se = err, zvalue = zvalue, p = p))
})

   fastLm <- compiler::cmpfun(function(x){
      Mod<- RcppEigen::fastLmPure(y = y, X = as.matrix(cbind(x, CVt), nrow=length(y)), method = 2)
      B2 <- Mod$coefficients[1]
      SE2 <- Mod$se[1]
      Z2 <- B2/SE2
      P2 <- 2*pnorm(-abs(Z2))
      if (length(ncol(X1)) == 0) {B1 = NULL} else {B1 = Mod$coefficients[2:(1+ncol(X1))]}
      if (length(ncol(X1)) == 0) {SE1 = NULL} else {SE1 = Mod$se[2:(1+ncol(X1))]}
      if (length(ncol(X1)) == 0) {Z1 = NULL} else {Z1 = B1/SE1}
      if (length(ncol(X1)) == 0) {P1 = NULL} else {P1 = 2*pnorm(-abs(Z1))}
      return(list(B2 = B2, SE2 = SE2, Z2 = Z2, P2 = P2, B1 = B1, SE1 = SE1, Z1 = Z1, P1 = P1))
})

Blink.LDRemove<-function(GDneo=NULL,LD=0.7,Porder=NULL,block=1000,LD.num =50){
#Objects: LD remove, especially length(Porder)>10000
#Authors: Yao Zhou
#Last update: 08/15/2016
#####################################
Blink.LDRemoveBlock<-function(GDneo=NULL,LD=NULL,Porder=NULL){

#Objects: Calculate LD and remove the correlated SNPs
#Authors: Yao Zhou
#Last Update:  03/03/16
#################################
    GDneo=as.matrix(GDneo)
    n=nrow(GDneo)
	corr=cor(GDneo)
	corr[is.na(corr)]=1
	corr[abs(corr)<=LD]=0
	corr[abs(corr)>LD]=1
	Psort=as.numeric(matrix(1,1,ncol(corr)))
	# print(ncol(corr))
	for(i in 2:ncol(corr)){
		p.a=Psort[1:(i-1)]
		p.b=as.numeric(corr[1:(i-1),i])
		index=(p.a==p.b)
		index[(p.a==0)&(p.b==0)]=FALSE
		if(sum(index)!=0) Psort[i]=0
	}
	seqQTN=Porder[Psort==1]
	#seqQTN=Porder[Psort==1,c(1,2),drop=F]
	return(seqQTN)
}

##############
    GDneo = as.matrix(GDneo)
    SNP.index = apply(GDneo, 2, sd) != 0
    GDneo = GDneo[, SNP.index]
    Porder = Porder[SNP.index]
    l = block
	seqQTN=NULL
	lp=length(Porder)
	k=ceiling(lp/l)
	GDneo=as.matrix(GDneo)
    n=nrow(GDneo)

	for(i in 1:k){
		bottom=(i-1)*l+1
		up=l*i
		if(up>lp) up = lp
		Porderb=Porder[bottom:up]

		index = seq(bottom:up)
		GDneob = GDneo[,index]
		# cat("i is ",i,"\n")
		# print(length(index))
		seqQTNs = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb)
		# print(seqQTN)
		seqQTN = append(seqQTN,seqQTNs)
		if(k >1){
		  index1 = which(Porder %in% seqQTN)
		  Porderb = Porder[index1]
		  GDneob = GDneo[,index1]
		  if(length(index1)>1){
		    seqQTN = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb)
		  }else{
		    seqQTN = Porderb
		  }

		}
		if(LD.num < length(seqQTN)) break
	}
	rm(GDneob,Porderb)
	return(seqQTN)
}

    if (is.null(maxStep)){maxStep=10}else{maxStep = maxStep}
    if (is.null(selector) ||selector == "L0"){requireNamespace(l0ara,quietly=TRUE)}
    if (selector == "L1"){requireNamespace(SIS,quietly=TRUE)}
    if (selector == "lars") {requireNamespace(lars,quietly=TRUE)}
    if (selector == "FBED") {requireNamespace(Rfast,quietly=TRUE)}
    if (selector == "vif") {requireNamespace(VIF,quietly=TRUE)}
    if (selector == "MARS") {requireNamespace(earth,quietly=TRUE)}
    if (selector == "MCP") {requireNamespace(ncvreg,quietly=TRUE)}
    if (selector == "SCAD") {requireNamespace(SIS,quietly=TRUE)}
    #if (selector == "SES") {requireNamespace(MXM,quietly=TRUE)}
    if (selector == "bess") {requireNamespace(BeSS,quietly=TRUE)}

   b.mat <- matrix(NA, maxStep, ncol(X))
   colnames(b.mat) <- colnames(X)
   rownames(b.mat)<- c(paste('step', 1:maxStep, sep=''))
   se.mat <- matrix(NA, maxStep, ncol(X))
   colnames(se.mat)<-colnames(X)
   rownames(se.mat)<- c(paste('step',1:maxStep,sep=''))
   z.mat <- matrix(NA, maxStep, ncol(X))
   colnames(z.mat)<-colnames(X)
   rownames(z.mat)<- c(paste('step',1:maxStep,sep=''))
   p.mat <- matrix(NA, maxStep, ncol(X))
   colnames(p.mat) <- colnames(X)
   rownames(p.mat) <- c(paste('step', 1:maxStep, sep=''))
   #QTNname<-vector()
   #QTNname<-list()
   yn = length(y)
   X2 <- X     ## initial setting for scanning SNPs
   X1 = NULL   ## initial setting for SNP covariates
   if (is.null(CV)) {CV=matrix(rep(1, yn), ncol=1); colnames(CV) <- c("mu")}else{
   if (!is.matrix(CV)) {CV = matrix(CV)}
   if (length(colnames(CV))==0) {colnames(CV) <- c(paste('CV', 1:ncol(CV), sep=''))}
   namet <- colnames(CV)
   CV = cbind(CV, rep(1, yn))
   colnames(CV) <- c(namet, "mu")
   rm(namet)
}


  test.name = NULL
  countor  = 0

  print("Memory used at begining of calculation")
  print(memory.size())

for(a in 1:maxStep) {
################## Beginning of maxStep loop ##################
   print(c("Iterative step at ", a), quote = FALSE)
   if (is.null(X1)){CVt = CV}else{
       CVt=as.matrix(cbind(X1, CV), nrow=yn)
       colnames(CVt)<-c(colnames(X1), colnames(CV))}

if (fun.selection == "fastLm") {
   if (is.null(chunk.num) || chunk.num == 1) {res <- apply(X2, 2, fastLm)} else {
     chunk.seq<-cut(seq_len(dim(X2)[2]), breaks= chunk.num, labels=FALSE)
     res <- vector(mode = "list", chunk.num)
     bili <- vector(mode="numeric", chunk.num)
     for (s in 1:chunk.num){
     index <- which(chunk.seq==s, arr.ind=TRUE)
     Xt<- X2[, index]
     m <- apply(Xt, 2, fastLm)
     res[[s]]<-m
     bili[s] <- length(index)
}
     rm(Xt); rm(m); rm(chunk.seq); rm(index)
}

#=================== Extract b-value  ========================#
 if (is.null(chunk.num) || chunk.num == 1) {
    if (is.null(X1)) {
    Bx1 = NULL
    Bx2 = unlist(lapply(res, function(x){c(x[[1]])}))
    names(Bx2)<-colnames(X2)
    B <-c(Bx1, Bx2)
    B <- B[c(colnames(X))]}else{
    tt <- lapply(res, function(x){c(x[[5]])})
    Bx1 <- matrix(unlist(tt), nrow = ncol(X1), ncol = ncol(X2))
    if (extraction == "min" || is.null(extraction)){Bx1<- Rfast::rowMins(Bx1, value = TRUE)}
    if (extraction == "mean" || extraction == "bagging-like"){Bx1<-Rfast::rowmeans(Bx1)}
    if (extraction == "1se") {Bx1 <- apply(Bx1, 1, function(x){x[which.min(x[x-min(x)>= sd(x)])]})}
    Bx1 <- as.numeric(Bx1)
    names(Bx1)  <- colnames(X1)
    Bx2 <- unlist(lapply(res, function(x){c(x[[1]])}))
    names(Bx2)<-colnames(X2)
    B <- c(Bx1, Bx2)
    B <- B[c(colnames(X))]
}
b.mat[a, ] <- B} else {
    if (is.null(X1)) {
    Bx1 = NULL
    Bx2 = NULL
    for (s in 1:chunk.num){Bx2 = c(Bx2, unlist(lapply(res[[s]], function(x){c(x[[1]])})))}
    names(Bx2)<-colnames(X2)
    B <-c(Bx1, Bx2)
    B <- B[c(colnames(X))]}else{
    Bx1 = matrix(NA, nrow = ncol(X1), ncol = ncol(X2))
    Bx2 = NULL
    aili<- cumsum(bili)
for (s in 1:chunk.num){
    tx1 <- lapply(res[[s]], function(x){c(x[[5]])})
    tx1 <- matrix(unlist(tx1), nrow = ncol(X1))
    s1 <- ifelse(s==1, 1, cumsum(bili[1:s])[s-1]+1)
    s2 <- aili[s]
    Bx1[, s1:s2] <- tx1}
if (extraction == "min" || is.null(extraction)){Bx1<- Rfast::rowMins(Bx1, value = TRUE)}
    if (extraction == "mean" || extraction == "bagging-like"){Bx1<-Rfast::rowmeans(Bx1)}
    if (extraction == "1se") {Bx1 <- apply(Bx1, 1, function(x){x[which.min(x[x-min(x)>= sd(x)])]})}
    Bx1 <- as.numeric(Bx1)
    names(Bx1)  <- colnames(X1)
    for (s in 1:chunk.num){Bx2 = c(Bx2, unlist(lapply(res[[s]], function(x){c(x[[1]])})))}
    names(Bx2)<-colnames(X2)
    B <- c(Bx1, Bx2)
    B <- B[c(colnames(X))]
}
    b.mat[a, ] <- B
}

#==================   Extract the se-value  ===================#
if (is.null(chunk.num) || chunk.num == 1) {
    if (is.null(X1)) {
    SEx1 = NULL
    SEx2 = unlist(lapply(res, function(x){c(x[[2]])}))
    names(SEx2)<-colnames(X2)
    SE <-c(SEx1, SEx2)
    SE <- SE[c(colnames(X))]}else{
    tt <- lapply(res, function(x){c(x[[6]])})
    SEx1 <- matrix(unlist(tt), nrow =ncol(X1), ncol = ncol(X2))
    if (extraction == "min" || is.null(extraction)){SEx1<- Rfast::rowMins(SEx1, value = TRUE)}
    if (extraction == "mean" || extraction == "bagging-like") {SEx1 <- Rfast::rowmeans(SEx1)}
    if (extraction == "1se") {SEx1 <- apply(SEx1, 1, function(x){x[which.min(x[x-min(x)>= sd(x)])]})}
    SEx1 <- as.numeric(SEx1)
    names(SEx1) <-colnames(X1)
    SEx2 <- unlist(lapply(res, function(x){c(x[[2]])}))
    names(SEx2)<-colnames(X2)
    SE <- c(SEx1, SEx2)
    SE <- SE[c(colnames(X))]
}
se.mat[a, ] <- SE} else {
    if (is.null(X1)) {
    SEx1 = NULL
    SEx2 = NULL
    for (s in 1:chunk.num){SEx2 = c(SEx2, unlist(lapply(res[[s]], function(x){c(x[[2]])})))}
    names(SEx2)<-colnames(X2)
    SE <-c(SEx1, SEx2)
    SE <- SE[c(colnames(X))]}else{
    SEx1 = matrix(NA, nrow = ncol(X1), ncol = ncol(X2))
    SEx2 = NULL
    aili<- cumsum(bili)
for (s in 1:chunk.num){
    tx1 <- lapply(res[[s]], function(x){c(x[[6]])})
    tx1 <- matrix(unlist(tx1), nrow = ncol(X1))
    s1 <- ifelse(s==1, 1, cumsum(bili[1:s])[s-1]+1)
    s2 <- aili[s]
    SEx1[, s1:s2] <- tx1}
    if (extraction == "min" || is.null(extraction)){SEx1<- Rfast::rowMins(SEx1, value = TRUE)}
    if (extraction == "mean" || extraction == "bagging-like") {SEx1 <- Rfast::rowmeans(SEx1)}
    if (extraction == "1se") {SEx1 <- apply(SEx1, 1, function(x){x[which.min(x[x-min(x)>= sd(x)])]})}
    SEx1 <- as.numeric(SEx1)
    names(SEx1) <-colnames(X1)
    for (s in 1:chunk.num){SEx2 = c(SEx2, unlist(lapply(res[[s]], function(x){c(x[[2]])})))}
    names(SEx2)<-colnames(X2)
    SE <- c(SEx1, SEx2)
    SE <- SE[c(colnames(X))]
}
    se.mat[a, ] <- SE
}

#===============  extract z-value  ====================#
if (is.null(chunk.num) || chunk.num == 1) {
    if (is.null(X1)) {
    Zx1 = NULL
    Zx2 = unlist(lapply(res, function(x){c(x[[3]])}))
    names(Zx2)<-colnames(X2)
    Z <-c(Zx1, Zx2)
    Z <- Z[c(colnames(X))]}else{
    tt <- lapply(res, function(x){c(x[[7]])})
    Zx1 <- matrix(unlist(tt), nrow =ncol(X1), ncol = ncol(X2))
    if (extraction == "min" || is.null(extraction)){Zx1<- Rfast::rowMins(Zx1, value = TRUE)}
    if (extraction == "mean" || extraction == "bagging-like") {Zx1<-Rfast::rowmeans(Zx1)}
    if (extraction == "1se") {Zx1 <- apply(Zx1, 1, function(x){x[which.min(x[x-min(x)>= sd(x)])]})}
    Zx1 <- as.numeric(Zx1)
    names(Zx1) <-colnames(X1)
    Zx2 <- unlist(lapply(res, function(x){c(x[[3]])}))
    names(Zx2)<-colnames(X2)
    Z <- c(Zx1, Zx2)
    Z <- Z[c(colnames(X))]
}
    z.mat[a, ] <- Z} else {
    if (is.null(X1)) {
       Zx1 = NULL
       Zx2 = NULL
       for (s in 1:chunk.num){Zx2 = c(Zx2, unlist(lapply(res[[s]], function(x){c(x[[3]])})))}
       names(Zx2)<-colnames(X2)
       Z <-c(Zx1, Zx2)
       Z <- Z[c(colnames(X))]}else{
   Zx1 = matrix(NA, nrow = ncol(X1), ncol = ncol(X2))
   Zx2 = NULL
   aili<- cumsum(bili)
for (s in 1:chunk.num){
   tx1 <- lapply(res[[s]], function(x){c(x[[7]])})
   tx1 <- matrix(unlist(tx1), nrow = ncol(X1))
   s1 <- ifelse(s==1, 1, cumsum(bili[1:s])[s-1]+1)
   s2 <- aili[s]
   Zx1[, s1:s2] <- tx1}
   if (extraction == "min" || is.null(extraction)){Zx1<- Rfast::rowMins(Zx1, value = TRUE)}
   if (extraction == "mean" || extraction == "bagging-like") {Zx1<-Rfast::rowmeans(Zx1)}
   if (extraction == "1se") {Zx1 <- apply(Zx1, 1, function(x){x[which.min(x[x-min(x)>= sd(x)])]})}
   Zx1 <- as.numeric(Zx1)
   names(Zx1) <-colnames(X1)
   for (s in 1:chunk.num){Zx2 = c(Zx2, unlist(lapply(res[[s]], function(x){c(x[[3]])})))}
   names(Zx2)<-colnames(X2)
   Z <- c(Zx1, Zx2)
   Z <- Z[c(colnames(X))]
}
   z.mat[a, ] <- Z
}

#================== extract p-value  =======================#
if (is.null(chunk.num) || chunk.num == 1) {
    if (is.null(X1)) {
    Px1 = NULL
    Px2 = unlist(lapply(res, function(x){c(x[[4]])}))
    names(Px2)<-colnames(X2)
    P <-c(Px1, Px2)
    P <- P[c(colnames(X))]}else{
    tt <- lapply(res, function(x){c(x[[8]])})
    Px1 <- matrix(unlist(tt), nrow =ncol(X1), ncol = ncol(X2))
    if (extraction == "min" || is.null(extraction)){Px1<- Rfast::rowMins(Px1, value = TRUE)}
    if (extraction == "mean") {Px1 <- Rfast::rowmeans(Px1)}
    if (extraction == "1se") {Px1 <- apply(Px1, 1, function(x){x[which.min(x[x-min(x)>= sd(x)])]})}
    if (extraction == "bagging-like"){Px1 <- 2*pnorm(-abs(Bx1/SEx1))}
    Px1 <- as.numeric(Px1)
    names(Px1) <-colnames(X1)
    Px2 <- unlist(lapply(res, function(x){c(x[[4]])}))
    names(Px2)<-colnames(X2)
    P <- c(Px1, Px2)
    P <- P[c(colnames(X))]
}
    p.mat[a, ] <- P} else {
    if (is.null(X1)) {
        Px1 = NULL
        Px2 = NULL
        for (s in 1:chunk.num){Px2 = c(Px2, unlist(lapply(res[[s]], function(x){c(x[[4]])})))}
       names(Px2)<-colnames(X2)
       P <-c(Px1, Px2)
       P <- P[c(colnames(X))]}else{
       Px1 = matrix(NA, nrow = ncol(X1), ncol = ncol(X2))
       Px2 = NULL
      aili<- cumsum(bili)
for (s in 1:chunk.num){
       tx1 <- lapply(res[[s]], function(x){c(x[[8]])})
       tx1 <- matrix(unlist(tx1), nrow = ncol(X1))
       s1 <- ifelse(s==1, 1, cumsum(bili[1:s])[s-1]+1)
       s2 <- aili[s]
       Px1[, s1:s2] <- tx1}
    if (extraction == "min" || is.null(extraction)){Px1<- Rfast::rowMins(Px1, value = TRUE)}
    if (extraction == "mean") {Px1 <- Rfast::rowmeans(Px1)}
    if (extraction == "1se") {Px1 <- apply(Px1, 1, function(x){x[which.min(x[x-min(x)>= sd(x)])]})}
    if (extraction == "bagging-like"){Px1 <- 2*pnorm(-abs(Bx1/SEx1))}
    Px1 <- as.numeric(Px1)
    names(Px1) <-colnames(X1)
    for (s in 1:chunk.num){Px2 = c(Px2, unlist(lapply(res[[s]], function(x){c(x[[4]])})))}
    names(Px2)<-colnames(X2)
    P <- c(Px1, Px2)
    P <- P[c(colnames(X))]
}
    p.mat[a, ] <- P

}

} else {## beginning of lm.vec
   if (is.null(chunk.num) || chunk.num == 1) {haha<-lm.vec(y, X2, CVt)}else{
     haha<- vector(mode="list", chunk.num)
     chunk.seq<-cut(seq(1, ncol(X2)), breaks= chunk.num, labels=FALSE)
     for (s in 1:chunk.num){
     index <- which(chunk.seq==s, arr.ind=TRUE)
     Xt<- X2[, index]
     m <- lm.vec(y, X2, CVt)
     haha[[s]]<-m
}
rm(Xt); rm(m); rm(chunk.seq); rm(index)
}

#=============== Extract b-value  =====================#
 if (length(haha) == 8) {
if (is.null(X1)) {
    Bx1 = NULL
    Bx2 = haha[[5]]
    names(Bx2)<-colnames(X2)
    B <-c(Bx1, Bx2)
    B <- B[c(colnames(X))]
}else{
    Bx1 <- haha[[1]][1:ncol(X1)]
    names(Bx1) <-colnames(X1)
    Bx2 <- haha[[5]]
    names(Bx2)<-colnames(X2)
    B <- c(Bx1, Bx2)
    B <- B[c(colnames(X))]
}
}else {
    Bx1 <- vector()
    Bx2 <- vector()
if (is.null(X1)) {
  Bx1 = NULL
  for (s in 1:chunk.num){Bx2 = c(Bx2, haha[[s]][[5]])}
  names(Bx2)<-colnames(X2)
  B = c(Bx1, Bx2)
  B <- B[c(colnames(X))]
}

if (!is.null(X1)) {
  for (s in 1:chunk.num){Bx2 = c(Bx2, haha[[s]][[5]])}
  Bx1 = c(haha[[1]][[1]][1:ncol(X1)])# same in each s
  names(Bx1)<-colnames(X1)
  names(Bx2)<-colnames(X2)
  B <- c(Bx1, Bx2)
  B <- B[c(colnames(X))]
}
}
b.mat[a, ] <- B

#============   Extract the se-value  =======================#
if (length(haha) == 8) {
if (is.null(X1)) {
    SEx1 = NULL
    SEx2 = haha[[6]]
    names(SEx2)<-colnames(X2)
    SE <-c(SEx1, SEx2)
    SE <- SE[c(colnames(X))]
}else{
    SEx1 <- haha[[2]][1:ncol(X1)]
    names(SEx1) <-colnames(X1)
    SEx2 <- haha[[6]]
    names(SEx2)<-colnames(X2)
    SE <- c(SEx1, SEx2)
    SE <- SE[c(colnames(X))]
}
}else {
  SEx1 <- vector()
  SEx2 <- vector()
if (is.null(X1)) {
  SEx1 = NULL
  for (s in 1:chunk.num){SEx2 = c(SEx2, haha[[s]][[6]])}
  names(SEx2)<-colnames(X2)
  SE = c(SEx1, SEx2)
  SE <- SE[c(colnames(X))]
}

if (!is.null(X1)) {
  for (s in 1:chunk.num){SEx2 = c(SEx2, haha[[s]][[6]])}
  SEx1 = c(haha[[1]][[2]][1:ncol(X1)]) #same in each s
  names(SEx1)<-colnames(X1)
  names(SEx2)<-colnames(X2)
  SE <- c(SEx1, SEx2)
  SE <- SE[c(colnames(X))]
}
}
se.mat[a, ] <- SE

#==================  Extract z-value  ======================#
if (length(haha) == 8) {
if (is.null(X1)) {
    Zx1 = NULL
    Zx2 = haha[[7]]
    names(Zx2)<-colnames(X2)
    Z <-c(Zx1, Zx2)
    Z <- Z[c(colnames(X))]
}else{
    Zx1 <- haha[[3]][1:ncol(X1)]
    names(Zx1) <-colnames(X1)
    Zx2 <- haha[[7]]
    names(Zx2)<-colnames(X2)
    Z <- c(Zx1, Zx2)
    Z <- Z[c(colnames(X))]
}
}else{
  Zx1 <- vector()
  Zx2 <- vector()
if (is.null(X1)) {
  Zx1 = NULL
  for (s in 1:chunk.num){Zx2 = c(Zx2, haha[[s]][[7]])}
  names(Zx2)<-colnames(X2)
  Z = c(Zx1, Zx2)
  Z <- Z[c(colnames(X))]
}

if (!is.null(X1)) {
  for (s in 1:chunk.num){Zx2 = c(Zx2, haha[[s]][[7]])}
  Zx1 = c(haha[[1]][[3]][1:ncol(X1)]) #same in each s
  names(Zx1)<-colnames(X1)
  names(Zx2)<-colnames(X2)
  Z<-c(Zx1, Zx2)
  Z <- Z[c(colnames(X))]
}
}
z.mat[a, ] <- Z

#================= extract p-value  ======================#
if (length(haha) == 8) {
if (is.null(X1)) {
    Px1 = NULL
    Px2 = haha[[8]]
    names(Px2)<-colnames(X2)
    P <-c(Px1, Px2)
    P <- P[c(colnames(X))]
}else{
    Px1 <- haha[[4]][1:ncol(X1)]
    names(Px1) <-colnames(X1)
    Px2 <- haha[[8]]
    names(Px2)<-colnames(X2)
    P <- c(Px1, Px2)
    P <- P[c(colnames(X))]
}
}else{
    Px1 <- vector()
    Px2 <- vector()
if (is.null(X1)) {
      Px1 = NULL
     for (s in 1:chunk.num){Px2 = c(Px2, haha[[s]][[8]])}
     names(Px2)<-colnames(X2)
     P = c(Px1, Px2)
     P <- P[c(colnames(X))]
}

if (!is.null(X1)) {
     for (s in 1:chunk.num){Px2 = c(Px2, haha[[s]][[8]])}
     Px1 = c(haha[[1]][[4]][1:ncol(X1)]) #same in each s
     names(Px1)<-colnames(X1)
     names(Px2)<-colnames(X2)
     P<-c(Px1, Px2)
     P <- P[c(colnames(X))]
}
}
p.mat[a, ] <- P
}

# End of lm.vec function
###############################################################
##=========== reset of pseudoQTNs by selector =================#
print(length(which(p.mat[a, ] < cutOff/ncol(X))))
 if (length(which(p.mat[a, ] < cutOff/ncol(X))) <= 1){
      X1 = matrix(X[, which.min(p.mat[a, ])], ncol = 1)
      colnames(X1)<-colnames(p.mat[, which.min(p.mat[a, ])])} else {
 if(length(which(p.mat[a, ] < cutOff/ncol(X))) <=2000){
      X1 = matrix(X[, which(p.mat[a, ] < cutOff/ncol(X))], nrow=nrow(X))
      colnames(X1)<-colnames(X[, which(p.mat[a, ] < cutOff/ncol(X))])
    }else{
      X1_PZ = matrix(X[, which(p.mat[a, ] < cutOff/ncol(X))], nrow=nrow(X))
      colnames(X1_PZ)<-colnames(X[, which(p.mat[a, ] < cutOff/ncol(X))])
##LD Removed

   X1_P<-p.mat[a,][colnames(X1_PZ)];X1_P<-as.matrix(X1_P);X1_P<-X1_P[order(X1_P[,1]),];X1_P<-as.matrix(X1_P)
   X1<-X1_PZ[,rownames(X1_P)]
   Psort=Blink.LDRemove(Porder=X1_P,GDneo=X1,LD=0.7)
   indexQTN<-X1_P[,1]%in%Psort; X1_newP<-as.matrix(X1_P[indexQTN,])
   X1<-X[,rownames(X1_newP)]
   rm(X1_P);rm(X1_newP);rm(indexQTN);rm(X1_PZ);rm(Psort)
print(dim(X1))
   #X1<-as.matrix(X1, nrow=nrow(X))
   name<-colnames(X1)
}
}
print(dim(X1))
   #X1<-as.matrix(X1, nrow=nrow(X))
   name<-colnames(X1)

   if (is.null(selector) ||selector == "L0"){
   if (ncol(X1)==1 || length(name)==1){name=name}else{
   fit <- l0ara::l0ara(X1, y, lam = log(yn))
   index<-which(fit$beta!=0)
   if (length(index)==0){name<-colnames(X[, which.min(p.mat[a, ])])}else{
   name<-colnames(X1[, index])}
   rm(index); rm(fit)}
   }

   # if (selector == "FBED"){
   # if (ncol(X1)==1 || length(name)==1){name=name}else{
   # fit <- Rfast::cor.fbed(y, X1, alpha = 0.01, K = 0)
   # index <- fit$res[,1]
   # remove(fit)
   # if (length(index)==0){name<- colnames(X[, which.min(p.mat[a, ])])}else{
   # name<-colnames(X1[, index])}
   # rm(index)}
   # }

   if (selector == "L1"){
   if (ncol(X1)==1 || length(name)==1){name=name}else{
   fit <- SIS::SIS(X1, y, family = "gaussian", penalty="lasso", tune = "bic", nfolds = 5)
   index<-fit$ix
   remove(fit)
   if (length(index)==0){name<-colnames(p.mat[a,][which.min(p.mat[a,])])}else{
   name<-colnames(X1[, index])}
   rm(index)}
   }

   if (selector == "lars"){
   if (ncol(X1)==1 || length(name)==1){name=name}else{
   fit <-lars::lars(X1, y, type = "lar", trace = FALSE)
   best_step <- fit$df[which.min(fit$Cp)]
   index<-coef(fit, s=best_step)
   remove(fit)
   index<-which(index != 0)
   if (length(index)==0){name<- colnames(X[, which.min(p.mat[a, ])])}else{
   name<- colnames(X1[, index])}
   rm(index)}
   }


   if (selector == "vif"){
   if (ncol(X1)==1 || length(name)==1){name=name}else{
   fit0 <- VIF::vif(y, X1, trace = FALSE)
   index<-fit0$select
   remove(fit0)
   if (length(index)==0){name<- colnames(X[, which.min(p.mat[a, ])])}else{
   name<-colnames(X1[, index])}
   rm(index)}
   }

if (selector == "MARS"){
   if (ncol(X1)==1 || length(name)==1){name=name}else{
   mydata<-as.data.frame(cbind(y, X1))
   colnames(mydata)<-c("y",colnames(X1))
   res <- earth::earth(y ~ . , data= mydata)
   index <- earth::evimp(res)
   remove(res);remove(mydata)
   if (length(index)==0){name<- colnames(X[, which.min(p.mat[a, ])])}else{
   name<-rownames(index)}
   rm(index)}
   }

   if (selector == "SCAD"){
   if (ncol(X1)==1 || length(name)==1){name=name}else{
   res <- SIS::SIS(X1, y, family = "gaussian", penalty="SCAD", tune = "bic", nfolds = 5)
   index<-res$ix
   remove(res)
   if (length(index)==0){name<- colnames(X[, which.min(p.mat[a, ])])}else{
   name<-colnames(X1)[index] }
   rm(index)}
   }

   if (selector == "MCP"){
   if (ncol(X1)==1 || length(name)==1){name=name}else{
   fit <- ncvreg::cv.ncvreg(X1, y, penalty = "MCP", gamma = 3)
   res<-coef(fit)[-1]
   remove(fit)
   index<-which(res != 0)
   remove(res)
   if (length(index)==0){name<- colnames(X[, which.min(p.mat[a, ])])}else{
   name<-names(index)}
   rm(index)}
   }

   if (selector == "stepwise"){
   if (ncol(X1)==1 || length(name)==1){name=name}else{
   #mydata<-as.data.frame(cbind(y, X1))
   #colnames(mydata)<-c("y",colnames(X1))
   #base.mod <- lm(y ~ 1 , data = mydata)
   #all.mod <- lm(y ~ . , data= mydata)
   #stepMod <- stats::step(base.mod, scope = list(lower = base.mod, upper = all.mod), direction = "both", trace = 0, steps = 1000)
   #shortlistedVars <- names(unlist(stepMod[[1]]))
   #index <- shortlistedVars[!shortlistedVars %in% "(Intercept)"]
   mydata<-cbind(y, X1)
   colnames(mydata)<-c("y",colnames(X1))
   res<-StepReg::stepwise(cbind(y,X1), "y", selection = "bidirection",sle = 0.01, sls = 0.01)
   index <- res$variate[-1]
   if (length(index)==0){name<-colnames(X[, which.min(p.mat[a, ])])}else{
   name<- index
   print(length(name))
   print(name)}
   rm(index);rm(mydata);rm(res)}
   }

# if (selector == "SES"){
   # if (ncol(X1)==1 || length(name)==1){name=name}else{
   # res <- MXM::SES(y, X1)
   # index <- res@selectedVars
   # remove(res)
   # if (length(index)==0){name<- colnames(X[, which.min(p.mat[a, ])])}else{
   # name<-colnames(X1[,index])}
   # rm(index)}
   # }

if (selector == "bess"){
   if (ncol(X1)==1 || length(name)==1){name=name}else{
   res <- BeSS::bess(X1,y)
   index <- colnames(res$bestmodel$model$xbest)
   remove(res)
   if (length(index)==0){name<- colnames(X[, which.min(p.mat[a, ])])}else{
   name<-index}
   rm(index)}
   }

#remove redundant SNPs in the same cluster cut by 0.7##
  if (length(name)>= 2){
  X1<-as.matrix(X[, c(name)],ncol=length(name))
   colnames(X1)<-name
   sigma = cov(X1)
   sigma.distance = as.dist(1 - abs(cov2cor(sigma)))
   fit = hclust(sigma.distance, method="single")
   clusters = cutree(fit, h=1-0.7)
   pvals.vector<- p.mat[a, colnames(X1)]
   names(pvals.vector)<- colnames(X1)
   top.selected = sapply(1:max(clusters), function(c) {
    cluster_elements = clusters==c
    top_within = which.min(pvals.vector[cluster_elements])
    if(length(top_within)==0) top_within = 1
    which(cluster_elements)[top_within]
    })
   name<-names(top.selected)

   QTNSelect<-name

   #print(length(name))
   rm(sigma);rm(sigma.distance);rm(fit);rm(clusters);rm(pvals.vector);rm(top.selected)
    }
print(length(name))
print(name)
###### end of redundant SNPs removing procedure #######

   if (identical(test.name, name) != TRUE) {test.name<-name} else {countor = countor +1; if (countor == 2) {print(c("Calculation ended at ", a), quote=FALSE); break}}
   X1<-as.matrix(X[, c(name)], nrow=nrow(X))
   colnames(X1)<-name

   index <- colnames(X)%in%name
   X2 <- matrix(X[, !index], nrow = nrow(X))
   colnames(X2)<- colnames(X)[!index]
   rm(index)
   gc(reset=TRUE)
}

############ End of the maxStep loop ##########################
   #set.seed(seed = NULL)             #because SIS's set.seed(0)
   rm(.Random.seed, envir=globalenv())#because SIS's set.seed(0)

   suggested.SNPs <- names(p.mat[a, which(p.mat[a, ]  < cutOff/ncol(X))])
   removeld.SNPs<-QTNSelect
   if (!is.null(GM)) {
       Pt <-p.mat[a, ]
       total.map <- cbind(GM, Pt)
       rownames(total.map) <- c(rownames(GM))
       colnames(total.map) <- c(colnames(GM),"p-value")} else {total.map <- NULL}
   total.value <- cbind(b.mat[a, ], se.mat[a, ], z.mat[a, ], p.mat[a, ])
   rownames(total.value) <- colnames(p.mat)
   colnames(total.value) <- c("coefficient", "se", "z-value", "p-value")
   if(file.output){
   # if(plot.style=="SEIR"){
  # SEIR.manhattan.plot(total.map,color="rainbow")
  # SEIR.QQ(total.map[,4])}
   if(plot.style=="CMplot"){
      CMplot::CMplot(total.map,plot.type=c("m","q"),LOG10=TRUE, ylim=NULL, threshold=cutOff,threshold.lty=c(2),cex.axis=1,
        threshold.lwd=1, threshold.col=c("red"), amplify=TRUE,bin.size=1e6,cex.lab=1.5,ylab="",
        chr.den.col=c("darkgreen", "yellow", "red"),signal.col=NULL,signal.cex=c(1),ylab.pos=3,
        signal.pch=c(19,19),file="jpg",memo="",dpi=600,file.output=TRUE,verbose=TRUE)}
   write.csv(total.map,paste("SEIR.", "name.of.trait",".GWAS.csv" ,sep = ""),quote=F,row.names=F)
  }
  return(list(suggested.SNPs = suggested.SNPs,removeld.SNPs=removeld.SNPs, total.map = total.map, total.value = total.value, p.mat = p.mat[1:a, ], se.mat = se.mat[1:a, ], z.mat = z.mat[1:a, ], b.mat = b.mat[1:a, ]))


}

################### End of SEIR function ######################
#################################################################





