#
#  growth : A Library of Normal Distribution Growth Curve Models
#  Copyright (C) 1998 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     carma(response, ccov=NULL, times=NULL, torder=0, interaction,
#	arma=c(0,0,0), parma=NULL, pre=NULL, position=NULL, iopt=T, resid=T,
#	transform="identity", delta=NULL, print.level=0, typsiz=abs(p),
#	ndigit=10, gradtol=0.00001, steptol=0.00001, iterlim=100, fscale=1,
#	stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    Function to fit the multivariate normal distribution with ARMA
#  and random intercept using Kalman filtering in continuous time

.First.lib <- function(lib, pkg) {
	library.dynam("growth", pkg, lib)
	provide(growth)
}
require(rmutil)

carma <- function(response, ccov=NULL, times=NULL, torder=0, interaction,
	arma=c(0,0,0), parma=NULL, pre=NULL, position=NULL, iopt=T, resid=T,
	transform="identity", delta=NULL, print.level=0, typsiz=abs(p),
	ndigit=10, gradtol=0.00001, steptol=0.00001, iterlim=100, fscale=1,
	stepmax=10*sqrt(p%*%p)){
kalman <- function(p){
	z <- .Fortran("kalman",
		np=as.integer(length(p)),
		par=as.double(p),
		like=double(1),
		xx=matrix(0,nrow=nlp+1,ncol=nlp+1),
		y=as.double(y),
		sse=double(1),
		nq=as.integer(nre),
		nlp=as.integer(nlp),
		ns=as.integer(ns),
		nt=as.integer(nt),
		model=as.integer(arma),
		t=as.double(times),
		nobs=as.integer(response$response$nobs),
		nod=as.integer(nod),
		as.integer(position),
		cv=as.double(ccov),
		ncv=as.integer(ncv),
		nxcv=as.integer(interaction),
		nx=as.integer(torder),
		p=double(length(p)),
		x=double(nlp+1),
		state=double(maxre*(nlp+1)),
		innov=double(nlp+1),
		cstate=complex(maxar*(nlp+1)),
		exx=complex(nlp+1),
		DUP=F)
	list(like=z$like/2,
		sse=z$sse,
		xx=z$xx)}
kalmanl <- function(p){
	z <- .Fortran("kalman",
		np=as.integer(length(p)),
		par=as.double(p),
		like=double(1),
		xx=matrix(0,nrow=nlp+1,ncol=nlp+1),
		y=as.double(y),
		sse=double(1),
		nq=as.integer(nre),
		nlp=as.integer(nlp),
		ns=as.integer(ns),
		nt=as.integer(nt),
		model=as.integer(arma),
		t=as.double(times),
		nobs=as.integer(response$response$nobs),
		nod=as.integer(nod),
		as.integer(position),
		cv=as.double(ccov),
		ncv=as.integer(ncv),
		nxcv=as.integer(interaction),
		nx=as.integer(torder),
		p=double(length(p)),
		x=double(nlp+1),
		state=double(maxre*(nlp+1)),
		innov=double(nlp+1),
		cstate=complex(maxar*(nlp+1)),
		exx=complex(nlp+1),
		DUP=F)
	z$like/2}
maxre <- 6
maxar <- 6
maxma <- 5
call <- sys.call()
if(!missing(transform))transform <- match.arg(transform,c("identity",
	"exp","square","sqrt","log"))
p <- parma
torder <- torder+1
if(!is.null(pre))for(i in 1:length(pre))if(pre[i]==0)pre[i] <- 0.01
linear <- NULL
if(!inherits(response,"repeated")){
	if(!inherits(response,"response"))response <- restovec(response,times,delta=delta)
	if(!is.null(ccov)){
		if(!inherits(ccov,"tccov")){
			ccname <- paste(deparse(substitute(ccov)))
			if((is.matrix(ccov)&&is.null(colnames(ccov)))){
				ccname <- paste(deparse(substitute(ccov)))
				if(ncol(ccov)>1){
					tmp <- NULL
					for(i in 1:ncol(ccov))tmp <- c(tmp,paste(ccname,i,sep=""))}
					ccname <- tmp}
			ccov <- tcctomat(ccov,names=ccname)}
		ncv <- ncol(ccov$ccov)}
	else ncv <- 0
	response <- rmna(response=response, ccov=ccov)
	if(!is.null(ccov))rm(ccov)}
else {
	if(is.null(ccov))response$ccov <- NULL
	else if(inherits(ccov,"formula"))
		response$ccov$ccov <- attr(finterp(ccov,envir=response,expand=F,name=paste(deparse(substitute(response)))),"model")[,-1,drop=F]
	else stop("ccov must be a W&R formula")
	ncv <- if(is.null(response$ccov$ccov)) 0
		 else  ncol(response$ccov$ccov)}
y <- response$response$y
ns <- length(response$response$nobs)
if(!is.null(response$ccov))ccov <- t(response$ccov$ccov)
if(!is.null(response$response$times)){
	ave <- mean(response$response$times)
	times <- response$response$times-ave}
else {
	if(any(arma>0)||!is.null(parma))stop("No times. ARMA cannot be fitted")
	else if(torder>1)stop("No times. Time trends cannot be fitted.")
	ave <- 0
	times <- 1:length(y)}
if(missing(interaction)) {
	if(ncv==0)interaction <- NULL
	else interaction <- rep(1,ncv)}
else if(length(interaction)!=ncv)
	stop(paste("The vector, interaction, must have length, ", ncv))
else interaction <- interaction+1
if(!is.null(interaction)&&any(interaction>torder))
	stop("Interactions cannot have higher order than torder")
nre <- length(pre)
if(nre>maxre)stop(paste("Only ",maxre," random effects allowed"))
if(is.null(position)){
	if(nre==1)position <- c(1,1)
	else if(nre==2)position <- rbind(c(1,1),c(2,2))
	else if(nre>2)stop("If random effects are required, position must be supplied")}
if(length(arma)!=3){
	if(length(arma)==1)arma <- c(arma,0,0)
	else stop("The model specification must have length 3")}
if(sum(arma)==0&!missing(parma)){
	if(length(parma)==1)arma[1] <- 1
	else stop("If an ARMA is required, arma must be supplied")}
if(arma[1]>maxar)stop(paste("Maximum order of AR is ",maxar))
if(arma[2]>maxma)stop(paste("Maximum order of MA is ",maxma))
if(arma[2]>=arma[1]&&arma[2]!=0)
	stop(paste("Order of MA(",arma[2],") >= order of AR(",arma[1],")"))
if(arma[1]==0&&arma[3]!=0)
	stop("Cannot have observation errors without AR structure")
np <- arma[1]+arma[2]+arma[3]
if(np>0&&length(parma)!=np)
	stop(paste(np,"initial arma estimates must be supplied"))
if(!is.null(position)){
	if(!is.matrix(position))position <- t(as.matrix(position))
	if(nrow(position)!=nre)
		stop(paste("Random effects position must have",nre,"rows"))
	if(ncol(position)!=2)stop("Position matrix must have two columns")
	if((any(position[,1]<1)|any(position[,1]>nre)|any(position[,2]<1)|
		any(position[,2]>nre))||position[,1]>position[,2])
		stop("Position for covariance matrix out of range")
	if(max(position[,2])>torder)
		warning("number of random effects greater than order of polynomial in time")
	nod <- nrow(position)
	if(nod>nre*(nre+1)/2)
		stop(paste("Only ",nre*(nre+1)/2," elements allowed in the random effects covariance matrix"))
	else if(nod<nre)stop("Not enough initial estimates of random effects")
	p <- c(p,pre)
	np <- length(p)}
else {
	position <- NULL
	nod <- 0}
if(!is.null(response$response$times)&&torder>length(unique(response$response$times)))
	stop("torder is too large for the number of distinct times")
nt <- sum(response$response$nobs)
nlp <- torder
if(ncv>0) nlp <- nlp+sum(interaction)
if(transform=="identity")jacob <- 0
else if(transform=="exp"){
	jacob <- -sum(y)
	y <- exp(y)}
else {
	if(any(y<0))stop("Negative response values: invalid transformation")
	else if(transform=="square"){
		jacob <- -sum(log(y[y>0]))
		y  <- y^2}
	else if(transform=="sqrt"){
		jacob <- sum(log(y[y>0]))/2
		y <- sqrt(y)}
	else if(any(y==0))stop("Zero response values: invalid transformation")
	else if(transform=="log"){
		jacob <- sum(log(y))
		y <- log(y)}}
if(!is.null(response$response$delta)){
	if(length(response$response$delta)==1)jacob <- jacob-nm*log(response$response$delta)
	else jacob <- jacob -sum(log(response$response$delta))}
if(iopt&np>0){
	if(fscale==1)fscale <- kalmanl(p)
	z0 <- nlm(kalmanl, p=p, hessian=T, print.level=print.level,
		typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
		steptol=steptol, iterlim=iterlim, fscale=fscale)
	p <- z0$estimate}
z2 <- kalman(p)
like <- z2$like+jacob
sse <- z2$sse
aic <- like+np+nlp+1
z1 <- .Fortran("back",
	xx=z2$xx,
	as.integer(nlp),
	DUP=F)
beta <- z1$xx[1:nlp,nlp+1]
z1 <- .Fortran("ttvert",
	xx=z1$xx,
	as.integer(nlp),
	DUP=F)
mse <- sse/(nt-nlp-np)
vbeta <- z1$xx[1:nlp,1:nlp]*mse
if(nlp>1){
	sd <- sqrt(diag(vbeta))
	corr <- vbeta/(sd%o%sd)}
else {
     sd <- sqrt(vbeta)
     corr <- vbeta/sd^2}
if(np>0){
	if(np==1){
		nlcov <- 1/z0$hessian
		nlse <- sqrt(nlcov)}
	else {
		if(any(is.na(z0$hessian)))a <- 0
		else a <- qr(z0$hessian)$rank
		if(a==np)nlcov <- solve(z0$hessian)
		else nlcov <- matrix(NA,ncol=np,nrow=np)
		nlse <- sqrt(diag(nlcov))}
	if(length(nlse)>1)nlcorr <- nlcov/(nlse%o%nlse)
	else nlcorr <- nlcov/nlse^2
	coef <- z0$estimate
	grad <- z0$gradient
	code <- z0$code
	iterations <- z0$iterations}
else nlcov <- nlse <- nlcorr <- coef <- grad <- code <- iterations <- NULL
if(nod>0){
	a <- matrix(0,nrow=max(position[,2]),ncol=max(position[,2]))
	ii <- sum(arma)
	for(i in 1:nod){
		ii <- ii+1
		a[position[i,1],position[i,2]] <- p[ii]}}
else a <- NULL
if(resid){
	coln <- NULL
	if(ncv>0)for(i in 1:ns)for(j in 1:response$response$nobs[i])
		coln <- rbind(coln,ccov[1:ncv,i])
	pred <- rep(beta[1],nt)
	cc <- matrix(rep(1,nt),ncol=1)
	if(torder>1)for(i in 1:(torder-1)) {
		cc <- cbind(cc,times^i)
		pred <- pred+beta[i+1]*times^i}
	if(ncv>0){
		jj <- torder
		for(i in 1:ncv)for(j in 1:interaction[i]){
			jj <- jj+1
			cc <- cbind(cc,cc[,j]*coln[,i])
			pred <- pred+beta[jj]*cc[,j]*coln[,i]}}
	z <- .Fortran("resid",
		np=as.integer(length(p)),
		par=as.double(coef),
		beta=as.double(beta),
		ave=as.double(ave),
		pred=double(nt),
		sdr=double(nt),
		res=double(nt),
		y=as.double(y),
		sse=as.double(sse),
		nq=as.integer(nre),
		nlp=as.integer(nlp),
		ns=as.integer(ns),
		nt=as.integer(nt),
		model=as.integer(arma),
		t=as.double(times),
		nobs=as.integer(response$response$nobs),
		nod=as.integer(nod),
		as.integer(position),
		cv=as.double(ccov),
		nxcv=as.integer(interaction),
		nx=as.integer(torder),
		ncv=as.integer(ncv),
		p=double(length(p)),
		x=double(nlp+1),
		DUP=F)
	rpred <- z$pred
	sdr <- z$sdr
	rres <- z$res}
else pred <- res <- rpred <- sdr <- rres <- NULL
z <- list(
	call=call,
	response=response$response,
	ccov=response$ccov,
	linear=linear,
	torder=torder-1,
	interaction=interaction-1,
	np=np,
	nlp=nlp,
	nre=nre,
	nod=nod,
	arma=arma,
	maxlike=like,
	aic=aic,
	df=nt-nlp-np-1,
	coefficients=coef,
	nlse=nlse,
	nlcov=nlcov,
	nlcorr=nlcorr,
	are=a,
	position=position,
	beta=beta,
	sse=sse,
	mse=mse,
	vbeta=vbeta,
	corr=corr,
	se=sd,
	transform=transform,
	pred=pred,
	rpred=rpred,
	sdrpred=sdr,
	rresiduals=rres,
	gradient=grad,
	iterations=iterations,
	code=code)
class(z) <- "carma"
if(resid)class(z) <- c(class(z),"recursive")
return(z)}

coefficients.carma <- function(z) list(beta=z$beta,coef=z$coefficients)
deviance.carma <- function(z) 2*z$maxlike
fitted.carma <- function(z, recursive=TRUE) if(recursive) z$rpred else z$pred
residuals.carma <- function(z, recursive=TRUE){
	if(recursive)return(z$rresiduals)
	else {
		if(z$transform=="exp")z$response$y <- exp(z$response$y)
		else if(z$transform=="square")z$response$y  <- z$response$y^2
		else if(z$transform=="sqrt")z$response$y <- sqrt(z$response$y)
		else if(z$transform=="log")z$response$y <- log(z$response$y)
		return((z$response$y-z$pred)/sqrt(z$mse))}}

print.carma <- function(z, digits = max(3, .Options$digits - 3)) {
	if(!is.null(z$ccov$ccov))nccov <- ncol(z$ccov$ccov)
	else nccov <- 0
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(!is.null(z$code)&&z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("Number of subjects    ",length(z$response$nobs),"\n")
	cat("Number of observations",length(z$response$y),"\n")
	if(!is.null(z$response$times))cat("Mean time             ",mean(z$response$times),"\n")
	cat("Transformation        ",z$trans,"\n\n")
	cat("-Log likelihood   ",z$maxlike,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n")
	if(!is.null(z$iterations))cat("Iterations        ",z$iterations,"\n\n")
	cat("Estimates of linear parameters\n")
	if(inherits(z$ccov$linear,"formula"))
		cat("Formula: ",deparse(z$ccov$linear[[2]]),"\n")
	coef.table <- cbind(z$beta, z$se)
	cname <- "(Intercept)"
	if(z$torder>0)cname <- c(cname,paste("t^",1:z$torder,sep=""))
	if(nccov>0)for(i in 1:nccov){
		cname <- c(cname,colnames(z$ccov$ccov)[i])
		if(!is.na(z$interaction)&z$interaction[i]>0){
			cname <- c(cname,paste(colnames(z$ccov$ccov)[i],".t^",1:z$interaction[i],sep=""))}}
	dimnames(coef.table) <- list(cname, c("estimate", "se"))
	print.default(coef.table, digits=digits, print.gap=2)
	if(z$nlp>1){
		cat("\nCorrelation matrix of linear parameters\n")
		dimnames(z$corr) <- list(seq(1,z$nlp),seq(1,z$nlp))
		print.default(z$corr, digits=digits)}
	cat("\n(REML) Variance ",z$mse,"\n")
	if(z$np>0){
		cat("\nEstimates of nonlinear parameters\n")
		n1 <- z$arma[1]+z$arma[2]
		if(n1>1) {
			z0 <- .Fortran("roots",
				as.integer(n1),
				as.double(z$coef[1:n1]),
				r=complex(n1),
				DUP=F)
			if(any(Im(z0$r)!=0))tmp <- z0$r
			else tmp <- Re(z0$r)
			title <- "Roots"}
		if(z$arma[1]>0){
			cat("Autoregression\n")
			n2 <- z$arma[1]
			if(n2==1&&z$arma[2]==0){
				tmp <- exp(-exp(z$coef[1]))
				title <- "AR"}
			coef.table <- cbind(z$coef[1:n2],z$nlse[1:n2],
				tmp[1:n2])
			dimnames(coef.table) <- list(seq(1,z$arma[1]),
				c("estimate","se",title))
			print.default(coef.table, digits=digits, print.gap=2)}
		if(z$arma[2]>0){
			cat("Moving average\n")
			n1 <- z$arma[1]+1
			n2 <- z$arma[1]+z$arma[2]
			coef.table <- cbind(z$coef[n1:n2],z$nlse[n1:n2],
				tmp[n1:n2])
			dimnames(coef.table) <- list(seq(1,z$arma[2]),
				c("estimate","se",title))
			print.default(coef.table, digits=digits, print.gap=2)}
		if(z$arma[3]>0){
			cat("Measurement error\n")
			n1 <- z$arma[1]+z$arma[2]+1
			coef.table <- cbind(z$coef[n1],z$nlse[n1],exp(z$coef[n1]))
			dimnames(coef.table) <- list("1",c("estimate","se",""))
			print.default(coef.table, digits=digits, print.gap=2)}}
	if(z$nod>0){
		cat("Estimated unscaled factored between subject covariance matrix\n")
		n1 <- z$arma[1]+z$arma[2]+z$arma[3]+1
		n2 <- n1+z$nre-1
		tmp <- vector(mode="numeric",z$nre)
		for (i in 1:z$nre) tmp[i] <- z$are[z$position[i,1],z$position[i,2]]
		coef.table <- cbind(z$position[1:z$nre,1],z$position[1:z$nre,2],
			tmp,z$nlse[n1:n2])
		dimnames(coef.table) <- list(rep("",z$nre),
			c("i1","i2","estimate","se"))
		print.default(coef.table, digits=digits, print.gap=2)
		cat("Estimated between subject covariance matrix\n")
		tmps <- max(z$position[i,2])
		wrk <- matrix(0,ncol=tmps,nrow=tmps)
		for(j in 1:tmps)for(i in 1:j){
			for(l in 1:i)wrk[i,j] <- wrk[i,j]+z$are[l,i]*
				z$are[l,j]*z$mse
			wrk[j,i] <- wrk[i,j]}
		dimnames(wrk) <- list(1:tmps,1:tmps)
		print.default(wrk, digits=digits)}
	if(z$np>1){
		cat("\nCorrelation matrix of nonlinear parameters\n")
		dimnames(z$nlcorr) <- list(seq(1,z$np),seq(1,z$np))
		print.default(z$nlcorr, digits=digits)}
	invisible(z)}

profile.carma <- function(z, times=NULL, ccov, plotse=T){
	nccov <- ncol(z$ccov$ccov)
	if(!is.null(z$ccov$ccov)){
		if(missing(ccov))stop("Covariate values must be supplied")
		else if(nccov!=length(ccov))
			stop("Incorrect number of covariate values")}
	z$ptimes <- if(is.null(times))seq(min(z$response$times),
		max(z$response$times),length.out=25) else times
	npts <- length(z$ptimes)
	z$pred <- rep(z$beta[1],npts)
	cc <- matrix(rep(1,npts),ncol=1)
	if(z$torder>0)for(i in 1:z$torder) {
		cc <- cbind(cc,(z$ptimes-mean(z$response$times))^i)
		z$pred <- z$pred+z$beta[i+1]*(z$ptimes-mean(z$response$times))^i}
	if(nccov>0){
		jj <- z$torder+1
		for(i in 1:nccov)for(j in 1:(z$interaction[i]+1)){
			jj <- jj+1
			cc <- cbind(cc,cc[,j]*ccov[i])
			z$pred <- z$pred+z$beta[jj]*cc[,j]*ccov[i]}}
	se <- NULL
	for(i in 1:npts) se <- c(se,sqrt(cc[i,]%*%z$vbeta%*%cc[i,]))
	if(z$transform=="exp"){
		if(plotse){
			se1 <- log(z$pred+2*se)
			se2 <- log(z$pred-2*se)}
		z$pred <- log(z$pred)}
	else if(z$transform=="square"){
		if(plotse){
			se1 <- sqrt(z$pred+2*se)
			se2 <- sqrt(z$pred-2*se)}
		z$pred  <- sqrt(z$pred)}
	else if(z$transform=="sqrt"){
		if(plotse){
			se1 <- (z$pred+2*se)^2
			se2 <- (z$pred-2*se)^2}
		z$pred <- z$pred^2}
	else if(z$transform=="log"){
		if(plotse){
			se1 <- exp(z$pred+2*se)
			se2 <- exp(z$pred-2*se)}
		z$pred <- exp(z$pred)}
	else {
		if(plotse){
			se1 <- z$pred+2*se
			se2 <- z$pred-2*se}}
	if(plotse)z$pse <- cbind(se1,se2)
	class(z) <- "profile"
	invisible(z)}
