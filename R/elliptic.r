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
#     elliptic(response, model="linear", distribution="elliptic", times=NULL,
#	dose=NULL, ccov=NULL, tvcov=NULL, nest=NULL, torder=0,
#	interaction=NULL, transform="identity", link="identity",
#	autocorr="exponential", pell=NULL, preg=NULL, pvar=var(y), varfn=NULL,
#	par=NULL, pre=NULL, delta=NULL, shfn=F, common=F,
#	envir=sys.frame(sys.parent()), print.level=0, ndigit=10,
#	gradtol=0.00001, steptol=0.00001, iterlim=100, fscale=1,
#	stepmax=10*sqrt(theta%*%theta), typsiz=abs(c(theta)))
#
#  DESCRIPTION
#
#    Function to fit the multivariate elliptical distribution with
# various autocorrelation functions, one or two levels of random
# effects, and nonlinear regression.

elliptic <- function(response, model="linear", distribution="elliptic",
	times=NULL, dose=NULL, ccov=NULL, tvcov=NULL, nest=NULL,
	torder=0, interaction=NULL, transform="identity",
	link="identity", autocorr="exponential", pell=NULL, preg=NULL,
	pvar=var(y), varfn=NULL, par=NULL, pre=NULL, delta=NULL,
	shfn=F, common=F, envir=sys.frame(sys.parent()),
	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
	iterlim=100, fscale=1, stepmax=10*sqrt(theta%*%theta),
	typsiz=abs(c(theta))){
plra <- function(theta){
	if(mdl==2)mu <- mu1(theta)
	if(cvar==1)varn <- sh1(theta)
	z <- .Fortran("plra",
		theta=as.double(theta),
		like=double(1),
		dist=as.integer(dst),
		as.double(rxl),
		x=as.double(times),
		as.double(y),
		tvcov=as.double(resp$tvcov$tvcov),
		ccov=as.double(resp$ccov$ccov),
		dose=as.double(dose),
		nobs=as.integer(resp$response$nobs),
		nest=as.integer(resp$response$nest),
		lnest=as.integer(lnest),
		dev=double(nm),
		nind=as.integer(nind),
		nld=as.integer(nld),
		nxrl=as.integer(nxrl),
		np=as.integer(np),
		npell=as.integer(npell),
		npv=as.integer(npv),
		npvl=as.integer(npvl),
		nccov=as.integer(nccov),
		npvar=as.integer(npvar),
		cvar=as.integer(cvar),
		npre=as.integer(npre),
		npar=as.integer(npar),
		link=as.integer(lnk),
		torder=as.integer(torder),
		inter=as.integer(interaction),
		model=as.integer(mdl),
		ar=as.integer(ar),
		tvc=as.integer(tvc),
		beta=double(npv2),
		betacov=double(npv2*npv2),
		v=double(nld*nld),
		sigsq=double(nld),
		ey=double(nld),
		tb=double(npvl),
		as.double(mu),
		as.double(varn),
		DUP=F)
	list(like=z$like,res=z$dev,beta=z$beta,betacov=z$betacov)}
plral <- function(theta){
	if(mdl==2)mu <- mu1(theta)
	if(cvar==1)varn <- sh1(theta)
	z <- .Fortran("plra",
		theta=as.double(theta),
		like=double(1),
		dist=as.integer(dst),
		as.double(rxl),
		x=as.double(times),
		as.double(y),
		tvcov=as.double(resp$tvcov$tvcov),
		ccov=as.double(resp$ccov$ccov),
		dose=as.double(dose),
		nobs=as.integer(resp$response$nobs),
		nest=as.integer(resp$response$nest),
		lnest=as.integer(lnest),
		dev=double(nm),
		nind=as.integer(nind),
		nld=as.integer(nld),
		nxrl=as.integer(nxrl),
		np=as.integer(np),
		npell=as.integer(npell),
		npv=as.integer(npv),
		npvl=as.integer(npvl),
		nccov=as.integer(nccov),
		npvar=as.integer(npvar),
		cvar=as.integer(cvar),
		npre=as.integer(npre),
		npar=as.integer(npar),
		link=as.integer(lnk),
		torder=as.integer(torder),
		inter=as.integer(interaction),
		model=as.integer(mdl),
		ar=as.integer(ar),
		tvc=as.integer(tvc),
		beta=double(npv2),
		betacov=double(npv2*npv2),
		v=double(nld*nld),
		sigsq=double(nld),
		ey=double(nld),
		tb=double(npvl),
		as.double(mu),
		as.double(varn),
		DUP=F)
	z$like}
call <- sys.call()
if(!is.function(model)&&!inherits(model,"formula")&&!is.null(model))
	model <- match.arg(model,c("linear","logistic","pkpd"))
tmp <- c("exponential","gaussian","cauchy","spherical","IOU")
ar <- match(autocorr <- match.arg(autocorr,tmp),tmp)
tmp <- c("elliptic","Student t")
dst <- match(distribution <- match.arg(distribution,tmp),tmp)
transform <- match.arg(transform,c("identity","exp","square","sqrt","log"))
tmp <- c("identity","exp","square","sqrt","log")
lnk <- match(link <- match.arg(link,tmp),tmp)
if(common){
	if(!is.function(model))stop("with common parameters, model must be a function")
	if(!is.function(varfn))stop("with common parameters, varfn must be a function")
	pvar <- NULL}
npar <- length(par)
tvc <- ifelse(!is.null(tvcov),1,0)
respenv <- inherits(response,"repeated")
envname <- if(respenv)paste(deparse(substitute(response)))
	else NULL
if(!respenv){
	resp <- if(inherits(response,"response")) response
		else restovec(response,times,nest=nest,delta=delta)
	if(is.null(ccov))nccov <- 0
	else {
		if(!inherits(ccov,"tccov")){
			ccname <- paste(deparse(substitute(ccov)))
			if((is.matrix(ccov)&&is.null(colnames(ccov)))){
				ccname <- paste(deparse(substitute(ccov)))
				if(ncol(ccov)>1){
					tmp <- NULL
					for(i in 1:ncol(ccov))tmp <- c(tmp,paste(ccname,i,sep=""))
					ccname <- tmp}}
			ccov <- tcctomat(ccov,names=ccname)}
		nccov <- ncol(ccov$ccov)}
	if(is.null(tvcov))tvc <- 0
	else {
		if(!inherits(tvcov,"tvcov")){
			tvcname <- paste(deparse(substitute(tvcov)))
			if(is.list(tvcov)&&ncol(tvcov[[1]])>1){
				if(is.null(colnames(tvcov[[1]]))){
					tvcname <- paste(deparse(substitute(tvcov)))
					tmp <- NULL
					for(i in 1:ncol(tvcov[[1]]))tmp <- c(tmp,paste(tvcname,i,sep=""))
					tvcname <- tmp}
				else tvcname <- colnames(tvcov[[1]])}
			tvcov <- tvctomat(tvcov,names=tvcname)}
		tvc <- ncol(tvcov$tvcov)}
	resp <- rmna(response=resp, tvcov=tvcov, ccov=ccov)
	if(!is.null(ccov))rm(ccov)
	if(!is.null(tvcov))rm(tvcov)}
else {
	if(is.character(model)){
		resp <- response
		if(is.null(ccov))resp$ccov <- NULL
		else if(inherits(ccov,"formula"))
			resp$ccov$ccov <- attr(finterp(ccov,envir=response,expand=F,name=paste(deparse(substitute(response)))),"model")[,-1,drop=F]
		else stop("ccov must be a W&R formula")
		if(is.null(tvcov))resp$tvcov <- NULL
		else if(inherits(tvcov,"formula"))
			resp$tvcov$tvcov <- attr(finterp(tvcov,envir=response,name=paste(deparse(substitute(response)))),"model")[,-1,drop=F]
		else stop("tvcov must be a W&R formula")}
	else resp <- rmna(response$response)
	nccov <- if(is.null(resp$ccov$ccov)) 0
		else  ncol(resp$ccov$ccov)
	tvc <- if(is.null(resp$tvcov$tvcov)) 0
		else  ncol(resp$tvcov$tvcov)}
if((inherits(envir,"repeated")&&
	(length(resp$response$nobs)!=length(envir$response$nobs)||
	any(resp$response$nobs!=envir$response$nobs)))||
	(inherits(envir,"tvcov")&&
	(length(resp$response$nobs)!=length(envir$tvcov$nobs)||
	any(resp$response$nobs!=envir$tvcov$nobs))))
	stop("response and envir objects are incompatible")
full <- !is.null(resp$ccov$ccov)
y <- resp$response$y
times <- resp$response$times
nind <- length(resp$response$nobs)
nld <- max(resp$response$nobs)
mu <- NULL
varn <- NULL
npell <- !is.null(pell)
if(!npell&&dst==2)stop("An estimate of the degrees of freedom must be supplied")
if(npell&&pell<=0){
	if(dst==1)stop("The elliptic power parameter must be positive.")
	else stop("The degrees of freedom parameter must be positive.")}
npre <- length(pre)
mu3 <- sh3 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")){
	type <- if(respenv||inherits(envir,"repeated"))"repeated"
		else if(inherits(envir,"tccov"))"tccov"
		else "tvcov"
	if(is.null(envname))envname <- paste(deparse(substitute(envir)))
	if(inherits(model,"formula")){
		mu3 <- if(respenv)finterp(model,envir=response,name=envname)
			else finterp(model,envir=envir,name=envname)
		class(model) <- c(class(model),type)}
	else if(is.function(model)){
		tmp <- parse(text=paste(deparse(model))[-1])
		class(model) <- type
		model <- if(respenv)fnenvir(model,envir=response,name=envname)
			else fnenvir(model,envir=envir,name=envname)
		mu3 <- model
		if(respenv)attr(mu3,"model") <- tmp}
	if(inherits(varfn,"formula")){
		sh3 <- if(respenv)finterp(varfn,envir=response,name=envname)
			else finterp(varfn,envir=envir,name=envname)
		class(varfn) <- c(class(varfn),type)}
	else if(is.function(varfn)){
		tmp <- parse(text=paste(deparse(varfn))[-1])
		class(varfn) <- type
		varfn <- if(respenv)fnenvir(varfn,envir=response,name=envname)
			else fnenvir(varfn,envir=envir,name=envname)
		sh3 <- varfn
		if(respenv)attr(sh3,"model") <- tmp}}
npr <- length(preg)
if(inherits(model,"formula")){
	mu2 <- if(respenv)finterp(model,envir=response,name=envname)
		else finterp(model,envir=envir,name=envname)
	npt1 <- length(attr(mu2,"parameters"))
	if(is.matrix(attr(mu2,"model"))){
		if(all(dim(attr(mu2,"model"))==1)){
			mu1 <- function(p) p[1]*rep(1,n)
			attributes(mu1) <- attributes(mu2)
			mu2 <- NULL}}
	else {
		if(npr!=npt1){
			cat("\nParameters are ")
			cat(attr(mu2,"parameters"),"\n")
			stop(paste("preg should have",npt1,"estimates"))}
		if(is.list(preg)){
			if(!is.null(names(preg))){
				o <- match(attr(mu2,"parameters"),names(preg))
				preg <- unlist(preg)[o]
				if(sum(!is.na(o))!=length(preg))stop("invalid estimates for mu - probably wrong names")}
			else preg <- unlist(preg)}}
	if(!is.null(mu2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			mu1 <- function(p) mu2(p)[cv]
			attributes(mu1) <- attributes(mu2)}
		else {
			mu1 <- mu2
			rm(mu2)}}}
else if(is.function(model))mu1 <- model
else mu1 <- NULL
if(!is.null(mu1)&&is.null(attributes(mu1))){
	attributes(mu1) <- if(is.function(model)){
		if(!inherits(model,"formulafn")){
			if(respenv)attributes(fnenvir(model,envir=response))
			else attributes(fnenvir(model,envir=envir))}
		else attributes(model)}
		else {
			if(respenv)attributes(fnenvir(mu1,envir=response))
			else attributes(fnenvir(mu1,envir=envir))}}
nlp <- if(is.function(mu1))length(attr(mu1,"parameters"))
	else if(is.null(mu1))NULL
	else npt1
if(!is.null(nlp)&&!common&&nlp!=npr)
	stop(paste("preg should have",nlp,"initial estimates"))
if(is.function(mu1))mdl <- 2
else if(model=="linear"){
	mdl <- 1
	torder <- torder+1}
else if(model=="logistic")mdl <- 3
else if(model=="pkpd")mdl <- 4
npvar <- length(pvar)
n1 <- if(common) 1 else length(preg)+1
n2 <- length(c(preg,pvar))
if(inherits(varfn,"formula")){
	sh2 <- if(respenv)finterp(varfn,envir=response,start=n1,name=envname)
		else finterp(varfn,envir=envir,start=n1,name=envname)
	npt2 <- length(attr(sh2,"parameters"))
	if(is.matrix(attr(sh2,"model"))){
		if(all(dim(attr(sh2,"model"))==1)){
			sh1 <- function(p) p[n1]*rep(1,n)
			attributes(sh1) <- attributes(sh2)
			sh2 <- NULL}}
	else {
		if(npvar!=npt2){
			cat("\nParameters are ")
			cat(attr(sh2,"parameters"),"\n")
			stop(paste("pvar should have",npt2,"estimates"))}
		if(is.list(pvar)){
			if(!is.null(names(pvar))){
				o <- match(attr(sh2,"parameters"),names(pvar))
				pvar <- unlist(pvar)[o]
				if(sum(!is.na(o))!=length(pvar))stop("invalid estimates for varfn - probably wrong names")}
			else pvar <- unlist(pvar)}}
	if(!is.null(sh2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			sh1 <- function(p) sh2(p)[cv]
			attributes(sh1) <- attributes(sh2)}
		else {
			sh1 <- sh2
			rm(sh2)}}}
else if(is.function(varfn))sh1 <- if(shfn) function(p) varfn(p[n1:n2],mu1(p))
		else function(p) varfn(p[n1:n2])
else sh1 <- NULL
if(!is.null(sh1)&&is.null(attributes(sh1)))
	attributes(sh1) <- if(is.function(varfn)){
		if(!inherits(varfn,"formulafn")){
			if(respenv)attributes(fnenvir(varfn,envir=response))
			else attributes(fnenvir(varfn,envir=envir))}
		else attributes(varfn)}
		else {
			if(respenv)attributes(fnenvir(sh1,envir=response))
			else attributes(fnenvir(sh1,envir=envir))}
nlp <- if(is.function(varfn))length(attr(sh1,"parameters"))-shfn
	else if(is.null(varfn))NULL
	else if(!inherits(varfn,"formula"))stop("varfn must be a function or formula")
	else npt2
if(!is.null(nlp)&&!common&&nlp!=npvar)
	stop(paste("pvar should have",nlp,"initial estimates"))
if(common){
	nlp <- length(unique(c(attr(mu1,"parameters"),attr(sh1,"parameters"))))
	if(nlp!=npr)stop(paste("with a common parameter model, preg should contain",nlp,"estimates"))}
if(mdl==1&&!is.null(interaction)){
	if(length(interaction)!=nccov)
		stop(paste(nccov,"interactions with time must be specified"))
	else if(any(interaction>torder-1))
		stop(paste("Interactions can be at most of order ",torder-1))
	else interaction <- interaction+1}
else interaction <- rep(1,nccov)
if(mdl==4&&tvc==0&&is.null(dose))stop("Doses required for PKPD model")
if(!is.null(resp$response$nest))lnest <- max(resp$response$nest)
else {
	lnest <- 0
	resp$response$nest <- rep(1,length(y))}
if(mdl==2&&length(mu1(preg))!=sum(resp$response$nobs))
	stop("The mean function must provide an estimate for each observation")
if(!is.null(times)){
	if(mdl==1){
		ave <- mean(times)
		times <- times-ave}
	else ave <- 0}
else {
	if(!is.null(par))stop("No times. AR cannot be fitted")
	if(torder>1)stop("No times. Time trends cannot be fitted.")
	ave <- 0}
if(full){
	rxl <- NULL
	for(i in 1:nind){
		tmp <- 1
		if(full)for(j in 1:ncol(resp$ccov$ccov))tmp <- tmp+resp$ccov$ccov[i,j]*2^(j-1)
		rxl <- c(rxl,tmp)}
	nxrl <- max(rxl)
	p <- preg}
else {
	nxrl <- 1
	rxl <- matrix(1,ncol=1,nrow=nind)
	p <- NULL}
nm <- sum(resp$response$nobs)
if(transform=="identity")jacob <- 0
else if(transform=="exp"){
	jacob <- -sum(y)
	y <- exp(y)}
else {
	if(any(y<0))stop("Nonpositive response values: invalid transformation")
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
if(!is.null(resp$response$delta)){
	if(length(resp$response$delta)==1)jacob <- jacob-nm*log(resp$response$delta)
	else jacob <- jacob -sum(log(resp$response$delta))}
if(mdl==1){
	if(lnk==1&&is.null(varfn)){
		npvl <- torder+sum(interaction)+tvc
		npv <- 0}
	else {
		npv <- torder+sum(interaction)+tvc
		npvl <- 0
		if(is.null(preg)){
			cat("Initial regression estimates must be supplied ")
			if(lnk!=1)stop(paste("for the linear model with ",link,"link"))
			else stop("for the linear model when there is a variance function")}}}
else {
	npvl <- 0
	if(mdl==2)npv <- length(preg)
	else if(mdl==3)npv <- if(tvc==0)4*nxrl else 4+nxrl-1
	else if(mdl==4){
		if(!is.function(varfn))npv <- 2+nxrl
		else npv <- 3}}
npv2 <- max(npv,npvl)
if(length(preg)!=npv)
	stop(paste(npv,"initial parameter estimates for the mean regression must be supplied"))
if(is.null(varfn))cvar <- 0
else if(is.function(sh1)){
	cvar <- 1
	if(length(sh1(if(common)preg else pvar))!=sum(resp$response$nobs))
		stop("The variance function or formula must provide an estimate for each observation")}
else if(varfn=="identity"){
	if(any(y<0))warning("identity variance function not recommended with negative responses")
	cvar <- 2}
else if(varfn=="square")cvar <- 3
else stop("Unknown variance function: choices are identity and square")
if(!(mdl==4&npvar==4)&&!cvar==1){
	if(any(pvar<=0))stop("All variance parameters must be positive")
	pvar <- log(pvar)}
theta <- c(preg,pvar)
if(npre>0){
	if(any(pre<=0))stop("All variance components must be positive")
	theta <- c(theta,log(pre))}
if(npar>0){
	if(par<=0|(par>=1&ar!=5))
		stop("Estimate of autocorrelation must lie between 0 and 1")
	theta <- if(ar!=5)c(theta,log(par/(1-par)))
		else c(theta,log(par))}
if(npell)theta <- c(theta,log(pell))
np <- npv+npvl*(mdl==1)+npvar+npre+npar+npell
if(fscale==1)fscale <- plral(theta)
z0 <- nlm(plral, p=theta, hessian=T, print.level=print.level,
	typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
p <- z0$estimate
z <- plra(p)
like <- z$like+nm*log(pi)/2+jacob
pred <- y-z$res
if(np-npvl==1){
	nlcov <- 1/z0$hessian
	nlse <- sqrt(nlcov)}
else {
	if(any(is.na(z0$hessian)))a <- 0
	else a <- qr(z0$hessian)$rank
	if(a==np-npvl)nlcov <- solve(z0$hessian)
	else nlcov <- matrix(NA,ncol=np-npvl,nrow=np-npvl)
	nlse <- sqrt(diag(nlcov))}
if(length(nlse)>1)nlcorr <- nlcov/(nlse%o%nlse)
else nlcorr <- as.matrix(nlcov/nlse^2)
dimnames(nlcorr) <- list(seq(1,np-npvl),seq(1,np-npvl))
betase <- betacorr <- NULL
if(mdl==1){
	if(lnk==1&&cvar==0){
		betacov <- matrix(z$betacov,ncol=npvl)
		bt <- exp(p[length(p)])
		if(npell>0)betacov <- betacov/(npvl*gamma(npvl/(2*bt)))*
			2^(1/bt)*gamma((npvl+2)/(2*bt))
		if(npvl==1)betase <- sqrt(betacov)
		else if(npvl>1)betase <- sqrt(diag(betacov))
		if(npvl>1){
			betacorr <- betacov/(betase%o%betase)
			dimnames(betacorr) <- list(seq(1,npvl),seq(1,npvl))}}
	else {
		betase <- nlse[1:npv]
		betacov <- nlcov[1:npv,1:npv]
		if(npvl>1){
			betacorr <- nlcorr[1:npv,1:npv]
			dimnames(betacorr) <- list(seq(1,npv),seq(1,npv))}}}
else betacov <- NULL
sigsq <- p[(npv+1):(npv+npvar)]
if(!(mdl==4&&npvar==4)&&!common)sigsq <- exp(sigsq)
if(npre>0)tausq <- exp(p[(npv+npvar+1):(npv+npvar+npre)])
else tausq <- 0
if(npar>0){
	rho <- exp(p[npv+npvar+npre+1])
	if(ar!=5)rho <- rho/(1+rho)}
else rho <- 0
if(is.null(pell)&&!common&&length(sigsq)==1&&((npar>0&&
	autocorr=="exponential")||(npar==0&&length(pre)==1))){
	nt <- sum(resp$response$nobs)
	mse <- sigsq*length(y)/(length(y)-npv-npvl*(mdl==1))
	coef <- NULL
	if(npar>0)coef <- c(coef,log(-log(rho)))
	if(length(pre)==1)coef <-
	c(coef,sqrt(exp(p[npv+npvar+1])/mse))
	z1 <- .Fortran("resid2",
		np=as.integer(npar+length(pre)),
		par=as.double(coef),
		ave=as.double(ave),
		pred=as.double(pred),
		rpred=double(nt),
		sdr=double(nt),
		res=double(nt),
		y=as.double(y),
		mse=as.double(mse),
		ns=as.integer(length(resp$response$nobs)),
		nt=as.integer(nt),
		model=as.integer(npar),
		t=as.double(times),
		nobs=as.integer(resp$response$nobs),
		nod=as.integer(length(pre)),
		p=double(length(p)),
		DUP=F)
	rpred <- z1$rpred
	sdr <- z1$sdr
	rres <- z1$res}
else rpred <- sdr <- rres <- NULL
if(!is.null(mu3))mu1 <- mu3
if(!is.null(sh3))sh1 <- sh3
z <- list(
	call=call,
	model=model,
	distribution=distribution,
	mu1=mu1,
	varfn=varfn,
	common=common,
	sh1=sh1,
	shfn=shfn,
	autocorr=autocorr,
	response=resp$response,
	transform=transform,
	torder=torder-1,
	interaction=interaction-1,
	ccov=resp$ccov,
	tvcov=resp$tvcov,
	full=full,
	link=link,
	maxlike=like,
	aic=like+np,
	df=nm-np,
	np=np,
	npell=npell,
	npv=npv,
	npvl=npvl,
	npvar=npvar,
	npar=npar,
	npre=npre,
	coefficients=p,
	beta=z$beta,
	betacov=betacov,
	betacorr=betacorr,
	betase=betase,
	nlse=nlse,
	nlcov=nlcov,
	nlcorr=nlcorr,
	sigsq=sigsq,
	tausq=tausq,
	rho=rho,
	residuals=z$res,
	pred=pred,
	rpred=rpred,
	sdrpred=sdr,
	rresiduals=rres,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z) <- "elliptic"
if(!is.null(rpred))class(z) <- c(class(z),"recursive")
return(z)}

coefficients.elliptic <- function(z) z$coefficients
deviance.elliptic <- function(z) 2*z$maxlike
residuals.elliptic <- function(z) z$residuals

print.elliptic <- function(z, digits = max(3, .Options$digits - 3)) {
	if(!is.null(z$ccov$ccov))nccov <- ncol(z$ccov$ccov)
	else nccov <- 0
	if(z$npell==0)cat("\nMultivariate normal distribution\n")
	else {
		if(z$distribution=="elliptical")cat("\nMultivariate elliptically-contoured distribution\n")
		else cat("\nMultivariate Student t distribution\n")}
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("Number of subjects    ",length(z$response$nobs),"\n")
	cat("Number of observations",length(z$response$y),"\n")
	if(!inherits(z$mu1,"formulafn")){
		if(z$model=="linear"){
			if(z$torder>0){
				cat("\nPolynomial model\n")
				cat("Times centred at  ",mean(z$response$times),"\n\n")}
			else cat("\nLinear model\n\n")}
		else if(z$model=="logistic")cat("\nGeneralized logistic model\n\n")
		else if(z$model=="pkpd")cat("\nPKPD model\n\n")}
	cat("Transformation:",z$trans,"\n")
	cat("Link function: ",z$link,"\n\n")
	cat("-Log likelihood   ",z$maxlike,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n\n")
	tvc <- !is.null(z$tvcov)
	if(!inherits(z$mu1,"formulafn")){
		cat("Location parameters\n")
		if(inherits(z$ccov$linear,"formula"))
			cat("Linear part: ",deparse(z$ccov$linear),sep="\n")}
	if(inherits(z$mu1,"formulafn")){
		if(z$common)cat("Location function\n")
		else cat("Parameters for location function\n")
		if(!is.null(attr(z$mu1,"formula")))cat(deparse(attr(z$mu1,"formula")),sep="\n")
		else if(!is.null(attr(z$mu1,"model"))){
			t <- deparse(attr(z$mu1,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		cname <- if(is.matrix(attr(z$mu1,"model")))
				colnames(attr(z$mu1,"model"))
			else attr(z$mu1,"parameters")
		coef.table <- cbind(z$coef[1:z$npv], z$nlse[1:z$npv])
		if(!z$common){
			dimnames(coef.table) <- list(cname, c("estimate", "se"))
			print.default(coef.table,digits=digits,print.gap=2)}}
	else if(z$model=="linear"){
		if(z$npell==1&&z$link=="identity"&&!is.function(z$varfn))
			cat("(Approximate s.e.)\n")
		tord <- z$torder+1+nccov
		if(tvc)tord <- tord+1
		coef.table <- cbind(z$beta,z$betase)
		cname <- "(Intercept)"
		if(z$torder>0)cname <- c(cname,paste("t^",1:z$torder,sep=""))
		if(nccov>0)for(i in 1:nccov){
			cname <- c(cname,colnames(z$ccov$ccov)[i])
			if(z$interaction[i]>0){
				cname <- c(cname,paste(colnames(z$ccov$ccov)[i],".t^",1:z$interaction[i],sep=""))}}
		if(tvc)cname <- c(cname,colnames(z$tvcov$tvcov))
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table, digits=digits, print.gap=2)
		if(z$npvl>1&&z$link=="identity"&&is.null(z$varfn)){
			cat("\nCorrelation matrix of linear parameters\n")
			print.default(z$betacorr, digits=digits)}}
	else if(z$model=="logistic"){
		coef.table <- cbind(z$coef[1:4], z$nlse[1:4])
		if(tvc)cname  <- c("kappa1","kappa3","kappa4","beta")
		else cname <- c("kappa1","kappa2","kappa3","kappa4")
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table, digits=digits, print.gap=2)
		if(z$full){
			if(tvc){
				coef.table <- cbind(z$coef[5:z$npv],z$nlse[5:z$npv])
				cname <- colnames(z$ccov$ccov)
				dimnames(coef.table) <- list(cname,c("estimate","se"))
				print.default(coef.table,digits=digits,print.gap=2)}
			else {
				for(i in 1:nccov){
					cat("   ",colnames(z$ccov$ccov)[i],"\n")
					coef.table <- cbind(z$coef[(i*4+1):((i+1)*4)],z$nlse[(i*4+1):((i+1)*4)])
					dimnames(coef.table) <- list(cname, c("estimate", "se"))
					print.default(coef.table, digits=digits, print.gap=2)}}}}
	else if(z$model=="pkpd"){
		coef.table <- cbind(z$coef[1:3], z$nlse[1:3])
		cname <- c("log k_a","log k_e","log V")
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table, digits=digits, print.gap=2)
		if(z$full){
			for(i in 1:nccov){
				cat("   ",colnames(z$ccov$ccov)[i],"\n")
				coef.table <- cbind(z$coef[i+3],z$nlse[i+3])
				dimnames(coef.table) <- list(cname[3], c("estimate", "se"))
				print.default(coef.table, digits=digits, print.gap=2)}}}
	if(inherits(z$sh1,"formulafn")){
		if(z$common)cat("\nVariance function\n")
		else {
			cat("\nVariance function parameters\n")
			cname <- NULL}
		if(!is.null(attr(z$sh1,"formula")))cat(deparse(attr(z$sh1,"formula")),sep="\n")
		else if(!is.null(attr(z$sh1,"model"))){
			t <- deparse(attr(z$sh1,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		cname <- c(cname,if(is.matrix(attr(z$sh1,"model")))
				colnames(attr(z$sh1,"model"))
			else attr(z$sh1,"parameters")[1:(length(attr(z$sh1,"parameters"))-z$shfn)])
		if(!z$common)coef.table <- cbind(z$coef[(z$npv+1):(z$npv+z$npvar)], z$nlse[(z$npv+1):(z$npv+z$npvar)])
		else cat("\nCommon parameters\n")
		dimnames(coef.table) <- list(unique(cname),c("estimate","se"))
		print.default(coef.table, digits=digits, print.gap=2)}
	else if(!is.function(z$model)&&!inherits(z$mu1,"formulafn")&&z$model=="pkpd"&&z$npvar==4){
		cat("\nVariance parameters\n")
		coef.table <- cbind(z$sigsq, z$nlse[(z$npv+1):(z$npv+4)])
		cname <- c("log k_a","log k_e","log V","power")
		dimnames(coef.table) <- list(cname, c("estimate", "se"))
		print.default(coef.table, digits=digits, print.gap=2)}
	else {
		if(z$npvar==1)cat("\nVariance\n")
		else cat("\nVariance parameters\n")
		if(!is.null(z$varfn)){
			cat(z$varfn,"function of the location parameter\n")
			vname <- "factor"
			if(z$npvar==2)vname <- c("constant",vname)}
		else if(z$npvar>1)
			vname <- c("(Intercept)",paste("t^",1:(z$npvar-1),sep=""))
		else vname <- ""
		coef.table <- cbind(z$coef[(z$npv+1):(z$npv+z$npvar)],
			z$nlse[(z$npv+1):(z$npv+z$npvar)],z$sigsq)
		dimnames(coef.table) <- list(vname,c("estimate","se","sigsq"))
		print.default(coef.table, digits=digits, print.gap=2)}
	if(z$npre>0){
		cat("\nVariance components\n")
		coef.table <- cbind(z$coef[(z$npv+z$npvar+1):(z$npv+z$npvar+z$npre)],z$nlse[(z$npv+z$npvar+1):(z$npv+z$npvar+z$npre)],z$tausq)
		if(z$npre==1)cname <- "tausq"
		else cname <- c("Level 1","Level 2")
		dimnames(coef.table) <- list(cname, c("estimate","se",""))
		print.default(coef.table, digits=digits, print.gap=2)}
	if(z$rho!=0){
		cat("\n",z$autocorr," autocorrelation\n",sep="")
		coef.table <- cbind(z$coef[z$npv+z$npvar+z$npre+1],z$nlse[z$npv+z$npvar+z$npre+1],z$rho)
		dimnames(coef.table) <- list("rho",c("estimate","se",""))
		print.default(coef.table, digits=digits, print.gap=2)}
	if(z$npell==1){
		coef.table <- cbind(z$coef[z$np-z$npvl],z$nlse[z$np-z$npvl],exp(z$coef[z$np-z$npvl]))
		if(z$distribution=="elliptic"){
			cat("\nElliptic distribution power parameter\n")
			dimnames(coef.table) <- list("",c("estimate","se","power"))}
		else {
			cat("\nDegrees of freedom parameter\n")
			dimnames(coef.table) <- list("",c("estimate","se","d.f."))}
		print.default(coef.table, digits=digits, print.gap=2)}
	if(z$np-z$npvl>1){
		cat("\nCorrelation matrix of nonlinear parameters\n")
		print.default(z$nlcorr, digits=digits)}
	invisible(z)}
