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
#     potthoff(response, x=NULL, ccov=NULL, times=NULL, torder=0, orthogonal=T)
#
#  DESCRIPTION
#
#    Function to fit the Potthoff and Roy (1964) growth curve model
# Y contains responses, X between unit design matrix, and Z within
# unit matrix

potthoff <- function(response, x=NULL, ccov=NULL, times=NULL, torder=0,
	orthogonal=T){
pcov <- function(z, s, x){
	p <- t(x)%x%z
	solve(p%*%(s%x%diag(n))%*%t(p))}
call <- sys.call()
if(is.data.frame(response))response <- as.matrix(response)
if(!is.matrix(response))stop("response must be a matrix or dataframe")
n <- nrow(response)
r <- ncol(response)
if(inherits(x,"formula")){
	formula <- x
	mt <- terms(x)
	if(is.numeric(mt[[2]]))x <- matrix(1,nrow=n)
	else {
		mf <- model.frame(mt,sys.frame(sys.parent()),na.action=na.fail)
		x <- model.matrix(mt, mf)}}
if(inherits(ccov,"formula")){
	formula <- ccov
	mt <- terms(ccov)
	if(is.numeric(mt[[2]]))ccov <- matrix(1,nrow=n)
	else {
		mf <- model.frame(mt,sys.frame(sys.parent()),na.action=na.fail)
		ccov <- model.matrix(mt, mf)}}
if(is.null(x)||nrow(x)!=n)stop(paste("x must be a matrix with",n,"rows"))
if(is.null(times))times <- 1:r
if(orthogonal)z <- rbind(rep(1,length(times)),t(contr.poly(times)))
else {
	z <- rep(1,r)
	tt <- times-sum(times)/r
	z <- rbind(z,tt)
	for(i in 2:(r-1))z <- rbind(z,tt^i)}
if(is.null(ccov)) xx <- matrix(rep(1,n),ncol=1)
else xx <- x
s <- t(response)%*%(diag(n)-xx%*%solve(t(xx)%*%xx)%*%t(xx))%*%response/n
ss <- solve(t(response)%*%(diag(n)-x%*%solve(t(x)%*%x)%*%t(x))%*%response/n)
b <- solve(t(xx)%*%xx)%*%t(xx)%*%response
if(!is.matrix(b))b <- matrix(b,nrow=1)
if(!is.numeric(torder)){
	zz <- diag(r)
	b1 <- b
	s1 <- s}
else {
	zz <- z[1:(torder+1),]
	if(!is.matrix(zz))zz <- matrix(zz,nrow=1)
	b1 <- b%*%ss%*%t(zz)%*%solve(zz%*%ss%*%t(zz))
	s1 <- response-xx%*%b1%*%zz
	s1 <- t(s1)%*%s1/n}
res <- response-xx%*%b1%*%zz
if(!is.matrix(b1))b1 <- matrix(b1,nrow=1)
like <- n*(r*(1+log(2*pi))+log(det(s1)))/2
aic <- like+length(b1)+r*(r+1)/2
pc <- pcov(zz,solve(s1),xx)
if(is.matrix(pc)){
	d <- sqrt(diag(pc))
	se <- matrix(d,ncol=ncol(b1),byrow=T)}
else se <- sqrt(pc)
corr <- pc
if(is.matrix(pc))
	for(i in 1:ncol(pc))for(j in 1:nrow(pc))corr[i,j] <- pc[i,j]/d[i]/d[j]
if(is.matrix(b1)&dim(b1)[1]>1){
	if(is.null(colnames(ccov))){
		tn <- "(Intercept)"
		if(ncol(ccov)>1)tn <- c(tn,paste("ccov",1:(ncol(ccov)-1),sep=""))
		colnames(ccov) <- tn}
	rownames(b1) <- colnames(ccov)}
else {
	b1 <- matrix(b1,nrow=1)
	rownames(b1) <- "Mean"}
if(is.matrix(se)&&dim(se)[1]>1)rownames(se) <- colnames(ccov)
else {
	se <- matrix(se,nrow=1)
	rownames(se) <- "Mean"}
if(!is.numeric(torder)){
	if(is.null(colnames(response)))colnames(b1) <- paste("t",1:ncol(b1),sep="")
	else colnames(b1) <- colnames(response)
	if(is.null(colnames(response)))colnames(se) <- paste("t",1:ncol(se),sep="")
	else colnames(se) <- colnames(response)}
d <- sqrt(diag(s1))
c1 <- s1
for(i in 2:r)for(j in 1:(i-1))c1[i,j] <- s1[i,j]/d[i]/d[j]
z <- list(
	call=call,
	y=response,
	x=x,
	time=z,
	torder=torder,
	ns=n,
	nt=n*r,
	df=n*r-(length(b1)+r*(r+1)/2),
	beta=b1,
	ccov=s1,
	maxlike=like,
	aic=aic,
	pcov=pc,
	pcorr=corr,
	se=se,
	corr=c1,
	residuals=res)
class(z) <- "potthoff"
return(z)}

coefficients.potthoff <- function(z) z$beta
deviance.potthoff <- function(z) 2*z$maxlike
residuals.potthoff <- function(z) z$residuals

print.potthoff <- function(z, digits = max(3, .Options$digits - 3)) {
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	cat("Number of subjects    ",z$ns,"\n")
	cat("Number of observations",z$nt,"\n")
	cat("-Log likelihood   ",z$maxlike,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n\n")
	if(is.null(colnames(z$beta))&&z$torder<=4){
		tn <- "(Intercept)"
		if(z$torder>0)tn <- c(tn,paste("t^",1:z$torder,sep=""))
		colnames(z$beta) <- tn
		colnames(z$se) <- tn}
	cat("Estimates of linear parameters\n")
	print(z$beta)
	cat("Standard errors\n")
	print(z$se)
	nlp <- length(z$beta)
	if(nlp>1){
		cat("\nCorrelation matrix of linear parameters\n")
		dimnames(z$pcorr) <- list(seq(1,nlp),seq(1,nlp))
		print.default(z$pcorr, digits=digits)}
	cat("\nEstimated covariance/correlation matrix\n")
	if(is.null(colnames(z$corr))){
		tn <- paste("t",1:ncol(z$corr),sep="")
		dimnames(z$corr) <- list(tn,tn)}
	print(z$corr)}
