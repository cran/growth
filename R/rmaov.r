#
#  growth : A Library of Normal Distribution Growth Curve Models
#  Copyright (C) 1999 J.K. Lindsey
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
#     rmaov(response, tvcov=NULL, ccov=NULL, analysis=T)
#
#  DESCRIPTION
#
#    Function to calculate the ANOVA table for a balanced split-plot design
#    Copyright (C) 1999 Ralf Goertz (ralf.goertz@uni-jena.de)

rmaov <- function(response, tvcov=NULL, ccov=NULL, analysis=T){
	if(!is.matrix(response)||!is.numeric(response))stop("response must be a matrix")
	if(is.null(tvcov)){
		nfac <- 1
		nrep <- nlev <- ncol(response)
		rfacname <- "repfac"}
	else if(is.vector(tvcov,mode="numeric")||is.factor(tvcov)){
		nfac <- 1
		nrep <- nlev <- length(unique(tvcov))
		rfacname <- paste(deparse(substitute(tvcov)))
		tvcov <- as.data.frame(tvcov)}
	else if(is.data.frame(tvcov)||is.matrix(tvcov)){
		nfac <- ncol(tvcov)
		if(!is.null(colnames(tvcov)))rfacname <- colnames(tvcov)
		else rfacname <- paste("repfac",1:ncol(tvcov),sep="")
		nlev <- NULL
		nrep <- 1
		for(i in 1:ncol(tvcov)){
			nlev <- c(nlev,length(unique(tvcov[,i])))
			nrep <- nrep*nlev[i]}}
	else stop("If specified, tvcov must be a vector, matrix, or dataframe containing the repeated factors")
	if(is.vector(ccov,mode="numeric")||is.factor(ccov)){
		tmp <- paste(deparse(substitute(ccov)))
		ccov <- list(ccov)
		names(ccov) <- tmp}
	else if(is.matrix(ccov))ccov <- data.frame(ccov)
	else if(is.list(ccov)&&!is.data.frame(ccov))for(i in 1:length(ccov))
		if(is.matrix(ccov[[i]]))ccov[[i]] <- as.vector(t(ccov[[i]]))
	else if(!is.null(ccov))print("ccov must be a vector, matrix, dataframe, or list of time-constant covariates")
	ncases <- ncol(response)*nrow(response)
	if(nrep!=ncol(response))stop("Number of columns in response must match product of levels")
	block <- gl(nrow(response),nrep)
	nc <- n <- ncol(response)
	response <- as.vector(t(response))
	for(i in 1:nfac){
		n <- n/nlev[i]
		com <- paste("fm",as.character(i),"<-gl(nlev[i],n,ncases",sep="")
		if(!is.null(levels(tvcov[,i])))
			com <- paste(com,",labels=levels(tvcov[,i])",sep="")
		com <- paste(com,")\n",sep="")
		eval(parse(text=com))}
	if(!is.null(ccov)){
		tccov <- ccov
		if(length(tccov[[1]])!=ncases)for(i in 1:length(tccov)){
			if(!is.factor(tccov[[i]]))
				tccov[[i]] <- as.factor(tccov[[i]])
			ccov[[i]] <- gl(nlevels(tccov[[i]]),1,ncases,labels=levels(tccov[[i]]))
			ccov[[i]][(rep(1:length(tccov[[i]]),rep(nc,length(tccov[[i]])))-1)*nc+rep(1:nc,length(tccov[[i]]))] <- tccov[[i]][rep(1:length(tccov[[i]]),rep(nc,length(tccov[[i]])))]}}
	com <- paste("b<-data.frame(response=response,block=block")
	for(i in 1:nfac)com <- paste(com,",",rfacname[i],"=fm",as.character(i),sep="")
	if(!is.null(ccov)){
		if(is.null(names(ccov)))
			names(ccov) <- paste("ccov",1:length(ccov),sep="")
		for(i in 1:length(ccov))
			com <- paste(com,",",names(ccov)[i],"=ccov$",names(ccov)[i],sep="")}
	com <- paste(com,")\n",sep="")
	eval(parse(text=com))
	if(!analysis)return(b)
	efacterm <- facterm <- paste(rfacname,collapse="*")
	if(!is.null(ccov))
		for(i in 1:length(ccov))
			facterm <- paste(facterm,names(ccov)[i],sep="*")
	com <- paste("res<-aov(response~",facterm,"+Error(block/(",efacterm,")),data=b)\n",sep="")
	eval(parse(text=com))
	res}


