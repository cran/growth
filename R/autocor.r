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
#     pergram(y)
#     corgram(y, wt=1, add=F, lty=1, ylim=NULL, xlab=NULL, ylab=NULL,
# main=NULL)
#
#  DESCRIPTION
#
#    Functions to compute and plot correlograms and periodograms

pergram <- function(y){
	ll <- length(y)
	len <- trunc(ll/2)
	fc <- fft(y)
	z <- cbind(2*pi*(1:len)/ll,
		(Re(fc[ll:(len+1)])^2+Im(fc[ll:(len+1)])^2)*ll/4/pi)
	class(z) <- "pergram"
	z}

#pergram <- function(y){
#	ll <- length(y)
#	len <- trunc(ll/2)
#	fc <- matrix(0,ncol=2,nrow=ll)
#	for(i in 1:ll){
#		a <- 2*pi*i*(0:(ll-1))/ll
#		fc[i,1] <- sum(y*sin(a))
#		fc[i,2] <- sum(y*cos(a))}
#	invisible(cbind(2*pi*(1:len)/ll,(fc[,1]^2+fc[,2]^2)*ll/4/pi))}

plot.pergram <- function(y, add=F, lty=1, xlab="Frequency", ylab="Periodogram",
	main="Periodogram", ylim=c(0,max(po[,2]))){
	if(inherits(y,"pergram"))po <- y
	else po <- pergram(y)
	if(add)lines(po[,1],po[,2],lty=lty)
	else plot(po[,1],po[,2],xlim=c(0,3.6),xlab=xlab,ylab=ylab,type="l",
		lty=lty,ylim=ylim)
	invisible(po)}

plot.cum <- function(z, ...) 
	UseMethod("plot.cum")

plot.cum.pergram <- function(y, xlab="Frequency", ylab="Periodogram",
	main="Cumulative periodogram",
	ylim=c(0,max(cpo+1.358/(a+0.12+0.11/a)))){
	if(inherits(y,"pergram")){
		len <- 2*nrow(y)
		po <- y}
	else {
		len <- length(y)
		po <- pergram(y)}
	cpo <- cumsum(po[,2])
	cpo <- cpo/max(cpo)
	a <- sqrt(len)
	pa <- 2*pi*(1:length(cpo))/len
	plot(po[,1],cpo,xlim=c(0,3.6),xlab=xlab,ylab=ylab,type="l",
		ylim=ylim,main=main)
	lines(po[,1],cpo+1.358/(a+0.12+0.11/a),lty=3)
	lines(po[,1],cpo-1.358/(a+0.12+0.11/a),lty=3)}

corgram <- function(y, wt=1, add=F, lty=1, ylim=NULL, xlab=NULL,
	ylab=NULL, main=NULL){
	len <- length(y)
	if(length(wt)==1)wt <- rep(1,len)
	wt <- wt>0
	num <- sum(wt)
	npt <- trunc(num/3)
	ave <- sum(y*wt)/num
	var <- sum((y-ave)^2*wt)/(num-1)
	co <- NULL
	for(lag in 1:npt){
		tmp <- sum((y-ave)[1:(len-lag)]*(y[(lag+1):len]-ave)*wt[1:(len-lag)])
		co <- c(co,tmp/var/sum(wt[1:(len-lag)]))}
	if(missing(ylim))ylim <- c(-(min(co)<0),1)
	if(missing(xlab))xlab <- "Lag"
	if(missing(ylab))ylab <- expression(rho)
	if(missing(main))main <- "Correlogram"
	if(add)lines(1:npt,co,lty=lty)
	else {
		plot(1:npt,co,type="l",lty=lty,ylim=ylim,xlab=xlab,
			ylab=ylab,main=main)
		limit <- rep(2/sqrt(num),2)
		lines(c(1,npt),limit,lty=3)
		if(min(co)<0)lines(c(1,npt),-limit,lty=3)}
	invisible(cbind(1:npt,co))}
