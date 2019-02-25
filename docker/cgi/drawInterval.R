#!/usr/bin/Rscript

# Copyright (C) 2018  Jochen Weile, Roth Lab
#
# This file is part of MaveClin.
#
# MaveClin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MaveClin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with MaveClin.  If not, see <https://www.gnu.org/licenses/>.

suppressMessages({
	library(cgir)
	library(RJSONIO)
	library(yogilog)
	library(maveclin)
})
setMessageSink("/dev/null")
log.dir <- Sys.getenv("MAVECLIN_LOGS",unset="/var/www/mavevis/logs/")
logger <- new.logger(paste0(log.dir,"cgi.log"),stdout=FALSE)

#define error handler
handler <- function() {
	tryCatch({
		msg <- paste0(
			geterrmessage(),"\n",
			paste(.traceback(),collapse="\n\t--> ")
		)
		cgir::respond500(paste(
			"Internal Server Error\n=============================\n",
			msg
		))
		logger$err(msg)
	},finally={
		quit(save="no",status=0)
	})
}
#register handler
options(
	error=handler
)

plotInterval <- function(lower,mid,upper,
		type=c("llr","probability"),
		resolution=50) {

	#check parameters
	stopifnot(is.numeric(lower), is.numeric(mid), is.numeric(upper))
	type <- match.arg(type,c("llr","probability"))
	if (type=="probability") {
		stopifnot(c(lower,mid,upper) >= 0, c(lower,mid,upper) <= 1)
	}

	#deterimine x-axis range
	range <- switch(type,
		llr=c(min(lower,-2),max(upper,2)),
		probability=c(0,1)
	)

	#require yogitools library for gradients
	library(yogitools)
	gradient <- yogitools::colmap(
		valStops=switch(type,llr=c(-2,0,2),probability=c(0,0.5,1)),
		colStops=c("darkolivegreen2","gold","firebrick3")
	)
	#determine gradient bins based on resolution
	xborders <- seq(range[[1]],range[[2]],length.out=resolution+1)
	xmids <- sapply(1:resolution,function(i) mean(xborders[i:(i+1)]))

	#draw plot
	op <- par(mar=c(4,0,0,0)+.1)
	plot(NA,type="n",
		xlim=range,ylim=c(0,1),
		xlab=switch(type,
			llr=expression(log[2]~"Likelihood Ratio"),
			probability="Posterior probability"
		),
		ylab="",
		axes=FALSE
	)
	axis(1)
	rect(
		xborders[-(resolution+1)],0.3,
		xborders[-1],0.7,
		border=NA,col=gradient(xmids)
	)
	segments(lower,0.5,upper,0.5,lty="dashed",lwd=2,)
	segments(c(lower,upper),0.1,c(lower,upper),0.9,lwd=2)
	segments(mid,0.2,mid,0.8,lwd=4)
	offset <- (upper-lower)/20
	segments(lower,c(0.1,0.9),lower+offset,c(0.1,0.9),lwd=2)
	segments(upper,c(0.1,0.9),upper-offset,c(0.1,0.9),lwd=2)
	par(op)
}


#Caching directory
cache.dir <- Sys.getenv("MAVECLIN_CACHE",unset="/var/www/maveclin/cache/")
templ.dir <- "/var/www/html/maveclin/httpdocs/"

inputGET <- readGET()

reqParam <- c("lower","mid","upper","type")
if (!all(reqParam %in% names(inputGET))) {
	respond400(paste(
		"Missing parameter(s):",
		paste(setdiff(reqParam,names(inputGET)),collapse=", ")
	))
	quit(save="no",status=0)
}

lower <- as.numeric(inputGET$lower)
mid <- as.numeric(inputGET$mid)
upper <- as.numeric(inputGET$upper)
if (is.na(lower) || is.na(mid) || is.na(upper)) {
	respond400(paste("lower/mid/upper must be numeric!"))
	quit(save="no",status=0)
}

if ("rebase" %in% names(inputGET) && as.logical(inputGET$rebase)) {
	lower <- lower/log(2)
	mid <- mid/log(2)
	upper <- upper/log(2)
}

outfile <- tempfile(fileext=".png")
png(outfile,320,80,bg="transparent")
plotInterval(lower,mid,upper,type=inputGET$type)
invisible(dev.off())
respondPNG(outfile)
file.remove(outfile)
