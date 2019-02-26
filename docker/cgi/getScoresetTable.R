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

#Caching directory
cache.dir <- Sys.getenv("MAVECLIN_CACHE",unset="/var/www/maveclin/cache/")
templ.dir <- "/var/www/html/maveclin/httpdocs/"

#turns an identifier into a hyperlink
linkify <- function(x) paste0("<a href=\"?q=",URLencode(x,reserved=TRUE),"\">",x,"</a>")

#turns a dataframe into an html table
df2html <- function(df) {
	stopifnot(inherits(df,"data.frame"))
	toRow <- function(x,tag="td") {
		paste0("<tr>",paste0("<",tag,">",x,"</",tag,">",collapse=""),"<tr/>")
	}
	header <- toRow(colnames(df),tag="th")
	body <- apply(df,1,toRow)
	paste("<table>",header,paste(body,collapse="\n"),"</table>",sep="\n")
}

inputPOST <- readPOST()

if (!c("urn") %in% names(inputPOST)) {
	respond400("Missing parameter(s)!")
	quit(save="no",status=0)
}

urn <- inputPOST$urn

logger$info("Querying table for",urn)

persist <- new.persistence.connection(paste0(cache.dir,"maveclin.db"))

if (!persist$isKnown(urn)) {
	respond400(paste("Unknown URN",urn))
	quit(save="no",status=0)
} 

vars <- persist$getVariants(urn)
#remove variants that are not part of a reference set (and drop the scoreset reference)
vars <- vars[!is.na(vars$refset),-2]
#linkify the accessions
vars$accession <- sapply(vars$accession,linkify)
#turn to HTML table
tableHTML <- df2html(format(vars,digits=2))

imgName <- paste0(urn,"_calibration.png")
imgFrom <- paste0(cache.dir,imgName)
imgTo <- paste0(templ.dir,"imgcache/",imgName)
invisible(file.copy(imgFrom,imgTo))

respondJSON(list(
	table=tableHTML,
	imgTarget=paste0("imgcache/",imgName)
))
