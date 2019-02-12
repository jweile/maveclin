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
logger <- new.logger(paste0(log.dir,"cgi.log"))

#define error handler
handler <- function() {
	tryCatch({
		msg <- paste0(
			geterrmessage(),"\n",
			paste(.traceback(),collapse="\n\t--> ")
		)
		cgir::respond500("Internal Server Error\n=============================\n")
		logger$err(msg)
	},finally={
		if (die) {
			quit(save="no",status=0)
		}
	})
}
#register handler
options(
	error=handler
)

#Caching directory
cache.dir <- Sys.getenv("MAVECLIN_CACHE",unset="/var/www/maveclin/cache/")
templ.dir <- "/var/www/html/maveclin/httpdocs/"

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

scoresetPage <- function(urn, persist) {

	setDetail <- persist$getScoreset(urn)
	if (setDetail$status == "calibrated") {

		vars <- persist$getVariants(urn)
		vars <- vars[!is.na(vars$refset),-2]
		varHtml <- df2html(format(vars,digits=2))

		imgName <- paste0(urn,"_calibration.png")
		imgFrom <- paste0(cache.dir,imgName)
		imgTo <- paste0(templ.dir,imgName)
		file.copy(imgFrom,imgTo)

		respondTemplateHTML(paste0(templ.dir,"caliScoreset.html"),c(setDetail,imgTarget=imgName,varTable=varHtml))

	} else {
		respondTemplateHTML(paste0(templ.dir,"newScoreset.html"),setDetail)
	}
}

variantPage <- function(acc, persist) {

	lo2post <- function(k) exp(k)/(1+exp(k))

	varDetail <- persist$getVariantDetail(acc)
	varDetail <- lapply(varDetail,function(x) {
		if (length(x) > 1) {
			paste(x,collapse=", ")
		} else if (length(x) < 1 || is.na(x)) {
			"N/A"
		} else { 
			x
		}
	})
	# varDetail$ci <- with(varDetail,sprintf("[ %.2f ; %.2f ]",llrCIleft/log(2),llrCIright/log(2)))
	varDetail$ciLeft <- varDetail$llrCIleft
	varDetail$ciRight <- varDetail$llrCIright
	# varDetail$posterior <- with(varDetail,sprintf(
	# 	"%.2f CI: [ %.2f ; %.2f ]",
	# 	lo2post(llr),lo2post(llrCIleft),lo2post(llrCIright)
	# ))
	# varDetail$llr <- sprintf("%.2f",varDetail$llr/log(2))

	respondTemplateHTML(paste0(templ.dir,"variant.html"),varDetail)	
}

searchResultPage <- function(query, persist) {
	
	linkify <- function(x) paste0("<a href=\"index.R?q=",URLencode(x,reserved=TRUE),"\">",x,"</a>")

	scoresets <- persist$searchScoresets(query)
	scoresets$urn <- sapply(scoresets$urn,linkify)
	variants <- persist$searchVariants(query)
	variants$accession <- sapply(variants$accession,linkify)
	ssTable <- if (nrow(scoresets) > 0) paste("<h4>Scoresets</h4>",df2html(scoresets)) else ""
	varTable <- if (nrow(variants) > 0) paste("<h4>Variants</h4>",df2html(variants)) else ""
	results <- paste(ssTable,varTable)

	respondTemplateHTML(paste0(templ.dir,"searchResult.html"),list(results=results))
}


#read data from HTTP GET and POST
inputGET <- readGET()
inputPOST <- readPOST()


#if there is a query, do a search
if ("q" %in% names(inputGET)) {

	persist <- new.persistence.connection(paste0(cache.dir,"maveclin.db"))

	#analyze query if it directly matches a type of entry
	#if it's a scoreset URN
	if (grepl("^urn:mavedb:\\d{8}-\\w{1}-\\d+$",inputGET$q) && 
			persist$isKnown(inputGET$q)) {

		scoresetPage(inputGET$q, persist)

	#or if it's a variant accession?
	} else if (grepl("^urn:mavedb:\\d{8}-\\w{1}-\\d+#\\d+$",inputGET$q) &&
			persist$hasVariant(inputGET$q)) {

		variantPage(inputGET$q, persist)

	#otherwise it's an open-ended search
	} else {
		searchResultPage(inputGET$q, persist)
	}

	#FIXME: This should probably be in a "finally" clause
	persist$close()
	quit(save="no",status=0)

#otherwise, if there is no query, it's just the main page
} else {
	respondTemplateHTML(paste0(templ.dir,"search.html"),list())	
	quit(save="no",status=0)
}