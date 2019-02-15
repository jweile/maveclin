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

##############################
# Synchronization script.
#
# Checks the contents of MaveDB via the API against a list of admin-provided parameters.
# If new data is available, a new flagged entry is added.
# If unprocessed data exists, it will be processed.

options(stringsAsFactors=FALSE)

library(rapimave)
library(hgvsParseR)
library(yogitools)
library(yogilog)
library(maveclin)

#Set up logger
logfile <- paste0(Sys.getenv("MAVECLIN_LOGS",unset="/var/www/maveclin/logs/"),"calibrationCron.log")
logger <- new.logger(logfile)
#Make sure any fatal errors will be logged before the script dies.
registerLogErrorHandler(logger)


# Before proceeding any further, check if the previous job is still running
disable <- FALSE
#When executed from cron, the script name shows up in both the /bin/sh shell wrapper and R processes!
# 2019-02-15_15:46:01 INFO: /bin/sh -c Rscript /setup/calibrationCron.R >> /var/www/maveclin/logs/cron.log 2>&1 
# 2019-02-15_15:46:01 INFO: /usr/lib/R/bin/exec/R --slave --no-restore --file=/setup/calibrationCron.R 
jobRX <- "'R --slave.*calibrationCron.R'"
njobs <- as.numeric(system(paste0("ps -eo command|grep ",jobRX,"|grep -cv grep"),intern=TRUE))
if (disable) {
	logger$info("calibration cancelled.")
	quit(save="no",status=0)
} else if (njobs > 1) {
	logger$info("Concurrent calibration detected. Aborting duplicate process.")
	detail <- system(paste0("ps -eo command|grep ",jobRX,"|grep -v grep"),intern=TRUE)
	logger$info(detail)
	quit(save="no",status=0)
}


#Performs score calibration on a given dataset
calibrate <- function(urn, persist, overrideCache=FALSE) {

	logger$info("Running calibration for ",urn)

	tryCatch({
		persist$setStatus(urn,"processing")

		params <- persist$getParameters(urn)
		scores <- persist$getVariants(urn)

		pngFile <- getCacheFile(paste0(urn,"_calibration.png"))
		dpi <- 100
		png(pngFile,7*dpi,5*dpi,res=dpi)

		caliScores <- map2bf(scores,
			ensembls=params$ensemblGeneID,symbols=params$symbol,
			minMaf=params$mafCutoff,flip=as.logical(params$flip),
			homozygous=as.logical(params$homozygous),
			drawPlot=TRUE,logger=logger
		)
		invisible(dev.off())

		logger$info("Updating database for",urn)

		persist$calibrateScores(caliScores)

		persist$setStatus(urn,"calibrated")

	},error=function(err) {
		logger$error("while calibrating dataset",urn,":\n",err)
		persist$setStatus(urn,"error")
	})

}


#Open Persistence connection
#FIXME: CRON is stupid and cannot see environment variables. Hence the hardcoding >:(
# dbfile <- getCacheFile("maveclin.db")
dbfile <- "/var/www/maveclin/cache/maveclin.db"
logger$info("Connecting to database",dbfile)
persist <- new.persistence.connection(dbfile)

pendingURNs <- persist$getPending()

lapply(pendingURNs, function(urn) {
	logger$info("Starting calibration for",urn)
	calibrate(urn,persist)
}) 

