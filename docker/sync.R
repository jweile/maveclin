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
# logfile <- getCacheFile("sync.log")
logfile <- paste0(Sys.getenv("MAVECLIN_LOGS",unset="/var/www/maveclin/logs/"),"sync.log")
logger <- new.logger(logfile)
#Make sure any fatal errors will be logged before the script dies.
registerLogErrorHandler(logger)

# Before proceeding any further, check if the previous sync is still running
disable <- FALSE
jobRX <- "'R --slave.*sync.R'"
njobs <- as.numeric(system(paste0("ps -eo command|grep ",jobRX,"|grep -cv grep"),intern=TRUE))
# njobs <- as.numeric(system("ps -eo command|grep sync.R|grep -cv grep",intern=TRUE))
if (disable) {
	logger$info("Synchronization cancelled.")
	quit(save="no",status=0)
} else if (njobs > 1) {
	logger$info("Concurrent synchronization detected. Aborting duplicate process.")
	quit(save="no",status=0)
}

logger$info("Starting new synchronization")


#Enters a new scoreset into the database
enterNewScoreset <- function(scoreset, persist, rmave) {

	urn <- scoreset$getURN()
	logger$info("Caching dataset ",urn)

	target <- scoreset$getTarget()
	#if it's not protein-coding, skip it
	if (target$getType() != "Protein coding") {
		logger$info("Skipping noncoding dataset ",urn)
		return()
	}

	#extract metadata for scoreset
	symbolRX <- "^[A-Z0-9-]+$|^C[0-9XY]+orf[0-9]+$"
	tname <- target$getName()
	if (!grepl(symbolRX,tname)) tname <- NA
	ensemblX <- target$getXrefEnsembl()
	ensembl <- if (!is.null(ensemblX)) ensemblX$getID() else NA	

	#load scores from MaveDB
	scores <- rmave$getScores(urn)

	#if it's not usable, skip it
	if (all(scores$hgvs_pro == "None")) {
		logger$info("Skipping dataset without amino acid changes",urn)
		return()
	}

	#insert data into database
	tryCatch({
		persist$newScoreset(urn,tname,ensembl)
		persist$newVariants(urn,scores)
	},error=function(err) {
		logger$error("Failed to cache",urn,"; Symbol:",tname,"; Ensembl:",ensembl)
	})

}


#Open API connection
rmave <- new.rapimave()
#Open Persistence connection
#FIXME: Cron is stupid and cannot see environment variables. Hence the hardcoding >:(
#Note: eventually Docker will need to set the varaibles in /etc/default/cron
# dbfile <- getCacheFile("maveclin.db")
dbfile <- "/var/www/maveclin/cache/maveclin.db"
logger$info("Connecting to database",dbfile)
persist <- new.persistence.connection(dbfile)

#Query list of scoresets
logger$info("Fetching scoreset list from MaveDB")
scoresets <- rmave$getAllScoreSets()

#iterate over scoresets and compare with our current maintained list
invisible(lapply(scoresets,function(scoreset) {

	urn <- scoreset$getURN()

	if (!persist$isKnown(urn)) {
		#it's unknown, so it'll need to be cached and tagged as "new"
		enterNewScoreset(scoreset,persist,rmave)

	} else {

		# switch(persist$getStatus(urn),
		# 	#"new" means fresh from MaveDB, not yet configured
		# 	new={
		# 		#can't do anything yet until curator provides parameters
		# 		logger$info(urn," still waiting for curation.")
		# 	},
		# 	#"configured" means the curator has provided parameters and it's ready for calibration
		# 	configured={
		# 		calibrate(urn, persist)
		# 	},
		# 	#"calibrated" means its ready to be used, but may need to be updated 
		# 	calibrated={
		# 		ageInDays <- difftime(Sys.time(),persist$getDate(urn),units="days")
		# 		if (ageInDays > 7) {
		# 			calibrate(urn, persist, overrideCache=TRUE)
		# 		}
		# 	}
		# )
		
	}

}))
