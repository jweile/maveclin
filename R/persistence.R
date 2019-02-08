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


#' Create a new persistence connection object
#' 
#' Constructor for creating a new object of type persistence.connection. 
#' 
#' Provides the following methods:
#' \begin{description}
#' \item{isKnown(urn)} checks whether a URN is known to the database
#' \item{getStatus(urn)} retrieves the status of the scoreset with the given URN
#' \item{setStatus(urn)} sets the status of the scoreset with the given URN to the given status. 
#'     Permissible values are "new", "configured", "calibrated".
#' \item{getParameters(urn)} retrieves the calibration parameters for the scoreset with the given URN
#' \item{newScoreset(urn,symbols,ensembl)} adds a new scoreset with the given URN to the database and initializes
#'     its gene symbols and ensembl gene names to the given values.
#' \item{close()} closes the database connection
#' \end{description}
#' 
#' @param dbfile the SQLite database file
#' @return a new persistence.connection object
#' @export
#' @examples
#' dbfile <- "maveclin.db"
#' persist <- new.persistence.connection(dbfile)
#' persist$isKnown("urn:mavedb:00000001-c-2")
#' persist$getStatus("urn:mavedb:00000001-c-2")
#' persist$setStatus("urn:mavedb:00000001-c-2","new")
#' persist$close()
#' 
new.persistence.connection <- function(dbfile) {

	library(DBI)
	# library(RMariaDB)
	library(RSQLite)

	# .con <- dbConnect(RMariaDB::MariaDB(), 
	# 	dbname="maveclin",
	# 	user=user, password=pwd, host=host
	# )
	.con <- dbConnect(RSQLite::SQLite(),dbfile)

	#check that the database has been properly initialized
	check <- dbGetQuery(.con,"SELECT name FROM sqlite_master WHERE type='table';")
	if (!all(c("scoresets","variants") %in% check$name)) {
		stop("Database was not initialized!")
	}

	#Scoreset URN validation regex
	urnRX <- "^urn:mavedb:\\d{8}-\\w{1}-\\d+$"
	#variant accession validation regex
	accRX <- "^urn:mavedb:\\d{8}-\\w{1}-\\d+#\\d+$"

	#checks whether dataset is known to database
	isKnown <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		res <- dbGetQuery(.con,
			sqlInterpolate(.con,
				"SELECT COUNT(*) FROM scoresets WHERE urn = ?urn;",
				urn=urn
			)
		)
		return(res[,1] > 0)
	}

	#returns the 'lastUpdated' date for the given scoreset URN
	getDate <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		res <- dbGetQuery(.con,
			sqlInterpolate(.con,
				"SELECT lastUpdated FROM scoresets WHERE urn = ?urn;",
				urn=urn
			)
		)
		if (nrow(res) < 1) {
			stop("Unknown URN ",urn)
		}
		return(as.Date(res$lastUpdated))
	}

	# returns the current dataset status for the given URN
	getStatus <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		res <- dbGetQuery(.con,
			sqlInterpolate(.con,
				"SELECT status FROM scoresets WHERE urn = ?urn;",
				urn=urn
			)
		)
		if (nrow(res) < 1) {
			stop("Unknown URN ",urn)
		}
		return(res$status)
	}

	#set the status for the given URN
	setStatus <- function(urn,status) {
		stopifnot(
			grepl(urnRX,urn),
			status %in% c("new","configured","calibrated")
		) 
		res <- dbSendStatement(.con,
			sqlInterpolate(.con,
				"UPDATE scoresets SET status = ?status WHERE urn = ?urn;",
				status=status,
				urn=urn
			)
		)
		dbClearResult(res)
	}

	#get the parameters for the scoreset
	getParameters <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		res <- dbGetQuery(.con,
			sqlInterpolate(.con,
				"SELECT symbol, ensemblGeneID, mafCutoff, flip, homozygous 
				FROM scoresets WHERE urn = ?urn;",
				urn=urn
			)
		)
		if (nrow(res) != 1) {
			stop("Unknown URN ",urn)
		}
		out <- res[1,,drop=TRUE]
		out$symbol <- strsplit(out$symbol,"\\|")[[1]]
		out$ensemblGeneID <- strsplit(out$ensemblGeneID,"\\|")[[1]]
		return(out)
	}


	#get the parameters for the scoreset
	getScoreset <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		res <- dbGetQuery(.con,
			sqlInterpolate(.con,
				"SELECT * FROM scoresets WHERE urn = ?urn;",
				urn=urn
			)
		)
		if (nrow(res) != 1) {
			stop("Unknown URN ",urn)
		}
		out <- res[1,,drop=TRUE]
		out$symbol <- strsplit(out$symbol,"\\|")[[1]]
		out$ensemblGeneID <- strsplit(out$ensemblGeneID,"\\|")[[1]]
		return(out)
	}

	#add a new scoreset 
	newScoreset <- function(urn,symbols,ensembl) {
		stopifnot(
			grepl(urnRX,urn)
			#TODO: validate symbols and ensembl
		)
		res <- dbSendStatement(.con,
			sqlInterpolate(.con,
				"INSERT INTO scoresets 
				(urn, symbol, ensemblGeneID, mafCutoff, flip, homozygous, lastUpdated, status) 
				VALUES (?urn, ?symbol, ?ensembl, 0, 0, 0, ?date, 'new');",
				urn=urn, 
				symbol=paste(symbols,collapse="|"),
				ensembl=paste(ensembl,collapse="|"),
				date=format(Sys.time(),"%Y-%m-%d")
			)
		)
		dbClearResult(res)
	}


	#enters a new set of variants from a scoreset
	newVariants <- function(urn,scores) {
		stopifnot(
			grepl(urnRX,urn),
			inherits(scores,"data.frame"),
			c("accession", "hgvs_pro", "score") %in% colnames(scores)
		)

		#determine which (if any) columns in the score table correspond to standard deviation or error
		sdCol <- which(colnames(scores) %in% c("sd","SD","stdev","STDEV","sigma","SIGMA"))
		seCol <- which(colnames(scores) %in% c("se","SE","stderr","STDERR","sem","SEM"))
		if (length(sdCol) > 1 || length(seCol) > 1) {
			stop("Ambiguous error columns!")
		}

		res <- dbSendStatement(.con,
			"INSERT INTO variants (accession, scoreset, hgvs_pro, score, sd, se) 
					VALUES ($accession, $scoreset, $hgvs_pro, $score, $sd, $se);"
		)
		dbBind(res, data.frame(
			accession = scores$accession,
			scoreset = urn,
			hgvs_pro = scores$hgvs_pro,
			score = scores$score,
			sd = if (length(sdCol)==1) scores[,sdCol] else NA,
			se = if (length(seCol)==1) scores[,seCol] else NA
		))
		if (dbGetRowsAffected(res) != nrow(scores)) {
			warning("Not all entries were updated!")
		}
		dbClearResult(res)
		#FIXME: Use dbBind instead of for-loop!
		# dbBegin(.con)
		# for (i in 1:nrow(scores)) {
		# 	res <- dbSendStatement(.con,
		# 		sqlInterpolate(.con,
					# "INSERT INTO variants (accession, scoreset, hgvs_pro, score, sd, se) 
					# VALUES (?acc, ?urn, ?hgvs, ?score, ?sd, ?se);",
		# 			acc=scores[i,"accession"],
		# 			urn=urn, 
		# 			hgvs=scores[i,"hgvs_pro"],
		# 			score=scores[i,"score"],
		# 			sd=if(length(sdCol) == 1) scores[i,sdCol] else NA,
		# 			se=if(length(seCol) == 1) scores[i,seCol] else NA
		# 		)
		# 	)
		# 	dbClearResult(res)
		# }
		# dbCommit(.con)
	}

	#retrieves the variant table for the given scoreset URN
	getVariants <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		res <- dbGetQuery(.con,
			sqlInterpolate(.con,
				"SELECT * FROM variants WHERE scoreset = ?urn",
				urn=urn
			)
		)
		return(res)
	}

	#checks whether a variant with the given accession exists
	hasVariant <- function(acc) {
		stopifnot(
			grepl(accRX,acc)
		)
		res <- dbGetQuery(.con,
			sqlInterpolate(.con,
				"SELECT COUNT(*) FROM variants WHERE accession = ?acc",
				acc=acc
			)
		)
		return(res[,1] > 0)
	}

	getVariantDetail <- function(acc) {
		stopifnot(
			grepl(accRX,acc)
		)
		res <- dbGetQuery(.con,
			sqlInterpolate(.con,
		 		"SELECT * FROM 
		 		variants INNER JOIN scoresets 
		 		ON variants.scoreset = scoresets.urn 
		 		AND variants.accession = ?acc;",
		 		acc=acc
			)
		)
		if (nrow(res) != 1) {
			stop("Unknown accession ",acc)
		}
		return(res[1,,drop=TRUE])
	}

	#applies the given calibrated scores to the variants in the database
	calibrateScores <- function(caliScores) {
		stopifnot(
			inherits(caliScores,"data.frame"),
			c("accession","llr","llrCIleft","llrCIright") %in% colnames(caliScores)
		)
		res <- dbSendStatement(.con,
			"UPDATE variants 
			SET flippedScore = $flippedScore, 
				clinsig = $clinsig, 
				maf = $maf, 
				hom = $hom, 
				refset = $set,
				llr = $llr, 
				llrCIleft = $llrCIleft, 
				llrCIright = $llrCIright 
			WHERE accession = $accession;"
		)
		dbBind(res, caliScores[c(
			"flippedScore","clinsig","maf","hom","set",
			"llr","llrCIleft","llrCIright","accession"
		)])
		if (dbGetRowsAffected(res) != nrow(caliScores)) {
			warning("Not all entries were updated!")
		}
		dbClearResult(res)

	}

	#close the db connection
	close <- function() {
		dbDisconnect(.con)
	}


	structure(list(
		isKnown=isKnown,
		getStatus=getStatus,
		getDate=getDate,
		setStatus=setStatus,
		getParameters=getParameters,
		getScoreset=getScoreset,
		newScoreset=newScoreset,
		newVariants=newVariants,
		getVariants=getVariants,
		hasVariant=hasVariant,
		getVariantDetail=getVariantDetail,
		calibrateScores=calibrateScores,
		close=close
	),class="persistence.connection")

}

