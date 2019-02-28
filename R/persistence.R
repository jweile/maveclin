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
#' Constructor for creating a new object of type persistence.connection. Each operation uses a mutex
#' lock to secure the database. 
#' 
#' Provides the following methods:
#' \begin{description}
#' \item{isKnown(urn)} checks whether a URN is known to the database
#' \item{isLocked()} checks whether the database is currently locked
#' \item{getDate(urn)} retrives the lastModified date for the scoreset with the given URN
#' \item{getStatus(urn)} retrieves the status of the scoreset with the given URN
#' \item{setStatus(urn)} sets the status of the scoreset with the given URN to the given status. 
#'     Permissible values are "new","pending","processing","calibrated",and "error".
#' \item{getParameters(urn)} retrieves the calibration parameters for the scoreset with the given URN
#' \item{setParameters(urn,symbol,ensemblGeneID,mafCutoff,flip,homozygous)} changes the parameters for 
#'     the dataset with the given URN to the given values.
#' \item{getScoreset(urn)} retrieves scoreset entry with the given URN
#' \item{newScoreset(urn,symbols,ensembl)} adds a new scoreset with the given URN to the database and initializes
#'     its gene symbols and ensembl gene names to the given values.
#' \item{getGoldStandard(urn)} retrives the set of gold standard variants for the given urn
#' \item{getVariants(urn)} retrives all variants for the given URN
#' \item{hasVariant(acc)} checks whether the variant with the given accession exists
#' \item{getVariantDetail(acc)} retrives all adata on the given variant
#' \item{calibrateScores(caliScores)} applies the given calibrated scores to the variants in the database
#'     caliScores should be a \code{data.frame} with at the columns: 
#'     accession, flippedScore, clinsig, maf, hom, set, llr, llrCIleft and llrCIright
#' \item{searchScoresets(query)} retrives a table of scoresets for which any field matches the query
#' \item{searchVariants(query)} retrieves a table of variants for which the hgvs string matches the query
#' \item{getPending()} retrieves a list of URNs of scoresets that are currently in the 'pending' state.
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
	library(flock)

	#Scoreset URN validation regex
	urnRX <- "^urn:mavedb:\\d{8}-\\w{1}-\\d+$"
	#variant accession validation regex
	accRX <- "^urn:mavedb:\\d{8}-\\w{1}-\\d+#\\d+$"
	#HGNC gene symbol validation regex
	symbolRX <- "^[A-Z0-9-]+$|^C[0-9XY]+orf[0-9]+$"
	#Ensembl Gene ID regex
	ensemblRX <- "^ENS[A-Z]+[0-9]{11}|[A-Z]{3}[0-9]{3}[A-Za-z](-[A-Za-z])?|CG[0-9]+|[A-Z0-9]+\\.[0-9]+|YM[A-Z][0-9]{3}[a-z][0-9]$"

	#internal fields for connection and lock objects
	.con <- NULL
	.lockfile <- sub(".db$",".lck",dbfile)
	.lock <- NULL

	#acquire mutex lock and open DB connection
	.connect <- function() {
		# .con <- dbConnect(RMariaDB::MariaDB(), 
		# 	dbname="maveclin",
		# 	user=user, password=pwd, host=host
		# )
		.lock <<- lock(.lockfile)
		.con <<- dbConnect(RSQLite::SQLite(),dbfile)
	}
	#disconnect and release mutex lock
	.disconnect <- function() {
		dbDisconnect(.con)
		unlock(.lock)
	}


	#check that the database has been properly initialized
	preflight <- function() {
		.connect()
		tryCatch({
			check <- dbGetQuery(.con,"SELECT name FROM sqlite_master WHERE type='table';")
		},finally={
			.disconnect()
		})	
		if (!all(c("scoresets","variants") %in% check$name)) {
			stop("Database was not initialized!")
		}
	}
	preflight()

	#checks whether the database is currently locked
	isLocked <- function() {
		is.locked(.lock)
	}

	#checks whether dataset is known to database
	isKnown <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		.connect()
		tryCatch({
			res <- dbGetQuery(.con,
				sqlInterpolate(.con,
					"SELECT COUNT(*) FROM scoresets WHERE urn = ?urn;",
					urn=urn
				)
			) 
			return(res[,1] > 0)
		},finally={
			.disconnect()
		})
		
	}

	#returns the 'lastUpdated' date for the given scoreset URN
	getDate <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		.connect()
		tryCatch({
			res <- dbGetQuery(.con,
				sqlInterpolate(.con,
					"SELECT lastUpdated FROM scoresets WHERE urn = ?urn;",
					urn=urn
				)
			)
		},finally={
			.disconnect()
		})
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
		.connect()
		tryCatch({
			res <- dbGetQuery(.con,
				sqlInterpolate(.con,
					"SELECT status FROM scoresets WHERE urn = ?urn;",
					urn=urn
				)
			)
		},finally={
			.disconnect()
		})
		if (nrow(res) < 1) {
			stop("Unknown URN ",urn)
		}
		return(res$status)
	}

	#set the status for the given URN
	setStatus <- function(urn,status) {
		stopifnot(
			grepl(urnRX,urn),
			status %in% c("new","pending","processing","calibrated","error")
		) 
		.connect()
		tryCatch({
			res <- dbSendStatement(.con,
				sqlInterpolate(.con,
					"UPDATE scoresets SET status = ?status WHERE urn = ?urn;",
					status=status,
					urn=urn
				)
			)
			dbClearResult(res)
		},finally={
			.disconnect()
		})
	}

	#get the parameters for the scoreset
	getParameters <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		.connect()
		tryCatch({
			res <- dbGetQuery(.con,
				sqlInterpolate(.con,
					"SELECT symbol, ensemblGeneID, mafCutoff, flip, homozygous 
					FROM scoresets WHERE urn = ?urn;",
					urn=urn
				)
			)
		},finally={
			.disconnect()
		})
		if (nrow(res) != 1) {
			stop("Unknown URN ",urn)
		}
		out <- res[1,,drop=TRUE]
		out$symbol <- strsplit(out$symbol,"\\|")[[1]]
		out$ensemblGeneID <- strsplit(out$ensemblGeneID,"\\|")[[1]]
		out$flip <- as.logical(out$flip)
		out$homozygous <- as.logical(out$homozygous)
		return(out)
	}

	#change the parameters for a given scoreset and set its state to "pending"
	setParameters <- function(urn,symbol,ensemblGeneID,mafCutoff,flip,homozygous) {
		#assert that all parameters are valid
		stopifnot(
			#urn must be valid URN
			length(urn) == 1, grepl(urnRX,urn),
			#symbol must be valid gene symbols
			length(symbol) > 0, grepl(symbolRX,symbol),
			#ensembl must be valid ensembl gene ids
			length(ensemblGeneID) > 0, grepl(ensemblRX,ensemblGeneID),
			#mafCutoff must be numeric
			length(mafCutoff) == 1, inherits(mafCutoff,"numeric"), !is.na(mafCutoff),
			#flip must be a boolean value
			length(flip) == 1, inherits(flip,"logical"), !is.na(flip),
			#homozygous must be a boolean value
			length(homozygous) == 1, inherits(homozygous,"logical"), !is.na(homozygous)
		) 
		if (length(symbol) > 1) {
			symbol <- paste(symbol,collapse="|")
		}
		if (length(ensemblGeneID) > 1) {
			ensemblGeneID <- paste(ensemblGeneID,collapse="|")
		}
		.connect()
		tryCatch({
			res <- dbSendStatement(.con,
				"UPDATE scoresets 
				SET symbol = $symbol, 
				ensemblGeneID = $ensemblGeneID,
				mafCutoff = $mafCutoff,
				flip = $flip,
				status = 'pending'
				WHERE urn = $urn;"
			)
			dbBind(res,list(
				symbol=symbol,
				ensemblGeneID=ensemblGeneID,
				mafCutoff=mafCutoff,
				flip=flip,
				urn=urn
			))
			rowsAffected <- dbGetRowsAffected(res)
			dbClearResult(res)
		},finally={
			.disconnect()
		})
		if (rowsAffected != 1) {
			stop("DB Update failed!")
		}
	}


	#get the parameters for the scoreset
	getScoreset <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		.connect()
		tryCatch({
			res <- dbGetQuery(.con,
				sqlInterpolate(.con,
					"SELECT * FROM scoresets WHERE urn = ?urn;",
					urn=urn
				)
			)
		},finally={
			.disconnect()
		})
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
			#urn must be a valid URN
			length(urn) == 1, grepl(urnRX,urn),
			#must provide same number of gene symbols and ensembl IDs
			length(symbols) == length(ensembl),
			#symbols must either be NA or valid gene symbols
			length(symbols) > 0, all(is.na(symbols)) || grepl(symbolRX,symbols),
			#ensembl must either be NA or valid ensembl gene IDs
			length(ensembl) > 0, all(is.na(ensembl)) || grepl(ensemblRX,ensembl)
		)
		if (length(symbols) > 1) {
			symbols <- paste(symbols,collapse="|")
		}
		if (length(ensembl) > 1) {
			ensembl <- paste(ensembl,collapse="|")
		}
		.connect()
		tryCatch({
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
		},finally={
			.disconnect()
		})
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

		.connect()
		tryCatch({
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
		},finally={
			.disconnect()
		})
	}

	#retrieves the table of gold standard variants for the given scoreset URN
	getGoldStandard <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		.connect()
		tryCatch({
			res <- dbGetQuery(.con,
				sqlInterpolate(.con,
					"SELECT * FROM variants WHERE scoreset = ?urn AND refset IS NOT NULL;",
					urn=urn
				)
			)
		},finally={
			.disconnect()
		})
		return(res)
	}

	#retrieves the variant table for the given scoreset URN
	getVariants <- function(urn) {
		stopifnot(
			grepl(urnRX,urn)
		)
		.connect()
		tryCatch({
			res <- dbGetQuery(.con,
				sqlInterpolate(.con,
					"SELECT * FROM variants WHERE scoreset = ?urn;",
					urn=urn
				)
			)
		},finally={
			.disconnect()
		})
		return(res)
	}

	#checks whether a variant with the given accession exists
	hasVariant <- function(acc) {
		stopifnot(
			grepl(accRX,acc)
		)
		.connect()
		tryCatch({
			res <- dbGetQuery(.con,
				sqlInterpolate(.con,
					"SELECT COUNT(*) FROM variants WHERE accession = ?acc;",
					acc=acc
				)
			)
		},finally={
			.disconnect()
		})
		return(res[,1] > 0)
	}

	#retrieves a list of all known data on the variant with the given accession
	getVariantDetail <- function(acc) {
		stopifnot(
			grepl(accRX,acc)
		)
		.connect()
		tryCatch({
			res <- dbGetQuery(.con,
				sqlInterpolate(.con,
			 		"SELECT * FROM 
			 		variants INNER JOIN scoresets 
			 		ON variants.scoreset = scoresets.urn 
			 		AND variants.accession = ?acc;",
			 		acc=acc
				)
			)
		},finally={
			.disconnect()
		})
		if (nrow(res) != 1) {
			stop("Unknown accession ",acc)
		}
		return(res[1,,drop=TRUE])
	}

	#applies the given calibrated scores to the variants in the database
	#caliScores data.frame with at least accession, flippedScore, clinsig, maf, hom, set, llr, llrCIleft and llrCIright columns
	calibrateScores <- function(caliScores) {
		stopifnot(
			inherits(caliScores,"data.frame"),
			c("accession","llr","llrCIleft","llrCIright") %in% colnames(caliScores)
		)
		.connect()
		tryCatch({
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
		},finally={
			.disconnect()
		})

	}

	#retrieves a table of scoreset for which any field matches the query
	#returns a dataframe with urn, symbol and ensemblGeneID columns
	searchScoresets <- function(query) {
		.connect()
		tryCatch({
			scoresets <- dbGetQuery(.con,"SELECT urn, symbol, ensemblGeneID FROM scoresets;")
		},finally={
			.disconnect()
		})
		hits <- apply(apply(scoresets,2,function(col) grepl(query,col,ignore.case=TRUE)),1,any)
		return(scoresets[hits,])
	}

	#retrieves a table of variants whose hgvs string matches the query
	#returns a dataframe with accession, symbol, hgvs_pro, clinsig and ensemblGeneID columns
	searchVariants <- function(query) {
		.connect()
		tryCatch({
			res <- dbSendStatement(.con,
		 		"SELECT accession, symbol, hgvs_pro, clinsig, ensemblGeneID FROM 
		 		variants INNER JOIN scoresets 
		 		ON variants.scoreset = scoresets.urn 
		 		AND hgvs_pro LIKE $pattern;"
		 	)
		 	#dbBind should automatically escape dangerous characters 
		 	#and prevent sql injection
		 	dbBind(res, list(pattern=paste0("%",query,"%")))
		 	variants <- dbFetch(res)
		 	dbClearResult(res)
		},finally={
			.disconnect()
		})
	 	return(variants)
		# variants <- dbGetQuery(.con,"SELECT accession, hgvs_pro, clinsig FROM variants;")
		# hits <- apply(apply(variants,2,function(col) grepl(query,col,ignore.case=TRUE)),1,any)
		# return(variants[hits,])
	}

	#returns a list of URNs of the scoresets that are currently in the 'pending' state.
	getPending <- function() {
		.connect()
		tryCatch({
			scoresets <- dbGetQuery(.con,"SELECT urn FROM scoresets WHERE status='pending';")
		},finally={
			.disconnect()
		})
		return(scoresets[,1])
	}


	structure(list(
		isKnown=isKnown,
		getStatus=getStatus,
		getDate=getDate,
		setStatus=setStatus,
		getParameters=getParameters,
		setParameters=setParameters,
		getScoreset=getScoreset,
		newScoreset=newScoreset,
		newVariants=newVariants,
		getVariants=getVariants,
		getGoldStandard=getGoldStandard,
		hasVariant=hasVariant,
		getVariantDetail=getVariantDetail,
		calibrateScores=calibrateScores,
		searchScoresets=searchScoresets,
		searchVariants=searchVariants,
		getPending=getPending,
		isLocked=isLocked
		# close=close
	),class="persistence.connection")

}

