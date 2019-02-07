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


#' Find cache file location by name
#' 
#' Finds the location for a cache file. The file does not necessary need to exist yet,
#' as this function is meant to be used determine to a location for both storage and retrieval.
#' 
#' Depending on the execution context, the storage location may differ. The cache location can 
#' be controlled with the environment variable \code{$MAVECLIN_CACHE}. This will be made use of within
#' the mavevis server docker container. If the variable is not set, a directory ".mavecache/"
#' will be created in the user's home directory to be used as the storage location.
#' 
#' @param name the name of the file, e.g. "P12456_alignment.fasta"
#' @return the full path to the file
#' @export
#' @examples
#' file <- getCacheFile("P12345_alignment.fasta")
#' 
getCacheFile <- function(name) {
	cache.loc <- Sys.getenv("MAVECLIN_CACHE",unset=NA)
	if (is.na(cache.loc)) {
		cache.loc <- paste0(Sys.getenv("HOME"),"/.mavecache/")
	}
	if (!file.exists(cache.loc)) {
		dir.create(cache.loc,showWarnings=FALSE,recursive=TRUE)
	}
	paste0(cache.loc,name)
}

#' Decode HTML strings
#' 
#' Decode HTML strings by converting special character escape sequences (e.g. \code{&gt;})
#' to their corresponding UTF-8 characters.
#' 
#' @param str input string
#' @return the decoded string
#' 
htmlDecode <- function(str) {
	str <- gsub("&gt;",">",str)
	#implement more as needed
	return(str)
}


#' Fetch missense variants in ClinVar for given gene
#' 
#' Makes a query to the NCBI webservice to fetch the ClinVar entries for the given gene
#' and filters them down to missense variants
#' 
#' @param gene gene name as in ClinVar (e.g. 'CALM1')
#' @param stagger logical determining whether HTTP requests should be staggered by
#'   a third of a second to avoid rejection by the server. Defaults to TRUE.
#' @param overrideCache logical determining whether local cache should be overridden
#' @param logger a yogilogger object to write log messages to. Defaults to NULL
#' @return a \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{hgvsc} The HGVS variant descriptor at the coding sequence (DNA) level
#'   \item \code{hgvsp} The HGVS variant descriptor at the amino acid sequence (Protein) level
#'   \item \code{clinsig} The ClinVar clinical significance string (e.g "Likely pathogenic")
#' }
#' @export
#' 
fetchClinvar <- function(gene,stagger=TRUE,overrideCache=FALSE,logger=NULL) {

	if (!is.null(logger)) {
		stopifnot(inherits(logger,"yogilogger"))
		library("yogilog")
	}

	cacheFile <- getCacheFile(paste0("clinvar_",gene,".csv"))

	if (!file.exists(cacheFile) || overrideCache) {

		library(httr)
		library(RJSONIO)
		library(yogitools)

		set_config(config(ssl_verifypeer = 0L))

		#This is a two-step process: First we have to query Clinvar for a list of matching
		# DB entries. Then we can make a query for the details of the matched entries.

		searchBase <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
		summaryBase <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

		clinvarIds <- character()

		if (stagger) {
			#stick to slightly less than 3 queries per second
			Sys.sleep(0.34)
		}

		if (!is.null(logger)) {
			logger$info("Querying Clinvar for ",gene)
		}
		#make HTTP GET request to find matching entries
		htr <- GET(searchBase,query=list(
			db="clinvar",
			term=paste0(gene,"[gene] AND single_gene[prop]"),
			retmax=1000,
			retmode="json"
		))
		if (http_status(htr)$category == "Success") {
			#parse returned content as JSON
			returnData <- fromJSON(content(htr,as="text",encoding="UTF-8"))

			if (is.null(returnData) || length(returnData) < 1) {
				stop("Server returned no data.")
			}

			if (as.numeric(returnData$esearchresult$count) > 1000) {
				warning("More than 1000 results, excess skipped.")
			}

			#extract the list of matches
			clinvarIds <- returnData$esearchresult$idlist

		} else {
			stop("server message: ",http_status(htr)$message)
		}

		if (length(clinvarIds) < 1) {
			stop("No results!")
		}

		if (stagger) {
			#stick to slightly less than 3 queries per second
			Sys.sleep(0.34)
		}

		if (!is.null(logger)) {
			logger$info("Fetching detailed data from Clinvar")
		}
		#make HTTP get request for details of the matches
		htr <- GET(summaryBase,query=list(
			db="clinvar",
			id=paste(clinvarIds,collapse=","),
			retmode="json"
		))
		if (http_status(htr)$category == "Success") {
			#parse JSON
			returnData <- fromJSON(content(htr,as="text",encoding="UTF-8"))

			if (is.null(returnData) || length(returnData) < 1) {
				stop("Server returned no data.")
			}

			#the first result is simply a repetition of the entries, so we discard it.
			#the rest are the actual details on the entries.
			results <- as.df(lapply(returnData$result[-1], function(vset) {
				#extract the variant description
				varStr <- htmlDecode(vset$variation_set[[1]]$cdna_change)
				#and the clinical significance statement
				clinsig <- vset$clinical_significance[["description"]]
				return(list(var=varStr,clinsig=clinsig))
			}))

		} else {
			stop("server message: ",http_status(htr)$message)
		}

		#filter the results down to only missense variants
		missense <- results[grepl("\\(p\\.\\w{3}\\d+\\w{3}\\)",results$var),]
		missense <- missense[!grepl("\\(p\\.\\w{3}\\d+Ter\\)",missense$var),]
		missense <- missense[!grepl("\\(p\\.\\w{3}\\d+\\w{3}fs\\)",missense$var),]
		#extract the HGVS descriptors at the coding and protein levels.
		missense$hgvsc <- extract.groups(missense$var,"(c\\.\\d+[ACGT]>[ACGT])")[,1]
		missense$hgvsp <- extract.groups(missense$var,"(p\\.\\w{3}\\d+\\w{3})")[,1]

		output <- missense[,c("hgvsc","hgvsp","clinsig")]

		if (!is.null(logger)) {
			logger$info("Caching Clinvar data for ",gene)
		}
		write.table(output,cacheFile,sep=",")

	} else {
		if (!is.null(logger)) {
			logger$info("Retrieving cached data for ",gene)
		}
		output <- read.csv(cacheFile)
	}

	return(output)
}




#' Fetch missense variants in GnomAD for given gene
#' 
#' Makes a query to the ExAC webservice to fetch the GnomAD entries for the given gene
#' and filters them down to missense variants.
#' 
#' @param ensemblID The Ensembl gene identifier (e.g. \code{ENSG00000198668})
#' @param overrideCache logical determining whether local cache should be overridden
#' @param logger a yogilogger object to write log messages to. Defaults to NULL
#' @return a \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{hgvsc} The HGVS variant descriptor at the coding sequence (DNA) level
#'   \item \code{hgvsp} The HGVS variant descriptor at the amino acid sequence (Protein) level
#'   \item \code{maf} The minor allele frequency.
#' }
#' @export
#' 
fetchGnomad <- function(ensemblID,overrideCache=FALSE,logger=NULL) {

	if (!is.null(logger)) {
		stopifnot(inherits(logger,"yogilogger"))
		library("yogilog")
	}

	cacheFile <- getCacheFile(paste0("gnomad_",ensemblID,".csv"))

	if (!file.exists(cacheFile) || overrideCache) {

		library(httr)
		library(RJSONIO)
		library(yogitools)

		exacURL <- "http://exac.hms.harvard.edu/rest/gene/variants_in_gene/"

		if (!is.null(logger)) {
			logger$info("Querying GnomAD for ",ensemblID)
		}
		#make HTTP GET request
		htr <- GET(paste0(exacURL,ensemblID))
		if (http_status(htr)$category == "Success") {
			#parse JSON
			returnData <- fromJSON(content(htr,as="text",encoding="UTF-8"))

			#filter down to only missense variants
			missense.idx <- which(sapply(returnData,`[[`,"category")=="missense_variant")

			missense.gnomad <- as.df(lapply(returnData[missense.idx],function(entry) {
				#extract hgvs and allele frequency data
				with(entry,list(
					hgvc = HGVSc,
					hgvsp = HGVSp,
					maf = if (is.null(allele_freq)) NA else allele_freq,
					hom = if (is.null(hom_count)) 0 else hom_count
				))
			}))

			#the annotation of missense is in gnomad is not correct, so we have to filter again!
			missense.gnomad <- missense.gnomad[grepl("^p\\.\\w{3}\\d+\\w{3}$",missense.gnomad$hgvsp),]
		} else {
			stop("server message: ",http_status(htr)$message)
		}

		#Write results to cache
		if (!is.null(logger)) {
			logger$info("Caching GnomAD data for ",ensemblID)
		}
		write.table(missense.gnomad,cacheFile,sep=",")


	} else {
		if (!is.null(logger)) {
			logger$info("Retrieving cached data for ",ensemblID)
		}
		missense.gnomad <- read.csv(cacheFile)
	}

	return(missense.gnomad)

}

#' Calculate posterior probability via Bayesian Inference
#' 
#' @param llrs Vector of log likelihood ratios (i.e. log(Bayes-Factor)) 
#' @param prior The prior probability
#' @return the posterior probability
calc.posterior <- function(llrs, prior) {
	#transform prior to log odds
	lpo <- log(prior/(1-prior))
	#calculate log odds posterior
	k <- sum(c(llrs,lpo))
	#transform to probability
	exp(k)/(1+exp(k))
}


#' Flip transformation
#' 
#' Applies the "flip" transformation to a numerical vector. I.e. all values 
#' greater than one are transformed by 1/x
#' 
#' @param xs numerical vector
#' @return the flipped scores as a numerical vector
flipScores <- function(xs) sapply(xs,function(x) if (x > 1) 1/x else x)


#' Build a table of gold standard variants and their associated scores
#' 
#' @param scores a dataframe with columns: hgvs_nt, hgvs_pro, score
#' @param ensembls a vector of ensembl IDs for the target gene
#' @param symbols a vector of gene symbols of the target gene
#' @param drawPlot boolean value indicating whether a plot of variant distributions
#'     should be drawn. Defaults to \code{TRUE}.
#' @param minMaf the minimum Minor Allele Frequency required for variants to be considered benign
#' @param flip logical; whether to apply the flip transformation to the scores
#' @param homozygous logical; filter benign variants to only those occuring homozygously.
#' @param overrideCache logical; whether to override local cache.
#' @param logger a yogilogger object to which to write log messages.
#' @return a \code{data.frame} containing the original \code{scores} input
#'   with an additional column \code{llr}, representing the log likelihood ratio
#'   (i.e. log Bayes Factor)
#' @export
#' 
goldStandardScores <- function(scores,ensembls,symbols,
		drawPlot=TRUE,minMaf=0,flip=FALSE,homozygous=FALSE,
		overrideCache=FALSE,logger=NULL) {

	if (!is.null(logger)) {
		stopifnot(inherits(logger,"yogilogger"))
		library("yogilog")
	}

	library(hash)
	library(yogitools)

	if (flip) {
		scores$score <- flipScores(scores$score)
	}

	score.idx <- hash(scores$hgvs_pro,scores$score)
	
	gnomad <- do.call(rbind,lapply(ensembls,fetchGnomad,overrideCache=overrideCache,logger=logger))
	clinvar <- do.call(rbind,lapply(symbols,fetchClinvar,overrideCache=overrideCache,logger=logger))

	if(!is.null(logger)) {
		logger$info("Building Gold Standard table")
	}

	gnomad <- gnomad[which(gnomad$maf > minMaf),]

	#Define pathogenic variants as those in Clinvar annotated as Pathogenic or Likely pathogenic
	patho <- clinvar[grepl("athogenic",clinvar$clinsig),]
	patho <- patho[!grepl("Conflict",patho$clinsig),]
	patho.vars <- unique(patho$hgvsp)

	#also extract benign and anything but benign lists for downstream filtering
	clin.benign <- clinvar[grepl("enign",clinvar$clinsig),]
	clin.nbenign <- clinvar[!grepl("enign",clinvar$clinsig),]

	#Define benign variants as those in Gnomad not labeled non-benign in clinvar
	benign.vars <- if (homozygous) {
		unique(c(setdiff(gnomad$hgvsp[which(gnomad$hom > 0)],clin.nbenign$hgvsp),clin.benign$hgvsp))
	} else {
		unique(c(setdiff(gnomad$hgvsp,clin.nbenign$hgvsp),clin.benign$hgvsp))
	}

	#Filter down to variants present in the map
	patho.vars <- patho.vars[has.key(patho.vars,score.idx)]
	benign.vars <- benign.vars[has.key(benign.vars,score.idx)]

	#Lookup the corresponding scores
	patho.scores <- values(score.idx,patho.vars)
	benign.scores <- values(score.idx,benign.vars)

	#Make a summary table of the gold standard variants
	summary.table <- rbind(
		if (length(patho.vars) == 0) NULL else as.df(lapply(patho.vars,function(v) {
			sig <- paste(unique(with(clinvar,clinsig[which(hgvsp==v)])),collapse="|")
			maf <- if (v %in% gnomad$hgvsp) {
				with(gnomad,maf[which(hgvsp==v)])
			} else NA
			hom <- if (v %in% gnomad$hgvsp) {
				sum(with(gnomad,hom[which(hgvsp==v)]))
			} else NA
			list(hgvsp=v,clinsig=sig,maf=max(maf),hom=hom,score=patho.scores[[v]],set="+")
		})),
		if (length(benign.vars) == 0) NULL else as.df(lapply(benign.vars,function(v) {
			maf <- if (v %in% gnomad$hgvsp) {
				with(gnomad,maf[which(hgvsp==v)])
			} else NA
			sig <- if (v %in% clinvar$hgvsp) {
				paste(unique(with(clinvar,clinsig[which(hgvsp==v)])),collapse="|")
			} else "GnomAD"
			hom <- if (v %in% gnomad$hgvsp) {
				sum(with(gnomad,hom[which(hgvsp==v)]))
			} else NA
			list(hgvsp=v,clinsig=sig,maf=max(maf),hom=hom,score=benign.scores[[v]],set="-")
		}))
	)

	if (length(patho.scores) == 0 || length(benign.scores) == 0) {
		warning("Insufficient reference set!")
		return(summary.table)
	}

	#draw a summary plot if desired
	if (drawPlot) {

		if(!is.null(logger)) {
			logger$info("Plotting score densities")
		}

		#Calculate distribution parameters and likelihood functions
		patho.m <- mean(patho.scores)
		patho.sd <- sd(patho.scores)
		dens.patho <- function(x) dnorm(x,patho.m,patho.sd)

		benign.m <- mean(benign.scores)
		benign.sd <- sd(benign.scores)		
		dens.benign <- function(x) dnorm(x,benign.m,benign.sd)

		xlim <- range(c(patho.scores,benign.scores))
		ylim <- c(0,max(dens.patho(patho.m),dens.benign(benign.m)))
		plot(
			0,type="n",
			xlim=xlim,ylim=ylim,
			xlab="score",ylab="density"
		)
		abline(v=benign.scores,col="darkolivegreen3")
		abline(v=patho.scores,col="firebrick3",lwd=2)
		plot(dens.patho,col="firebrick3",lwd=2,lty="dashed",xlim=xlim,add=TRUE)
		plot(dens.benign,col="darkolivegreen3",lwd=2,lty="dashed",xlim=xlim,add=TRUE)
	}

	return(summary.table)
}

#' Calibrate a variant map to Bayes Factors
#' 
#' @param scores a dataframe with columns: hgvs_nt, hgvs_pro, score
#' @param ensembls a vector of ensembl IDs for the target gene
#' @param symbols a vector of gene symbols of the target gene
#' @param drawPlot boolean value indicating whether a plot of variant distributions
#'     should be drawn. Defaults to \code{TRUE}.
#' @param minMaf the minimum Minor Allele Frequency required for variants to be considered benign
#' @param flip logical; whether to apply the flip transformation to the scores
#' @param homozygous logical; filter benign variants to only those occuring homozygously.
#' @param overrideCache logical; whether to override local cache.
#' @param logger a yogilogger object to which to write log messages.
#' @return a \code{data.frame} containing the original \code{scores} input
#'   with an additional column \code{llr}, representing the log likelihood ratio
#'   (i.e. log Bayes Factor)
#' @export
#' 
map2bf <- function(scores,ensembls,symbols,drawPlot=TRUE,minMaf=0,flip=FALSE,
			homozygous=FALSE,overrideCache=FALSE,logger=NULL) {

	if (!is.null(logger)) {
		stopifnot(inherits(logger,"yogilogger"))
		library("yogilog")
	}

	#We use drawPlot=FALSE here, because we draw a better plot below if desired
	summary.table <- goldStandardScores(scores,ensembls,symbols,
		drawPlot=FALSE,minMaf=minMaf,flip=flip,homozygous=homozygous,
		overrideCache=overrideCache,logger=logger
	)

	if (length(unique(summary.table$set)) != 2) {
		stop("Insufficient reference set!")
	}

	if (!is.null(logger)) {
		logger$info("Calculating Bayes Factors")
	}

	if (flip) {
		scores$flippedScore <- flipScores(scores$score)
	}

	rownames(summary.table) <- summary.table$hgvsp
	scores$clinsig <- summary.table[scores$hgvs_pro,"clinsig"]
	scores$maf <- summary.table[scores$hgvs_pro,"maf"]
	scores$hom <- summary.table[scores$hgvs_pro,"hom"]
	scores$set <- summary.table[scores$hgvs_pro,"set"]

	#TODO: Apply test to check if distributions are significantly different

	patho.scores <- with(summary.table,score[which(set == "+")])
	benign.scores <- with(summary.table,score[which(set == "-")])

	#Calculate distribution parameters
	patho.m <- mean(patho.scores)
	patho.sd <- sd(patho.scores)
	benign.m <- mean(benign.scores)
	benign.sd <- sd(benign.scores)

	#likelihood functions
	likely.patho <- function(x) pnorm(x,patho.m,patho.sd,lower.tail=FALSE)
	likely.benign <- function(x) pnorm(x,benign.m,benign.sd,lower.tail=TRUE)
	#log likelihood ratio function
	llr <- function(x) log(likely.patho(x)) - log(likely.benign(x))

	#draw a summary plot if desired
	if (drawPlot) {
		xlim <- range(c(patho.scores,benign.scores))
		plot(
			0,type="n",
			xlim=xlim,ylim=c(0,1),
			xlab="score",ylab="probability"
		)
		abline(v=patho.scores,col="firebrick3")
		abline(v=benign.scores,col="darkolivegreen3")
		plot(likely.patho,col="firebrick3",lwd=2,xlim=xlim,add=TRUE)
		plot(likely.benign,col="darkolivegreen3",lwd=2,xlim=xlim,add=TRUE)
	}

	#calculate log likelihood ratio (log Bayes Factor) for each variant in the table
	scores$llr <- sapply(scores$score,llr)

	#build confidence interval preferably on stderr, but if not available, use stdev
	error <- with(scores, if (all(is.na(se))) sd else se)

	#calculate 90% confidence interval for LLR
	scores$llrCIleft=sapply(with(scores,qnorm(p=0.95,mean=score,sd=error)),llr)
	scores$llrCIright=sapply(with(scores,qnorm(p=0.05,mean=score,sd=error)),llr)
	# scores$llrCI <- apply(scores[,c("llrCIleft","llrCIright")],1,function(vs)sprintf("[%.02f;%.02f]",vs[[1]],vs[[2]]))

	#return the updated score table
	return(scores)

}


#' Calibrate a variant map from MaveDB to Bayes Factors
#' 
#' @param ssid the scoreset id (URN) from MaveDB
#' @param ensembls a vector of ensembl IDs for the target gene
#' @param symbols a vector of gene symbols of the target gene
#' @param drawPlot boolean value indicating whether a plot of variant distributions
#'     should be drawn. Defaults to \code{TRUE}.
#' @param minMaf the minimum Minor Allele Frequency required for variants to be considered benign
#' @param flip logical; whether to apply the flip transformation to the scores
#' @param homozygous logical; filter benign variants to only those occuring homozygously.
#' @param overrideCache logical; whether to override local cache.
#' @param logger a yogilogger object to which to write log messages.
#' @return a \code{data.frame} containing the original \code{scores} input
#'   with an additional column \code{llr}, representing the log likelihood ratio
#'   (i.e. log Bayes Factor)
#' @export
#' @examples
#' 
#' maveUrn <- "urn:mavedb:00000001-c-2"
#' ensembls <- c("ENSG00000198668","ENSG00000143933","ENSG00000160014")
#' symbols <- c("CALM1","CALM2","CALM3")
#' calmBFs <- map2bf.mavedb(maveUrn,ensembls,symbols)
#' 
map2bf.mavedb <- function(ssid,ensembls,symbols,
			drawPlot=TRUE,minMaf=0,flip=FALSE,homozygous=FALSE,
			overrideCache=FALSE,logger=NULL) {

	if (!is.null(logger)) {
		stopifnot(inherits(logger,"yogilogger"))
		library("yogilog")
	}

	nonstandardCall <- minMaf != 0 || flip || homozygous

	#FIXME: Should cache file name reflect gene names?
	cacheFile <- getCacheFile(paste0("bayesFactors_",ssid,".csv"))

	if (!file.exists(cacheFile) || overrideCache || nonstandardCall) {

		scoreCacheFile <- getCacheFile(paste0(ssid,".csv"))

		if (!file.exists(scoreCacheFile) || overrideCache) {

			library(rapimave)
			mave <- new.rapimave()
			if (!is.null(logger)) {
				logger$info("Querying MaveDB for ",ssid)
			}
			scores <- mave$getScores(ssid)
			write.table(scores,scoreCacheFile,sep=",")

		} else {

			if (!is.null(logger)) {
				logger$info("Retrieving cached scores for ",ssid)
			}
			scores <- read.csv(scoreCacheFile,stringsAsFactors=FALSE)
		}

		out <- map2bf(scores,ensembls,symbols,minMaf=minMaf,flip=flip,
			homozygous=homozygous,overrideCache=overrideCache,logger=logger
		)

		if (!is.null(logger)) {
			logger$info("Caching Bayes Factor results")
		}
		if (!nonstandardCall) {
			write.table(out,cacheFile,sep=",")
		}

	} else {

		if (!is.null(logger)) {
			logger$info("Retrieving cached Bayes Factors for ",ssid)
		}
		out <- read.csv(cacheFile,stringsAsFactors=FALSE)
	}

	return(out)
}


#' Calibrate a variant map from a file to Bayes Factors
#' 
#' @param csvfile the scoreset file in csv format
#' @param ensembls a vector of ensembl IDs for the target gene
#' @param symbols a vector of gene symbols of the target gene
#' @param drawPlot boolean value indicating whether a plot of variant distributions
#'     should be drawn. Defaults to \code{TRUE}.
#' @param minMaf the minimum Minor Allele Frequency required for variants to be considered benign
#' @param flip logical; whether to apply the flip transformation to the scores
#' @param homozygous logical; filter benign variants to only those occuring homozygously.
#' @param overrideCache logical; whether to override local cache.
#' @param logger a yogilogger object to which to write log messages.
#' @return a \code{data.frame} containing the original \code{scores} input
#'   with an additional column \code{llr}, representing the log likelihood ratio
#'   (i.e. log Bayes Factor)
#' @export
map2bf.file <- function(csvfile,ensembls,symbols,drawPlot=TRUE,minMaf=0,flip=FALSE,
		homozygous=FALSE,overrideCache=FALSE,logger=NULL) {

	scores <- read.csv(csvfile,stringsAsFactors=FALSE)

	if (!all(c("hgvs_pro","score") %in% colnames(scores))) {
		stop("File must contain at least columns 'hgvs_pro' and 'score'")
	}

	map2bf(scores,ensembls,symbols,minMaf=minMaf,flip=flip,
		homozygous=homozygous,overrideCache=overrideCache,logger=logger
	)
}
