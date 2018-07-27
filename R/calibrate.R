
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
#' @return a \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{hgvsc} The HGVS variant descriptor at the coding sequence (DNA) level
#'   \item \code{hgvsp} The HGVS variant descriptor at the amino acid sequence (Protein) level
#'   \item \code{clinsig} The ClinVar clinical significance string (e.g "Likely pathogenic")
#' }
#' @export
#' 
fetchClinvar <- function(gene) {

	library(httr)
	library(RJSONIO)
	library(yogitools)

	set_config(config(ssl_verifypeer = 0L))

	#This is a two-step process: First we have to query Clinvar for a list of matching
	# DB entries. Then we can make a query for the details of the matched entries.

	searchBase <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
	summaryBase <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

	clinvarIds <- character()

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

	return(missense[,c("hgvsc","hgvsp","clinsig")])
}




#' Fetch missense variants in GnomAD for given gene
#' 
#' Makes a query to the ExAC webservice to fetch the GnomAD entries for the given gene
#' and filters them down to missense variants.
#' 
#' @param ensemblID The Ensembl gene identifier (e.g. \code{ENSG00000198668})
#' @return a \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{hgvsc} The HGVS variant descriptor at the coding sequence (DNA) level
#'   \item \code{hgvsp} The HGVS variant descriptor at the amino acid sequence (Protein) level
#'   \item \code{maf} The minor allele frequency.
#' }
#' @export
#' 
fetchGnomad <- function(ensemblID) {

	library(httr)
	library(RJSONIO)
	library(yogitools)

	exacURL <- "http://exac.hms.harvard.edu/rest/gene/variants_in_gene/"

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
	k <- sum(c(llrs,prior))
	#transform to probability
	exp(k)/(1+exp(k))
}


flipScores <- function(xs) sapply(xs,function(x) if (x > 1) 1/x else x)


#' Build a table of gold standard variants and their associated scores
#' 
#' @param scores a dataframe with columns: hgvs_nt, hgvs_pro, score
#' @param ensembls a vector of ensembl IDs for the target gene
#' @param symbols a vector of gene symbols of the target gene
#' @param drawPlot boolean value indicating whether a plot of variant distributions
#'     should be drawn. Defaults to \code{TRUE}.
#' @return a \code{data.frame} containing the original \code{scores} input
#'   with an additional column \code{llr}, representing the log likelihood ratio
#'   (i.e. log Bayes Factor)
#' @export
#' 
goldStandardScores <- function(scores,ensembls,symbols,
		drawPlot=TRUE,minMaf=0,flip=FALSE,homozygous=FALSE) {

	library(hash)
	library(yogitools)

	score.idx <- hash(scores$hgvs_pro,scores$score)
	
	gnomad <- do.call(rbind,lapply(ensembls,fetchGnomad))
	clinvar <- do.call(rbind,lapply(symbols,fetchClinvar))

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
		as.df(lapply(patho.vars,function(v) {
			sig <- paste(with(clinvar,clinsig[which(hgvsp==v)]),collapse="|")
			maf <- if (v %in% gnomad$hgvsp) {
				with(gnomad,maf[which(hgvsp==v)])
			} else NA
			hom <- if (v %in% gnomad$hgvsp) {
				sum(with(gnomad,hom[which(hgvsp==v)]))
			} else NA
			list(hgvsp=v,clinsig=sig,maf=maf,hom=hom,score=patho.scores[[v]])
		})),
		as.df(lapply(benign.vars,function(v) {
			maf <- if (v %in% gnomad$hgvsp) {
				with(gnomad,maf[which(hgvsp==v)])
			} else NA
			sig <- if (v %in% clinvar$hgvsp) {
				paste(with(clinvar,clinsig[which(hgvsp==v)]),collapse="|")
			} else "GnomAD"
			hom <- if (v %in% gnomad$hgvsp) {
				sum(with(gnomad,hom[which(hgvsp==v)]))
			} else NA
			list(hgvsp=v,clinsig=sig,maf=maf,hom=hom,score=benign.scores[[v]])
		}))
	)

	#draw a summary plot if desired
	if (drawPlot) {
		#Calculate distribution parameters
		patho.m <- mean(patho.scores)
		patho.sd <- sd(patho.scores)
		benign.m <- mean(benign.scores)
		benign.sd <- sd(benign.scores)

		#likelihood functions
		dens.patho <- function(x) dnorm(x,patho.m,patho.sd)
		dens.benign <- function(x) dnorm(x,benign.m,benign.sd)

		xlim <- range(c(patho.scores,benign.scores))
		plot(
			0,type="n",
			xlim=xlim,ylim=c(0,max(dens.patho(patho.m),dens.benign(benign.m))),
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
#' @return a \code{data.frame} containing the original \code{scores} input
#'   with an additional column \code{llr}, representing the log likelihood ratio
#'   (i.e. log Bayes Factor)
#' @export
#' 
map2bf <- function(scores,ensembls,symbols,drawPlot=TRUE,minMaf=0,flip=FALSE) {

	# library(hash)

	# score.idx <- hash(scores$hgvs_pro,scores$score)
	
	# gnomad <- do.call(rbind,lapply(ensembls,fetchGnomad))
	# clinvar <- do.call(rbind,lapply(symbols,fetchClinvar))

	# #Define pathogenic variants as those in Clinvar annotated as Pathogenic or Likely pathogenic
	# patho <- clinvar[grepl("athogenic",clinvar$clinsig),]
	# patho <- patho[!grepl("Conflict",patho$clinsig),]
	# patho.vars <- unique(patho$hgvsp)

	# #Define benign variants as those in Gnomad not labeled pathogenic in clinvar
	# benign.vars <- setdiff(unique(gnomad$hgvsp),patho.vars)

	# #Filter down to variants present in the map
	# patho.vars <- patho.vars[has.key(patho.vars,score.idx)]
	# benign.vars <- benign.vars[has.key(benign.vars,score.idx)]

	# #Lookup the corresponding scores
	# patho.scores <- values(score.idx,patho.vars)
	# benign.scores <- values(score.idx,benign.vars)

	summary.table <- goldStandardScores(scores,ensembls,symbols,FALSE)

	#TODO: Apply test to check if distributinos are significantly different

	patho.scores <- with(summary.table,score[clinsig != "GnomAD"])
	benign.scores <- with(summary.table,score[clinsig == "GnomAD"])

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
#' @return a \code{data.frame} containing the original \code{scores} input
#'   with an additional column \code{llr}, representing the log likelihood ratio
#'   (i.e. log Bayes Factor)
#' @export
#' 
map2bf.mavedb <- function(ssid,ensembls,symbols,drawPlot=TRUE,minMaf=0,flip=FALSE) {

	library(rapimave)
	mave <- new.rapimave()

	scores <- mave$getScores(ssid)

	map2bf(scores,ensembls,symbols,minMaf,flip)
}


#' Calibrate a variant map from a file to Bayes Factors
#' 
#' @param csvfile the scoreset file in csv format
#' @param ensembls a vector of ensembl IDs for the target gene
#' @param symbols a vector of gene symbols of the target gene
#' @param drawPlot boolean value indicating whether a plot of variant distributions
#'     should be drawn. Defaults to \code{TRUE}.
#' @return a \code{data.frame} containing the original \code{scores} input
#'   with an additional column \code{llr}, representing the log likelihood ratio
#'   (i.e. log Bayes Factor)
#' @export
#' 
map2bf.file <- function(csvfile,ensembls,symbols,drawPlot=TRUE,minMaf=0,flip=FALSE) {

	scores <- read.csv(csvfile,stringsAsFactors=FALSE)

	if (!all(c("hgvs_pro","score") %in% colnames(scores))) {
		stop("File must contain at least columns 'hgvs_pro' and 'score'")
	}

	map2bf(scores,ensembls,symbols,minMaf,flip)
}
