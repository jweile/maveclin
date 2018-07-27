

#' New Family Object
#'
#' Creates a new family object, allowing for calculation of a Bayes factor. The Bayes Factor
#' Represents the likelyhood ratio of a phenotype being caused by an allele in question.
#'
#' A family object offers the following functions:
#' \itemize{
#'   \item \code{add.parent(pheno)} Adds a parent to the family, with boolean parameter \code{pheno}
#'      indicating whether the parent has the phenotype/disease.
#'   \item \code{add.sibling(pheno)} Adds a sibling to the family, with boolean parameter \code{pheno}
#'      indicating whether the sibling has the phenotype/disease.
#'   \item \code{add.child(pheno)} Adds a child to the family, with boolean parameter \code{pheno}
#'      indicating whether the child has the phenotype/disease.
#'   \item \code{log.bayes.factor()} Returns the overall log likelihood ratio (i.e. log Bayes-Factor) 
#'      for this family.
#' }
#'
#' @param maf Minor Allele Frequency of allele in question
#' @param etiologic.fraction The fraction of symptomatic patients that carry the genotype
#' @param penetrance The fraction of genotype carriers that are symptomatic
#' @param self The genotype of the patient him/herself. "heterozyogus" or "homozygous"
#' @param modus The mode of inheritance. "dominant" or "recessive"
#' @export
#' @return a new Family object
#' @examples
#' family <- new.family(
#'   maf=0.05,
#'   etiologic.fraction=0.8,
#'   penetrance=0.8,
#'   self="heterozygous",
#'   modus="dominant"
#' )
#' family$add.parent(pheno=TRUE)
#' family$add.parent(pheno=FALSE)
#' family$add.sibling(pheno=TRUE)
#' family$summary()
new.family <- function(maf,etiologic.fraction,penetrance,
	self=c("heterozygous","homozygous"),
	modus=c("dominant","recessive")) {

	#If patient is heterozygous for recessive disease, they should be unaffected!
	if (modus=="recessive" && self=="heterozygous") {
		stop("No implementation for heterozygous recessive!")
	}

	#set up object's fields
	.maf <- maf
	.self <- match.arg(self,choices=c("heterozygous","homozygous"))
	.modus <- match.arg(modus,choices=c("recessive","dominant"))
	.members <- list()

	# adds a parent to the family object
	# @param pheno A boolean value indicating whether parent is symptomatic
	add.parent <- function(pheno) {
		if (is.null(pheno) || is.na(pheno) || !is.logical(pheno)) {
			stop("Parameter phenotype must be boolean!")
		}
		if (sum(sapply(.members,`[[`,"relation")=="parent") >= 2) {
			stop("No more than two parents allowed!")
		}
		.members <<- c(.members,list(list(relation="parent",pheno=pheno)))
	}

	# adds a sibling to the family object
	# @param pheno A boolean value indicating whether sibling is symptomatic
	add.sibling <- function(pheno) {
		if (is.null(pheno) || is.na(pheno) || !is.logical(pheno)) {
			stop("Parameter phenotype must be boolean!")
		}
		.members <<- c(.members,list(list(relation="sibling",pheno=pheno)))
	}

	# adds a child to the family object
	# @param pheno A boolean value indicating whether child is symptomatic
	add.child <- function(pheno) {
		if (is.null(pheno) || is.na(pheno) || !is.logical(pheno)) {
			stop("Parameter phenotype must be boolean!")
		}
		.members <<- c(.members,list(list(relation="child",pheno=pheno)))
	}

	# print a summary for this family
	summary <- function() {
		lcaus <- likelihoods.causal()
		lncaus <- likelihoods.notcausal()
		for (i in 1:length(.members)) {
			member <- .members[[i]]
			geno <- genoProb(member)
			cat(member$relation,":\n")
			cat("\tPhenotype: ",member$pheno,"\n")
			cat("\tP(Genotype): ",geno,"\n")
			cat("\tL(causal): ",lcaus[[i]],"\n")
			cat("\tL(not causal): ",lncaus[[i]],"\n")
		}
		cat("================================
			Likelihood ratio: ",exp(log.bayes.factor()),"\n")
	}

	#calculate the log likelihood ratio, i.e. log(Bayes Factor)
	log.bayes.factor <- function() {
		loglikelyFor <- sum(log(likelihoods.causal()))
		loglikelyAgainst <- sum(log(likelihoods.notcausal()))
		return(loglikelyFor-loglikelyAgainst)
	}

	#calculate the individual family members' likelihoods under the causal hypothesis
	likelihoods.causal <- function() {
		sapply(.members,function(member) {
			g <- genoProb(member)
			if (member$pheno) {
				1 - (1 - g * penetrance) * etiologic.fraction
			} else {
				1 - (g * (1 - (g * (1 - penetrance))))
			}
		})
	}

	#calculate the individual family members' likelihoods under the non-causal hypothesis
	likelihoods.notcausal <- function() {
		sapply(.members,function(member) {
			if (member$pheno) {
				1-etiologic.fraction
			} else {
				1
			}
		})
	}

	#calculate the probability of the relevant genotype for a given family member
	genoProb <- function(member) {

		#alpha = frequency of heterozygousity
		a <- sqrt(.maf+0.25) - 0.5

		switch(member$relation,
			parent={
				switch(.modus,
					recessive={
						#assume homozygous
						(2*a^3+a^4)/(a^2+2*a^3+a^4)
					},
					dominant={
						switch(.self,
							heterozygous={
								(2+(a/(1-a^2)))/4
							},
							homozygous={
								1
							}
						)
					}
				)
			},
			sibling={
				switch(.modus,
					recessive={
						#assume homozygous
						(0.25*a^2+a^3+a^4)/(a^2+2*a^3+a^4)
					},
					dominant={
						switch(.self,
							heterozygous={
								(-a^4-0.5*a^3+1.25*a^2+0.5*a)/(-a^4-a^3+a^2+a)
							},
							homozygous={
								(a^4+2*a^3+0.75*a^2)/(a^4+2*a^3+a^2)
							}
						)
						
					}
				)
			},
			child={
				switch(.modus,
					recessive={
						#assume homozygous
						0.5*a+a^2
					},
					dominant={
						switch(.self,
							heterozygous={
								0.5*a^2-0.25*a+0.5
							},
							homozygous={
								1
							}
						)
					}
				)
			}
		)

	}

	structure(list(
		add.parent=add.parent,
		add.sibling=add.sibling,
		add.child=add.child,
		summary=summary,
		log.bayes.factor=log.bayes.factor
	),class="family")
}
