\name{Gene_del}
\alias{gene deletion}
\alias{genetic mutants }
\alias{Gene_del}
\title{Gene_del, a function to create Gene Deletion mutants}
\description{Gene_del interprets the boolean rules inherent to the model and returns a mutant model when supplied with the Wild-type model and a vector of gene names to be deleted. }

\usage{
Gene_del(query_genes,fba_object,return_reactions)
}

\arguments{

\item{query_genes}{ query_genes should be a vector the elements of which are the genes that are to be deleted from the model.} 

\item{fba_object}{A flux balance model containing the necessary elements for simulation} 

\item{return_reactions}{a boolean specifying if the reaction numbers deleted as a consequence of the gene deletion should be returned or if the mutant model itself be returned by the function} 

}

\examples{
#To create a mutant model .
{data(Ec_iAF1260_flux1)
	ho<-Gene_del(c("b3040"),Ec_iAF1260_flux1,return_reactions=TRUE)
}
}
