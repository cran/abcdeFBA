\name{Gene_del}
\alias{gene deletion}
\alias{genetic mutants }
\alias{Gene_del}
\title{Gene_del, a function to create Gene Deletion mutants}
\description{Gene-protein are not one-to-one relationships due to multi-meric enzymes 
and iso-enzymes, Gene_del interprets the boolean rules inherent to the model and returns 
a mutant model when supplied with the Wild-type model and a vector of gene names to be deleted. }
\usage{
Gene_del(fba_object,gene_set)
}
\arguments{
\item{fba_object}{Is a list containing the data required to perform flux
	balance analysis. The elements of the list are mat which is the 
	stoichiometric matrix, dir which gives the direction of the equality 
	constraints, obj specifies the objective function for the simulation, 
	bounds specifies the lower and upper inequality constraints, rhs is the 
	right hand side of the steady state expression, types refers to the 
	numeric nature of the variables which in case of FBA happens to be 
	"Continuous", max is a Boolean specifying the type of optimization,
	"Maximization" by default, all_genes is all the genes present in the 
	model, gpr contains boolean expressions of gene essentiality for the 
	corresponding reactions in the model, metabolite_name contains list of 
	all the metabolites, reaction_list contains all the reactions present in
	the model, compartment is a numeric identifier for each reaction the key
	for which is in comp_name.}
\item{gene_set}{a character vector containing the gene names, frequently as bnums eg."b0001"}
}
\examples{
	#To create a mutant model .
	#data(Ecoli_core)
	#Ec_mutant<-Gene_del(Ecoli_core,c("b0001","b0002"))
	#Ec_mutant will have reactions corresponding to the genes specified deleted. 
	#FBA_solve(Ec_mutant)
}
