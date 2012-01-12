\name{FLUX_VAR_ANALYSIS}
\alias{Flux Variability Analysis}
\alias{FLUX_VAR_ANALYSIS}
\title{FLUX_VAR_ANALYSIS, a function to perform a flux variability analysis.}
\description{FVA can indicate the decrease in network robustness caused by non-lethal deletions.}
\usage{
    FLUX_VAR_ANALYSIS(fba_object,reactions,filename)
    }
\arguments{
	\item{fba_object}{
	Is a list containing the data required to perform flux
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
	\item{reactions}{
	subset of reaction numbers on which to perform FVA, if left empty, all
	reactions of the network will be subject to FVA
	}	
	\item{filename}{
	A string which will be the name of the file containing the 
	output of the FVA simulation}
}
\examples{
	#Simple flux variability Analysis
	data(Ecoli_core)
	FLUX_VAR_ANALYSIS(Ecoli_core, filename="Wt_FVA.xls")
}
