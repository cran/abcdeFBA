\name{BYPASS_REACTIONS_SUBSTRATE}
\alias{BYPASS_REACTIONS_SUBSTRATE}
\alias{bypass reactions for metabolites}
\alias{alternative pathways}
\title{BYPASS_REACTIONS_SUBSTRATE, to find the influx and efflux reactions available to a particular metabolite/substrate}
\description{this function computes the production and 
	consumption fluxes available to a particular 
	metabolite from the fba_object and returns a 
	list of with elements Production and Consumption}
\usage{
BYPASS_REACTIONS_SUBSTRATE(substrate_number,fba_object)
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
	\item{substrate_number}{this is the metabolite number of interest 
	which may be found by using the SEARCH_metabolite function }
}
\examples{
	#To find the Bypass Reactions for a Substrate
	data(Ecoli_core)	
	ATP_prod_consump<-BYPASS_REACTIONS_SUBSTRATE(17,Ecoli_core)
	
}
