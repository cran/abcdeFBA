\name{Exhaustive_double_deletion}
\alias{double knockout} 
\alias{parallel multi-core simulation}
\alias{double synthetic lethality}
\alias{Exhaustive_double_deletion}    
\title{Exhaustive_double_deletion, a function for computing synthetic double knockouts.}
\description{
	Exhaustive_double_deletion enables parallel computing of double knockouts by 
	splitting and running the simulation as different instances on the same multi-core 
	machine. The results of the simulation are stored into files named result1, 
	result2 etc.which may be reassembled using the BFDProcessor function
}
\usage{
    Exhaustive_double_deletion(fba_object,thread_no,core_number)
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
	\item{thread_no}{This specifies the thread number, it is useful to 
	think of the thread number as a chunk of the complete combination of 
	all pairs of reactions that can be formed. This depends on the number of 
	the cores you are employing which must obviously be a fixed number. 
	For instance if you have a 40 core computer you can divide the double-knockout 
	simulation into 40 chunks. In that case the thread number would refer to 
	the chunks of simulation as 0-39, with thread 0 performing the first 1/40th 
	part of the simulation}
	\item{core_number}{core_number specifies the number of cores that are available 
	to you for performing double knockout simulations. It should be a fixed number 
	for a particular instance of simulations on one model}
}
\examples{
	#Performing double knockouts in a 2 core computer 
	data(Ecoli_core)	
	Exhaustive_double_deletion(Ecoli_core,0,2)
	#New instance of R
	data(Ecoli_core)
	Exhaustive_double_deletion(Ecoli_core,1,2)
}
