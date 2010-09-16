\name{SEARCH_metabolite}
\alias{SEARCH_metabolite}
\alias{Search-metabolite}
\title{SEARCH_metabolite, a function to search for metabolites in a model using a simple text query}
\description{this function helps to search for metabolites in a model 
	using a simple text query, it returns a list of possible hits 
	along with their metabolite numbers and locations in the 
	compartments of the model}
    \usage{
    SEARCH_metabolite(metabolite_name,fba_object)
    }
\arguments{
      \item{fba_object}{is a list containing the data required to 
	perform flux balance analysis. The elements of the list 
	are mat which is the stoichiometric matrix, dir which gives 
	the direction of the equality constraints, obj specifies the 
	objective function for the simulation, bounds specifies the 
	lower and upper inequality constraints, rhs is the right hand 
	side of the steady state expression, types refers to the numeric 
	nature of the variables which in case of FBA happens to be "Continuous", 
	max is a Boolean specifying the type of optimization,
	"Maximization" by default".}
      \item{metabolite_name}{is a character string containing a part 
	or the complete name of the metabolite to be searched for in
	 the fba object}
	}
\examples{
	#To search for the metabolite in the model
	data(Ecoli_core)	
	SEARCH_metabolite("ATP",Ecoli_core)
}