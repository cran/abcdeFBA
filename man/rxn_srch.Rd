\name{SEARCH_reaction}
\alias{Search-reaction}
\alias{SEARCH_reaction}
\title{SEARCH_reaction, a function to search for reactions in a model using a simple text query}
\description{this function helps to search for reactions in a model using a simple text query, it returns a list of possible hits along with their reaction numbers, these reaction numbers are to be used with other perturbation/knockout/optimality/robustness analysis functions}
\usage{
SEARCH_reaction(react_name,fba_object)
}
\arguments{
	\item{fba_object}{is a list containing the data required to 
	perform flux balance analysis. The elements of the list are 
	mat which is the stoichiometric matrix, dir which gives the 
	direction of the equality constraints, obj specifies the 
	objective function for the simulation, bounds specifies the 
	lower and upper inequality constraints, rhs is the right hand
	side of the steady state expression, types refers to the numeric 
	nature of the variables which in case of FBA happens to be "Continuous",
	max is a Boolean specifying the type of optimization,"Maximization" by default".}
	\item{react_name}{is a character string containing a part or the 
	complete name of the reaction to be searched for in the fba object}
	}
\examples{
	#Search for reactions involving glucose
	data(Ecoli_core)	
	SEARCH_reaction("glucose",Ecoli_core)
	}
