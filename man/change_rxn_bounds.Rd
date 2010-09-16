\name{CHANGE_RXN_BOUNDS}
\alias{Change Reaction Bounds}
\alias{CHANGE_RXN_BOUNDS}    
\title{CHANGE_RXN_BOUNDS, Change Reaction Bounds}
\description{this function helps to modify the bounds on the fba_object, 
		it returns an object of the same type as the model, inspired 
		by the COBRA-function of the same name}
\usage{
CHANGE_RXN_BOUNDS(reaction_number,fba_object,lower_bound,upper_bound)
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
	max is a Boolean specifying the type of optimization,"Maximization" by default".
	}
	\item{lower_bound}{the new value of the lower bound}
	\item{upper_bound}{the new value of the upper bound}
	\item{reaction_number}{the reaction number of the reaction the bounds of which have to be changed}
	}
\examples{
	#Changing Reaction Bounds, to simulate the deletion/inhibition of a reaction.
	data(Ecoli_core)	
	Ec_mutant<-CHANGE_RXN_BOUNDS(reaction_number=36,fba_object=Ecoli_core,
	lower_bound=0,upper_bound=0)	
	}
