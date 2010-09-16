\name{FBA_solve}
\alias{solve FBA problem}
\alias{FBA_solve}
\title{FBA_solve, a function to solve CBM problems}
\description{
	This function sugar-coats Rglpk_solve_LP 
	which is a function provided by the Rglpk package, 
	FBA_solve solves FBA problems using Rglpk_solve_LP 
	and gives the solver output in a CBM context including 
	solver error messages and graceful degradation.}
    \usage{
    FBA_solve(fba_object, precision,verbosity,maximize)
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
	\item{precision}{A number indicating the precision of the flux solution 
	after the decimal point, defaults to 6}
	\item{verbosity}{is a Boolean indicating if the verbose output of the LP solver 
	should be displayed during simulation, defaults to FALSE}
	\item{maximize}{is a Boolean that can over-ride the default mode of optimization 
	(maximization) and minimize if FALSE}
}

\examples{# Flux Balance Analysis performed on a core-metabolism model of E.coli
	data(Ecoli_core)
	FBA_solve(fba_object=Ecoli_core,precision=6,verbosity=FALSE,maximize=TRUE)
	}
