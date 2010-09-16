\name{CHANGE_OBJ_FUNCTION}
\alias{change objective}
\alias{objective function}
\alias{CHANGE_OBJ_FUNCTION}
\title{CHANGE_OBJ_FUNCTION, a function to change the objective for optimization}
\description{a function akin to the COBRA function to change the objective function for FBA }
\usage{
CHANGE_OBJ_FUNCTION(objective_reaction,fba_object,new_obj_weight,old_obj_weight)
}
\arguments{
\item{fba_object}{is a list containing the data required to perform 
		flux balance analysis. The elements of the list are 
		mat which is the stoichiometric matrix, dir which 
		gives the direction of the equality constraints, 
		obj specifies the objective function for the simulation, 
		bounds specifies the lower and upper inequality constraints, 
		rhs is the right hand side of the steady state expression, 
		types refers to the numeric nature of the variables which 
		in case of FBA happens to be "Continuous", max is a Boolean 
		specifying the type of optimization,"Maximization" by default".}
\item{objective_reaction}{a reaction number which is to be made the new objective
		function; retrieved using the SEARCH reaction function}
\item{new_obj_weight}{the weight of the new objective, defaults to 1 but 
			can be any number from 0~1}
\item{old_obj_weight}{the weight of the old objective, defaults to 0 but 
			can be any number from 0~(1-new_obj_weight)}
}
\examples{
	#To change the objective function of the model.
	data(Ecoli_core)
	ec_new_obj<-CHANGE_OBJ_FUNCTION(11,Ecoli_core,0.5,0.5)
	#ec_new_obj will be identical to the Ecoli_core model except that 
	#the objective function would change
	FBA_solve(ec_new_obj)
}
