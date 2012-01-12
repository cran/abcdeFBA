\name{View_objective}
\alias{Objective function}
\alias{View_objective}    
\title{
	View_objective, a function that shows the components of the  
	objective function.}
\description{
	This function takes the argument of type fba_object and prints out the 
	components of the objective function, i.e substrates and products involved
	and their corresponding pseudo-stoichiometry }
\usage{
    View_objective(fba_object)
    }
\arguments{
      \item{fba_object}{
Is a list containing the data required to perform flux balance analysis.
The elements of the list are mat which is the stoichiometric matrix, dir which 
gives the direction of the equality constraints, obj specifies the objective 
function for the simulation, bounds specifies the lower and upper inequality 
constraints, rhs is the right hand side of the steady state expression, 
types refers to the numeric nature of the variables which in case of FBA happens
to be "Continuous", max is a Boolean specifying the type of optimization,
"Maximization" by default, all_genes is all the genes present in the model, 
gpr contains boolean expressions of gene essentiality for the corresponding 
reactions in the model, metabolite_name contains list of all the metabolites, 
reaction_list contains all the reactions present in the model, compartment 
is a numeric identifier for each reaction the key for which is in comp_name}
}

\examples{
	#Viewing the objective function of the E.coli core metabolism
	data(Ecoli_core)
	View_objective(Ecoli_core)
}

