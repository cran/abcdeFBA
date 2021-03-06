\name{CHANGE_OBJ_FUNCTION}
\alias{change objective}
\alias{objective function}
\alias{CHANGE_OBJ_FUNCTION}
\title{CHANGE_OBJ_FUNCTION, a function to change the objective for optimization}
\description{a function akin to the COBRA function to change the objective 
function for FBA }

\usage{
CHANGE_OBJ_FUNCTION(obj_reaction,fba_object,new_wt,old_wt)
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

\item{obj_reaction}{a reaction number which is to be made the new 
objective function; retrieved using the SEARCH reaction function}

\item{new_wt}{the weight of the new objective, defaults to 1 but 
can be any number from 0~1}

\item{old_wt}{the weight of the old objective, defaults to 0 but 
can be any number from 0~(1-new_obj_weight) or any
other co-efficient if you wish for a customized 
objective function.}

}

\examples{
#To change the objective function of the model.
data(Ecoli_core)
ec_new_obj<-CHANGE_OBJ_FUNCTION(11,Ecoli_core,0.5,0.5)
#ec_new_obj will be identical to the Ecoli_core model except that 
#the objective function would change
FBA_solve(ec_new_obj)
}
