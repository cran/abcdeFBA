\name{FVA_robustness}
\alias{Robustness analysis using FVA}
\alias{FVA_robustness}    
\title{FVA_robustness, a function to determine the change in robustness 
of the network caused by a mutation.}

\description{
FVA_robustness builds upon FLUX_VAR_ANALYSIS to 
give the user a convenient function to examine 
the change in network robustness caused by a user 
supplied mutation. Graphical results for absolute 
flux span comparisons may be expected in the working directory}
\usage{
FVA_robustness(fba_object, mutation)
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

\item{mutation}{
A reaction number obtained using the SEARCH_reaction,or in the non-lethal
reaction "list" result of an Exhaustive_single_deletion run. 
This reaction will be deleted during the FVA analysis}
}

\examples{
#Mutant network robustness, removing Formate Exchange which has no effect 
#on fluxes 
data(Ecoli_core)	
#FVA_robustness(Ecoli_core,25)	
}
