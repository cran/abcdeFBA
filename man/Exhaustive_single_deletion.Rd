\name{Exhaustive_single_deletion}
\alias{single knockouts}
\alias{Exhaustive_single_deletion}    
\title{
Exhaustive_single_deletion, a function that deletes each reaction 
in the network one at a time and returns results describing reaction lethality}
\description{
This function takes the argument of type fba_object which would be a FBA object 
and performs an exhaustive deletion of all the reactions in the reaction network 
and returns a list of biomass generated for all the deletions, the sub-optimal 
deletions, the super-optimal deletions and non-lethal deletions, also generates
a pdf containing a distribution of the fatal reactions according 
to their sub-systems and a histogram of the biomass distribution 
for each deletion.}

\usage{
Exhaustive_single_deletion(fba_object,reactions,plot_to_file)
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
	
\item{reactions}{
is a vector containing the reaction numbers to be deleted,if none are supplied all the
reactions in the model are deleted iteratively.}

\item{plot_to_file}{
is a boolean indicating if a PDF output of the results of the single deletion 
should be made.}
}

\examples{
#Performing an exhaustive single reaction deletion
data(Ecoli_core)
#Results<-Exhaustive_single_deletion(Ecoli_core)
}

