\name{PHPP}
\alias{PHPP}    
\alias{phenotypic phase plane analysis}
\title{PHPP, phenotypic phase plane analysis}
\description{this function helps in performing a phenotypic 
	phase plane analysis, a visualization of the effect 
	of two input fluxes on the value of the objective function.}
\usage{
	PHPP(reaction_number,fba_object,PCS,flux_range,
	ret_OBJ_mat,surf_col,divs,dimension,animate,objective)
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

\item{PCS}{a string to search for the primary carbon source 
of the model, for example PCS="D glucose" will display all 
the reactions with that string in them, you need to choose 
the appropriate reaction number of the Carbon source before 
continuing. In case you are performing a PhPP for alternate 
carbon sources then the primary carbon source is supposed to 
be shut-down to get a correct picture of the flux cone.}

\item{reaction_number}{two reaction numbers that specify the two input 
fluxes that make the x and y axes of the PhPP}

\item{flux_range}{the range between which the input fluxes are to be varied}

\item{ret_OBJ_mat}{boolean indicating if the PHPP matrix should be returned}

\item{surf_col}{character string for surface color}

\item{divs}{number of divisions of ramp- note this increases computation 
time of PHPP by n-squared}

\item{dimension}{characters- "2" or "3" specifying if the PHPP should be 
viewed in 2-D or 3-D}

\item{animate}{boolean indicating if the 3-D plot should be animated;
spins it around once on each axes}

\item{objective}{reaction_number to be used as the objective parameter for PhPP}
}

\examples{
#Performing a phenotypic phase plane analysis of 
#glucose and oxygen in Core E.coli Metabolism
data(Ecoli_core)
#PHPP(reaction_number=c(28,36),fba_object=Ecoli_core,
#PCS="glucose",flux_range=c(1,15),ret_OBJ_mat=FALSE,surf_col="red")
# a menu pops up asking to select the primary carbon source, 
#select D glucose for Ecoli_core	
}
