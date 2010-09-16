\name{PHPP}
\alias{PHPP}    
\alias{phenotypic phase plane analysis}
\title{PHPP, phenotypic phase plane analysis}
\description{this function helps in performing a phenotypic 
	phase plane analysis, a visualization of the effect 
	of two input fluxes on the value of the objective function.}
\usage{
PHPP(reaction_number,fba_object,PCS,flux_range)
}
\arguments{
	\item{fba_object}{is a list containing the data required to 
		perform flux balance analysis. The elements 
		of the list are mat which is the stoichiometric 
		matrix, dir which gives the direction of the 
		equality constraints, obj specifies the objective 
		function for the simulation, bounds specifies the 
		lower and upper inequality constraints, rhs is the 
		right hand side of the steady state expression, 
		types refers to the numeric nature of the variables 
		which in case of FBA happens to be "Continuous", 
		max is a Boolean specifying the type of optimization,
		"Maximization" by default".}
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
}
\examples{
	#Performing a phenotypic phase plane analysis of 
	#glucose and oxygen in Core E.coli Metabolism
	data(Ecoli_core)
	PHPP(reaction_number=c(28,36),fba_object=Ecoli_core,
	PCS="glucose",flux_range=c(1,15))
	# a menu pops up asking to select the primary carbon source, select D glucose for Ecoli_core	
	}
