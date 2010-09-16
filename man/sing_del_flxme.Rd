\name{SINGLE_DEL_FLUXOME}
\alias{SINGLE_DEL_FLUXOME}
\alias{fluxome graphs}
\alias{comparison of single deletion flux distributions}
\title{SINGLE_DEL_FLUXOME, a function to generate comparative fluxome graphs}

\description{To analyse the effect of reaction deletions 
	on the fluxome requires a contextual visualization, simple 
	plots give little insight on what the results of the simulation 
	mean. This function uses the annotations inherent to SBML 
	models and generates comparative overlapping fluxome bar 
	graphs depicting the overlap/change in flux based on 
	Sub-system wise classification and generates PDF's 
	for the same. The color scheme is green for wild-type 
	and red for mutant. Overlaps of red and green generate 
	brown while overshot mutant fluxes show up as magenta-pink,
	also separate PDF's are generated for increased and decreased fluxes.}
	\usage{
	    SINGLE_DEL_FLUXOME(fba_object, deletion_number)
		}
\arguments{
	\item{fba_object}{is a list containing the data required 
	to perform flux balance analysis. The elements of the list 
	are mat which is the stoichiometric matrix, dir which gives 
	the direction of the equality constraints, obj specifies the 
	objective function for the simulation, bounds specifies the 
	lower and upper inequality constraints, rhs is the right hand 
	side of the steady state expression, types refers to the numeric 
	nature of the variables which in case of FBA happens to be 
	"Continuous", max is a Boolean specifying the type of optimization,
	"Maximization" by default".}
      \item{deletion_number}{Reaction number of the reaction to be 
	deleted. This number may be obtained by using the function 
	SEARCH_reaction for a text search through the fba_object-reaction list.}
	}
\examples{
	#Compare the flux distributions of the mutant with the wild-type
	data(Ecoli_core)	
	SINGLE_DEL_FLUXOME(Ecoli_core,36)	
	}
