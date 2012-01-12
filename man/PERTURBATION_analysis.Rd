\name{PERTURBATION_analysis}
\alias{PERTURBATION_analysis}
\alias{ramp-up}
\alias{ramp-down} 
\alias{robustness analysis}
\title{PERTURBATION_analysis, a function for robustness analysis.}
\description{ Robustness analysis is described procedurally in the COBRA-ToolBox manual. This function encodes the basic principle of the procedure and returns the result as an X versus Y list and generates a plot on successful completion.}
\usage{
    PERTURBATION_analysis(reaction_number,fba_object,y_axis_rxn=NULL,plot_to_file=FALSE,write_FLD_file=FALSE,ret_FLD_matrix=FALSE)
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
\item{reaction_number}{Reaction number of the reaction to be ramped down. 
	This number may be obtained by using the function SEARCH_reaction 
	for a text search through the FBA_obj-reaction list.}
\item{y_axis_rxn}{Reaction to be plotted on the y axis in case it is not 
	the objective fuction.}
\item{plot_to_file}{Boolean indicating if the graph to be printed to a file.}
\item{write_FLD_file}{The flux vector for each point of perturbation is 
	stored into the columns of a matrix, this boolean indicates if 
	that matrix should be written out to a tab-delimited file. }
\item{ret_FLD_matrix}{The flux vector for each point of perturbation is 
	stored into the columns of a matrix, this boolean indicates if 
	that matrix should be added to the list returned by this function. }
	}


\examples{
	#Perturbation analysis of fluxes
	data(Ecoli_core)
	Ec_xy_O2<-PERTURBATION_analysis(reaction_number=36,fba_object=Ecoli_core)
}

