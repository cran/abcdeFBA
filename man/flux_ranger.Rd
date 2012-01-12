\name{Flux_Ranger}
\alias{Flux_Ranger}
\description{a function to create a linear ramp for the specified flux
automatically/intelligently.}
\title{Flux_Ranger - a function to create a flux ramp.}
\usage{Flux_Ranger(reaction_number,fba_object,divs,art_limit_range)}

\arguments{
\item{reaction_number}{reaction of which the ramp is to be created}
\item{fba_object}{list containing the necessary elements making up
the flux balance model.}
\item{divs}{the number of divisions for the ramp}
\item{art_limit_range}{in case the "intelligent ramp" misbehaves
you can force the ramp into a particular numeric range by concatenating the 
lower and upper numeric limits into this variable.}
}

\examples{
#Creating a ramp for any reaction
data(Ecoli_core)
flux_range<-Flux_Ranger(reaction_number=12,fba_object=Ecoli_core,divs=10)	
}
