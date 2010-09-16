###- PHPP is a linear optimization procedure which can be used to study the value of the objective function (a desired phenotype) as two variables (external substrates) vary simultaneously.###
 
PHPP<-function(reaction_number=c(28,36),fba_object,PCS="D glucose",flux_range=c(1,15)){
candidates<-grep(PCS,fba_object$reaction_list,ignore.case=TRUE) 
message("\nProbable Carbon Sources:") 
print(fba_object$reaction_list[candidates])
C_S<-as.numeric(readline("\nSelect the Carbon source from the list above by serial number"))

	if(length(which(reaction_number==candidates[C_S]))==0)
	{
	message("Alternate carbon source PhPP")
	fba_object$bounds$lower$val[candidates[C_S]]=0
	OBJ_MATRIX<-matrix(0,flux_range[2],flux_range[2])
	for(i in flux_range[1]:flux_range[2])
		{
		print(i)
		for(j in flux_range[1]:flux_range[2])
			{
			fba_object$bounds$lower$val[reaction_number[1]]=-1*i
			fba_object$bounds$lower$val[reaction_number[2]]=-1*j
			FBA_SOL_TEMP<-FBA_solve(fba_object,7)
			OBJ_MATRIX[i,j]=FBA_SOL_TEMP$objective
			}
		}
			zlim <- range(flux_range)
			zlen <- zlim[2] - zlim[1] + 1
			colorlut<-rainbow(zlen)
			col <- colorlut[ OBJ_MATRIX*(zlim[1]+5) ] # assign colors to heights for each point
			persp3d(c(flux_range[1]:flux_range[2]),c(flux_range[1]:flux_range[2]),
			OBJ_MATRIX,xlim=range(flux_range),ylim=range(flux_range),zlim=range(OBJ_MATRIX),
			xlab=fba_object$reaction_list[reaction_number[1]],ylab=fba_object$reaction_list[reaction_number[2]],
			zlab=fba_object$reaction_list[which(fba_object$obj==1)],color=col, alpha=0.95, back="lines")
	}
###############################################################################################################################################
	if(length(which(reaction_number==candidates[C_S]))>0)
	{
	message("Primary carbon source PhPP")
	OBJ_MATRIX<-matrix(0,flux_range[2],flux_range[2])
	for(i in flux_range[1]:flux_range[2])
		{
		print(i)
		for(j in flux_range[1]:flux_range[2])
			{
			fba_object$bounds$lower$val[reaction_number[1]]=i*(-1)
			fba_object$bounds$lower$val[reaction_number[2]]=j*(-1)
			FBA_SOL_TEMP<-FBA_solve(fba_object,7)
			OBJ_MATRIX[i,j]=FBA_SOL_TEMP$objective
			}
		}
			zlim <- range(flux_range)
			zlen <- zlim[2] - zlim[1] + 1
			colorlut<-rainbow(zlen)
			col <- colorlut[ OBJ_MATRIX*(zlim[1]+5) ] # assign colors to heights for each point
			persp3d(c(flux_range[1]:flux_range[2]),c(flux_range[1]:flux_range[2]),OBJ_MATRIX,
			xlim=range(flux_range),ylim=range(flux_range),zlim=range(OBJ_MATRIX),
			xlab=fba_object$reaction_list[reaction_number[1]],ylab=fba_object$reaction_list[reaction_number[2]],
			zlab=fba_object$reaction_list[which(fba_object$obj==1)],color=col, alpha=0.95, back="lines")	
	}
}
