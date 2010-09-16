FLUXOME_DIFFERENTIALS<-function(flux_vector,reaction_number,fba_object){
mean_flux=0
non_linear_fluxes=rep(0,dim(flux_vector)[1])

for(i in 1:dim(flux_vector)[1])
	{
	mean_flux[i]=mean(abs(flux_vector[i,2:16]))
		if(mean_flux[i]>flux_vector[i,1])
		{
		non_linear_fluxes[i]=1
		}
	}

	comparison_vector=cbind(flux_vector[,1],mean_flux,non_linear_fluxes)
	write.table(comparison_vector,file=paste(fba_object$reaction_list[reaction_number],
	"_DiffAnalysis",".xls",sep=""),row.names=fba_object$reaction_list,col.names=FALSE,sep="\t")
}

########################################################################################################################################
PERTURBATION_analysis<-function(reaction_number,fba_object)
{
	message(fba_object$reaction_list[reaction_number])
	if(fba_object$bounds$lower$val[reaction_number]<0 && fba_object$bounds$upper$val[reaction_number]>0){flag=0}#reversible reactions
	if(fba_object$bounds$lower$val[reaction_number]<0 && fba_object$bounds$upper$val[reaction_number]==0){flag=-1}#exchange reaction
	if(fba_object$bounds$lower$val[reaction_number]==0 && fba_object$bounds$upper$val[reaction_number]>0){flag=1}#closed exchange/transport/internal/secretion

	if(flag==0)
	{
		FBA_SOL_WT<-FBA_solve(fba_object)
		#---------------------------------------#	it is best to use the bound from the WT solution for reversible rxns
		if(FBA_SOL_WT$fluxes[reaction_number]>0)
		{
		flux_range=seq(0,round(FBA_SOL_WT$fluxes[reaction_number],9),length.out=20)
			biomass_vector=0		
			flux_vector=FBA_SOL_WT$fluxes			
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$upper$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective
			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}
		}
		#--------------------------------------#
		if(FBA_SOL_WT$fluxes[reaction_number]<0)	
		{
		flux_range=seq(round(FBA_SOL_WT$fluxes[reaction_number],9),0,length.out=20)
			biomass_vector=0		
			flux_vector=FBA_SOL_WT$fluxes			
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$lower$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective
			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}
		}
		#--------------------------------------#	
	}
	
	if(flag==-1)
	{	#---------------------------------------# reactions constrained for secretion but not intake are generally exchange reactions	
		FBA_SOL_WT<-FBA_solve(fba_object)
		flux_range=seq(round(FBA_SOL_WT$fluxes[reaction_number],9),0,length.out=20)		
			biomass_vector=0
			flux_vector=FBA_SOL_WT$fluxes		
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$lower$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective
			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}		
	}

	if(flag==1)	
	{	#-------------------------------------# positively fluxed reactions are probably closed exchange reactions	
		FBA_SOL_WT<-FBA_solve(fba_object)
		if(FBA_SOL_WT$fluxes[reaction_number]==0)		
		{
		flux_range=seq(-20,0,length.out=20)		
			biomass_vector=0
			flux_vector=FBA_SOL_WT$fluxes		
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$lower$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective
			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}
		}

		if(FBA_SOL_WT$fluxes[reaction_number]>0)		
		{
		flux_range=seq(0,round(FBA_SOL_WT$fluxes[reaction_number],9),length.out=20)		
			biomass_vector=0
			flux_vector=FBA_SOL_WT$fluxes		
			for(i in 1:length(flux_range))
			{
			fba_object$bounds$upper$val[reaction_number]=flux_range[i]
			FBA_SOL_TEMP<-FBA_solve(fba_object)
			biomass_vector[i]=FBA_SOL_TEMP$objective
			flux_vector=cbind(flux_vector,FBA_SOL_TEMP$fluxes)		
			}
		}	
	}

	perturbation_result=list()
	perturbation_result$x=flux_range
	perturbation_result$y=biomass_vector

	plot(perturbation_result,type="b",pch="*",col="blue",xlab="mmol/gDW/hr",
	main="Metabolic Perturbation Analysis",sub=fba_object$reaction_list[reaction_number],
	ylab=fba_object$reaction_list[which(fba_object$obj==1)])

	png(paste(fba_object$reaction_list[reaction_number],".png",sep=""))
	plot(perturbation_result,type="b",pch="*",col="blue",xlab="mmol/gDW/hr",
	main="Metabolic Perturbation Analysis",sub=fba_object$reaction_list[reaction_number],
	ylab=fba_object$reaction_list[which(fba_object$obj==1)])
	dev.off()

	unique_list<-sort(unique(fba_object$sub_system))
	vec1=vector()
	for(i in 1:length(unique_list))
	{vec1<-c(vec1,which(unique_list[i]==fba_object$sub_system))}
	
	flux_SS<-cbind(fba_object$reaction_list[vec1],flux_vector[vec1,],fba_object$reaction_list[vec1])
	write.table(flux_SS,file=paste("Ramp",fba_object$reaction_list[reaction_number],
	".xls"),sep="\t",quote=F,row.names=fba_object$sub_system[vec1],col.names=F)

	FLUXOME_DIFFERENTIALS(flux_vector,reaction_number,fba_object)
	return(perturbation_result)
}
