Exhaustive_double_deletion<-function(fba_object,thread_no=0,core_number=1){

message("Generating Reaction Combinations...")
reacn_combos<-combn(1:dim(fba_object$mat)[2],2)
message("Done!")

j=0
ko_1x=0
ko_2x=0
ko_stat=0
message("Starting Double Knockout Simulation...")
	if(core_number==1)
	{	
		for(i in 1:dim(reacn_combos)[2])
			{
			print(paste(i,dim(reacn_combos)[2],sep=" "))
			temp_lb1=fba_object$bounds$lower$val[reacn_combos[,i][1]]
			temp_lb2=fba_object$bounds$lower$val[reacn_combos[,i][2]]
			temp_ub1=fba_object$bounds$upper$val[reacn_combos[,i][1]]
			temp_ub2=fba_object$bounds$upper$val[reacn_combos[,i][2]]
			
			fba_object$bounds$lower$val[reacn_combos[,i][1]]=0
			fba_object$bounds$lower$val[reacn_combos[,i][2]]=0
			fba_object$bounds$upper$val[reacn_combos[,i][1]]=0
			fba_object$bounds$upper$val[reacn_combos[,i][2]]=0
			
			FBA_MUTANT<-FBA_solve(fba_object)
		
			if(FBA_MUTANT$objective==0)
			{
			j=j+1
			ko_1x[j]=reacn_combos[,i][1]
			ko_2x[j]=reacn_combos[,i][2]
			ko_stat[j]=FBA_MUTANT$status
			}	
		
			fba_object$bounds$lower$val[reacn_combos[,i][1]]=temp_lb1
			fba_object$bounds$lower$val[reacn_combos[,i][2]]=temp_lb2
			fba_object$bounds$upper$val[reacn_combos[,i][1]]=temp_ub1
			fba_object$bounds$upper$val[reacn_combos[,i][2]]=temp_ub2
			}
	message("End of simulation.")
	flux_pairs<-cbind(ko_1x,ko_2x,ko_stat)
	message("Writing output to file...")
	write.table(flux_pairs,file=paste("results",thread_no+1,sep=""),sep="\t",row.names=TRUE, col.names=FALSE,quote=FALSE)
	message("Complete!")
}

	if(core_number>1)
	{
	OP_size<-round(dim(reacn_combos)[2]/core_number)
	I_1=1+(OP_size*thread_no)
	I_2=I_1+OP_size-1

		if(I_2>dim(reacn_combos)[2]){I_2=dim(reacn_combos)[2]}
			for(i in I_1:I_2)
			{
			print(paste(i,I_2,sep=" "))
			temp_lb1=fba_object$bounds$lower$val[reacn_combos[,i][1]]
			temp_lb2=fba_object$bounds$lower$val[reacn_combos[,i][2]]
			temp_ub1=fba_object$bounds$upper$val[reacn_combos[,i][1]]
			temp_ub2=fba_object$bounds$upper$val[reacn_combos[,i][2]]
			
			fba_object$bounds$lower$val[reacn_combos[,i][1]]=0
			fba_object$bounds$lower$val[reacn_combos[,i][2]]=0
			fba_object$bounds$upper$val[reacn_combos[,i][1]]=0
			fba_object$bounds$upper$val[reacn_combos[,i][2]]=0
			
			FBA_MUTANT<-FBA_solve(fba_object)
		
			if(FBA_MUTANT$objective==0)
			{
			j=j+1
			ko_1x[j]=reacn_combos[,i][1]
			ko_2x[j]=reacn_combos[,i][2]
			ko_stat[j]=FBA_MUTANT$status
			}	
		
			fba_object$bounds$lower$val[reacn_combos[,i][1]]=temp_lb1
			fba_object$bounds$lower$val[reacn_combos[,i][2]]=temp_lb2
			fba_object$bounds$upper$val[reacn_combos[,i][1]]=temp_ub1
			fba_object$bounds$upper$val[reacn_combos[,i][2]]=temp_ub2
			}
	message("End of Simulation")
	flux_pairs<-cbind(ko_1x,ko_2x,ko_stat)
	message("Writing output to file...")
	write.table(flux_pairs,file=paste("results",thread_no+1,sep=""),sep="\t",row.names=FALSE, col.names=FALSE)
	message("Complete!")
	}

}

