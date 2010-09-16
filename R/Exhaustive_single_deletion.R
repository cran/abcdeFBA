### Code to delete fluxes in Rglpk for all the reactions and doing an FBA and storing the biomass changes for all the simulations 
Exhaustive_single_deletion<-function(fba_object){
dim_mat<-dim(fba_object$mat)
deletion_number=0
sub_optimal_deletion=0
super_optimal_deletion=0
no_effect_deletion=0
index1=0
index2=0
index3=0
index4=0

status=0
Exhaus_sol_list=0
init_fba_sol<-FBA_solve(fba_object,6)
optimum=init_fba_sol$objective
for(i in 1:dim_mat[2])
	{
	temp_lb<-fba_object$bounds$lower$val[i]
	temp_ub<-fba_object$bounds$upper$val[i]
	
	fba_object$bounds$lower$val[i]=0
	fba_object$bounds$upper$val[i]=0

	fba_sol<-FBA_solve(fba_object,6)
	Exhaus_sol_list[i]<-fba_sol$objective

	if(fba_sol$objective==0 && fba_sol$status==0)
		{
		index1=index1+1
		deletion_number[index1]=i
		}

if(fba_sol$objective>0 && fba_sol$objective<optimum && fba_sol$status==0)
{
index2=index2+1
sub_optimal_deletion[index2]=i
}

if(fba_sol$objective>optimum && fba_sol$status==0)
{
index3=index3+1
super_optimal_deletion[index3]=i
}

if(fba_sol$objective==optimum && fba_sol$status==0)
{
index4=index4+1
no_effect_deletion[index4]=i
}

fba_object$bounds$lower$val[i]=temp_lb
fba_object$bounds$upper$val[i]=temp_ub

}
flux_sings<-list(deletion_number,status)
lethal_dels<-fba_object$reaction_list[deletion_number]
systemic_lethality<-fba_object$sub_system[deletion_number]

single_fatality<-cbind(fba_object$reaction_list[deletion_number],fba_object$sub_system[deletion_number])


write.table(single_fatality,file="single_fatals",sep="\t",row.names=FALSE,col.names=FALSE)
uniq_subsys<-unique(fba_object$sub_system[deletion_number])
total_subsys<-length(fba_object$sub_system[deletion_number])
freq_subsys<-rbind(c(1,length(which(uniq_subsys[1]==fba_object$sub_system[deletion_number]))))

for(i in 2:length(uniq_subsys))
{
freq_subsys<-rbind(freq_subsys,c(i,length(which(uniq_subsys[i]==fba_object$sub_system[deletion_number]))))
}

graph_legend=cbind(freq_subsys,uniq_subsys[freq_subsys[,1]])
write.table(graph_legend,file="Graph_Legend.xls",sep="\t",quote=F,row.names=F,col.names=F)
barplot(freq_subsys[,2],col="blue",xlab="Reaction Subsystem",ylab="Number of Lethal Knockouts")
pdf("singleKOresults.pdf")
barplot(rev(sort(freq_subsys[,2])),col="red",xlab="Reaction Subsystem",ylab="Number of Lethal Knockouts")
hist(Exhaus_sol_list,col="purple",main="Objective F(x) Dist for Single Knockouts",xlab=fba_object$reaction_list[which(fba_object$obj==1)])
dev.off()
result<-list(lethal_dels,systemic_lethality)

Results=list()
Results$biomass_all=Exhaus_sol_list
Results$lethal_dels=deletion_number
Results$sub_optimal_dels=sub_optimal_deletion
Results$super_optimal_dels=super_optimal_deletion
Results$no_effect_deletions=no_effect_deletion
return(Results)
}
