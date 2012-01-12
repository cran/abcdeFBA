# take input of query genes
# check against model and remove the genes from query which do not exist in the model
# retrieve all gpr's related to the query genes
# cycle through all gpr's
# obtained all the gene names in the gpr
# match genes and convert them to 0, take setdiff and convert to 1 the rest
# change AND and OR to & and |
# interpret the boolean expression and solve for status using parse eval

Gene_del<-function(query_genes=NULL,fba_object,return_reactions=FALSE)
{
require(abcdeFBA)
	#----------- Check genes -----------#
genes_NA_model<-setdiff(query_genes,fba_object$all_genes)
if(length(genes_NA_model)>0)
	{
	message("These genes are not present in the model")
	message(genes_NA_model)
	}
query_genes<-setdiff(query_genes,genes_NA_model)
if(length(query_genes)>0)
{
	gpr_ix=vector()

	for(i in 1:length(query_genes))
		{
		gpr_ix<-c(gpr_ix,grep(query_genes[i],fba_object$gpr))
		}
		gpr_ix<-unique(gpr_ix)		
		enlisted_gprs<-fba_object$gpr[gpr_ix]

	# go through the gpr, replace "(" with "( " and ")" with " )"
	for(i in 1:length(enlisted_gprs))
	{
	enlisted_gprs[i]<-gsub("\\(","( ",enlisted_gprs[i])
	enlisted_gprs[i]<-gsub("\\)"," )",enlisted_gprs[i])
	enlisted_gprs[i]<-gsub("and","&",enlisted_gprs[i])
	enlisted_gprs[i]<-gsub("or","|",enlisted_gprs[i])
	
	gpr_split<-strsplit(enlisted_gprs[i]," ")[[1]]
	allgenes_gpr<-gpr_split[grep("[0-9]",gpr_split)]
	not_in_query<-setdiff(allgenes_gpr,query_genes)
	
		if(length(not_in_query!=0))	
		{	
			for(j in 1:length(not_in_query))
				{
				enlisted_gprs[i]<-gsub(not_in_query[j],"1",enlisted_gprs[i])	
				}
		}
	
		for(j in 1:length(query_genes))
			{
			enlisted_gprs[i]<-gsub(query_genes[j],"0",enlisted_gprs[i])
			}
	}
	
	Effect=vector()
	
		for(i in 1:length(enlisted_gprs))
		{
			Effect<-c(Effect,eval(parse(text=enlisted_gprs[i])))	
		
		}
	
		KillSwitch<-rep(TRUE,length(fba_object$gpr))
		KillSwitch[gpr_ix]<-Effect
		Switch_off<-which(KillSwitch==0)
	
		if(return_reactions==TRUE)
			{return(Switch_off)}
		
		if(return_reactions==FALSE)
			{return(CHANGE_RXN_BOUNDS(Switch_off,fba_object,0,0))}
}else(message("the query list is empty"))
}
