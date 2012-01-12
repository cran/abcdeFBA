# convert_B takes in a complex boolean gene protein relation and converts it into a mathematical form suitable for evaluation
# query genes contains the genes which are to be turned off in the boolean gpr by setting them to zero.

convert_B<-function(gpr,query_genes) # converts the GPR relation from SBML to a format that can be evaluated for TRUTH---- "(" "0","+","1" ")"-
{

	boo_exp=vector()
	gpr_split<-strsplit(gpr," ")[[1]]

	if(length(gpr_split)==1){ # if the boolean exp is only one gene name
		for(i in 1:length(query_genes))
		{
			if(length(grep(query_genes[i],gpr_split))!=0)
			{boo_exp[i]="0"}else{boo_exp[i]="1"}
		}
		return(c("(",boo_exp,"+",boo_exp,")")) # padding to avoid error in con_to_dec
		}

	if(length(gpr_split)>1){ # functions if boolean exp is multi-gene logical
	for(i in 1:length(gpr_split))
		{
		if(gpr_split[i]=="(")
		{boo_exp[i]="("}

		if(gpr_split[i]=="")
		{boo_exp[i]="_"}

		if(length(grep("b",gpr_split[i]))>0)
		{
			for(j in 1:length(query_genes))
			{
				if(length(grep(query_genes[j],gpr_split[i]))!=0)
				{
				boo_exp[i]="0"
				}
			}
		}

		if(gpr_split[i]==")")
		{boo_exp[i]=")"}
	
		if(gpr_split[i]=="and")
		{boo_exp[i]="*"}

		if(gpr_split[i]=="or")
		{boo_exp[i]="+"}
	
		}
		
	boo_exp[which(is.na(boo_exp))]="1"
	boo_exp<-boo_exp[which(boo_exp!="_")]
	return(boo_exp)
	}
}

eval_N_bit<-function(enbit) # evaluates an N-bit boolean expression of the form "(", "(", "0", "*", "1",")", "+", "(", "1", "*", "1" ")"")"
		{
		bits=vector()
		operand=vector()

		for(i in 2:(length(enbit)-1))
			{
			if(i%%2==0)
				{operand=c(operand,enbit[i])}
			if(i%%2==1)
				{bits=c(bits,enbit[i])}
			}
		if(unique(bits)=="+")
		{return(max(as.numeric(operand)))}
		if(unique(bits)=="*")
		{return(prod(as.numeric(operand)))}
		}

pop_ix<-function(lb)
		{
		out=list()
		out$pop=lb$lb[length(lb$lb)]
		out$lb=rev(rev(lb$lb)[-1])	
		return(out)	
		}

solve_boolean<-function(boo_exp)
		{
		lb=list()
		lb$pop=vector()
		lb$lb=vector()
	
		i=0
		while(length(which(boo_exp=="("))!=1)
			{
			i=i+1
			if(boo_exp[i]=="("){lb$lb<-c(lb$lb,i)}
			if(boo_exp[i]==")")
				{
				beg=pop_ix(lb)$pop
				end=i
				lb=pop_ix(lb)
				etb<-as.character(eval_N_bit(boo_exp[beg:end]))
				boo_exp=c(boo_exp[1:(beg-1)],etb,boo_exp[(end+1):length(boo_exp)])		
				i=1
				lb$lb=vector()
				}
			}
	
		return(as.numeric(eval_N_bit(boo_exp)))
		}
# first take query gene_set and parse through the fba_object$gpr relations and store the indices
#  

Gene_del<-function(fba_object,gene_set){

	gpr_set=vector()
	print(gpr_set)	
	for(i in 1:length(gene_set))
	{
	gpr_set<-c(gpr_set,grep(gene_set[i],fba_object$gpr))
	}
	kill_switch_vector=rep(1,length(fba_object$gpr))

	for(i in 1:length(gpr_set))
	{
	boo_exp<-convert_B(fba_object$gpr[gpr_set[i]],gene_set)
	kill_switch_vector[gpr_set[i]]<-solve_boolean(boo_exp)
	}
	g_del_mutant<-CHANGE_RXN_BOUNDS(which(kill_switch_vector!=1),fba_object,0,0)	
	return(g_del_mutant)
}

