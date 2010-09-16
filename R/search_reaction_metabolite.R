SEARCH_reaction<-function(react_name,fba_object)
{
	print(fba_object$reaction_list[grep(react_name,fba_object$reaction_list,ignore.case=TRUE)])
	print(grep(react_name,fba_object$reaction_list,ignore.case=TRUE))
}

SEARCH_metabolite<-function(metabolite_name,fba_object)
{
	print(cbind(fba_object$metabolite_name[grep(metabolite_name,fba_object$metabolite_name,
	ignore.case=TRUE)],grep(metabolite_name,fba_object$metabolite_name,ignore.case=TRUE),
	fba_object$comp_name[fba_object$compartment[grep(metabolite_name,fba_object$metabolite_name,
	ignore.case=TRUE)]]))
}
