add1summary <- function(fit, scope, alpha = 0.05, adjust.method = "fdr", sort.by = "p-value")UseMethod("add1summary")
add1summary.default <- function(fit, scope, alpha = 0.05, adjust.method = "fdr", sort.by = "p-value"){
    
    #####errors messages for wrong inputs#####
    if(alpha>1 || alpha <0){
        stop("\nSignificance level alpha should be smaller than 1 and larger than 0\n")
    }
    
    if ((adjust.method != "fdr") && (adjust.method != "holm") && (adjust.method != "hochberg") && (adjust.method != "hommel") && (adjust.method != "bonferroni") && (adjust.method != "BH") && (adjust.method != "BY") && (adjust.method != "none"))
      stop("adjust.method must be a valid method for p.adjust()")








################# start here ################################################
################# start here ################################################
################# start here ################################################
################# start here ################################################
################# start here ################################################
################# start here ################################################

    if (missing(scope)) 
		stop("no terms in scope")
    else {
        if (is.list(scope)) {
            var_lower <- if (!is.null(var_lower <- scope$lower)) 
                attr(terms(update.formula(fit, var_lower)), "factors")
            else numeric()
            var_upper <- if (!is.null(var_upper <- scope$upper)) 
                attr(terms(update.formula(fit, var_upper)), "factors")
        }
        else {
            var_upper <- if (!is.null(var_upper <- scope)) 
                attr(terms(update.formula(fit, scope)), "factors")
            var_lower <- numeric()
        }
    }
    
    var_current <- attr(terms(fit), "factors")
    scope <- factor.scope(var_current, list(add = var_upper, drop = var_lower))    	

	if (length(scope$add) > 0)
	{
		################# Outcome Table #################
	    out_tab <- data.frame("Resid Dev" = rep(NA,length(scope$add)+1), "AIC" = rep(NA,length(scope$add)+1), "BIC" = rep(NA,length(scope$add)+1), "adj.rsq" = rep(NA,length(scope$add)+1), "PRESS" = rep(NA,length(scope$add)+1), "max_pvalue" = rep(NA,length(scope$add)+1), "max VIF" = rep(NA,length(scope$add)+1), "correction" = rep(NA,length(scope$add)+1))
	    colnames(out_tab)[8] = paste("pass", adjust.method, "correction")

	    rownames(out_tab) <- c("<none>", paste("+",scope$add))
	    
	    ########### Current Model #############################
	    out_tab[1,1:5] <- c(deviance(fit), 
	                        AIC(fit), 
	                        BIC(fit), 
	                        ifelse(is.null(summary(fit)$"adj.r.squared"), NA, summary(fit)$"adj.r.squared"), 
	                        sum((residuals(fit)/(1-lm.influence(fit)$hat))^2)) ## last one is PRESS residuals
	    out_tab[1,7] = ifelse((length(attr(terms(fit),"term.labels")) > 1), max(vif(fit)), NA) 
	
	    fit_pval <- drop1(fit, test="F")$"Pr(>F)"[-1] #save p-values of the current model
	    if(length(fit_pval) == 0){ #for the null model
	        out_tab[1,6] <- out_tab[1,8]  <- NA
	    }else{
	        out_tab[1,6] <- max(fit_pval) #the max p-value of the current model
	        out_tab[1,8] <- ifelse(sum(p.adjust(fit_pval, method = adjust.method) <= alpha) < length(fit_pval), FALSE, TRUE) #pvalue cut-off (no correction)
	    }
	    
		########### Prospective Models #############################
	    for (n in 1:length(scope$add))
	    {
	        fit2 <- update(fit, paste(".~. +", scope$add[n]))
	        out_tab[n+1,1:5] <- c(deviance(fit2), 
	                              AIC(fit2), 
	                              BIC(fit2), 
	                              ifelse(is.null(summary(fit2)$"adj.r.squared"), NA, summary(fit2)$"adj.r.squared"), 
	                              sum((residuals(fit2)/(1-lm.influence(fit2)$hat))^2))
		    out_tab[n+1,7] = ifelse((length(attr(terms(fit2),"term.labels")) > 1), max(vif(fit2)), NA) 
	        
	        fit2_pval <- drop1(fit2, test="F")$"Pr(>F)"[-1] #save p-values of each model
	        out_tab[n+1,6] <- max(fit2_pval) #the max p-value of each model
	        out_tab[n+1,8] <- ifelse(sum(p.adjust(fit2_pval, method = adjust.method) <= alpha) < length(fit2_pval), FALSE, TRUE)  #pvalue cut-off (no correction)
	    }
	    
		##### sort table ##################################################
		if (sort.by == "p-value")
			sort_var = 6
		else if (sort.by == "AIC")
			sort_var = 2
		else if (sort.by == "BIC")
			sort_var = 3
		else if (sort.by == "r-adj")
			sort_var = 4
		else if (sort.by == "PRESS")
			sort_var = 5
	    out_tab = out_tab[order(as.numeric(out_tab[,sort_var])),]

	    out_tab[,1:7] <- data.frame(lapply(out_tab[,1:7], round, 5))
	} else
	{
		out_tab = NULL ### there is nothing in scope to add
		print("no terms in scope for adding to object.")
	}
	
    return(out_tab)





}
