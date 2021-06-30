add1SignifReg <-
function(fit, scope, alpha = 0.05, criterion = "p-value", adjust.method = "fdr", override = FALSE, print.step = FALSE)UseMethod("add1SignifReg")

add1SignifReg.default <- function(fit, scope, alpha = 0.05, criterion = "p-value", adjust.method = "fdr", override = FALSE, print.step = FALSE)
{
    
    #####errors messages for wrong inputs#####
    if(criterion != "p-value" && criterion != "AIC" && criterion != "BIC" && criterion != "r-adj" && criterion != "PRESS"){
        stop("\ncriterion should be one of the following: p-value, AIC, BIC, r-adj, or PRESS\n")
    }
    if(criterion == "r-adj" && is.null(summary(fit)$"adj.r.squared")){
      stop("\nr-adj not a valid criterion for a glm model\n")
    }
    if ((adjust.method != "fdr") && (adjust.method != "holm") && (adjust.method != "hochberg") && (adjust.method != "hommel") && (adjust.method != "bonferroni") && (adjust.method != "BH") && (adjust.method != "BY") && (adjust.method != "none"))
       stop("adjust.method must be a valid method for p.adjust()")
  
    if(alpha>1 || alpha <0){
        stop("\nSignificance level alpha should be smaller than 1 and larger than 0\n")
    }
    if(override != FALSE && override != TRUE){
        stop("\noverride should be either TRUE or FALSE\n")
    }
    
    
   
   
   
   
   
   
  deviance_current_model = deviance(fit)
  
   
   

################# start here ################################################
################# start here ################################################
################# start here ################################################
################# start here ################################################
################# start here ################################################
################# start here ################################################
   
   
	out_tab = add1summary(fit = fit, scope = scope, alpha = alpha, adjust.method = adjust.method, sort.by = criterion)   
	if ((print.step == TRUE) && (!is.null(out_tab)))
    print(out_tab)   
  
	if (!is.null(out_tab)) 
	{
	  correction_value = out_tab[,8] ## p-value correction True of False for each prospective model
	  vars_pass_correction = setdiff(which(correction_value == TRUE), which(rownames(out_tab) == "<none>"))
	
    var_index = NULL
    if(criterion == "p-value")
    {
       	crit_value <- as.numeric(out_tab[,6]) 
		    if (length(vars_pass_correction) > 0)
    		{
	        var_index <- which.min(crit_value[vars_pass_correction])[1] #select a model having the smallest AIC, BIC. index 0 is the current model
	        fit <- update(fit, paste(".~. ", rownames(out_tab)[vars_pass_correction[var_index]]))			
	     	} else if (override == TRUE)
	     	{
	     	  var_index <- which.min(crit_value)[1] #select a model having the smallest AIC, BIC. index 0 is the current model
	     	  fit <- update(fit, paste(".~. ", rownames(out_tab)[var_index]))			
	     	}
    	
    } else 
    {
    	if(criterion == "AIC")
        	crit_value <- as.numeric(out_tab[,2]) 
	    else if(criterion == "BIC") 
    	    crit_value <- as.numeric(out_tab[,3]) 
	    else if(criterion == "r-adj")
    	    crit_value <- -as.numeric(out_tab[,4]) 
	    else if(criterion == "PRESS")
    	    crit_value <- as.numeric(out_tab[,5]) 

	   	### For AIC the new model needs to have smaller AIC and pass correction
	   	if ( (length(vars_pass_correction) > 0) && (min(crit_value[vars_pass_correction]) < (crit_value[which(rownames(out_tab) == "<none>")])))
	  	{
	        var_index <- which.min(crit_value[vars_pass_correction])[1] #select a model having the smallest AIC, BIC. index 0 is the current model
	        fit <- update(fit, paste(".~. ", rownames(out_tab)[vars_pass_correction[var_index]]))			
	  	} else if (override == TRUE)
	  	{
	  	  var_index <- which.min(crit_value)[1] #select a model having the smallest AIC, BIC. index 0 is the current model
	  	  fit <- update(fit, paste(".~. ", rownames(out_tab)[var_index]))			
	  	}

    } 

    
    if (!is.null(var_index))
    {
       fit$step.info = data.frame("Step" = rownames(out_tab)[vars_pass_correction[var_index]], "Df" = 1, Deviance = round((deviance(fit) - deviance_current_model),5), "Resid Df" = fit$df.residual, "Resid Dev" = as.numeric(out_tab[vars_pass_correction[var_index],1]), AIC = as.numeric(out_tab[vars_pass_correction[var_index],2]), BIC = as.numeric(out_tab[vars_pass_correction[var_index],3]), "adj.rsq" = as.numeric(out_tab[vars_pass_correction[var_index],4]), "PRESS" = as.numeric(out_tab[vars_pass_correction[var_index],5]), "max_pvalue" = as.numeric(out_tab[vars_pass_correction[var_index],6]), "max VIF" = as.numeric(out_tab[vars_pass_correction[var_index],7]), "correction" = out_tab[vars_pass_correction[var_index],8])
       colnames(fit$step.info)[12] = paste("pass", adjust.method, "correction")
       fit$step.info[,5:11] <- data.frame(lapply(fit$step.info[,5:11], round, 5)) 
    }
  }
   
   return(fit)
   
}   
   
   
   
   
   
   