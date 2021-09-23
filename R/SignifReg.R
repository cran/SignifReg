SignifReg <-
function(fit,scope, alpha = 0.05,direction = "forward",criterion = "p-value",adjust.method = "fdr",trace=FALSE)UseMethod("SignifReg")
SignifReg.default <-
function(fit,scope, alpha = 0.05,direction = "forward",criterion = "p-value",adjust.method = "fdr",trace=FALSE){
    

  
  #####errors messages for wrong inputs#####
    
    #if(criterion=="p-value"){
    #  if(direction=="both"){
    #    stop("\ncriterion p-value cannot be used for both\n")
    #  }
    #}
    
    if(direction != "forward" && direction != "backward" && direction != "both"){
        stop("\ndirection should be one of the following: forward, backward, and both\n")
    }
    if ((adjust.method != "fdr") && (adjust.method != "holm") && (adjust.method != "hochberg") && (adjust.method != "hommel") && (adjust.method != "bonferroni") && (adjust.method != "BH") && (adjust.method != "BY") && (adjust.method != "none"))
        stop("adjust.method must be a valid method for p.adjust()")
    if(criterion != "p-value" && criterion != "AIC" && criterion != "BIC" && criterion != "r-adj" && criterion != "PRESS"){
        stop("\ncriterion should be one of the following: p-value, AIC, BIC, r-adj, or PRESS\n")
    }
    if(alpha>1 || alpha <0){
        stop("\nSignificance level alpha should be smaller than 1 and larger than 0\n")
    }
    if(trace != FALSE && trace != TRUE){
        stop("\ntrace should be either TRUE or FALSE\n")
    }
    if (criterion == "r-adj")
    {
      if (is.null(summary(fit)$"adj.r.squared"))
         stop("\nr-adj not a valid criterion for glm\n")
    }
       
    
    
    
    
    ################################# start  here######################################
    ################################# start  here######################################
    ################################# start  here######################################
    ################################# start  here######################################
    ################################# start  here######################################
   


    fit_pval <- drop1(fit, test="F")$"Pr(>F)"[-1]
    if(length(fit_pval) == 0){ #for the null model
      passed.correction  <- max_pval <- max_vif <- NA
    }else{
      max_pval <- max(fit_pval) #the max p-value of the current model
      passed.correction <- ifelse(sum(p.adjust(fit_pval, method = adjust.method) <= alpha) < length(fit_pval), FALSE, TRUE) #pvalue cut-off (no correction)
      max_vif = ifelse((length(attr(terms(fit),"term.labels")) > 1), max(vif(fit)), NA) 
    }
    adj.rsq = ifelse(is.null(summary(fit)$"adj.r.squared"), NA, summary(fit)$"adj.r.squared")
    steps.info = data.frame("Step" = "", "Df" = NA, Deviance = "", "Resid Df" = fit$df.residual, "Resid Dev" = deviance(fit), AIC = AIC(fit), BIC = BIC(fit), "adj.rsq" = adj.rsq, "PRESS" = sum((residuals(fit)/(1-lm.influence(fit)$hat))^2), "max_pvalue" = max_pval, "max VIF" = max_vif, "correction" = passed.correction)
    colnames(steps.info)[12] = paste("pass", adjust.method, "correction")
    steps.info[,5:11] <- data.frame(lapply(steps.info[,5:11], round, 5)) 
    
    
    ################################# forward process ######################################
    if (direction == "forward")
    {       
       	if (!missing(scope)) 	
	    	repeat
    		{
	          if(trace==TRUE) 
    	      	print(fit)
    	    	
		      	new_fit = add1SignifReg(fit, scope = scope, alpha = alpha, criterion = criterion, adjust.method = adjust.method, print.step = trace)
		      	if (identical(new_fit, fit))
		      		break
			

		       	fit = new_fit        	
		       	steps.info = rbind(steps.info, fit$step.info)
			
        }	
     }   	
        	
    ################################# backward process ######################################
    if (direction == "backward")
    {       
		repeat
		{
	        if(trace==TRUE) 
    	    	print(fit)
    	    	
			new_fit = drop1SignifReg(fit, scope = scope, alpha = alpha, criterion = criterion, adjust.method = adjust.method, print.step = trace)
			if (identical(new_fit, fit))
				break
			fit = new_fit        	
			steps.info = rbind(steps.info, fit$step.info)
		}	
    }   	
        	        	
        	
        	
        	
    ################################# stepwise process ######################################
    if (direction == "both")
    {       
    	
      if (missing(scope)) 
      {
          var_lower = numeric()
          var_upper = attr(terms(formula(fit)), "factors")
          scope = formula(fit)
      }
      else 
      {
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
      scope_var <- factor.scope(var_current, list(add = var_upper, drop = var_lower))    
      if ((length(scope_var$add) == 0) && (length(scope_var$drop) == 0))
          stop("nothing to add or drop")
      
      
      if(trace==TRUE) 
        print(fit)    	    	
      
      old_fit = fit				
    	repeat
    	{

    	  if (length(scope_var$drop) > 0)
    	  {
    	      new_fit_drop = drop1SignifReg(fit, scope = scope, alpha = alpha, criterion = criterion, adjust.method = adjust.method, print.step = FALSE)
    	      new_fit = new_fit_drop
    	      print_steps = drop1summary(fit = fit, scope = scope, alpha = alpha, adjust.method = adjust.method, sort.by = criterion)   
    	  }

    	  if (length(scope_var$add) > 0)
    	  {
    	      new_fit_add = add1SignifReg(fit, scope = scope, alpha = alpha, criterion = criterion, adjust.method = adjust.method, print.step = FALSE)
    	      new_fit = new_fit_add
    	      print_steps = add1summary(fit = fit, scope = scope, alpha = alpha, adjust.method = adjust.method, sort.by = criterion)   
    	  }
    	
    	  
    	  if ((length(scope_var$add) > 0) && (length(scope_var$drop) > 0)) ## choose between step forward and step backward
    	  {
    	    print_steps = rbind(drop1summary(fit = fit, scope = scope, alpha = alpha, adjust.method = adjust.method, sort.by = criterion),print_steps[-which(rownames(print_steps) == "<none>"),])
    	    #### check best model out of adding and dropping
    	    if (criterion == "AIC"){
    	      if (AIC(new_fit_add) < AIC(new_fit_drop))
    	        new_fit = new_fit_add
    	      else 
    	        new_fit = new_fit_drop
    	    } else
    	    if (criterion == "BIC"){
    	      if (BIC(new_fit_add) < BIC(new_fit_drop))
    	        new_fit = new_fit_add
    	      else 
    	        new_fit = new_fit_drop
    	    } else 	    
    	    if (criterion == "r-adj"){
    	      if (summary(new_fit_add)$"adj.r.squared" > summary(new_fit_drop)$"adj.r.squared")
    	        new_fit = new_fit_add
    	      else 
    	        new_fit = new_fit_drop
    	    } else 	    
    	    if (criterion == "PRESS"){
    	      if (sum((residuals(new_fit_add)/(1-lm.influence(new_fit_add)$hat))^2) < sum((residuals(new_fit_drop)/(1-lm.influence(new_fit_drop)$hat))^2))
    	        new_fit = new_fit_add
    	      else 
    	        new_fit = new_fit_drop
    	    } else
    	    if (criterion == "p-value"){
    	      if (identical(fit, new_fit_drop))
    	        new_fit = new_fit_add
    	      else if (identical(fit, new_fit_add))
    	        new_fit = new_fit_drop
    	      else if (max(drop1(new_fit_add, test="F")$"Pr(>F)"[-1]) <  max(drop1(new_fit_drop, test="F")$"Pr(>F)"[-1]))
    	        new_fit = new_fit_add
    	      else 
    	        new_fit = new_fit_drop
    	    }
    	  }
    	  
    	  
			if(trace == TRUE)
			{
			  if (criterion == "p-value")
			    sort_var = 6
			  else if (criterion == "AIC")
			    sort_var = 2
			  else if (criterion == "BIC")
			    sort_var = 3
			  else if (criterion == "r-adj")
			    sort_var = 4
			  else if (criterion == "PRESS")
			    sort_var = 5
			  print(print_steps[order(as.numeric(print_steps[,sort_var])),])
			}
  			if (identical(fit, new_fit))  # cant add nor drop
	  			break

		  	if (identical(new_fit, old_fit))  # stepped back
			  	break

  			old_fit = fit				
	  		fit = new_fit        	
	  		steps.info = rbind(steps.info, fit$step.info)
	  		var_current <- attr(terms(fit), "factors")
	  		scope_var <- factor.scope(var_current, list(add = var_upper, drop = var_lower))    
	  		
	  		
	  		if(trace==TRUE) 
	  		  print(fit)    	    	
    	}	
    }   	
        	        	
        	
 	
    fit$steps.info = steps.info    	
    return(fit)
}
