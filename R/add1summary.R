add1summary <- function(fit, scope=eval(fit$call$data), alpha = 0.05)UseMethod("add1summary")
add1summary.default <- function(fit, scope=eval(fit$call$data), alpha = 0.05){
    
    #####errors messages for wrong inputs#####
    if(alpha>1 || alpha <0){
        stop("\nSignificance level alpha should be smaller than 1 and larger than 0\n")
    }
    
    #######################correction functions########################
    
    #correction function for Bonferroni, FDR, and default
    none <- function(pvalues){ #function for no correction
        max_p <- max(pvalues) #save maximum p-value
        result <- ifelse(max_p <= alpha, TRUE, FALSE)
        return(result) #when satisfy corrected criterion return result=TRUE, otherwise result=FALSE
    }
    bonferroni <- function(pvalues){
        max_p <- max(pvalues) #save maximum p-value
        result <- ifelse(max_p <= alpha/length(pvalues), TRUE, FALSE)#alpha/number of tests
        return(result) #when satisfy corrected criterion return result=TRUE, otherwise result=FALSE
    }
    fdr <- function(pvalues){
        pi <- sort(pvalues) #sorting pvalues in ascending order
        j <- c(1:length(pi))
        t <- 0
        for(n in j){
            if(pi[n] <= (n/length(pvalues))*(alpha/sum(j^(-1)))) t <- 1+t
        }
        result <- ifelse (t == length(pi), TRUE, FALSE) #if true, that means every variable satisfies the FDR criterion.
        return(result) #when satisfy corrected criterion return result=TRUE, otherwise result=FALSE
    }
    

    #######################PRESS residuals########################
	PRESS <- function(linear.model) {
	  #' calculate the predictive residuals
	  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
	  #' calculate the PRESS
	  PRESS <- sum(pr^2)
	  return(PRESS)
	}

    ####################### VIF ########################
	vif <- function(mod, ...) 
	{
	    if (any(is.na(coef(mod)))) 
	        stop ("there are aliased coefficients in the model")
	    v <- vcov(mod)
	    assign <- attr(model.matrix(mod), "assign")
	    if (names(coefficients(mod)[1]) == "(Intercept)") {
	        v <- v[-1, -1]
	        assign <- assign[-1]
	    }
	    else warning("No intercept: vifs may not be sensible.")
	    terms <- labels(terms(mod))
	    n.terms <- length(terms)
	    if (n.terms < 2) stop("model contains fewer than 2 terms")
	    R <- cov2cor(v)
	    detR <- det(R)
	    result <- matrix(0, n.terms, 3)
	    rownames(result) <- terms
	    colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
	    for (term in 1:n.terms) {
	        subs <- which(assign == term)
	        result[term, 1] <- det(as.matrix(R[subs, subs])) *
	            det(as.matrix(R[-subs, -subs])) / detR
	        result[term, 2] <- length(subs)
	    }
	    if (all(result[, 2] == 1)) result <- result[, 1]
	    else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
	    result
	}

		crit_column_to_sort = 6 ## this is used to sort and present the out_tab (ouput table) sorted by this criterion
    
    #########################add1.summary proccess#############################
    
    data <- eval(fit$call$data) #save data
    response <- fit$terms[[2]] #reponse variable
    reg_var <- variable.names(fit)[-1] #save all the dependent vars of 'fit' model
    
    if(is.data.frame(scope)){ #when scope is data.frame
        var <- colnames(scope)
        res_index <- which(var==response) #index of a response variable
        var <- var[-res_index] #save predictor variables by excluding a response variable
        if (length(reg_var!=0)){ #in the case of not null model
            for(i in 1:length(reg_var)){ #selecting predictors not included in the current model (candidate predictors)
                var_index <- which(var==reg_var[i])
                var <- var[-var_index]
            }
        }
    }else if(class(scope)=="formula"){ #when scope is formula
        if(all.vars(scope)[1]=="."){
            var <- all.vars(scope)[-1]
        }else if(all.vars(scope)[1] == response){
            var <- all.vars(scope)[c(-1,-2)]
        }else stop("\nwrong scope\n")
    }else stop("\nwrong scope\n")
    
    if(length(var) == 0) stop("\nwrong scope\n")
    
    out_tab <- matrix(nrow=length(var)+1, ncol=10) #outcome table
    colnames(out_tab) <- c("RSS", "AIC", "BIC", "adj.rsq", "PRESS", "max_pvalue", "max VIF", "alpha_cut_off", "Bonferroni", "FDR")
    rownames(out_tab) <- c("<none>", paste("+",var))
    out_tab[1,1] <- deviance(fit) #RSS of the current model
    out_tab[1,2] <- AIC(fit) #AIC of the current model
    out_tab[1,3] <- BIC(fit) #BIC of the current model
    out_tab[1,4] <- summary(fit)$"adj.r.squared" #adjusted r-square of the current model
    out_tab[1,5] <- PRESS(fit) #Press of the current model
    
    if (length(fit$coefficients) > 2) ## one of them is beta0, so need 3 or more coef
        out_tab[1,7] <- max(vif(fit)) 
    else
        out_tab[1,7] <- NA

    fit_pval <- drop1(fit, test="F")$"Pr(>F)"[-1] #save p-values of the current model
    if(length(fit_pval) == 0){ #for the null model
        out_tab[1,6] <- out_tab[1,8] <- out_tab[1,9] <-out_tab[1,10] <- NA
    }else{
        out_tab[1,6] <- max(fit_pval) #the max p-value of the current model
        out_tab[1,8] <- none(fit_pval) #pvalue cut-off (no correction)
        out_tab[1,9] <- bonferroni(fit_pval) #Bonfferoni cut-off
        out_tab[1,10] <- fdr(fit_pval) #FDR cut-off
    }
    
    reg <- paste(response,"~",paste(reg_var,collapse="+"))
    for (n in 1:length(var)){
        fit2 <- update(fit, paste(".~. +", var[n]))
        out_tab[n+1,1] <- deviance(fit2) #RSS of each model
        out_tab[n+1,2] <- AIC(fit2) #save value of AIC value for each model ##########
        out_tab[n+1,3] <- BIC(fit2) #save value of BIC value for each model ##########
        out_tab[n+1,4] <- summary(fit2)$"adj.r.squared" #save adjusted r-square
        out_tab[n+1,5] <- PRESS(fit2)
        
        if (length(fit2$coefficients) > 2) ## one of them is beta0, so need 3 or more coef
            out_tab[n+1,7] <- max(vif(fit2))
        else
            out_tab[n+1,7] <- NA

        fit2_pval <- drop1(fit2, test="F")$"Pr(>F)"[-1] #save p-values of each model
        out_tab[n+1,6] <- max(fit2_pval) #the max p-value of each model
        out_tab[n+1,8] <- none(fit2_pval) #pvalue cut-off (no correction)
        out_tab[n+1,9] <- bonferroni(fit2_pval) #Bonfferoni cut-off
        out_tab[n+1,10] <- fdr(fit2_pval) #FDR cut-off
    }
    
    out_tab = out_tab[order(as.numeric(out_tab[,crit_column_to_sort])),]
    out_tab <- round(as.table(out_tab, dimnames),5)
    out_tab[,8:10] <- as.character(as.logical(out_tab[,8:10]))
    return(out_tab)
}
