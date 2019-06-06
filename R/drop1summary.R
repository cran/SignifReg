drop1summary <- function(fit, scope=eval(fit$call$data), alpha = 0.05)UseMethod("drop1summary")
drop1summary.default <- function(fit, scope=eval(fit$call$data), alpha = 0.05){
    
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
    
    #########################drop1.summary proccess############################
    
    data <- eval(fit$call$data) #save data
    response <- fit$terms[[2]] #reponse variable
    reg_var <- variable.names(fit)[-1] #save all the dependent vars of 'fit' model
    
    if(is.data.frame(scope)){ #when scope is data.frame
        var <- colnames(scope)
        res_index <- which(var==response) #index of a response variable
        var <- var[-res_index] #save predictor variables by excluding a response variable
    }else if(class(scope)=="formula"){ #when scope is formula
        if(all.vars(scope)[1]=="."){
            var <- all.vars(scope)[-1]
        }else if(all.vars(scope)[1] == response){
            var <- all.vars(scope)[c(-1,-2)]
        }else stop("\nwrong scope\n")
    }else stop("\nwrong scope\n")
    
    out_tab <- matrix(nrow=length(var)+1, ncol=8) #outcome table
    colnames(out_tab) <- c("RSS", "AIC", "BIC", "adj.rsq","max_pvalue", "alpha_cut-off", "Bonferroni", "FDR")
    rownames(out_tab) <- c("<none>", paste("-",var))
    out_tab[1,1] <- deviance(fit) #RSS of the current model
    out_tab[1,2] <- AIC(fit) #AIC of the current model
    out_tab[1,3] <- BIC(fit) #BIC of the current model
    out_tab[1,4] <- summary(fit)$"adj.r.squared" #adjusted r-square of the current model
    
    fit_pval <- drop1(fit, test="F")$"Pr(>F)"[-1] #save p-values of the current model
    if(length(fit_pval)==0) stop("\nwrong scope\n")
    out_tab[1,5] <- max(fit_pval) #the max p-value of the current model
    out_tab[1,6] <- none(fit_pval) #pvalue cut-off (no correction)
    out_tab[1,7] <- bonferroni(fit_pval) #Bonfferoni cut-off
    out_tab[1,8] <- fdr(fit_pval) #FDR cut-off
    
    reg <- paste(response,"~",paste(reg_var,collapse="+"))
    for (n in 1:length(var)){
        #fit2 <- lm(paste(reg,"-",var[n]),data)
        fit2 <- update(fit, paste(".~. -", var[n]))
        out_tab[n+1,1] <- deviance(fit2) #RSS of each model
        out_tab[n+1,2] <- AIC(fit2) #save value of AIC value for each model ##########
        out_tab[n+1,3] <- BIC(fit2) #save value of BIC value for each model ##########
        out_tab[n+1,4] <- summary(fit2)$"adj.r.squared" #save adjusted r-square
        
        fit2_pval <- drop1(fit2, test="F")$"Pr(>F)"[-1] #save p-values of each model
        if(length(fit2_pval)==0)  stop("\nwrong scope\n")
        out_tab[n+1,5] <- max(fit2_pval) #the max p-value of each model
        out_tab[n+1,6] <- none(fit2_pval) #pvalue cut-off (no correction)
        out_tab[n+1,7] <- bonferroni(fit2_pval) #Bonfferoni cut-off
        out_tab[n+1,8] <- fdr(fit2_pval) #FDR cut-off
    }
    
    out_tab <- round(as.table(out_tab, dimnames),5)
    out_tab[,6:8] <- as.character(as.logical(out_tab[,6:8]))
    return(out_tab)
}
