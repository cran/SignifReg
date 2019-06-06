drop1SignifReg <-
function(fit, scope=eval(fit$call$data), alpha = 0.05, criterion = "p-value", correction = "FDR", override = FALSE)UseMethod("drop1SignifReg")

drop1SignifReg.default <- function(fit, scope=eval(fit$call$data), alpha = 0.05, criterion = "p-value", correction = "FDR", override = FALSE){
    
    #####errors messages for wrong inputs#####
    if(criterion != "p-value" && criterion != "AIC" && criterion != "BIC" && criterion != "r-adj"){
        stop("\ncriterion should be one of the following: p-value, AIC, BIC, and r-adj\n")
    }
    if(correction != "FDR" && correction != "Bonferroni" && correction !="Bonf" && correction != "None"){
        stop("\ncorrection should be one of the following: FDR, Bonferroni, Bonf, None\n")
    }
    if(alpha>1 || alpha <0){
        stop("\nSignificance level alpha should be smaller than 1 and larger than 0\n")
    }
    if(override != FALSE && override != TRUE){
        stop("\noverride should be either TRUE or FALSE\n")
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
    
    
    #########################drop1 proccess############################
    data <- eval(fit$call$data) #save data
    response <- fit$terms[[2]] #reponse variable
    reg_var <- variable.names(fit)[-1] #save all the dependent vars of 'fit' model
    
    if(is.data.frame(scope)){ #when scope is data.frame
        var <- colnames(scope)
        res_index <- which(var==response) #index of a response variable
        var <- var[-res_index] #save predictor variables by excluding a response variable
        var <- var[var %in% reg_var] #remove candidate predictor(s) variable which is not in the fit model
    }else if(class(scope)=="formula"){ #when scope is formula
        if(all.vars(scope)[1]=="."){
            var <- all.vars(scope)[-1]
            var <- var[var %in% reg_var] #remove candidate predictor(s) variable which is not in the fit model
            
        }else if(all.vars(scope)[1] == response){
            var <- all.vars(scope)[c(-1,-2)]
            var <- var[var %in% reg_var] #remove candidate predictor(s) variable which is not in the fit model
        }else stop("\nwrong scope\n")
    }else stop("\nwrong scope\n")
    
    out_tab <- matrix(nrow=length(var)+1, ncol=8) #outcome table
    colnames(out_tab) <- c("RSS", "AIC", "BIC", "adj.rsq","max_pvalue", "alpha_cut-off", "Bonferroni", "FDR")
    rownames(out_tab) <- c("<none>", paste("-", var))
    out_tab[1,1] <- deviance(fit) #RSS of the current model
    out_tab[1,2] <- AIC(fit) #AIC of the current model
    out_tab[1,3] <- BIC(fit) #BIC of the current model
    out_tab[1,4] <- summary(fit)$"adj.r.squared" #adjusted r-square of the current model
    
    fit_pval <- drop1(fit, test="F")$"Pr(>F)"[-1] #save p-values of the current model
    out_tab[1,5] <- max(fit_pval) #the max p-value of the current model
    out_tab[1,6] <- none(fit_pval) #pvalue cut-off (no correction)
    out_tab[1,7] <- bonferroni(fit_pval) #Bonfferoni cut-off
    out_tab[1,8] <- fdr(fit_pval) #FDR cut-off
    
    crit_value <- matrix() #criterion values
    if(criterion == "AIC"){
        crit_value[1] <-  out_tab[1,2] #save AIC as criterion value #########
    } else if(criterion == "BIC") {
        crit_value[1] <- out_tab[1,3] #save BIC as criterion value #########
    } else if(criterion == "r-adj"){
        crit_value[1] <- -out_tab[1,4] #save -adjusted r-square
    }
    
    #reg <- paste(response,"~",paste(reg_var,collapse="+"))
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
        
        if (criterion == "p-value"){
            crit_value[n+1] <- out_tab[n+1,5] #save the biggest p-value of each model
        } else if (criterion == "AIC"){
            crit_value[n+1] <- out_tab[n+1,2] #save AIC value of each model ##########
        } else if (criterion == "BIC"){
            crit_value[n+1] <- out_tab[n+1,3] #save BIC value of each model ##########
        } else if (criterion == "r-adj"){
            crit_value[n+1] <- -out_tab[n+1,4] #save -adjusted r-square
        }
    }
    
    out_tab <- round(as.table(out_tab, dimnames),5)
    out_tab[,6:8] <- as.character(as.logical(out_tab[,6:8]))
    print(out_tab)
    
    crit_value[which(is.na(crit_value=="NA"))]  <- 99 #remove the model of which p-value is NA
    ########################
    var_index <- which(crit_value==min(crit_value[-1]))-1
    if(length(var_index)>1) var_index <- var_index[1] #when p-value is the same(models for tie), select the first one
    
    fit2 <- update(fit, paste(".~. -", var[var_index]))
    pvalues <- drop1(fit2,test="F")$"Pr(>F)"[-1] #p-values
    
    if(override == FALSE){
        if(correction=="FDR"){
            result <- fdr(pvalues)
        }else if(correction=="Bonfferoni" | correction=="Bonf"){
            result <- bonferroni(pvalues) #Bonfferoni cut-off
        }else if(correction=="None"){
            result <- none(pvalues)
        }
        if(criterion == "p-value"){
            if(result == FALSE){
                fit <- fit2
            }else if(result == TRUE){
                
                fit <- fit
            }
        }else if(criterion == "AIC" | criterion =="BIC" | criterion =="r-adj"){
            if(result == TRUE & min(crit_value[-1]) > crit_value[1]){
                fit <- fit
            }else fit <- fit2
        }
    }else if(override == TRUE){ #do not check correction)
        if(criterion == "p-value"){
            fit <- fit2
        }else if(criterion == "AIC" | criterion =="BIC" | criterion =="r-adj"){
            if(min(crit_value[-1]) > crit_value[1]){
                fit <- fit
            }else fit <- fit2
        }
    }
    return(fit)
}
