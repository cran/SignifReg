SignifReg <-
function(fit,scope=eval(fit$call$data), alpha = 0.05,direction = "forward",criterion = "p-value",correction = "FDR",trace=FALSE)UseMethod("SignifReg")
SignifReg.default <-
function(fit,scope=eval(fit$call$data), alpha = 0.05,direction = "forward",criterion = "p-value",correction = "FDR",trace=FALSE){
    
    #####errors messages for wrong inputs#####
    
    #if(criterion=="p-value"){
    #  if(direction=="both"){
    #    stop("\ncriterion p-value cannot be used for both\n")
    #  }
    #}
    
    if(direction != "forward" && direction != "backward" && direction != "both"){
        stop("\ndirection should be one of the following: forward, backward, and both\n")
    }
    if(criterion != "p-value" && criterion != "AIC" && criterion != "BIC" && criterion != "r-adj"){
        stop("\ncriterion should be one of the following: p-value, AIC, BIC, and r-adj\n")
    }
    if(correction != "FDR" && correction != "Bonferroni" && correction !="Bonf" && correction != "None"){
        stop("\ncorrection should be one of the following: FDR, Bonferroni, Bonf, None \n")
    }
    if(alpha>1 || alpha <0){
        stop("\nSignificance level alpha should be smaller than 1 and larger than 0\n")
    }
    if(trace != FALSE && trace != TRUE){
        stop("\ntrace should be either TRUE or FALSE\n")
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
    
    #########################forward1 function#############################
    
    forward1 <- function(fit,override=FALSE){
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
            
            if(length(reg_var) != 0){
                var <- var[-which(var %in% reg_var)] #remove candidate predictor(s) variable which already in the fit model
            }
        }else stop("\nwrong scope\n")
        
        if(length(var)==0){ #no more candidate variable left. stop the while loop
            #fit <- fit
            switch <- FALSE
            out_tab <- NA
            
        }else{
            out_tab <- matrix(nrow=length(var)+1, ncol=8) #outcome table
            colnames(out_tab) <- c("RSS", "AIC", "BIC", "adj.rsq","max_pvalue", "alpha_cut_off", "Bonferroni", "FDR")
            rownames(out_tab) <- c("<none>", paste("+", var))
            out_tab[1,1] <- deviance(fit) #RSS of the current model
            out_tab[1,2] <- AIC(fit) #AIC of the current model
            out_tab[1,3] <- BIC(fit) #BIC of the current model
            out_tab[1,4] <- summary(fit)$"adj.r.squared" #adjusted r-square of the current model
            
            fit_pval <- drop1(fit, test="F")$"Pr(>F)"[-1] #save p-values of the current model
            if(length(fit_pval) == 0){ #for the null model
                out_tab[1,5] <- out_tab[1,6] <- out_tab[1,7] <-out_tab[1,8] <- NA
            }else{
                out_tab[1,5] <- max(fit_pval) #the max p-value of the current model
                out_tab[1,6] <- none(fit_pval) #pvalue cut-off (no correction)
                out_tab[1,7] <- bonferroni(fit_pval) #Bonferroni cut-off
                out_tab[1,8] <- fdr(fit_pval) #FDR cut-off
            }
            crit_value <- matrix() #criterion values
            if(criterion == "p-value"){
                crit_value[1] <- 999 #don't compare with the current model if p-value is a criterion (exclude the current model in the candidate models)
            }else if(criterion == "AIC"){
                crit_value[1] <-  out_tab[1,2] #save AIC as criterion value #########
            } else if(criterion == "BIC") {
                crit_value[1] <- out_tab[1,3] #save BIC as criterion value #########
            } else if(criterion == "r-adj"){
                crit_value[1] <- -out_tab[1,4] #save -adjusted r-square
            }
            
            for (n in 1:length(var)){
                #fit2 <- lm(paste(reg,"+",var[n]),data)
                fit2 <- update(fit, paste(".~. +", var[n]))
                out_tab[n+1,1] <- deviance(fit2) #RSS of each model
                out_tab[n+1,2] <- AIC(fit2) #save value of AIC value for each model ##########
                out_tab[n+1,3] <- BIC(fit2) #save value of BIC value for each model ##########
                out_tab[n+1,4] <- summary(fit2)$"adj.r.squared" #save adjusted r-square
                
                fit2_pval <- drop1(fit2, test="F")$"Pr(>F)"[-1] #save p-values of each model
                out_tab[n+1,5] <- max(fit2_pval) #the max p-value of each model
                out_tab[n+1,6] <- none(fit2_pval) #pvalue cut-off (no correction)
                out_tab[n+1,7] <- bonferroni(fit2_pval) #Bonferroni cut-off
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
            #if(trace==TRUE){cat("\n"); print(out_tab)}
            
            crit_value[which(is.na(crit_value=="NA"))]  <- 99 #remove the model of which p-value is NA
            ##################################
            
            var_index <- which(crit_value==min(crit_value[-1]))-1
            if(length(var_index)>1) var_index <- var_index[1] #when p-value is the same(models for tie), select the first one
            
            fit2 <- update(fit, paste(".~. +", var[var_index]))
            pvalues <- drop1(fit2,test="F")$"Pr(>F)"[-1] #p-values
            
            if(override == FALSE){
                if(correction=="FDR"){
                    result <- fdr(pvalues)
                }else if(correction=="Bonferroni" | correction=="Bonf"){
                    result <- bonferroni(pvalues) #Bonferroni cut-off
                }else if(correction=="None"){
                    result <- none(pvalues)
                }
                if(criterion == "p-value"){
                    if(direction == "forward"){
                        if(result == TRUE){
                            fit <- fit2
                            switch <- TRUE
                        }else if(result == FALSE){
                            #fit <- fit
                            switch <- FALSE
                        }
                    }else if(direction == "both"){
                        fit <- fit2
                        switch <- result #save correction1 result
                    }
                }else if(criterion == "AIC" | criterion =="BIC" | criterion =="r-adj"){
                    if(direction == "forward"){
                        if(result == TRUE & min(crit_value[-1]) < crit_value[1]){
                            fit <- fit2
                            switch <- TRUE
                        }else{
                            #fit <- fit
                            switch <- FALSE
                        }
                    }else if(direction == "both"){
                        fit <- fit2
                        switch <- result #save correction1 result
                    }
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
        }
        
        list1 <- list(fit, switch, out_tab)
        return(list1)
    }
    
    ######################################################################
    #########################backward1 function############################
    backward1 <- function(fit, override=FALSE){
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
        
        if(length(var) == 0){ #no left variable in scope
            #fit <- fit
            switch <- FALSE
            out_tab <- NA
        }else{
            if(length(reg_var) !=0){
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
                out_tab[1,7] <- bonferroni(fit_pval) #Bonferroni cut-off
                out_tab[1,8] <- fdr(fit_pval) #FDR cut-off
                
                crit_value <- matrix() #criterion values
                if(criterion == "p-value"){
                    crit_value[1] <- 999 #do not compare with the current model if criterion is p-value
                } else if(criterion == "AIC"){
                    crit_value[1] <-  out_tab[1,2] #save AIC as criterion value #########
                } else if(criterion == "BIC") {
                    crit_value[1] <- out_tab[1,3] #save BIC as criterion value #########
                } else if(criterion == "r-adj"){
                    crit_value[1] <- -out_tab[1,4] #save -adjusted r-square
                }
                
                for (n in 1:length(var)){
                    fit2 <- update(fit, paste(".~. -", var[n]))
                    
                    out_tab[n+1,1] <- deviance(fit2) #RSS of each model
                    out_tab[n+1,2] <- AIC(fit2) #save value of AIC value for each model ##########
                    out_tab[n+1,3] <- BIC(fit2) #save value of BIC value for each model ##########
                    out_tab[n+1,4] <- summary(fit2)$"adj.r.squared" #save adjusted r-square
                    
                    fit2_pval <- drop1(fit2, test="F")$"Pr(>F)"[-1] #save p-values of each model
                    #if(length(fit2_pval)==0)  stop("\nwrong scope or outcome is a null model\n")
                    if(length(fit2_pval)==0){ #if fit2 is a null model
                        out_tab[n+1,5] <- NA
                        out_tab[n+1,6] <- NA
                        out_tab[n+1,7] <- NA
                        out_tab[n+1,8] <- NA
                    }else{
                        out_tab[n+1,5] <- max(fit2_pval) #the max p-value of each model
                        out_tab[n+1,6] <- none(fit2_pval) #pvalue cut-off (no correction)
                        out_tab[n+1,7] <- bonferroni(fit2_pval) #Bonferroni cut-off
                        out_tab[n+1,8] <- fdr(fit2_pval) #FDR cut-off
                    }
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
                #if(trace==TRUE){cat("\n"); print(out_tab)}
                
                crit_value[which(is.na(crit_value=="NA"))]  <- 99 #remove the model of which p-value is NA
                ########################
                var_index <- which(crit_value==min(crit_value[-1]))-1
                #if(length(var_index)>1) var_index <- var_index[1] #when p-value is the same(models for tie), select the first one
                
                if(var_index == 0){ #if the current model has one predictor and criterion is p-value
                    fit2 <- update(fit, paste(".~. -", var[var_index+1])) #make a null model
                }else fit2 <- update(fit, paste(".~. -", var[var_index]))
                pvalues <- drop1(fit2,test="F")$"Pr(>F)"[-1] #p-values
                
                if(override == FALSE){
                    if(length(pvalues)==0){ #if fit2 is a null model, donot check the correction
                        result <- TRUE
                    }else{ #if fit2 is not a null model, check the correction
                        if(correction=="FDR"){
                            result <- fdr(pvalues)
                        }else if(correction=="Bonferroni" | correction=="Bonf"){
                            result <- bonferroni(pvalues) #Bonferroni cut-off
                        }else if(correction=="None"){
                            result <- none(pvalues)
                        }
                    }
                    if(criterion == "p-value"){
                        if (direction == "backward"){
                            if (length(reg_var)>1){ # if 'fit' model has more than one predictor
                                if(result == FALSE){
                                    fit <- fit2
                                    switch <- TRUE #when a correction result is FALSE, run the while loop
                                }else if(result == TRUE){
                                    
                                    if(correction=="FDR"){ #check the current model is significant with correction
                                        if(out_tab[1,8] == TRUE){ #if the current model is already significant even at the beginning
                                            switch <- FALSE
                                        }else{fit<-fit2; switch <- TRUE}  # if the current model is not significant, but fit2 is significant
                                    }else if(correction=="Bonferroni" | correction=="Bonf"){
                                        if(out_tab[1,7] == TRUE){
                                            switch <- FALSE
                                        }else{fit<-fit2; switch <- TRUE} # if the current model is not significant, but fit2 is significant
                                    }else if(correction=="None"){ #if the current model is already significant
                                        if(out_tab[1,6] == FALSE){
                                            switch <- FALSE
                                        }else{fit<-fit2; switch <- TRUE}
                                    }
                                    
                                }
                            }else if (length(reg_var)==1){ #fit model has only one variable
                                if (result == TRUE){ #when the model has only one predictor
                                    if(correction=="FDR"){ #check the only predictor is significant with correction
                                        if(out_tab[1,8] == FALSE){
                                            fit <- fit2; switch <- FALSE #########################
                                        }else switch <- FALSE # if only predictor is significant, stop the while loop
                                    }else if(correction=="Bonferroni" | correction=="Bonf"){
                                        if(out_tab[1,7] == FALSE){
                                            fit <- fit2; switch <- FALSE
                                        }else switch <- FALSE # if only predictor is significant, stop the while loop
                                    }else if(correction=="None"){
                                        if(out_tab[1,6] == FALSE){
                                            fit <- fit2; switch <- FALSE
                                        }else switch <- FALSE # if only predictor is significant, stop the while loop
                                    }
                                    #switch <- FALSE # if only predictor is significant, stop the while loop
                                }else if(length(var == 1)) { #when result is false and scope has one variable, make a null model
                                    fit2 <- update(fit, paste(".~. -", var[n]))
                                    switch <- FALSE
                                    #stop("\noutcome is a null model\n")
                                }else{ #should keep a variable not in the scope even if it is not significant
                                    fit2 <- fit
                                    switch <- FALSE
                                }
                            }else switch <- FALSE #when no predictor is left, stop the while loop
                        }else if(direction == "both"){
                            fit <- fit2
                            switch <- result #send correction1 result as switch
                        }
                    }else if(criterion == "AIC" | criterion =="BIC" | criterion =="r-adj"){
                        if(direction=="backward"){
                            if(result == FALSE | min(crit_value[-1]) < crit_value[1]){
                                fit <- fit2
                                switch <- TRUE #run the while loop
                            }else if(length(variable.names(fit2)[-1])==0){ #if fit2 is a null model and fit is better than fit2
                                #fit <- fit
                                switch <- FALSE
                            }else{
                                fit <- fit2
                                switch <- FALSE #escape the while loop when there is no variable decreasing AIC when dropped
                            }
                        }else if(direction == "both"){
                            fit <- fit2
                            switch <- result #send correction1 result as switch
                        }
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
                
            }else{ #if length(reg_var) == 0, no variable in the model
                #fit <- fit
                switch <- FALSE
            }
        }
        list1 <- list(fit, switch, out_tab)
        return(list1)
    }
    
    #########################stepwise function#########################
    
    stepwise <- function(fit){
        list_forward <- forward1(fit) #forward model
        fit_forward <- list_forward[[1]]
        switch_forward <- list_forward[[2]]
        out_tab_forward <- list_forward[[3]]
        
        list_backward <- backward1(fit) #backward model
        fit_backward <- list_backward[[1]]
        switch_backward <- list_backward[[2]]
        out_tab_backward <- list_backward[[3]]
        
        fit_current <- fit #current model
        pvalues <- drop1(fit_current,test="F")$"Pr(>F)"[-1] #p-values
        if(length(pvalues) == 0){ #if the current model is a null model
            switch_current <- TRUE
        }else{
            if(correction=="FDR"){
                switch_current <- fdr(pvalues)
            }else if(correction=="Bonferroni" | correction=="Bonf"){
                switch_current <- bonferroni(pvalues) #Bonferroni cut-off
            }else if(correction=="None"){
                switch_current <- none(pvalues)
            }
        }
        if (criterion=="AIC") {
            crit_forward <- AIC(fit_forward)
            crit_backward <- AIC(fit_backward)
            crit_current <- AIC(fit_current)
        } else if (criterion == "BIC"){
            crit_forward <- BIC(fit_forward)
            crit_backward <- BIC(fit_backward)
            crit_current <- BIC(fit_current)
        } else if (criterion == "r-adj"){
            crit_forward <- -summary(fit_forward)$"adj.r.squared" #save -adjusted r-square
            crit_backward <- -summary(fit_backward)$"adj.r.squared" #save -adjusted r-square
            crit_current <- -summary(fit_current)$"adj.r.squared" #save -adjusted r-square
        } else if (criterion == "p-value"){
            crit_forward <- max(drop1(fit_forward, test="F")$"Pr(>F)"[-1]) #save the max p-values of each model
            
            if(length(drop1(fit_backward, test="F")$"Pr(>F)"[-1])==0){
                crit_backward <- 99 #if backward model is a null model, don't choose it
            }else crit_backward <- max(drop1(fit_backward, test="F")$"Pr(>F)"[-1]) #in the case that backward model is not a null model
            
            if(length(drop1(fit_current, test="F")$"Pr(>F)"[-1])==0){
                crit_current <- 99 #if the current model is a null model, don't choose it
            }else crit_current <- max(drop1(fit_current, test="F")$"Pr(>F)"[-1]) #in the case that current model is not a null model
            
        }
        models <- data.frame(model=c("current","forward","backward"),crit_value=c(crit_current,crit_forward,crit_backward),result=c(switch_current, switch_forward, switch_backward))
        switch_index <- which(models$result==TRUE) #model number having value 'TRUE' for correction1
        if(length(switch_index)==0 & models[1,2]==models[2,2] & models[1,2]==models[3,2] & models[2,2]==models[3,2]){ #when every model has values, 'FALSE', for correction results and identical critical values
            #sig_models <- models #because every model includes some insignificant predictors after checking it through correction1
            selected_model <- models[3,]
            fit <- fit_backward #select a backward model when every model includes some insignificant predictors
            switch <- FALSE
        }else if(length(switch_index)==0){#when every model has values, 'FALSE', for correction results
            #sig_models <- models #because every model includes some insignificant predictors after checking it through correction1
            selected_model <- models[3,]
            fit <- fit_backward #select a backward model when every model includes some insignificant predictors
            switch <- TRUE
        }else{
            sig_models <- models[switch_index,] #save models having only significant predictors
            crit_index <- which(sig_models$crit_value==min(sig_models$crit_value)) #model having minimum criterion
            if(length(crit_index) !=1) crit_index <- crit_index[1] #if there is the same critical values among forward, backward, and current model, choose the current model
            selected_model <- sig_models[crit_index,] #model having minimum criterion
            
            if(selected_model$model == "current"){
                fit <- fit_current
                switch <- FALSE
            }else if (selected_model$model == "forward"){
                fit <- fit_forward
                switch <- TRUE
            }else if (selected_model$model == "backward"){
                fit <- fit_backward
                switch <- TRUE
            }
        }
        list1 <- list(fit, switch, out_tab_forward, out_tab_backward)
        
        return (list1)
    }
    
    #################################forward process######################################
    if (direction == "forward"){
        if(trace==TRUE) print(fit)
        switch <- TRUE #turn on the while loop
        
        while (switch==TRUE){
            list1 <- forward1(fit)
            fit <- list1[[1]]
            switch <- list1[[2]]
            out_tab <- list1[[3]]
            if(trace==TRUE){cat("\n"); print(out_tab); print(fit); cat("\n\n")}
        }
    }
    #################################backward process######################################
    else if(direction == "backward"){
        if(trace==TRUE) print(fit)
        #test to check the case the current model is the best from beginning.
        list1 <- backward1(fit)
        fit <- list1[[1]]
        switch <- list1[[2]]
        out_tab <- list1[[3]]
        if(trace==TRUE){cat("\n"); print(out_tab); print(fit); cat("\n\n")}
        
        while (switch == TRUE){ #for backward the while loop works if the criterion is not satisfied (switch=False)
            list1 <- backward1(fit)
            fit <- list1[[1]]
            switch <- list1[[2]]
            out_tab <- list1[[3]]
            if(trace==TRUE & switch == TRUE){cat("\n"); print(out_tab); print(fit); cat("\n\n")}
        }
    }
    #############stepwise process####################
    else if(direction == "both"){
        if(trace==TRUE) print(fit)
        reg_var <- variable.names(fit)[-1] #save all the dependent vars of 'fit' model
        switch <- TRUE
        if(length(reg_var) == 0){ #if the first model is a null model
            list1 <- forward1(fit)
            fit <- list1[[1]]
            switch <- list1[[2]]
            out_tab <- list1[[3]]
            if(trace==TRUE){cat("\n"); print(out_tab); print(fit); cat("\n\n")}
        }
        list1 <- stepwise(fit) #when the current model is the best from the beginning
        fit <- list1[[1]]
        switch <- list1[[2]]
        forward_out_tab <- list1[[3]]
        backward_out_tab <- list1[[4]]
        if(trace==TRUE){cat("\n"); if(is.table(forward_out_tab)==TRUE){print(forward_out_tab)}; print(backward_out_tab); print(fit); cat("\n\n")}
        
        while (switch==TRUE){
            list1 <- stepwise(fit)
            fit <- list1[[1]]
            switch <- list1[[2]]
            forward_out_tab <- list1[[3]]
            backward_out_tab <- list1[[4]]
            if(trace==TRUE){cat("\n"); print(forward_out_tab); print(backward_out_tab); print(fit); cat("\n\n")}
        }
        fit <- list1[[1]]
    }
    return(fit)
}
