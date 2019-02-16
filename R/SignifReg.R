SignifReg <-
function(scope,data,alpha = 0.05,direction = "forward",criterion = "p-value",correction = "FDR",trace=TRUE)UseMethod("SignifReg")
SignifReg.default <-
function(scope,data,alpha = 0.05,direction = "forward",criterion = "p-value",correction = "FDR",trace=TRUE){
    
    #####errors messages for wrong inputs#####
    if(criterion=="p-value"){
        if(direction=="step_null" || direction=="step_full"){
            stop("\ncriterion p-value cannot be used for step_null or step_full\n")
        }
    }
    if(direction != "forward" && direction != "backward" && direction != "step_null" && direction != "step_full"){
        stop("\ndirection should be one of the following: forward, backward, step_null, and step_full\n")
    }
    if(criterion != "p-value" && criterion != "AIC" && criterion != "BIC" && criterion != "r-adj" && criterion != "Cp"){
        stop("\ncriterion should be one of the following: p-value, AIC, BIC, r-adj, and Cp\n")
    }
    if(correction != "FDR" && correction != "Bonferroni" && correction !="Bonf" && correction != "None"){
        stop("\ncorrection should be one of the following: FDR, Bonferroni, Bonf, None F\n")
    }
    if(alpha>1 || alpha <0){
        stop("\nSignificance level alpha should be smaller than 1 and larger than 0\n")
    }
    if(length(scope)!=3){
        stop("\nscope should include response variable and predictor variables\n")
    }
    
    ###################################################################
    ###################################################################
    ########################Mallows' Cp function#######################
    ###################################################################
    ###################################################################
    Cp <- function(fit){
        response <- all.vars(scope)[1] #reponse variable
        var <- all.vars(scope)[-1] #predictor variables
        if(any(var==".")){
            var <- colnames(data)
            res_index <- which(var==response) #index of a response variable
            var<-var[-res_index] #save predictor variables by excluding a response variable
        }
        reg <- paste(var, collapse=" + ")
        reg <- paste(response,"~",reg)
        full_model <- lm(reg,data)
        MSE_full <- sum(summary(full_model)$"residuals"^2)/full_model$"df.residual"
        SSE_p <- sum(summary(fit)$"residuals"^2)
        n <- nrow(data) #number of obs
        p <- summary(fit)$"df"[1] #number of predictors + 1
        Cp <- SSE_p/MSE_full+2*p-n
        return(Cp)
    }
    
    ###################################################################
    ###################################################################
    #######################correction function#########################
    ###################################################################
    ###################################################################
    #correction function for Bonferroni, FDR, and default
    correction1 <- function(pvalues){
        max_p <- max(pvalues) #save maximum p-value
        if (correction=="None"){ #without any correction
            result <- ifelse(max_p < alpha, TRUE, FALSE)
        }else if (correction=="Bonf"|correction=="Bonferroni"){ #bonferroni test
            result <- ifelse(max_p < alpha/length(pvalues), TRUE, FALSE)#alpha/number of tests
        }else if(correction=="FDR"){ #computing FDR
            pi <- sort(pvalues) #sorting pvalues in ascending order
            j <- c(1:length(pi))
            t <- 0
            for(n in j){
                if(pi[n] <= (n/length(pvalues))*(alpha/sum(j^(-1)))) t <- 1+t
            }
            result <- ifelse (t == length(pi), TRUE, FALSE) #if true, that means every variable satisfies the FDR criterion.
        } else result <- NULL
        return(result) #when satisfy corrected criterion return result=TRUE, otherwise result=FALSE
    }
    
    ###################################################################
    ###################################################################
    ########################forward function###########################
    ###################################################################
    ###################################################################
    forward1 <- function(reg){
        fit <- lm(reg,data) #creating a regression model
        reg_var <- attr(fit$terms,"term.labels") #save all the dependent vars of 'fit' model
        var <- all.vars(scope)[-1] #predictor variable list which is wanted to check
        if(any(var==".")){ #in the case when user uses '.' in the scope command
            response <- all.vars(scope)[1]
            var <- colnames(data)
            res_index <- which(var==response) #index of a response variable
            var <- var[-res_index] #save predictor variables by excluding a response variable
        }
        if (length(reg_var!=0)){ #in the case of not null model
            for(i in 1:length(reg_var)){ #selecting predictors not included in the current model (candidate predictors)
                var_index <- which(var==reg_var[i])
                var <- var[-var_index]
            }
        }
        if(criterion == "AIC"){
            crit_value <- AIC(fit) #save AIC as criterion value #########
        } else if(criterion == "BIC") {
            crit_value <- BIC(fit) #save BIC as criterion value #########
        } else if(criterion == "Cp"){
            crit_value <- Cp(fit)
        } else if(criterion == "r-adj"){
            crit_value <- -summary(fit)$"adj.r.squared" #save -adjusted r-square
        }
        i <- c(1:length(var))
        x <- matrix()
        for (n in i){
            fit2 <- lm(paste(reg,"+",var[n]),data)
            if (criterion == "p-value"){
                x[n] <- max(drop1(fit2, test="F")$"Pr(>F)"[-1]) #save the biggest p-value of each model
            } else if (criterion == "AIC"){
                x[n] <- AIC(fit2) #save AIC value of each model ##########
            } else if (criterion == "BIC"){
                x[n] <- BIC(fit2) #save BIC value of each model ##########
            } else if (criterion == "Cp"){
                x[n] <- Cp(fit2)
            } else if (criterion == "r-adj"){
                x[n] <- -summary(fit2)$"adj.r.squared" #save -adjusted r-square
            }
        }
        
        x[which(is.na(x=="NA"))]  <- 99 #remove the model of which p-value is NA
        
        var_index <- which(x==min(x)) #select a model having the smallest p-value, AIC, BIC ###########
        if(length(var_index)>1) var_index <- var_index[1] #when p-value is the same(models for tie), select the first one
        if (length(reg_var) == 0){ # in the case of a null model
            response <- all.vars(scope)[1] #reponse variable
            reg2 <- paste(response,"~",var[var_index]) #model having one predictor
        }else reg2 <- paste(reg,"+",var[var_index]) #non-null model
        
        fit2 <- lm(reg2,data) # a model having the smallest p-value, AIC, BIC
        pvalues <- drop1(fit2,test="F")$"Pr(>F)"[-1] #p-values
        result <- correction1(pvalues) #save correction result as a switch to keep the while loop
        
        if (criterion == "p-value" & direction == "forward"){
            
            if (result==TRUE){
                switch <- TRUE #keep runing the while loop when correction is true
                reg <- reg2 #overide the model
            }else switch <- FALSE
            
            list1 <- list(reg,switch)
            
        } else if (criterion == "AIC" | criterion =="BIC" | criterion =="Cp" | criterion =="r-adj"){
            if(criterion == "AIC"){
                crit_value2 <- AIC(fit2) #save AIC #######
            } else if (criterion == "BIC"){
                crit_value2 <- BIC(fit2)
            } else if (criterion == "Cp"){
                crit_value2 <- Cp(fit2)
            } else if (criterion == "r-adj"){
                crit_value2 <- -summary(fit2)$"adj.r.squared" #save -adjusted r-square
            }
            #compare AIC, BIC, Mallows Cp, adjusted R-square with the previous model##############
            if(direction == "forward"){
                if (crit_value2 < crit_value & result == TRUE){ # in the case the new model has smaller criterion and also all significant predictors
                    reg <- reg2
                    switch <- TRUE #run the while loop
                } else switch <- FALSE #escape from the while loop
            }else if(direction == "step_null" | direction =="step_full"){
                reg <- reg2
                switch <- result #save correction1 result
            }
            list1 <- list(reg,switch)
        }
        return(list1)
    }
    
    ###################################################################
    ###################################################################
    #######################backward function###########################
    ###################################################################
    ###################################################################
    backward1 <- function(reg){
        fit <- lm(reg,data)
        var <- attr(fit$terms,"term.labels") #save all the independent variables in 'fit' model
        
        if(criterion =="p-value" & direction =="backward"){
            if (length(var)>1){ # if 'fit' model has more than one predictor
                pvalues <- drop1(fit,test="F")$"Pr(>F)"[-1] #p-values
                max_p <- max(pvalues) #maximum p-value
                var_index <- which(pvalues == max_p) #index number of variable having the biggest p-value
                var <- var[-var_index]
                reg2 <- paste(var,collapse = " + ") #to make a new model removing a variable having the biggest p-value
                reg2 <- paste(scope[2],scope[1],reg2) #to make a new model removing a variable having the biggest p-value
                fit2 <- lm(reg2,data)
                pvalues <- drop1(fit2,test="F")$"Pr(>F)"[-1] #p-values
                result <- correction1(pvalues) #correction
                
                if(result==FALSE){
                    switch <- TRUE #Only when a correction result is FALSE, run the while loop
                    reg <- reg2 #overide the model
                }else{
                    switch <- FALSE #stop when every predictor is significant (result = TRUE)
                    reg <- reg2 #overide the model
                }
            }else if (length(var)==1){
                pvalues <- drop1(fit,test="F")$"Pr(>F)"[-1] #p-values
                result <- correction1(pvalues)
                if (result == TRUE){ #when the model has only one predictor
                    switch <- FALSE # if only predictor is significant, stop the while loop
                }else{
                    response <- all.vars(scope)[1] #response variable
                    reg <- paste(response,"~ 1") #make a null model
                    switch <- FALSE # if only predictor is insignificant, return a null model and stop the loop
                }
            }else switch <- FALSE #when no predictor is left, stop the while loop
            
        } else if(criterion == "AIC" | criterion =="BIC" | criterion == "Cp" | criterion == "r-adj"){
            if(criterion == "AIC"){
                crit_value <- AIC(fit) #save AIC #######
            } else if (criterion == "BIC"){
                crit_value <- BIC(fit) #save BIC
            } else if (criterion == "Cp"){
                crit_value <- Cp(fit) #save Mallow's Cp
            } else if (criterion == "r-adj"){
                crit_value <- -summary(fit)$"adj.r.squared" #save -adjusted r-square
            }
            crit_matrix <- matrix()
            reg_matrix <- matrix()
            
            for(i in 1:length(var)){
                var2 <- var[-i]
                reg2 <- var2[1]
                
                if (length(var2) > 1){
                    for(j in 2:length(var2)){
                        reg2 <- paste(reg2, "+",var2[j])
                    }
                } else if (length(var2)==1){
                    reg2 <- var2
                } else reg2 <- 1
                
                reg2 <- paste(scope[2],scope[1],reg2)
                fit2 <- lm(reg2,data)
                reg_matrix[i] <- reg2
                if(criterion=="AIC"){
                    crit_matrix[i] <- AIC(fit2)
                }else if (criterion =="BIC"){
                    crit_matrix[i] <- BIC(fit2)
                }else if (criterion =="Cp"){
                    crit_matrix[i] <- Cp(fit2)
                }else if (criterion =="r-adj"){
                    crit_matrix[i] <- -summary(fit2)$"adj.r.squared" #save -adjusted r-square
                }
            }
            crit_value2 <- min(crit_matrix)
            var_index <- which(crit_matrix==min(crit_matrix)) #variable num making the minimum AIC when dropped
            reg2 <- reg_matrix[var_index] #regression equation droping the variable # save reg <- reg2
            fit2 <- lm(reg2, data)
            pvalues <- drop1(fit2,test="F")$"Pr(>F)"[-1] #p-values
            if (length(pvalues) != 0){
                result <- correction1(pvalues) #save correction result as a switch to keep the while loop
            } else result <- TRUE #result is true, when a null model is.
            #compare criterion with the previous model##############
            if(direction == "backward"){
                if (crit_value2 < crit_value | result == FALSE){ #in the case a new model has smaller criterion and all siginificant predictors
                    reg <- reg2
                    switch <- TRUE #run the while loop
                }else{
                    reg <- reg2
                    switch <- FALSE #escape the while loop when there is no variable decreasing AIC when dropped
                }
            }else if(direction == "step_null" | direction == "step_full"){
                reg <- reg2
                switch <- result #send correction1 result as switch
            }
        }
        list1 <- list(reg,switch)
        return(list1)
    }
    ###################################################################
    ###################################################################
    #########################stepwise function#########################
    ###################################################################
    ###################################################################
    stepwise <- function(reg){
        list_forward <- forward1(reg) #forward model
        reg_forward <- as.character(list_forward[1])
        switch_forward <- unlist(list_forward[2])
        fit_forward <- lm(reg_forward,data)
        
        list_backward <- backward1(reg) #backward model
        reg_backward <- as.character(list_backward[1])
        switch_backward <- unlist(list_backward[2])
        fit_backward <- lm(reg_backward,data)
        
        fit_current <- lm(reg,data) #current model
        reg_current <- reg
        pvalues <- drop1(fit_current,test="F")$"Pr(>F)"[-1] #p-values
        switch_current <- correction1(pvalues)
        
        if (criterion=="AIC") {
            crit_forward <- AIC(fit_forward)
            crit_backward <- AIC(fit_backward)
            crit_current <- AIC(fit_current)
        } else if (criterion == "BIC"){
            crit_forward <- BIC(fit_forward)
            crit_backward <- BIC(fit_backward)
            crit_current <- BIC(fit_current)
        } else if (criterion == "Cp"){
            crit_forward <- Cp(fit_forward)
            crit_backward <- Cp(fit_backward)
            crit_current <- Cp(fit_current)
        } else if (criterion == "r-adj"){
            crit_forward <- -summary(fit_forward)$"adj.r.squared" #save -adjusted r-square
            crit_backward <- -summary(fit_backward)$"adj.r.squared" #save -adjusted r-square
            crit_current <- -summary(fit_current)$"adj.r.squared" #save -adjusted r-square
        }
        models <- data.frame(model=c("current","forward","backward"),crit_value=c(crit_current,crit_forward,crit_backward),result=c(switch_current, switch_forward, switch_backward))
        switch_index <- which(models$result==TRUE) #model number having value 'TRUE' for correction1
        if(length(switch_index)==0){ #when every model has values, 'FALSE', for correction1 results
            #sig_models <- models #because every model includes some insignificant predictors after checking it through correction1
            selected_model <- models[3,]
            reg <- reg_backward #select a backward model when every model includes some insignificant predictors
            switch <- TRUE
        }else{
            sig_models <- models[switch_index,] #save models having only significant predictors
            
            crit_index <- which(sig_models$crit_value==min(sig_models$crit_value)) #model having minimum criterion
            selected_model <- sig_models[crit_index,] #model having minimum criterion
            
            if(selected_model$model == "current"){
                reg <- reg_current
                switch <- FALSE
            }else if (selected_model$model == "forward"){
                reg <- reg_forward
                switch <- TRUE
            }else if (selected_model$model == "backward"){
                reg <- reg_backward
                switch <- TRUE
            }
        }
        list1 <- list(reg,switch)
        
        return (list1)
    }
    
    ##################################################
    ##################################################
    #############forward process######################
    ##################################################
    if (direction == "forward"){
        response <- all.vars(scope)[1] #reponse variable
        var <- all.vars(scope)[-1] #predictor variables
        reg <- paste(response,"~","1") #null model
        switch <- TRUE #turn on the while loop
        
        while (switch==TRUE){
            if (trace == TRUE){
                cat("\n")
                print(summary(lm(reg,data))$coeff)
            }
            list1 <- forward1(reg)
            reg <- as.character(list1[1])
            switch <- unlist(list1[2])
        }
        fit <- lm(reg, data)
        print(summary(fit))
    }
    
    ##################################################
    ##################################################
    ############backward process######################
    ##################################################
    else if (direction == "backward"){
        reg <- scope
        #test p-value to check the case every variables are significant from beginning.
        fit <- lm(reg,data) #start from a full model
        pvalues <- drop1(fit,test="F")$"Pr(>F)"[-1] #p-values
        result <- correction1(pvalues) #result of a test by corrected criterion
        switch <- ifelse(result==TRUE, FALSE, TRUE) #run the while loop only when the result of correction is not satisfied
        
        while (switch == TRUE){ #for backward the while loop works if the criterion is not satisfied (switch=False)
            if (trace == TRUE){
                cat("\n")
                print(summary(lm(reg,data))$coeff)
            }
            list1 <- backward1(reg)
            reg <- as.character(list1[1])
            switch <- unlist(list1[2])
        }
        fit <- lm(reg,data)
        print(summary(fit))
    }
    
    ##################################################
    ##################################################
    #############step-null process####################
    ##################################################
    else if (direction == "step_null"){
        response <- all.vars(scope)[1] #reponse variable
        var <- all.vars(scope)[-1] #predictor variables
        reg <- paste(response,"~","1") #null model
        if (trace == TRUE){
            cat("\n")
            print(summary(lm(reg,data))$coeff)
        }
        list1 <- forward1(reg) #add the first variable
        switch <- unlist(list1[2]) #switch to decide whether running the while loop
        if (switch == TRUE){ #if switch is not true after the 1st foward, choose a null model
            reg <- as.character(list1[1])
        }
        while (switch==TRUE){
            if (trace == TRUE){
                cat("\n")
                print(summary(lm(reg,data))$coeff)
            }
            list1 <- stepwise(reg)
            reg <- as.character(list1[1])
            switch <- unlist(list1[2])
        }
        reg <- as.character(list1[1])
        fit <- lm(reg,data)
        print(summary(fit))
    }
    ##################################################
    ##################################################
    #############step-full process####################
    ##################################################
    else if(direction=="step_full"){
        reg <- scope
        fit <- lm(reg,data) #start from a full model
        pvalues <- drop1(fit,test="F")$"Pr(>F)"[-1] #p-values
        result <- correction1(pvalues) #result of a test by corrected criterion
        switch <- ifelse(result==TRUE, FALSE, TRUE) #run the while loop only when the result of correction is not satisfied
        
        if (switch == TRUE){
            list1 <- backward1(reg)
            reg <- as.character(list1[1])
            result <- unlist(list1[2]) #save correction1 result
            switch <- ifelse(result==TRUE,FALSE,TRUE) #run the while loop if correction1 == FALSE as a result of backward
        }
        
        while(switch == TRUE){
            if (trace == TRUE){
                cat("\n")
                print(summary(lm(reg,data))$coeff)
            }
            list1 <- stepwise(reg)
            reg <- as.character(list1[1])
            switch <- unlist(list1[2])
        }
        reg <- as.character(list1[1])
        fit <- lm(reg,data)
        print(summary(fit))
    }
    return(fit)
    }
