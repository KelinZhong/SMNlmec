

#' @title Get a summary of results from SMNlmec.est.
#' @import rstan
#' @param object An object of class \code{SMNlmecfit}.
#' @examples 
#' See the example in SMNlmec.est.
#' @exportMethod SMNlmec.summary
setMethod("SMNlmec.summary","SMNlmecfit",function(object){
  temp_dist <- object@dist_set
  temp_struc <- object@struc_set
  temp_model <- object@stan_object
  temp_criteria <- object@model_criteria
  
  if(temp_dist == "Normal") {
    
    if(temp_struc == "UNC") {
      print(temp_model,par=c("beta","sigma2","sigmab2"), digits = 3)
      print(temp_criteria)
    }
    
    if(temp_struc == "DEC") {
      print(temp_model,par=c("beta","sigma2","phi1","phi2","D1"), digits = 3)
      print(temp_criteria)
    }
    
    if(temp_struc == "CAR") {
      print(temp_model,par=c("beta","sigma2","phi1","D1"), digits = 3)
      print(temp_criteria)
    }
  }
  
  if(temp_dist == "Student") {
    
    if(temp_struc == "UNC") {
      print(temp_model,par=c("beta","sigma2","sigmab2","nu"), digits = 3)
      print(temp_criteria)
    }
    
    if(temp_struc == "DEC") {
      print(temp_model,par=c("beta","sigma2","phi1","phi2","D1","nu"), digits = 3)
      print(temp_criteria)
    }
    
    if(temp_struc == "CAR") {
      print(temp_model,par=c("beta","sigma2","phi1","D1","nu"), digits = 3)
      print(temp_criteria)
    }
  }
  
  if(temp_dist == "Slash") {
    
    if(temp_struc == "UNC") {
      print(temp_model,par=c("beta","sigma2","sigmab2","nu"), digits = 3)
      print(temp_criteria)
    }
    
    if(temp_struc == "DEC") {
      print(temp_model,par=c("beta","sigma2","phi1","phi2","D1","nu"), digits = 3)
      print(temp_criteria)
    }
    
    if(temp_struc == "CAR") {
      print(temp_model,par=c("beta","sigma2","phi1","D1","nu"), digits = 3)
      print(temp_criteria)
    }
  }
})


