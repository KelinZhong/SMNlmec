
#' @import rstan
#' @import methods



#' @title SMNlmecfit
#' @slot stan_object stanfit object from rstan.
#' @slot model_criteria list, model selection criteria.
#' @slot dist_set character, the name of distribution.
#' @slot struc_set character, the name of correlation structure.
#' @exportClass SMNlmecfit
setClass(
  Class = "SMNlmecfit",
  representation(
    stan_object = "stanfit",
    model_criteria = "list",
    dist_set = "character",
    struc_set = "character"
  )
)

#' @title Creator Function for SMNlmecfit Class
#' @description A function to create objects of class \code{SMNlmecfit}.
#' @name SMNlmecfit.creator
#' @param stan_object stanfit object from rstan.
#' @param model_criteria list, model selection criteria.
#' @param dist_set character, the name of distribution.
#' @param struc_set character, the name of correlation structure.
#' @export
SMNlmecfit.creator <- function(stan_object, model_criteria, dist_set, struc_set){

  if(!inherits(stan_object, "stanfit")){
    stop("stan_object must be a stanfit from rstan package.")
  }

  if(!is.list(model_criteria)){
    stop("model_criteria must be a list.")
  }

  if(!is.character(dist_set) || length(dist_set) != 1) {
    stop("dist must be a single character string.")
  }

  if(!is.character(struc_set) || length(struc_set) != 1) {
    stop("structure must be a single character string.")
  }

  methods::new("SMNlmecfit", stan_object = stan_object, model_criteria = model_criteria,
      dist_set = dist_set, struc_set = struc_set)
}



#' @title Summary Method for SMNlmecfit Class
#' @description A generic function to provide a summary for objects of class \code{SMNlmecfit}.
#' @param object An object of class \code{SMNlmecfit}.
#' @return A summary of model estimations, R-hats, standard errors, and criteria.
#' @export
#' @rdname SMNlmec.summary
setGeneric("SMNlmec.summary",function(object){
  standardGeneric("SMNlmec.summary")
})


#' @rdname SMNlmec.summary
#' @aliases SMNlmec.summary,SMNlmecfit-method
#' @param object An object of class \code{SMNlmecfit}.
#' @return A printed summary of the SMNlmecfit object.
#' @export
setMethod("SMNlmec.summary","SMNlmecfit",function(object){
  temp_dist <- object@dist_set
  temp_struc <- object@struc_set
  temp_model <- object@stan_object
  temp_criteria <- object@model_criteria
  crit <- as.data.frame(cbind(temp_criteria$LPML,temp_criteria$DIC,temp_criteria$EAIC,temp_criteria$EBIC))
  colnames(crit) <- c("LPML","DIC","EAIC","EBIC")

  if(temp_dist == "Normal") {

    if(temp_struc == "UNC") {
      print(temp_model,par=c("beta","sigma2","sigmab2"), digits = 3)
      print(crit)
    }

    if(temp_struc == "DEC") {
      print(temp_model,par=c("beta","sigma2","phi1","phi2","D1"), digits = 3)
      print(crit)
    }

    if(temp_struc == "CAR") {
      print(temp_model,par=c("beta","sigma2","phi1","D1"), digits = 3)
      print(crit)
    }
  }

  if(temp_dist == "Student") {

    if(temp_struc == "UNC") {
      print(temp_model,par=c("beta","sigma2","sigmab2","nu"), digits = 3)
      print(crit)
    }

    if(temp_struc == "DEC") {
      print(temp_model,par=c("beta","sigma2","phi1","phi2","D1","nu"), digits = 3)
      print(crit)
    }

    if(temp_struc == "CAR") {
      print(temp_model,par=c("beta","sigma2","phi1","D1","nu"), digits = 3)
      print(crit)
    }
  }

  if(temp_dist == "Slash") {

    if(temp_struc == "UNC") {
      print(temp_model,par=c("beta","sigma2","sigmab2","nu"), digits = 3)
      print(crit)
    }

    if(temp_struc == "DEC") {
      print(temp_model,par=c("beta","sigma2","phi1","phi2","D1","nu"), digits = 3)
      print(crit)
    }

    if(temp_struc == "CAR") {
      print(temp_model,par=c("beta","sigma2","phi1","D1","nu"), digits = 3)
      print(crit)
    }
  }
})




