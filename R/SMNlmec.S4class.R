
#' @import rstan

#' @export
setClass(
  Class = "SMNlmecfit",
  representation(
    stan_object = "stanfit",
    model_criteria = "data.frame",
    dist_set = "character",
    struc_set = "character"
  )
)


#' @export
SMNlmecfit.creator <- function(stan_object, model_criteria, dist_set, struc_set){

  if(!inherits(stan_object, "stanfit")){
    stop("stan_object must be a stanfit from rstan package.")
  }

  if(!is.data.frame(model_criteria)){
    stop("model_criteria must be a data frame.")
  }

  if(!is.character(dist_set) || length(dist_set) != 1) {
    stop("dist must be a single character string.")
  }

  if(!is.character(struc_set) || length(struc_set) != 1) {
    stop("structure must be a single character string.")
  }

  new("SMNlmecfit", stan_object = stan_object, model_criteria = model_criteria,
      dist_set = dist_set, struc_set = struc_set)
}

#' @export
setGeneric("SMNlmec.summary",function(object){
 standardGeneric("SMNlmec.summary")
})


#' @export
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

