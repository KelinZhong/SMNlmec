
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




