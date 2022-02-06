#' Append List
#'
#' The function appends objects into given list.
#'
#' @param lst a list to append object to
#' @param obj object to append
#'
#' @return a list with appended object
#' @export
#'
#' @examples
#' myList <- list()
#' myList <- lappend(myList, mtcars)
lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

#' Saves RData object
#'
#' The function saves an object into specified directory naming with the name of the object.
#'
#' @param object an object to save
#' @param DATA_PATH a directory path to save object
#' @param suffix (optional, Default = "") adds suffix string to the saved object
#'
#' @return
#' @export
#'
#' @examples
#' save_object(mtcars, DATA_PATH="./", suffix="myresults")
save_object <- function(object, DATA_PATH, suffix = "") {
  object_name = deparse(substitute(object))
  file_name = paste(object_name, suffix, "RData", sep = ".")
  file_path = paste(DATA_PATH, file_name, sep="/")
  save(list = eval(object_name), file = file_path)
}



