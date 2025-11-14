#' shinymle
#'
#' @returns a shinyapp
#' @export
#'
#' @examples
#' \dontrun{shinymle()}
shinymle<-function(){
  shiny::runApp(system.file("SHINY", package = "MATH4753F25fabian"
                            ), launch.browser = TRUE)
}
