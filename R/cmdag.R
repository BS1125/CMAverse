#' DAG For Causal Mediation Analysis
#'
#' Create the directed acyclic graph(DAG) for causal mediation analysis.
#'
#' @param outcome the variable name of the outcome
#' @param exposure the variable name of the exposure
#' @param mediator a character or vector of variable name(s) of the mediator(s)
#' @param prec (optional) a character or vector of variable name(s) of the pre-exposure confounder(s)
#' @param postc (optional) a character or vector of variable name(s) of the post-exposure confounder(s)
#' @param x.outcome x coordinate of the outcome. Default is \code{4}.
#' @param x.exposure x coordinate of the exposure. Default is \code{0}.
#' @param x.mediator x coordinate of the mediator. Default is \code{2}.
#' @param x.prec x coordinate of the pre-exposure confounder. Default is \code{2}.
#' @param x.postc x coordinate of the post-exposure confounder. Default is \code{2}.
#' @param y.outcome y coordinate of the outcome. Default is \code{0}.
#' @param y.exposure y coordinate of the exposure. Default is \code{0}.
#' @param y.mediator y coordinate of the mediator. Default is \code{1}.
#' @param y.prec y coordinate of the pre-exposure confounder. Default is \code{2}.
#' @param y.postc y coordinate of the post-exposure confounder. Default is \code{-0.5}.
#' @param ... additional arguments passed to \code{ggdag()}. See \link[ggdag]{ggdag} for details.
#'
#' @export
#' @examples
#' ## prec and postc are empty
#' cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
#'      prec = NULL, postc = NULL, node = TRUE, text_col = "white")
#'
#' ## postc is empty
#' cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
#'      prec = c("C1", "c2", "c3"), postc = NULL, node = FALSE, text_col = "black")
#'
#' ## prec and postc aren't empty
#' cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
#'      prec = c("C1", "c2", "c3"), postc = c("L1", "L2"), node = FALSE, text_col = "black")

cmdag <- function(outcome = NULL, exposure = NULL, mediator = NULL,
                  prec = NULL, postc = NULL,
                  x.outcome = 4, x.exposure = 0, x.mediator = 2, x.prec = 2, x.postc = 2,
                  y.outcome = 0, y.exposure = 0, y.mediator = 1, y.prec = 2, y.postc = -0.5,
                  ...) {

  usethis::use_package("ggdag")
  usethis::use_package("gridExtra")
  usethis::use_package("grid")

  require(ggdag)
  require(gridExtra)
  require(grid)

  if (is.null(outcome) | is.null(exposure) | is.null(mediator)) stop("Unspecified outcome, exposure or mediator")

  if (is.null(prec) && is.null(postc)) {

    dag <- dagify(Y ~ A + M, M ~ A, coords = list(x = c(A = x.exposure, Y = x.outcome, M = x.mediator),
                                                  y = c(A = y.exposure, Y = y.outcome, M = y.mediator)))

    grid.arrange(ggdag(dag, ...) + theme_dag_blank(),
                 bottom = textGrob(paste0("A(Exposure):", exposure, "\n M(Mediator):",
                                          paste(mediator, collapse = ", "),
                                          "\n Y(Outcome):", outcome), x = 1,
                                   hjust = 1, gp = gpar(fontface = 3L, fontsize = 8)))

  } else if (!is.null(prec) && is.null(postc)) {

    dag <- dagify(Y ~ A + M + C, M ~ A + C, A ~ C,
                  coords = list(x = c(A = x.exposure, Y = x.outcome, M = x.mediator, C = x.prec),
                                y = c(A = y.exposure, Y = y.outcome, M = y.mediator, C = y.prec)))
 
    grid.arrange(ggdag(dag, ...) + theme_dag_blank(),
                 bottom = textGrob(paste0("A(Exposure):", exposure, "\n M(Mediator):", paste(mediator, collapse = ", "),
                                          "\n Y(Outcome):", outcome, "\n C(Pre-exposure confounders): prec"), x = 1,
                                   hjust = 1, gp = gpar(fontface = 3L, fontsize = 8)))

  } else if (is.null(prec) && !is.null(postc)) {

    dag <- dagify(Y ~ A + M + L, M ~ A + L, L ~ A,
                  coords = list(x = c(A = x.exposure, Y = x.outcome, M = x.mediator, L = x.postc),
                                y = c(A = y.exposure, Y = y.outcome, M = y.mediator, L = y.postc)))

    grid.arrange(ggdag(dag, ...) + theme_dag_blank(),
                 bottom = textGrob(paste0("A(Exposure):", exposure, "\n M(Mediator):",
                                          paste(mediator, collapse = ", "), "\n Y(Outcome):",
                                          outcome, "\n L(Post-exposure confounders): postc"),
                                   x = 1,
                                   hjust = 1, gp = gpar(fontface = 3L, fontsize = 8)))

  } else if (!is.null(prec) && !is.null(postc)) {

    dag <- dagify(Y ~ A + M + C + L, M ~ A + C + L, A ~ C, L ~ A,
                  coords = list(x = c(A = x.exposure, Y = x.outcome, M = x.mediator, C = x.prec, L = x.postc),
                                y = c(A = y.exposure, Y = y.outcome, M = y.mediator, C = y.prec, L = y.postc)))

    grid.arrange(ggdag(dag, ...) + theme_dag_blank(),
                 bottom = textGrob(paste0("A(Exposure):", exposure, "\n M(Mediator):",
                                          paste(mediator, collapse = ", "), "\n Y(Outcome):",
                                          outcome, "\n C(Pre-exposure confounders): prec \n L(Post-exposure confounders): postc"), x = 1,
                                   hjust = 1, gp = gpar(fontface = 3L, fontsize = 8)))

  }

}
