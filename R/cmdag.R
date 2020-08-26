#' DAG For Causal Mediation Analysis
#'
#' Plot the directed acyclic graph (DAG) for causal mediation analysis.
#'
#' @param outcome the variable name of the outcome
#' @param exposure the variable name of the exposure
#' @param mediator a character or vector of variable name(s) of the mediator(s)
#' @param basec (optional) a character or vector of variable name(s) of the baseline confounder(s)
#' @param postc (optional) a character or vector of variable name(s) of the post-exposure 
#' confounder(s) affected by the exposure
#' @param x.outcome x coordinate of the outcome. Default is \code{4}.
#' @param x.exposure x coordinate of the exposure. Default is \code{0}.
#' @param x.mediator x coordinate of the mediator. Default is \code{2}.
#' @param x.basec x coordinate of the pre-exposure confounder. Default is \code{2}.
#' @param x.postc x coordinate of the post-exposure confounder. Default is \code{2}.
#' @param y.outcome y coordinate of the outcome. Default is \code{0}.
#' @param y.exposure y coordinate of the exposure. Default is \code{0}.
#' @param y.mediator y coordinate of the mediator. Default is \code{1}.
#' @param y.basec y coordinate of the baseline confounder. Default is \code{2}.
#' @param y.postc y coordinate of the post-exposure confounder. Default is \code{-0.5}.
#' @param caption.width line width in characters for the caption. Default is \code{50}.
#' @param caption.size text size in pts for the caption. Default is \code{10}.
#' @param ... additional arguments passed to \code{ggdag()}. See \link[ggdag]{ggdag} for details.
#'
#' @examples
#' ## basec and postc are empty
#' cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
#'      basec = NULL, postc = NULL, node = TRUE, text_col = "white")
#'
#'#' ## postc is empty
#' cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
#'      basec = c("C1", "C2", "C3"), postc = NULL, node = FALSE, text_col = "black")
#'      
#' ## basec is empty
#' cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
#'      basec = NULL, postc = c("L1", "L2"), node = FALSE, text_col = "black")
#'
#' ## basec and postc aren't empty
#' cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
#'      basec = c("C1", "C2", "C3"), postc = c("L1", "L2"), node = FALSE, text_col = "black")
#' 
#' @importFrom ggplot2 ggproto theme labs element_text
#' @importFrom ggdag dagify ggdag theme_dag_blank
#' @importFrom stringr str_wrap
#'
#' @export

cmdag <- function(outcome = NULL, exposure = NULL, mediator = NULL,
                  basec = NULL, postc = NULL,
                  x.outcome = 4, x.exposure = 0, x.mediator = 2, x.basec = 2, x.postc = 2,
                  y.outcome = 0, y.exposure = 0, y.mediator = 1, y.basec = 2, y.postc = -0.5,
                  caption.width = 50, caption.size = 10, ...) {
  if (is.null(outcome) | is.null(exposure) | is.null(mediator)) stop("Unspecified outcome, exposure or mediator")
  assign2glob <- function(key, val, pos) assign(key, val, envir=as.environment(pos))
  assign2glob("ggproto", ggplot2::ggproto, 1L)
  if (is.null(basec) && is.null(postc)) {
    dag <- dagify(Y ~ A + M, M ~ A, coords = list(x = c(A = x.exposure, Y = x.outcome, M = x.mediator),
                                                  y = c(A = y.exposure, Y = y.outcome, M = y.mediator)))
    p <- ggdag(dag, ...) + theme_dag_blank() + 
      labs(caption = paste0("A(Exposure): ", str_wrap(exposure, width = caption.width), 
                            "\nM(Mediator): ", 
                            str_wrap(paste(mediator, collapse = ", "), width = caption.width),
                            "\nY(Outcome): ", str_wrap(outcome, width = caption.width))) +
      theme(plot.caption = element_text(hjust = 0, face = "italic", size = caption.size))
    print(p)
  } else if (!is.null(basec) && is.null(postc)) {
    dag <- dagify(Y ~ A + M + C, M ~ A + C, A ~ C,
                  coords = list(x = c(A = x.exposure, Y = x.outcome, M = x.mediator, C = x.basec),
                                y = c(A = y.exposure, Y = y.outcome, M = y.mediator, C = y.basec)))
    p <- ggdag(dag, ...) + theme_dag_blank() + 
      labs(caption = paste0("A(Exposure): ", str_wrap(exposure, width = caption.width), 
                            "\nM(Mediator): ", 
                            str_wrap(paste(mediator, collapse = ", "), width = caption.width),
                            "\nY(Outcome): ", str_wrap(outcome, width = caption.width), 
                            "\nC(Baseline confounders): ", 
                            str_wrap(paste(basec, collapse = ", "), width = caption.width))) +
      theme(plot.caption = element_text(hjust = 0, face= "italic", size = caption.size))
    print(p)
  } else if (is.null(basec) && !is.null(postc)) {
    dag <- dagify(Y ~ A + M + L, M ~ A + L, L ~ A,
                  coords = list(x = c(A = x.exposure, Y = x.outcome, M = x.mediator, L = x.postc),
                                y = c(A = y.exposure, Y = y.outcome, M = y.mediator, L = y.postc)))
    p <- ggdag(dag, ...) + theme_dag_blank() + 
      labs(caption = paste0("A(Exposure): ", str_wrap(exposure, width = caption.width), 
                            "\nM(Mediator): ", 
                            str_wrap(paste(mediator, collapse = ", "), width = caption.width),
                            "\nY(Outcome): ", str_wrap(outcome, width = caption.width), 
                            "\nL(Post-exposure confounders): ", 
                            str_wrap(paste(postc, collapse = ", "), width = caption.width))) +
      theme(plot.caption = element_text(hjust = 0, face = "italic", size = caption.size))
    print(p)
  } else if (!is.null(basec) && !is.null(postc)) {
    dag <- dagify(Y ~ A + M + C + L, M ~ A + C + L, A ~ C, L ~ A,
                  coords = list(x = c(A = x.exposure, Y = x.outcome, M = x.mediator, C = x.basec, L = x.postc),
                                y = c(A = y.exposure, Y = y.outcome, M = y.mediator, C = y.basec, L = y.postc)))
    p <- ggdag(dag, ...) + theme_dag_blank() + 
      labs(caption = paste0("A(Exposure): ", str_wrap(exposure, width = caption.width), 
                            "\nM(Mediator): ", 
                            str_wrap(paste(mediator, collapse = ", "), width = caption.width),
                            "\nY(Outcome): ", str_wrap(outcome, width = caption.width), 
                            "\nC(Baseline confounders): ", 
                            str_wrap(paste(basec, collapse = ", "), width = caption.width),
                            "\nL(Post-exposure confounders): ", 
                            str_wrap(paste(postc, collapse = ", "), width = caption.width))) +
      theme(plot.caption = element_text(hjust = 0, face= "italic", size = caption.size))
    print(p)
  }
  rm(ggproto, envir = .GlobalEnv)
}
