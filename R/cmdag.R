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
                                          outcome, "\n L(Post-exposure confounders):",
                                          paste(postc, collapse = ", ")),
                                   x = 1,
                                   hjust = 1, gp = gpar(fontface = 3L, fontsize = 8)))

  } else if (!is.null(prec) && !is.null(postc)) {

    dag <- dagify(Y ~ A + M + C + L, M ~ A + C + L, A ~ C, L ~ A,
                  coords = list(x = c(A = x.exposure, Y = x.outcome, M = x.mediator, C = x.prec, L = x.postc),
                                y = c(A = y.exposure, Y = y.outcome, M = y.mediator, C = y.prec, L = y.postc)))

    grid.arrange(ggdag(dag, ...) + theme_dag_blank(),
                 bottom = textGrob(paste0("A(Exposure):", exposure, "\n M(Mediator):",
                                          paste(mediator, collapse = ", "), "\n Y(Outcome):",
                                          outcome, "\n C(Pre-exposure confounders): prec \n L(Post-exposure confounders):",
                                          paste(postc, collapse = ", ")), x = 1,
                                   hjust = 1, gp = gpar(fontface = 3L, fontsize = 8)))

  }

}
