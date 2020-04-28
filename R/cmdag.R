cmdag <- function(outcome = NULL, exposure = NULL, mediator = NULL,
                  prec = NULL, postc = NULL) {

  require(ggdag)
  require(gridExtra)
  require(grid)

  if (is.null(postc)) {

  dag <- dagify(Y~A+M+C,M~A+C,A~C,coords=list(x=c(A=0,Y=4,M=2,C=2),y=c(A=0,Y=0,M=1,C=2)))

  grid.arrange(ggdag(dag,node=FALSE,text_col = "black") + theme_dag_blank(),
               bottom = textGrob(paste0("A(Exposure):", exposure, "\n M(Mediator):", paste(mediator, collapse = ", "),
                                        "\n Y(Outcome):", outcome, "\n C(Pre-exposure confounders)"), x = 1,
                                 hjust = 1, gp = gpar(fontface = 3L, fontsize = 8)))

  } else {

    dag <- dagify(Y~A+M+C+L,M~A+C+L,A~C,L~A,coords=list(x=c(A=0,Y=4,M=2,C=2,L=2),y=c(A=0,Y=0,M=1,C=2,L=-0.5)))

    grid.arrange(ggdag(dag,node=FALSE,text_col = "black") + theme_dag_blank(),
                 bottom = textGrob(paste0("A(Exposure):", exposure, "\n M(Mediator):",
                                          paste(mediator, collapse = ", "), "\n Y(Outcome):",
                                          outcome, "\n C(Pre-exposure confounders) \n L(Post-exposure confounders):",
                                          paste(postc, collapse = ", ")),
                                   x = 1,
                                   hjust = 1, gp = gpar(fontface = 3L, fontsize = 8)))


  }

}
