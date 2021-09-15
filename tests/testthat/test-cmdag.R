context("cmdag works correctly")

test_that("cmdag works correctly", {

  p1 <- cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
          basec = NULL, postc = NULL, node = TRUE, text_col = "white")
  p2 <- cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
          basec = c("C1", "C2", "C3"), postc = NULL, node = FALSE, text_col = "black")
  p3 <- cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
          basec = NULL, postc = c("L1", "L2"), node = FALSE, text_col = "black")
  p4 <- cmdag(outcome = "Y", exposure = "A", mediator = c("M1", "M2"),
        basec = c("C1", "C2", "C3"), postc = c("L1", "L2"), node = FALSE, text_col = "black")

  # test
  expect_equal(p1, NULL)
  expect_equal(p2, NULL)
  expect_equal(p3, NULL)
  expect_equal(p4, NULL)

})