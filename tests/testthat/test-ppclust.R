

test_that("PPCLUST", {
  data("diflogadenoma")
  data <- diflogadenoma[,10:13]
  cl <- ppclust(data, alpha = 1e-5)
  data("diflogcarcinoma")
  data <- diflogcarcinoma[,38:50]
  cl <- ppclust(data, alpha = 1e-5)
  data <- diflogcarcinoma[, 38:39]
  expect_error(ppclust(data, alpha = 1e-5))
  data <- diflogcarcinoma[1, 38:50]
  expect_error(ppclust(data, alpha = 1e-5))
  expect_error(ppclust("string", alpha = 1e-5))
  expect_error(ppclust(1:10, alpha = 1e-5))
  expect_error(ppclust(diflogadenoma, alpha = 1e-5))
  expect_error(ppclust(diflogadenoma[,-1], alpha = -1))
  expect_error(ppclust(diflogadenoma[,-1], alpha = 2))
})

test_that("PPCLUST-H", {
  data("diflogcarcinoma")
  data <- diflogcarcinoma[,38:50]
  cl <- ppclust_h(data, alpha = 1e-5)
  data <- diflogcarcinoma[, 38:39]
  expect_error(ppclust_h(data, alpha = 1e-5))
  data <- diflogcarcinoma[1, 38:50]
  expect_error(ppclust_h(data, alpha = 1e-5))
  expect_error(ppclust_h("string", alpha = 1e-5))
  expect_error(ppclust_h(1:10, alpha = 1e-5))
  expect_error(ppclust_h(diflogadenoma, alpha = 1e-5))
  expect_error(ppclust_h(diflogadenoma[,-1], alpha = -1))
  expect_error(ppclust_h(diflogadenoma[,-1], alpha = 2))
})


