context("mixed model tools")

test_that("get_reTrms matches lme4::mkReTrms", {

  skip_if(!requireNamespace("lme4", quietly=TRUE), 
          "lme4 package unavailable")
  
  dat <- data.frame(x = rnorm(20), y = rnorm(20), z = factor(sample(letters[1:3], 20, replace=T)),
                    group = factor(sample(letters[4:6], 20, replace=T)),
                    id = factor(sample(letters[7:9], 20, replace=T)))

  form1 <- ~x + (1|group)
  r1 <- lme4::mkReTrms(lme4::findbars(form1), dat)
  r2 <- get_reTrms(form1, dat)
  expect_identical(r2$Z, Matrix::t(r1$Zt))
  expect_identical(r1$cnms, r2$cnms)
  attributes(r1$flist)$assign <- NULL
  expect_identical(r1$flist, r2$flist)

  form1 <- ~x + (x||group)
  r1 <- lme4::mkReTrms(lme4::findbars(form1), dat)
  r2 <- get_reTrms(form1, dat)
  expect_identical(r2$Z, Matrix::t(r1$Zt))
  expect_identical(r1$cnms, r2$cnms)
  attributes(r1$flist)$assign <- NULL
  expect_identical(r1$flist, r2$flist)

  form1 <- ~x + (x||group) + (1|id)
  r1 <- lme4::mkReTrms(lme4::findbars(form1), dat)
  r2 <- get_reTrms(form1, dat)
  expect_identical(r2$Z, Matrix::t(r1$Zt))
  expect_identical(r1$cnms, r2$cnms)
  attributes(r1$flist)$assign <- NULL
  expect_identical(r1$flist, r2$flist)

  form1 <- ~x + (x||group) + (y||id)
  r1 <- lme4::mkReTrms(lme4::findbars(form1), dat)
  r2 <- get_reTrms(form1, dat)
  expect_identical(r2$Z, Matrix::t(r1$Zt))
  expect_identical(r1$cnms, r2$cnms)
  attributes(r1$flist)$assign <- NULL
  expect_identical(r1$flist, r2$flist)

  form1 <- ~x + (x*y||group) + (y||id)
  r1 <- lme4::mkReTrms(lme4::findbars(form1), dat)
  r2 <- get_reTrms(form1, dat)
  expect_identical(r2$Z, Matrix::t(r1$Zt))
  expect_identical(r1$cnms, r2$cnms)
  attributes(r1$flist)$assign <- NULL
  expect_identical(r1$flist, r2$flist)

  form1 <- ~(1|group)
  r1 <- lme4::mkReTrms(lme4::findbars(form1), dat)
  r2 <- get_reTrms(form1, dat)
  expect_identical(r2$Z, Matrix::t(r1$Zt))
  expect_identical(r1$cnms, r2$cnms)
  attributes(r1$flist)$assign <- NULL
  expect_identical(r1$flist, r2$flist)

  form1 <- ~(1|group) + x
  r1 <- lme4::mkReTrms(lme4::findbars(form1), dat)
  r2 <- get_reTrms(form1, dat)
  expect_identical(r2$Z, Matrix::t(r1$Zt))
  expect_identical(r1$cnms, r2$cnms)
  attributes(r1$flist)$assign <- NULL
  expect_identical(r1$flist, r2$flist)
})

test_that("get_reTrms errors correctly", {
  form1 <- ~x + (x+z||group) + (y||id)
  expect_error(get_reTrms(form1, dat))

  form1 <- ~x + (x|group:id)
  expect_error(get_reTrms(form1, dat))

  form1 <- ~x + (x|group/id)
  expect_error(get_reTrms(form1, dat))

  expect_identical(find_bars(NULL), NULL)
})
