library(orthGS)


## ---------------------------------------------- ##
#                 Testing orthG                    #
## ---------------------------------------------- ##
test_that("orthG() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- orthG(sp = "Pp")
  a <- orthG(set = c("Pp", "Psy", "Psm", "Ap"))

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 3)
  expect_equal(ncol(a), 19)

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 9)
  expect_equal(ncol(b), 19)

})

## ---------------------------------------------- ##
#                 Testing orthP                    #
## ---------------------------------------------- ##
test_that("orthP() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- orthP(phylo_id = "Zm_GS1b_1", set = c("Arabidopsis thaliana", "Oryza sativa"))
  b <- orthP(phylo_id = "Zm_GS1b_1", set = "Ath")
  c <- orthP(phylo_id = "Zm_GS1b_1", set = "all")

  expect_is(a, 'list')
  expect_equal(length(a), 3)
  expect_is(a[[2]], 'character')
  expect_is(a[[3]], 'character')
  expect_true("Zm_GS1b_1" %in% a[[3]])

  expect_is(b, 'character')

  expect_is(c, 'list')
  expect_equal(length(c), 3)
  expect_is(c[[2]], 'character')
  expect_is(c[[3]], 'character')
  expect_true("Zm_GS1b_1" %in% a[[3]])

})

## ---------------------------------------------- ##
#               Testing getseqGS                   #
## ---------------------------------------------- ##
test_that("getseqGS() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- getseqGS(id = "Pp_GS1b_2", molecule = "Prot")
  b <- getseqGS(id = "Pp_GS1b_2", molecule = "CDS")

  expect_is(a, 'character')
  expect_equal(nchar(a), 357)

  expect_is(b, 'character')
  expect_equal(nchar(b), 1074)

})

## ---------------------------------------------- ##
#                Testing subsetGS                  #
## ---------------------------------------------- ##
test_that("subsetqGS() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- subsetGS(sp = c("Arabidopsis thaliana", "Oryza sativa"))
  b <- subsetGS(sp = "Ath")

  expect_is(a, "data.frame")
  expect_equal(dim(a), c(10,19))
  expect_is(b, "data.frame")
  expect_equal(dim(b), c(6,19))

})
