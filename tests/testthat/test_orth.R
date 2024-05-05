library(orthGS)

## ---------------------------------------------- ##
#                 Testing orthG                    #
## ---------------------------------------------- ##
test_that("orthG() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- orthG(set = "Pp")
  b <- orthG(set = c("Pp", "Psy", "Psm", "Ap"))
  c <- orthG(set = "all")

  expect_is(a, 'list')
  expect_is(a[[1]], 'matrix')
  expect_is(a[[2]], 'igraph')
  expect_equal(nrow(a[[1]]), 3)

  expect_is(b, 'list')
  expect_is(b[[1]], 'matrix')
  expect_is(b[[2]], 'igraph')
  expect_equal(nrow(b[[1]]), 11)
  expect_equal(ncol(b[[1]]), 11)

  expect_is(c, 'list')
  expect_is(c[[1]], 'matrix')
  expect_is(c[[2]], 'igraph')
  expect_equal(nrow(c[[1]]), 155)
  expect_equal(ncol(c[[1]]), 155)

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

  a <- getseqGS(phylo_id = "Pp_GS1b_2", molecule = "Prot")
  b <- getseqGS(phylo_id = "Pp_GS1b_2", molecule = "CDS")

  expect_is(a, 'character')
  expect_equal(nchar(a), 357)

  expect_is(b, 'character')
  expect_equal(nchar(b), 1074)

})

## ---------------------------------------------- ##
#                Testing subsetGS                  #
## ---------------------------------------------- ##
test_that("subsetGS() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- subsetGS(sp = c("Arabidopsis thaliana", "Oryza sativa"))
  b <- subsetGS(sp = "Ath")

  expect_is(a, "data.frame")
  expect_equal(dim(a), c(10,23))
  expect_is(b, "data.frame")
  expect_equal(dim(b), c(6,23))

})


## ---------------------------------------------- ##
#                Testing speciesGS                 #
## ---------------------------------------------- ##
test_that("subsetGS() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- speciesGS(c("Pinus pinaster", "Ath"))

  expect_is(a, "NULL")

})
