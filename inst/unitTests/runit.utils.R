test.formatLong <- function() {
  df <- read.csv(system.file("csv","frog2001pcru.csv", package = "unmarked"))
  umf <- formatLong(df, type = "unmarkedFrameOccu")
  ## Add some assertions...
}





test.rowProds <- function()
{
    m <- matrix(1:9/10, 3)
    rowProds <- unmarked:::rowProds
    checkEquals(rowProds(m), apply(m, 1, prod))    
}





test.tranProbs <- function()
{
    N <- 0:100
    omega <- 0.5
    gamma <- 2
    delta <- 5
    tranProbs <- unmarked:::tranProbs
    tranProbsR <- unmarked:::tranProbsR
    tp1 <- tranProbs(N, omega, gamma, delta, "constant")
    tp1r <- tranProbsR(N, omega, gamma, delta, "constant")
    checkEquals(tp1, tp1r)
    
    checkTrue(all.equal(colSums(tp1), rep(1, 101)))
    
}
    
    
