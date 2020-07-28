test.occuRN.fit <- function() {
  
  data(birds)
  woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin)

  # survey occasion-specific detection probabilities
  fm_C <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="C")   
  fm_R <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="R")

  # check that output matches
  checkEqualsNumeric(coef(fm_C),coef(fm_R),tol=1e-5)

  # check output is correct
  checkEqualsNumeric(coef(fm_C),
    c(0.7921122,-1.8328867,0.4268205,-0.1442194,0.4634105,0.7787513,
      0.8008794,1.0569827,0.8048578,0.8779660,0.9374874,0.7064848),tol=1e-5)

}

test.occuRN.na <- function() {

  data(birds)
  woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin)
  
  #Remove one observation
  woodthrushUMF@y[1,1] <- NA

  fm_C <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="C")   
  fm_R <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="R")

  # check that output matches
  checkEqualsNumeric(coef(fm_C),coef(fm_R),tol=1e-5)

  # check output is correct
  checkEqualsNumeric(coef(fm_C),
    c(0.793042, -1.902789, 0.494098, -0.074573, 0.53074, 0.845903,
    0.867936, 1.123959, 0.871912, 0.944917, 1.004499, 0.773679), tol=1e-5)

  #Remove entire site
  woodthrush.bin_na <- woodthrush.bin
  woodthrush.bin_na[1,] <- NA
  woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin_na)

  fm_C <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="C")   
  fm_R <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="R")

  # check that site was removed
  checkEqualsNumeric(fm_C@sitesRemoved,1)

  # check that output matches
  checkEqualsNumeric(coef(fm_C),coef(fm_R),tol=1e-5)

  # check output is correct
  checkEqualsNumeric(coef(fm_C),
    c(0.783066, -1.920232, 0.448369, -0.009701, 0.490085, 0.814767,
    0.837669, 1.097903, 0.842467, 0.916831, 0.976707, 0.740672), tol=1e-5)
}

#Test that parboot works
test.occuRN.parboot <- function(){
  
  data(birds)
  woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin)
  fm <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="C")  
  
  chisq2 <- function(fm) {
   observed <- getY(fm)
   expected <- fitted(fm)
   sum((observed - expected)^2/expected, na.rm=T)  
  }
  
  set.seed(123)
  pb <- parboot(fm, statistic=chisq2, nsim=2)
  checkEqualsNumeric(as.numeric(pb@t.star), c(342.2285, 318.0965), tol=1e-5)
}
