test.powerAnalysis <- function(){
  forms <- list(state=~elev, det=~1)
  coefs <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))
  design <- list(M=300, J=8) # 300 sites, 8 occasions per site
  occu_umf <- simulate("occu", formulas=forms, coefs=coefs, design=design)

  template_model <- occu(~1~elev, occu_umf)
  checkException(powerAnalysis(template_model))

  effect_sizes <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))

  set.seed(123)
  pa <- powerAnalysis(template_model, coefs=effect_sizes, alpha=0.05)
  checkTrue(inherits(pa, "unmarkedPower"))
  s <- summary(pa)$Power
  checkTrue(s[2]>0.7)

  # fewer sites
  set.seed(123)
  pa2 <- powerAnalysis(template_model, effect_sizes, design=list(M=50, J=3))
  checkTrue(all(summary(pa2)$Power < 0.3))

  # more sites
  set.seed(123)
  pa3 <- powerAnalysis(template_model, effect_sizes, design=list(M=400, J=4))
  checkTrue(summary(pa3)$Power[2] >0.7)

  # set null
  set.seed(123)
  nul <- list(state=c(intercept=5, elev=0), det=c(intercept=0))
  pa4 <- powerAnalysis(template_model, effect_sizes, nulls=nul)
  checkTrue(summary(pa4)$Power[1]==1)
  checkEqualsNumeric(summary(pa4)$Null, c(5,0,0))

  # list
  pl <- unmarkedPowerList(list(pa, pa2, pa3, pa4))
  checkTrue(inherits(pl, "unmarkedPowerList"))

  # generate list
  scenarios <- expand.grid(M=c(50,100), J=c(2,3))
  pl <- unmarkedPowerList(template_model, effect_sizes, design=scenarios)
  checkTrue(inherits(pl, "unmarkedPowerList"))
}
