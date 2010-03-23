test.occu.fit <- function() {
	M <- 50
	J <- 4
	
	# simulation with factors
	detParms <- c(-0.5, 0.5)
	occParms <- c(0, 1, 0, 0.5, 1, 1)
	
	abundance.cov1 <- rnorm(M)
	abundance.cov2 <- as.factor(rep(1:5, each = M/5))
	detection.cov <- rnorm(M*J)
	X <- model.matrix(~ abundance.cov1 + abundance.cov2)
	V <- cbind(rep(1, M*J), detection.cov)
	
	psi.probs <- plogis(X %*% occParms)
	occ <- rbinom(M, 1, psi.probs)
	
	p <- plogis(V %*% detParms)
	p <- matrix(p, M, J, byrow = TRUE)
	
	y <- matrix(NA, M, J)
	for(i in 1:J) {
		y[,i] <- rbinom(M, 1, p[,i])*occ
	}
	
	#abundance.cov1[sample(1:M, 100)] <- NA
	#detection.cov[sample(1:M, 100)] <- NA
	umf <- unmarkedFrameOccu(y,
			siteCovs = data.frame(abundance.cov1, abundance.cov2),
			obsCovs = data.frame(detection.cov))
	
	fm <- occu(~detection.cov ~ abundance.cov1 + abundance.cov2 , umf)
	
	## check that detection is w/in 3 SE of truth
	checkTrue(all(detParms - 3*SE(fm['det']) < coef(fm,'det')) &
			all(coef(fm,'det') < detParms + 3*SE(fm['det'])))

	## check that occupancy params are w/in 3 SE of truth
	checkTrue(all(occParms - 3*SE(fm['state']) < coef(fm,'state')) &
				all(coef(fm,'state') < occParms + 3*SE(fm['state'])))
}
