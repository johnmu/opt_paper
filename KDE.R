Hel100 = rep(0, 10)
Hel1000 = rep(0, 10)
Hel10000 = rep(0, 10)
mass = pgamma(1, shape = 2, scale = 0.1) - pgamma(0, shape = 2, scale = 0.1)
for(k in 1:10) {
	sample = as.numeric(as.matrix(read.table(paste("Gamma_", as.character(k), sep = ""), head = F)))
	for(i in 2:4) {
		KDE = density(sample[1:10^i], bw = "nrd0", adjust = 1, kernel = "gaussian", from = 0, to = 1, n = 2^10)
		trueVal = dgamma(KDE$x, shape = 2, scale = 0.1)
		Hel = 0;
		for(j in 1:length(KDE$x)) {
			Hel = Hel + (sqrt(trueVal[j]) - sqrt(KDE$y[j]))^2
		}
		if(i == 2) {
			Hel100[k] = sqrt(Hel / 2 / length(KDE$x))
		}
		else if(i == 3) {
			Hel1000[k] = sqrt(Hel / 2 / length(KDE$x))
		}
		else {
			Hel10000[k] = sqrt(Hel / 2 / length(KDE$x))
		}
	}
}
