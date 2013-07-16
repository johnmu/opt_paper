args <- commandArgs(T)
dim = as.numeric(args[1])
num = as.numeric(args[2])
me = rep(0, dim)
sig = rep(0, dim)
mat = matrix(0, dim, dim)
for(i in 1:dim) {
    me[i] = as.numeric(args[2 * i + 1])
    sig[i] = as.numeric(args[2 * i + 2])
}
ind = 2 * dim + 3
for(i in 1:dim) {
    for(j in 1:dim) {
        mat[i, j] = as.numeric(args[ind])
        ind = ind + 1
    }
}
cat("", file = "dat.txt", append = F)
for(i in 1:num) {
    tmp = rep(0, dim)
    for(j in 1:dim) {
        tmp[j] = rnorm(1, me[j], sig[j])
    }
    onesam = rep(0, dim)
    for(j in 1:dim) {
        for(k in 1:dim) {
            onesam[j] = onesam[j] + mat[j, k] * tmp[k]
        }
    }
    cat(onesam, "\n", file = "dat.txt", append = T)
}
