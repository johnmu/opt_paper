library(rgl)
POINTDIR = "output/points.txt"
VALDIR = "output/val.txt"
SIMDIR = "output/simpleces.txt"
points = as.matrix(read.table(POINTDIR, header = F))
val = as.matrix(read.table(VALDIR, header = F))
simpleces = as.matrix(read.table(SIMDIR, header = F))
simpleces = simpleces + 1
open3d()
rgl.viewpoint(scale = c(1, 1, 1 / max(val)))
axes3d(xat = seq(0, 1, length.out = 6), yat = seq(0, 1, length.out = 6), zat = seq(0, floor(max(val) + 1), length.out = 6))
cols = cm.colors(floor(max(val) + 1))
for(i in 1:nrow(simpleces)) {
    x = c(points[simpleces[i, 1], 1], points[simpleces[i, 2], 1], points[simpleces[i, 3], 1])
    y = c(points[simpleces[i, 1], 2], points[simpleces[i, 2], 2], points[simpleces[i, 3], 2])
    #z = c(val[simpleces[i, 1], 1], val[simpleces[i, 2], 1], val[simpleces[i, 3], 1])
    z = rep(0, 3)
    triangles3d(x, y, z, col = cols[floor(z) + 1], front = "fill", back = "fill", lit = F)
    triangles3d(x, y, z, front = "line", back = "line", lit = F)
}
