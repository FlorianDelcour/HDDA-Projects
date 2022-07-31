setwd("C:/Unif/Master 1/Q1/HDDA/Project 3")

library(Thermimage)
library(imager)
library(magrittr)
image1 <- rotate270.matrix(t(as.matrix(load.image("DelcourMakedonsky_Mix1.jpg"))))
image2 <- rotate270.matrix(t(as.matrix(load.image("DelcourMakedonsky_Mix2.jpg"))))
image3 <- rotate270.matrix(t(as.matrix(load.image("DelcourMakedonsky_Mix3.jpg"))))
image_tot <- cbind(as.vector(image1), as.vector(image2), as.vector(image3))




library(fastICA)
set.seed(1)
ICA <- fastICA(image_tot, 3, alg.typ = "deflation", fun = "logcosh", alpha = 1, 
                  method = "C", row.norm = FALSE, maxit = 200, 
                  tol = 0.0000001, verbose = TRUE)

library(Matrix)
est_image_1 <- matrix(ICA$S[,1], nrow = 300, ncol = 300)
est_image_2 <- matrix(ICA$S[,2], nrow = 300, ncol = 300)
est_image_3 <- matrix(ICA$S[,3], nrow = 300, ncol = 300)

mix_matrix <- ICA$A
require("grDevices")

image(est_image_1, col = grey(seq(0, 1, length = 256)), axes = FALSE)
dev.print(device = png, filename = "image1.png", width = 300, height = 300)
image(est_image_2, col = grey(seq(0, 1, length = 256)), axes = FALSE)
dev.print(device = png, filename = "image2.png", width = 300, height = 300)
image(est_image_3, col = grey(seq(0, 1, length = 256)), axes = FALSE)
dev.print(device = png, filename = "image3.png", width = 300, height = 300)


library(corrplot)
corrplot(cor(ICA$S))
cov_images <- cov(ICA$S)

hist(est_image_1)
hist(est_image_2)
hist(est_image_3)

