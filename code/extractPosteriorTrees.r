library("ape")

trees1 <- ape::read.nexus('output/sample30.run1.t')
trees2 <- ape::read.nexus('output/sample30.run2.t')

l1 <- length(trees1)
l2 <- length(trees2)

trees1 <- trees1[(l1 %/% 2):l1]
trees2 <- trees2[(l1 %/% 2):l2]

trees <- c(trees1, trees2)
idc <- (1:1000) * floor(length(trees)/1000)

write.tree(trees[idc], 'data/posterior.tree')



