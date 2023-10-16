
xplxplDistNum <- function(opt.path, pars, n.init, n.tot) {
  # compute distance matrix to previuos points
  dist = dist(opt.path[pars], method = "euclidean")
  dist.matrix = as.matrix(dist)
  dist.matrix.low = dist.matrix * lower.tri(dist.matrix) # lower triangle distance matrix
  
  # compute pairwise distances (take 0 as first distance, because there is no pairwise 
  # distance of the very first point)
  pairwise.dist = c(0, diag(dist.matrix[,-1])) 
  # ignore the initial points, i.e. focus on the proposed points (0 is inserted 
  # because first proposed point has no distance to previosly proposed points, 
  # only to initial points that we are not interested in)
  pairwise.dist.prop = c(0, pairwise.dist[(n.init + 2):n.tot])
  
  # compute minimal distance
  # take minimum of lower triangle distance matrix row by row, because only distances 
  # to previous points matter not to future ones 
  min.dist = apply(dist.matrix.low, MARGIN = 1,  
                      FUN = function(x) min(x[x > 0]))
  # ignore the initial points, i.e. focus on the proposed points
  min.dist.prop = min.dist[(n.init + 1):n.tot]
  
  #compute maximal distance
  max.dist = apply(dist.matrix.low, MARGIN = 1, 
                   FUN = function(x) max(x[x > 0]))
  # ignore the initial points, i.e. focus on the proposed points
  max.dist.prop = max.dist[(n.init + 1):n.tot]
  
  # compute mean distance vector whose element i contains mean distance from 
  # proposed point i to all previous points untill iteration i
  acc.dist = apply(dist.matrix.low, MARGIN = 1, FUN = sum)
  mean.dist = acc.dist / 1:n.tot
  # ignore the initial points, i.e. focus on the proposed points
  mean.dist.prop = mean.dist[(n.init + 1):n.tot]
  
  
  xpl.xpl.distance = data.frame("Pairwise Euclidean Distance" = pairwise.dist.prop, 
                                "Minimal Euclidean Distance" = min.dist.prop,
                                "Maximal Euclidean Distance" = max.dist.prop, 
                                "Mean Euclidean Distance" = mean.dist.prop)
  xpl.xpl.distance

}


