function d = NormL1_dual(x,weights,params)

d = norm(x./weights,inf);
