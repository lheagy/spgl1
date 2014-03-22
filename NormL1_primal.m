function p = NormL1_primal(x,weights,params)

p = norm(x.*weights,1);
