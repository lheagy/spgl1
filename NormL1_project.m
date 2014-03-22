function x = NormL1_project(x,weights,tau,params)

if isreal(x)
   x = oneProjector(x,weights,tau);
else
   xa  = abs(x);
   idx = xa < eps;
   xc  = oneProjector(xa,weights,tau);
   xc  = xc ./ xa; xc(idx) = 0;
   x   = x .* xc;
end
