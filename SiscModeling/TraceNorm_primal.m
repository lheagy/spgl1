% Authors: Rajiv Kumar, Hassan Mansour, and Aleksandr Aravkin. 


function p = TraceNorm_primal(x, weights, params)

p = 0.5 * norm(x.*weights)^2; 

end