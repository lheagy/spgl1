function [ok] = spgl1SamplingTest()

m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));

params.A = opGaussian(m, n); % initial operator
params.xTrue = x0;  % need this to generate data


% first run: don't resample
opts = spgSetParms('optTol',1e-4, ...
                   'project', @NormL1_project, ...
                   'primal_norm', @NormL1_primal, ...
                   'dual_norm', @NormL1_dual, ...
                   'subspaceMin', 0, ...
                   'resample',   0, ...
                   'linear',     1 ...
                   );

sigma = 1e-3;  
%% Run three separate root finding problems


% first: no sampling 

[xLSnosamp, r, g, info] = spgl1(@funForward, @b, [], sigma, [], opts, params);


% second: resampling when we hit Pareto curve
opts.resample = 1;
[xLSsampPareto, r, g, info] = spgl1(@funForward, @b, [], sigma, [], opts, params);


% third: resampling at every iteration  DOESNT WORK
%opts.resample = 2;
%[xLSsampAll, r, g, info] = spgl1(@funForward, @b, [], sigma, [], opts, params);


%[xLS, r, g, info] = spgl1Classic(A, b, [], sigma, 0*x0, opts);

%%


%sigma = 1e0;
%[xMY, r, g, info] = nlspgl1(@funForward, @funMY, params, b, [], sigma, 0*x0, opts );


%% Plot the results

 figure(1)
 plot(1:n, x0 + 2, 1:n, xLSnosamp(1:n), 1:n, xLSsampPareto - 2, '-k', 'Linewidth', 1.5);
 legend('true', 'Regular', 'redrawn');
% 
 %figure(1)
 %plot(1:n, x0 + 2, 1:n, xLSnosamp(1:n), 1:n, xLSsampPareto - 2, 1:n, xLSsampAll - 4, '-k', 'Linewidth', 1.5);

 
 
% figure(1)
% plot(1:n, x0 + 2, 1:n, xLSnosamp(1:n), 1:n, xHUB(1:n) - 2, 1:n, xST(1:n) -4,  '-k', 'Linewidth', 1.5);
% 
% legend('true', 'LS', 'Hub', 'ST');
% 
% figure(2);
% plot(1:m, b - A*x0 + 1.5, 1:m, b - A*xLS(1:n), 1:m, b - A*xHUB(1:n) - 1.5,  1:m, b-A*xST(1:n) - 3, '-k', 'Linewidth', 2); 
% legend('true', 'LS', 'Hub', 'ST');

    function [data params] = b(flag, params)
        
     if(~flag) % use current operator
        data = params.A*params.xTrue;
     else      % renew the operator, and then use it
        params.A = opGaussian(size(params.A, 1), size(params.A, 2));
        data = params.A*params.xTrue;
     end
    end

    function f = funForward(x, g, params)
    if isempty(g)
        f = params.A*x;
    else
        f = params.A'*g;
    end
end
 


end

