% SISC script for small Movielens Example
% Authors: Rajiv Kumar, Hassan Mansour, and Aleksandr Aravkin. 

% Instructions: 
% This script depends on fast Omega on LR' operators, which can be obtained
% from outside libraries. You can find either package in the folder
% 'Outside Libraries', unzip it, and place it in the same directory as
% spgl1General. 

% This script then installs mex files for whichever library you want. 





clear all;clc;close all;
% add everything in the path
addpath(genpath('../../'))


% Which forward model to use? 

model = 'partXY';       % Operator from Vanderecken's code; seems to be fastest
%model = 'XonOmega';    % operator from SVT package


switch(model)
    case{'XonOmega'}
        fprintf('Installing mex files for SVT package...\n');
        run ../../SVT_MostRecent/utility/install_mex.m
        
    case{'partXY'}
        fprintf('Installing mex files for Riemmanian package...\n');
        run ../../RiemannianMatrixCompletion_6Jan2014/startup.m
        run ../../RiemannianMatrixCompletion_6Jan2014/Install_mex.m

    
    otherwise{'unknown model selected'};
end



% Read in data
fprintf('Loading data...\n');
D=dlmread('ratings.dat');
userid=6040;
movie=3952;
%Define the intial incomplete matrix
Dinit=zeros(userid,movie);
l=length(D(:,1));
% fill the matrix with the given data
for Z=1:l
    k=D(Z,1);
    m=D(Z,3);
    Dinit(k,m)=D(Z,5);
end


%% Set experimental parameters

fprintf('Running solver...\n');


Dtest=vec(Dinit);
%rank=[10 20 30 50];
rank = 10; % 

Sigmafact=[5e-1 3e-1 2e-1];
%Sigmafact = 5e-1;


initfact = 1e-4;
opts.tol = 1e-5;
options = spgSetParms('optTol',1e-6, ...
                    'bpTol', 1e-6,...
                    'decTol',1e-4,...
                    'project', @TraceNorm_project, ...
                    'primal_norm', @TraceNorm_primal, ...
                    'dual_norm', @TraceNorm_dual, ...
                    'proxy', 1, ...
                    'ignorePErr', 1, ...
                    'iterations', 1000);
%
ind=find(Dtest); % keep track of index where data is zero
for J = 1
    s = RandStream('mt19937ar','Seed',J);
    RandStream.setGlobalStream(s);
    ind1=randperm(length(ind));
    % data Initilization
    ind2=ind1(1:floor(length(ind1)/2));
    prob.Omega=ind(ind2);
    ind3 = setdiff(ind1,ind2);
    ind3 = ind(ind3);
    b=Dtest(prob.Omega);
    prob.n1 = userid;
    prob.n2 = movie;
    [prob.Omega_i, prob.Omega_j] = ind2sub([prob.n1,prob.n2], prob.Omega);
    prob.m = length(prob.Omega);
   
    clear vars ind1 ind2
    switch(model)
        case{'partXY'}
            params.afun = @(L,R)partXY(L', R',prob.Omega_i, prob.Omega_j,prob.m)';
        case{'XonOmega'}
            params.afun = @(L,R)XonOmega(L, R, prob.Omega_i,prob.Omega_j);
        otherwise
            error('Unknown model selected');
    end
    params.afunT = @(x)sparse(prob.Omega_i, prob.Omega_j,x*1,prob.n1,prob.n2,prob.m);
    
    % Function Handle
    params.numr = userid;
    params.numc = movie;
    params.ls = 1;
    params.mode=1;
    params.funForward = @NLfunForward_partXY;
    for k=1:length(Sigmafact)
        sigmafact=Sigmafact(k);
        for i=1:length(rank)
            params.nr=rank(i);
            [U,S,V] = svds(reshape(params.afunT(b),prob.n1,prob.n2), params.nr, 'L', opts);
            Linit=U*sqrt(S);
            Rinit=V*sqrt(S);
            xinit  = initfact*[vec(Linit);vec(Rinit)]; % Initial guess
            clearvars Linit Rinit U S V
            tau = norm(xinit,1);
            sigma=sigmafact*norm(b,2);
            opts.funPenalty = @funLS;
            tic
            [xLS,r1,g1,info] = spgl1General(@NLfunForward_partXY,b,tau,sigma,xinit,options,params);
            toc

            e = params.numr*params.nr;
            L1 = xLS(1:e);
            R1 = xLS(e+1:end);
            L1 = reshape(L1,params.numr,params.nr);
            R1 = reshape(R1,params.numc,params.nr);
            prob1.n1 = userid;
            prob1.n2 = movie;
            [prob1.Omega_i, prob1.Omega_j] = ind2sub([prob1.n1,prob1.n2], ind3);
            prob1.m = length(ind3);
            params1.afun = @(L,R)partXY(L', R',prob1.Omega_i, prob1.Omega_j,prob1.m)';
            Drec = params1.afun(L1,R1);
            Dorig=vec(Dinit);
            Dorigtest=Dorig(ind3);
            NNorm(J,k,i) = 0.5 * (norm(L1,'fro')^2 + norm(R1,'fro')^2);
            SNR(J,k,i) = -20*log10(norm(Dorigtest-Drec,'fro')/norm(Dorigtest,'fro'))
            RMSE(J,k,i)=sqrt(sum((Dorigtest-Drec).^2)/length(ind3));
            Tau(J,k,i)= info.tau;
            iter(J,k,i)=info.iter;
            time(J,k,i) = info.timeTotal;
            clearvars L1 R1 xls Drec Dorig Dorigtest
        end
    end
end
save('Movielens_Cleverinit.mat','SNR','RMSE','Tau','iter','time','NNorm');
