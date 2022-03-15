function [Xfinal,A,r] = get_GLV_tab(n,m,sigma,C,minmax)

% This function genrates a GLV model and integrates alternative m steady states 
% Yogev Yonatan, yogev.yn@gmail.com
% 
% Inputs:
% n - number of species in the model
% m - number of alternative steady states
% sigma - the characteristic interactions' stregth
% C - interaction density matrix
% minimax - array of two numbers, represint the minimal and maximal
% fraction of species in every steadt states
% 
% Outputs:
% Xfinal - Table of m alternative steady states
% A and r - The generated interaction matrix (A) and growth rates (r)


%% NANsSearch ensures no NANs in the final correlation matrix.
NANsSearch = true(1);
% while NANsSearch

r = rand(n,1);

% [A,G] = getER(species_num,k);
A= rand(n);
A = A<C;
A = sigma*sprandn(A);
A(eye(n)==1)=-1;

%% set initial condions:

min_species = round(min(minmax)*n);
max_species = round(max(minmax)*n);
% exam_species = randi([min_species,max_species],1,exams_num);
X0 = zeros(n,m);

    for n= 1:m
        
        idxs= randperm(n);
        stop= randi([min_species,max_species]);
        X0(idxs(1:stop),n) = rand(1,stop);
    end
    
    %% create Xfinal:
    Xfinal=zeros(n,m);
    for i = 1:m
        
        lotka_voltera = @(t,x,A,r) x.*(r+A*x);
        f =@(t,x) lotka_voltera(t,x,A,r);
        
        [~,X] = ode45(f,[0 100],X0(:,i));
%         histogram(X(end,:),10)
        Xfinal(:,i) = X(end,:);
    end
    Xfinal(Xfinal<10^-10)=0;
%     NANsSearch = sum(sum(isnan(corr(Xfinal'))));
    
end