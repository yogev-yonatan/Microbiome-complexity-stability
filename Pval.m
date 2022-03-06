function [p,beta,gamma]=Pval(n,Ds,shuf_num,gamma_i,beta_i)
%Calculate the statisctical significance of complexity-stabilty pattern
%Yogev Yonatan, yogev.yn@gmail.com

%This function takes data points n and D^2, stacked in the vectors n (number of species) and
%Ds (squared DOC slope), respectively.
%gamma_i and beta_i are the parameters space that the algorithm checks.

%The funcrion returns:
%p - the empirical p value associated with the pattern.
%beta, - the beta and gamma values, assiciated with the best fitting curve.


min_ratio=.9;
l=length(n);
[n_area,D_area]=meshgrid(n,Ds);
n_area=n_area(:);
D_area=D_area(:);

R=nan(length(gamma_i),length(beta_i));
for i=1:length(beta_i)
    A=beta_i(i);
    for j=1:length(gamma_i)
        a=-gamma_i(j);
        r=sum(Ds<A*n.^a);
        if r/l>min_ratio
        area=sum(D_area<A*n_area.^a);
        R(j,i)=r/area;
        end
    end
end
[Rt,j]=max(R);
[statistic,i]=max(Rt);
beta=beta_i(i);
gamma=-gamma_i(j(i));

statistic_sh=nan(1,shuf_num);
for k=1:shuf_num
    Ds=Ds;
    ns=n(randperm(l));
    Rs=nan(length(gamma_i),length(beta_i));
    for i=1:length(beta_i)
        A=beta_i(i);
        for j=1:length(gamma_i)
            a=-gamma_i(j);
            r=sum(Ds<A*ns.^-a);
            if r/l>min_ratio
            area=sum(D_area<A*n_area.^a);
            Rs(j,i)=r/area;
            end            
        end
    end
    Rs=max(Rs);
    statistic_sh(k)=max(Rs);
end
p=sum(statistic_sh>statistic)./shuf_num;
end