function [r_tilta,test_outlier]=test_residual(r_cap,sigma2_cap_r_cap,nu,alpha)

t = tinv((1-alpha/2),nu-1); % t is a value of the t-distribution
kesi_tau = t.*sqrt(nu)./sqrt(nu-1+t.^2);

% kesi_tau=0.8*(1.9538-1.9527)+1.9527;
m=size(r_cap,2);
test_outlier=zeros(m,1);
r_tilta=r_cap./sqrt(sigma2_cap_r_cap);
for k=1:m
    if r_tilta(k)>-kesi_tau & r_tilta(k)<kesi_tau
        test_outlier(k)=1;
    end    
    end
end


