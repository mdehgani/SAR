function [coeff_all,e_cap,sigma0_2_cap,sigma2_cap_e_cap,sigma2_x_cap]=adjustment_LS(t,y,pl,sigma02)

% Model: 'a*sin(2*pi*x+b)+c*x_d'  u = 4
 u=4; %Number of unknowns
Qy=sigma02*inv(pl);

%*************************************************************
A=[sin(2*pi*t) cos(2*pi*t) t ones(length(t),1)];
X=(inv(A'*inv(Qy)*A))*(A'*inv(Qy)*y);
%*************************************************************
coeff_all(1)=abs(sqrt(X(1)^2+X(2)^2));
if X(1)>0 && X(2)>0
    coeff_all(2)=atan(abs((X(2)/X(1))));
end
if X(2)>0 && X(1)<0
        coeff_all(2)=2*pi-atan(abs((X(2)/X(1))));
end
if X(2)<0 && X(1)>0
        coeff_all(2)=pi-atan(abs((X(2)/X(1))));
end
if X(2)<0 && X(1)<0
            coeff_all(2)=pi+atan(abs((X(2)/X(1))));
end
% X(2)
% X(1)
coeff_all(3)=X(3);
coeff_all(4)=X(4);

y_model_positive=coeff_all(1)*sin(2*pi*t+coeff_all(2))+coeff_all(3)*t+coeff_all(4);
res_positive=y-y_model_positive;
res_abs_positive=res_positive'*res_positive;

y_model_negative=-coeff_all(1)*sin(2*pi*t+coeff_all(2))+coeff_all(3)*t+coeff_all(4);
res_negative=y-y_model_negative;
res_abs_negative=res_negative'*res_negative;

if res_abs_negative<res_abs_positive
    coeff_all(1)=-coeff_all(1);
end
%*************************************************************
e_cap=y-A*X;
df=length(t)-u;
% sigma0_2_cap=e_cap'*(inv(Qy))*e_cap/(df);
sigma0_2_cap=e_cap'*(pl)*e_cap/(df);


Qx_cap=sigma0_2_cap*inv(A'*inv(Qy)*A);
sigma2_x_cap=diag(Qx_cap);

y_cap=A*X;
Qy_cap=A*Qx_cap*A';

Q_cap_y=sigma0_2_cap*inv(pl);
Qe_cap=Q_cap_y-Qy_cap;
% Qe_cap=Qy-Qy_cap;

sigma2_cap_e_cap=diag(Qe_cap);


    
