function y=Df_ext(x)

abs_tol=5e-4;

Dn=[1/2 1/6 -1/180 1/5040 -1/151200];
Df1=@(t) Dn(1)+Dn(2).*t;
% f2=@(x) Bn(1)+Bn(2).*x+Bn(3).*x.^2;
% f4=@(x) Bn(1)+Bn(2).*x+Bn(3).*x.^2+Bn(4).*x.^4;
% f6=@(x) Bn(1)+Bn(2).*x+Bn(3).*x.^2+Bn(4).*x.^4+Bn(5).*x.^6;
% f8=@(x) Bn(1)+Bn(2).*x+Bn(3).*x.^2+Bn(4).*x.^4+Bn(5).*x.^6+Bn(6).*x.^8;
Df_o=@(t) (1-(1+t).*exp(-t))./(1-exp(-t)).^2;

y=Df_o(x);
x_id=abs(x)<abs_tol;
x_sp=x(x_id);
if sum(x_id)~=0
    y(x_id)=Df1(x_sp);
end
end