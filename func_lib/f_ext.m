function y=f_ext(x)

abs_tol=5e-7;

Bn=[1 1/2 1/12 -1/720 1/30240 -1/1209600];
f1=@(t) Bn(1)+Bn(2).*t;
% f2=@(x) Bn(1)+Bn(2).*x+Bn(3).*x.^2;
% f4=@(x) Bn(1)+Bn(2).*x+Bn(3).*x.^2+Bn(4).*x.^4;
% f6=@(x) Bn(1)+Bn(2).*x+Bn(3).*x.^2+Bn(4).*x.^4+Bn(5).*x.^6;
% f8=@(x) Bn(1)+Bn(2).*x+Bn(3).*x.^2+Bn(4).*x.^4+Bn(5).*x.^6+Bn(6).*x.^8;
f_o=@(t) t./(1-exp(-t));

y=f_o(x);
x_id=abs(x)<abs_tol;
x_sp=x(x_id);
if sum(x_id)~=0
    y(x_id)=f1(x_sp);
end
end