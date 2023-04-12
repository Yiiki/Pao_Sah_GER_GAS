function y=gausf(x)
% for BP-device in qflps_model_v1 paper

beta =[7.3534
    5.1961
    1.2351
    2.5424];
H=beta(1);
Vgs0=beta(2);
bet=beta(3);
eta0=beta(4);
% Gauss
y=eta0+H.*exp(-(x-Vgs0).^2.*bet.^-2);
end