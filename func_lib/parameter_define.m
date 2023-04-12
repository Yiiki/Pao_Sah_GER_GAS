function par=parameter_define()
% device parameters
par.Vt=0.0258549;% kT/q, 1 for T=300K, 0.1 for 30K, etc.
par.g0=0.00258031*1.0;% Cox/Cq, 1 for 20nm, 0.1 for 200nm, etc.
par.vfb=0;% flat band voltage

% SRH recombination
par.qL2m=1.60217662/9.10938356;% L= 1um, q*L^2/m0, mo free electrons mass

% materials parameters
par.mue=50e-4;% mobility, unit m^2/(V.s)
par.muh=10e-3;

% BP, armchair

% par.me=0.15;% electrons mass
% par.mh=0.14;% holes mass
% par.eg=0.69;% eg=psb+psc
% par.psb=0.4;% fermi potential for electrons
% par.psc=0.29;% fermi potential for holes

% WSe2

par.me=0.35;% electrons mass
par.mh=0.46;% holes mass
par.eg=1.62;% eg=psb+psc
par.psb=1.62*0.6;% fermi potential for electrons
par.psc=1.62-1.62*0.6;% fermi potential for holes

end