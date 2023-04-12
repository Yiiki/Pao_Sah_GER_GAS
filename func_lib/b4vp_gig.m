function sol=b4vp_gig(vgs0,vds0,ini,N,num)
if isempty(gcp('nocreate'))
    parpool(num)
end
vgs=vgs0;
vds=vds0;

spot=linspace(0,vgs,N);

tag=ones(1,length(spot));

sol_ini(length(tag))=main_ode(0,1e-3);
% Set a couple of warnings to temporarily issue errors (exceptions)
s = warning('error','MATLAB:bvp4c:RelTolNotMet');
pctRunOnAll warning('error','MATLAB:bvp4c:RelTolNotMet')
parfor i=1:length(tag)
   fprintf('\n at i=%d th loop\n',i) 
   try
       warning on verbose
       sol_ini(i)=main_ode(spot(i),vds,ini);
   catch ME
   % Exception-handling part
    tag(i)=0;
    switch ME.identifier
      case 'MATLAB:bvp4c:RelTolNotMet'
         fprintf('warning, RelTol Can''t Met\n');

      case 'MATLAB:bvp4c:SingJac'
         fprintf('error, Singular Jac\n');

      otherwise
         rethrow(ME)
    end
   end
       
end
warning(s);
sol.tag=tag;
sol.ini=sol_ini(tag==1);
sol.vgs=spot(tag==1);

end