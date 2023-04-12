function sol=b4vp_dig(vds0,N,num)
if isempty(gcp('nocreate'))
    parpool(num)
end
vds=vds0;

vds_spot=linspace(0,vds,N);
spot=vds_spot(2:end);

tag=ones(1,length(spot));
tim=tag;

sol_ini(length(tag))=main_ode(0,1e-3);
% Set a couple of warnings to temporarily issue errors (exceptions)
s = warning('error','MATLAB:bvp4c:RelTolNotMet');
pctRunOnAll warning('error','MATLAB:bvp4c:RelTolNotMet')
parfor i=1:length(tag)
   fprintf('\n at i=%d th loop\n',i) 
   try
       warning on verbose
       tstart=tic;
       sol_ini(i)=main_ode(0,spot(i));
       tim(i)=toc(tstart);
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
sol.vds=spot(tag==1);
sol.tim=tim(tag==1);
end