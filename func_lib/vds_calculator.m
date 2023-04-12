function sol=vds_calculator(vds1,vds0,ini)

% calculate the drain current nearly at given vgs0 and vds0
% rely on the iteration initial solution strategy
% when vds < 0.1, dv = 0.1.*vds
% when vds > 0.1, dv = 0.01

if vds1~=vds0
vds=vds1;

if abs(vds1-vds0) >= 0.1
    dv=0.01;
else
    dv=0.1.*abs(vds1-vds0);
end

imax=ceil(abs(vds1-vds0)./dv)+1;

vds_list=linspace(vds0,vds1,imax);

dv=vds_list(2)-vds_list(1);

vd=vds0.*ones(1,imax);

% Set a couple of warnings to temporarily issue errors (exceptions)
s = warning('error','MATLAB:bvp4c:RelTolNotMet');

if nargin <3
    ini_tag=1;
    steta=0;

    while steta~=1
   steta=1;
   try
   % Regular processing part
   
        % initial solution for vgs=0 at vds = vds0
        if ini_tag==1
            ini_tag=0;
            if vds<=4
                vds_a=0;
                vds_b=vds;
                % vds is not strange
                sol_vec(1)=main_ode(vd(1),vds);
            elseif vds>4
                vds_a=4;
                vds_b=vds;
                % % vds need guide
                sol_aid=main_ode(vd(1),4);
                sol_vec(1)=main_ode(vd(1),vds,sol_aid);
            end
        else
            if vds_a==0
                sol_aid=main_ode(vd(1),vds_b);
                vds_a=vds_b;
                vds_b=vds;
                sol_vec(1)=main_ode(vd(1),vds,sol_aid);
            else
                sol_aid=main_ode(vd(1),vds_b,sol_aid);
                vds_a=vds_b;
                vds_b=vds;
                sol_vec(1)=main_ode(vd(1),vds,sol_aid);
            end
        end
   %
   catch ME
        % Exception-handling part
        switch ME.identifier
            case 'MATLAB:bvp4c:RelTolNotMet'
                fprintf('warning, RelTol Can''t Met at initial solution\n');
                steta=0;
            case 'MATLAB:bvp4c:SingJac'
                fprintf('error, Singular Jac at initial solution\n');
                steta=0;
            otherwise
                rethrow(ME)
        end
        % exception-handling
        if steta==0
            vds_b=0.5.*(vds_a+vds_b);
        end
   end
    end
else
    sol_vec(1)=ini;
end

if vds1>vds0
    bfwd=1;
else
    bfwd=-1;
end

% preallocating structure array
sol_vec(imax)=sol_vec(1);


start_tag=2;
final_tag=imax;

while bfwd*vd(end)<bfwd*vds

for i=start_tag:final_tag
    fprintf('\n cal: vds0=%f, --> %f --> vds1=%f, at %d / %d th loop with dv=%f\n',...
        vds0,vd(i-1),vds1,i,imax,dv)
    % exception state
    state=1;
while state==1
   state=0;
try
   % Regular processing part
   sol_gen=main_ode(0,vd(i-1)+dv,sol_vec(i-1));
%    sol3_ex=main_ode(0.7,1,sol2);
catch ME
   % Exception-handling part
   switch ME.identifier
      case 'MATLAB:bvp4c:RelTolNotMet'
         fprintf('warning, RelTol Can''t Met\n');
         state=1;
      case 'MATLAB:bvp4c:SingJac'
         fprintf('error, Singular Jac\n');
         state=1;
      otherwise
         rethrow(ME)
   end
   % exception-handling
   if state==1
       dv=0.5.*dv;
   end
end
end
sol_vec(i)=sol_gen;
vd(i)=vd(i-1)+dv;
end

if bfwd*vd(end)<bfwd*vds
    
    del_max=ceil(abs(vds-vd(end))./abs(dv));
    
    nor_list=linspace(vd(end),vds,del_max+1);
    
    dv=nor_list(2)-nor_list(1);
    
    imax=imax+del_max;
    
    vd(imax)=0;
    sol_vec(imax)=sol_vec(1);
    start_tag=final_tag+1;
    final_tag=imax;
end

end

% prevent from vg surmount the vgs preallocated
% it is redundant-design.
if bfwd*vd(end)>bfwd*vds
    sol_vec(end)=main_ode(0,vds,sol_vec(end));
end

warning(s);

sol=sol_vec(end);
else
    sol=ini;
end
end