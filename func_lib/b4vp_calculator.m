function sol=b4vp_calculator(vds0,vgs0,ini)

% calculate the drain current nearly at given vgs0 and vds0
% rely on the iteration initial solution strategy
% when vds < 0.1, dv = 0.1.*vds
% when vds > 0.1, dv = 0.01

vds=vds0;
vgs=vgs0;

if abs(vds) > 0.1
    dv0=0.01;
else
    dv0=0.1.*abs(vds);
end

if abs(vgs) > dv0
    dv=dv0.*sign(vgs);
else 
    dv=vgs;
end

imax=ceil(vgs./dv)+1;

vgs_list=linspace(0,vgs,imax);
dv=vgs_list(2)-vgs_list(1);

vg=zeros(1,imax);
% vg(1)=0;

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
                sol_vec(1)=main_ode(vg(1),vds);
            elseif vds>4
                vds_a=4;
                vds_b=vds;
                % % vds need guide
                sol_aid=main_ode(vg(1),4);
                sol_vec(1)=main_ode(vg(1),vds,sol_aid);
            end
        else
            if vds_a==0
                sol_aid=main_ode(vg(1),vds_b);
                vds_a=vds_b;
                vds_b=vds;
                sol_vec(1)=main_ode(vg(1),vds,sol_aid);
            else
                sol_aid=main_ode(vg(1),vds_b,sol_aid);
                vds_a=vds_b;
                vds_b=vds;
                sol_vec(1)=main_ode(vg(1),vds,sol_aid);
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

% preallocating structure array
sol_vec(imax)=sol_vec(1);


start_tag=2;
final_tag=imax;

while abs(vg(end))<abs(vgs)

for i=start_tag:final_tag
    fprintf('At %d / %d th loop, vgs=%f, dv=%f\n',i,imax,vg(i-1),dv)
    % exception state
    state=1;
while state==1
   state=0;
try
   % Regular processing part
   sol_gen=main_ode(vg(i-1)+dv,vds,sol_vec(i-1));
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
vg(i)=vg(i-1)+dv;
end

if abs(vg(end))<abs(vgs)
    
    del_max=ceil(abs(vgs-vg(end))./abs(dv));
    
    nor_list=linspace(vg(end),vgs,del_max+1);
    
    dv=nor_list(2)-nor_list(1);
    
    imax=imax+del_max;
    
    vg(imax)=0;
    sol_vec(imax)=sol_vec(1);
    start_tag=final_tag+1;
    final_tag=imax;
end

end

% prevent from vg surmount the vgs preallocated
% it is redundant-design.
if abs(vg(end))>abs(vgs)
    sol_vec(end)=main_ode(vgs,vds,sol_vec(end));
end

warning(s);
% sol.vg=vg(end);
solp_set=sol_post(sol_vec(end),vg(end),vds);
% sol=mean(solp_set.in+solp_set.ip);
sol.cur=solp_set.in+solp_set.ip;
sol.pos=solp_set;
sol.ini=sol_vec(end);
end