function x_sol=interp_sol(a_vgs,a_vds,a_sol, ...
                          b_vgs,b_vds,b_sol, ...
                          x_vgs,x_vds)
% initial
x_sol=b_sol;
% re-write point a
a_sol_x=b_sol.x;
a_sol_y=interp1(a_sol.x,a_sol.y',a_sol_x)';
a_sol_yp=interp1(a_sol.x,a_sol.yp',a_sol_x)';
% calculate point x
x_sol_y=x_sol.y;
x_sol_yp=x_sol.yp;

% decide predict dim
err_cond=abs((a_vds-x_vds).*(b_vgs-x_vgs)-(a_vgs-x_vgs).*(b_vds-x_vds))<1e-12;
ovr_cond=(sum(sum((a_sol_y-b_sol.y).^2)+sum((a_sol_yp-b_sol.yp).^2))==0);
if (a_vds==b_vds)&&(a_vgs~=b_vgs)
    dim=1;
elseif (a_vds~=b_vds)&&((a_vgs==b_vgs)||err_cond)
    dim=2;% assume vds-distinguished by default
elseif (a_vds~=b_vds)&&(a_vgs~=b_vgs)
    error('non-colinear initials provided, a(vgs,vds), b(vgs,vds), x(vgs,vds)')
elseif (a_vds==b_vds)&&(a_vgs==b_vgs)&&ovr_cond
    dim=3;% two are same.
else
    error('inconsistent initials provided, vds and vgs same while solutions not.')
end

switch dim
    case 1
    
    for i=1:x_sol.stats.nmeshpoints
        x_sol_y(:,i)=interp1([a_vgs;b_vgs],[a_sol_y(:,i),b_sol.y(:,i)]',x_vgs, 'linear', 'extrap')';
        x_sol_yp(:,i)=interp1([a_vgs;b_vgs],[a_sol_yp(:,i),b_sol.yp(:,i)]',x_vgs, 'linear', 'extrap')';
    end

    case 2
    
    for i=1:x_sol.stats.nmeshpoints
        x_sol_y(:,i)=interp1([a_vds;b_vds],[a_sol_y(:,i),b_sol.y(:,i)]',x_vds, 'linear', 'extrap')';
        x_sol_yp(:,i)=interp1([a_vds;b_vds],[a_sol_yp(:,i),b_sol.yp(:,i)]',x_vds, 'linear', 'extrap')';
    end

    case 3
    
    for i=1:x_sol.stats.nmeshpoints
        x_sol_y(:,i)=b_sol.y(:,i);
        x_sol_yp(:,i)=b_sol.yp(:,i);
    end

end

x_sol.y=x_sol_y;
x_sol.yp=x_sol_yp;
end