function xi=sol_bnd(a,b,xi_n,xi_p)
    xi=fzero(@(x)bdeq(a,b,xi_n,xi_p,x),ap_sol(a,b,xi_n,xi_p));
end