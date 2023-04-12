function xi=ap_sol(a,b,xi_n,xi_p)
    xi=((max(0,xi_p)>xi_n).*a.*xi_n+(min(0,xi_n)<xi_p).*b.*xi_p)...
        ./(1+(max(0,xi_p)>xi_n).*a+(min(0,xi_n)<xi_p).*b);
end