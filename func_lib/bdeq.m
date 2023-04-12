function res=bdeq(a,b,xi_n,xi_p,xi)
   res=a.*log_exp_plus(xi-xi_n)-b.*log_exp_plus(xi_p-xi)+xi;
end