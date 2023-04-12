function y=log_exp_minus(x)
for i=1:length(x)
    if x(i)>1
        % original form
        y(i)=log(exp(x(i))-1);
    else
        % Generalized Puiseux Series
        y(i)=log(x(i))+x(i)./2+x(i).^2./24-x(i).^4./2880+x(i).^6./181440-x(i).^8./9676800+x(i).^10./479001600;
    end
end
end