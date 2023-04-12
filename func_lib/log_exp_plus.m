function y=log_exp_plus(x)
% % old version code
% y=x-x;
% for i=1:length(x)
%     if x(i)>0
%         if isinf(exp(x(i)))
%             y(i)=x(i);
%         else
%         % original form
%         y(i)=log(exp(x(i))+1);
%         end
%     else
%         % Generalized Puiseux Series
%         y(i)=exp(x(i))-exp(2.*x(i))./2+exp(3.*x(i))./3-exp(4.*x(i))./4+exp(5.*x(i))./5-exp(6.*x(i))./6+exp(7.*x(i))./7-exp(8.*x(i))./8;
%     end
% end

% numerical limitation is set by the exp(-inf)
bdn=-6;
y=(x<bdn).*log_exp_plus_expand(x)+(x>=bdn).*log(exp(x)+1);
y(isnan(y))=x(isnan(y));

    function y=log_exp_plus_expand(x)
        y=exp(x)-exp(2.*x)./2+exp(3.*x)./3-exp(4.*x)./4+exp(5.*x)./5-exp(6.*x)./6+exp(7.*x)./7-exp(8.*x)./8;
    end
end