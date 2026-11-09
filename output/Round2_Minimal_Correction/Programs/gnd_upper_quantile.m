function q = gnd_upper_quantile(mu,variance,kurtosis,alpha_tail)
% q satisfies P(X>q)=alpha_tail. Inverse upper gamma avoids 1-CDF subtraction.
validateattributes(alpha_tail,{'numeric'},{'scalar','>',0,'<',1});
[a,b]=gnd_parameters([],mu,variance,kurtosis);
if alpha_tail==0.5
    q=mu;
elseif alpha_tail<0.5
    q=mu+a*gammaincinv(2*alpha_tail,1/b,'upper')^(1/b);
else
    q=mu-a*gammaincinv(2*(1-alpha_tail),1/b,'upper')^(1/b);
end
end
