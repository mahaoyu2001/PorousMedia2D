function probability = gnd_survival(s,mu,variance,kurtosis)
% Stable upper-tail probability; Gaussian is the kurtosis=3 case.
if nargin<4, kurtosis=3; end
[a,b]=gnd_parameters([],mu,variance,kurtosis);
q=gammainc(abs((s-mu)/a).^b,1/b,'upper');
probability=0.5*q;
left=s<mu;
probability(left)=1-0.5*q(left);
end
