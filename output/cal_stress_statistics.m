function [mu, variance, skewness, kurtosis] = cal_stress_statistics(Sigma,f)
% Revised Eqs. (21)-(23), with the dilute identification rho*A0 ~= f.
% Rows: S11, S22, S12. Columns: porosities. Zero-f skewness/kurtosis undefined.
Sigma=Sigma(:); f=f(:).';
assert(numel(Sigma)==3 && all(isfinite(Sigma)) && Sigma(3)==0,...
    'Sigma must be [Sigma1; Sigma2; 0] for the present biaxial model.');
assert(all(isfinite(f)) && all(f>=0 & f<1),'Require 0 <= f < 1.');
kn=cal_stress_coef(Sigma(1),Sigma(2));
K2=kn(:,2); K3=kn(:,3); K4=kn(:,4);
mu=Sigma./(1-f);
variance=K2.*f./(1-f).^2;
skewness=nan(3,numel(f)); kurtosis=skewness;
valid=f>0;
for i=1:3
    if K2(i)>0
        skewness(i,valid)=(K3(i)/K2(i)^1.5)./sqrt(f(valid));
        kurtosis(i,valid)=3+(K4(i)/K2(i)^2)./f(valid);
    end
end
end
