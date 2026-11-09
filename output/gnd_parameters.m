function [alpha,beta,pdf,cdf] = gnd_parameters(x,mu,variance,kurtosis)
% Same symmetric GND family; solve ordinary kurtosis, not excess kurtosis.
validateattributes(mu,{'numeric'},{'scalar','real','finite'});
validateattributes(variance,{'numeric'},{'scalar','real','finite','positive'});
if nargin==3
    beta=2;
else
    validateattributes(kurtosis,{'numeric'},{'scalar','real','finite','>=',3});
    if abs(kurtosis-3)<1e-12
        beta=2;
    else
        eq=@(b) gammaln(5./b)+gammaln(1./b)-2*gammaln(3./b)-log(kurtosis);
        lo=1;
        while eq(lo)<0, lo=lo/2; end
        beta=fzero(eq,[lo,2]);
    end
end
alpha=exp(0.5*(log(variance)+gammaln(1/beta)-gammaln(3/beta)));
pdfh=@(s) exp(log(beta)-log(2)-log(alpha)-gammaln(1/beta) ...
    -(abs((s-mu)/alpha)).^beta);
cdfh=@(s) 0.5+0.5*sign(s-mu).*gammainc((abs((s-mu)/alpha)).^beta,1/beta);
if isempty(x), pdf=pdfh; cdf=cdfh; else, pdf=pdfh(x); cdf=cdfh(x); end
end
