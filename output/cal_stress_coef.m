function [kn,corr,corr2] = cal_stress_coef(S1,S2)
% Same single-pore coefficients as manuscript Eq. (17); no phase factors here.
% Numeric polynomials avoid loading a coefficient file from the wrong folder.
validateattributes(S1,{'numeric'},{'scalar','real','finite'});
validateattributes(S2,{'numeric'},{'scalar','real','finite'});
kn=zeros(3,4);
kn(:,2)=[5*S1^2/4-S1*S2+S2^2/4;...
    S1^2/4-S1*S2+5*S2^2/4;(S1^2+S2^2)/4];
kn(1,4)=335*S1^4/448-323*S1^3*S2/280+797*S1^2*S2^2/1120 ...
    -57*S1*S2^3/280+51*S2^4/2240;
kn(2,4)=51*S1^4/2240-57*S1^3*S2/280+797*S1^2*S2^2/1120 ...
    -323*S1*S2^3/280+335*S2^4/448;
kn(3,4)=51*S1^4/2240+3*S1^3*S2/140+41*S1^2*S2^2/1120 ...
    +3*S1*S2^3/140+51*S2^4/2240;
if nargout>1
    S=(S1+S2)/2; D=(S1-S2)/2;
    corr=diag(kn(:,2));
    corr(1,2)=(-S^2+3*D^2)/2; corr(2,1)=corr(1,2);
end
if nargout>2
    % h_i^2*h_j^2/u^2 is a polynomial of degree <=6 in u; four-point
    % Gauss-Legendre and 32 angular nodes integrate it exactly to roundoff.
    nodes=[0.0694318442029737,0.330009478207572,0.669990521792428,0.930568155797026];
    weights=[0.173927422568727,0.326072577431273,0.326072577431273,0.173927422568727];
    theta=(0:31)*2*pi/32; corr2=zeros(3);
    for k=1:4
        u=nodes(k); B=D*u*(3*u-2);
        h=[-(S+2*D)*u*cos(2*theta)+B*cos(4*theta);...
            (S-2*D)*u*cos(2*theta)-B*cos(4*theta);...
            -S*u*sin(2*theta)+B*sin(4*theta)];
        h2=h.^2;
        corr2=corr2+weights(k)*(h2*h2.')/(32*u^2);
    end
end
end
