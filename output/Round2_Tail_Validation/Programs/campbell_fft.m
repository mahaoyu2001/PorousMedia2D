function out=campbell_fft(component,f,cfg)
% Full compensated Poisson exponent evaluated by Levy-measure quadrature + FFT.
% Only small-amplitude expm1 evaluation uses a converged series (not a global
% fourth-cumulant truncation). The strong-field exponent retains exp(i*t*h).
% u=a^2/r^2; eta=rho*pi*a^2=f, consistent with revised dilute moment formulas.
arguments
 component (1,1) double {mustBeMember(component,[1 2 3])}
 f (1,1) double {mustBePositive}
 cfg struct
end
p=tail_paths();name=sprintf('%s_S%d_f%04d',cfg.Name,component,round(f*10000));
cache=fullfile(p.cacheDir,[name '.mat']);
if isfile(cache)
 saved=load(cache,'out','savedCfg');
 if isequal(saved.savedCfg,cfg),out=saved.out;return,end
end
timer=tic;N=cfg.N;L=cfg.L;dx=L/N;eta=f;scale=1/(1-f);
theta=((0:cfg.Ntheta-1)+.5)*(2*pi/cfg.Ntheta);
[v,w]=composite_gauss(0,log(1/cfg.Ucut),cfg.Nradial,8);
u=exp(-v);nu=zeros(N,1);S=.5;D=-.5;
for start=1:64:numel(u)
 ix=start:min(start+63,numel(u));uu=u(ix);bb=D*uu.*(3*uu-2);
 if component==1
  h=-(S+2*D)*uu*cos(2*theta)+bb*cos(4*theta);
 elseif component==2
  h=(S-2*D)*uu*cos(2*theta)-bb*cos(4*theta);
 else
  h=-S*uu*sin(2*theta)+bb*sin(4*theta);
 end
 h=h*scale;weights=repmat((w(ix)./uu)/cfg.Ntheta,1,cfg.Ntheta);
 z=h(:)/dx;left=floor(z);fraction=z-left;
 nu=nu+accumarray(mod(left,N)+1,weights(:).*(1-fraction),[N 1],@sum,0)...
       +accumarray(mod(left+1,N)+1,weights(:).*fraction,[N 1],@sum,0);
end
nu=eta*nu;
x=(0:N-1)'*dx;x(x>=L/2)=x(x>=L/2)-L;
t=[0:N/2-1 -N/2:-1]'*(2*pi/L);
firstMoment=sum(nu.*x);
% Weak-field moments integrate the original exponent on 0<u<Ucut.
% Choose enough points to integrate every retained polynomial moment exactly
% up to floating-point error, including the 40-term high-frequency check.
weakRadialOrder=max(24,cfg.SeriesOrder+4);
weakAngularOrder=2^nextpow2(4*cfg.SeriesOrder+1);
[un,uw]=gauss_nodes(weakRadialOrder);un=un*cfg.Ucut;uw=uw*cfg.Ucut;
th=((0:weakAngularOrder-1)+.5)*2*pi/weakAngularOrder;b=D*un.*(3*un-2);
if component==1,hh=-(S+2*D)*un*cos(2*th)+b*cos(4*th);
elseif component==2,hh=(S-2*D)*un*cos(2*th)-b*cos(4*th);
else,hh=-S*un*sin(2*th)+b*sin(4*th);end
hh=hh*scale;weakMom=zeros(1,cfg.SeriesOrder);
weak=zeros(N,1);term=ones(N,1);
for n=1:cfg.SeriesOrder
 term=term.*(1i*t/n);
 if n>=2
  weakMom(n)=eta*sum(uw.*mean(hh.^n,2)./un.^2);
  weak=weak+term*weakMom(n);
 end
end
exponent=N*ifft(nu)-sum(nu)-1i*t*firstMoment+weak;
assert(max(real(exponent))<1e-7,'Characteristic modulus exceeds one');
phi=exp(exponent);phi(1)=1;
mass=real(fft(phi))/N;negativeMass=-sum(mass(mass<0));
assert(negativeMass<1e-7 && abs(sum(mass)-1)<1e-9,'Invalid inversion probability mass');
mass=max(mass,0);mass=mass/sum(mass);xx=fftshift(x);pp=fftshift(mass);
cc=cumsum(pp)-.5*pp;
meanModel=[0 1 0];mu=meanModel(component)/(1-f);
thresholds=.5;if component==2,thresholds=[1.5 2 2.5];end
survival=zeros(size(thresholds));gp=survival;
positive=2:N/2;tt=t(positive);ph=phi(positive);
for j=1:numel(thresholds)
 z=thresholds(j)-mu;
 survival(j)=1-interp1(xx,cc,z,'linear');
 gp(j)=.5+(2/L)*(-.5*z+sum(imag(exp(-1i*tt*z).*ph)./tt));
end
iq=find(cc>=.99,1);q=interp1(cc(iq-1:iq),xx(iq-1:iq),.99)+mu;
z=q-mu;qCheck=.5+(2/L)*(-.5*z+sum(imag(exp(-1i*tt*z).*ph)./tt));
kn=cal_stress_coef(0,1);
k2=sum(nu.*x.^2)+weakMom(2);k3=sum(nu.*x.^3)+weakMom(3);k4=sum(nu.*x.^4)+weakMom(4);
target2=f*kn(component,2)/(1-f)^2;target4=f*kn(component,4)/(1-f)^4;
boundaryMass=sum(pp(abs(xx)>L/2-1));
phaseBound=max(abs(t))*scale*2.5*cfg.Ucut;
% Bound exp-series remainder integrated over u, using |h| <= 2.5*u/(1-f).
order=cfg.SeriesOrder+1;
remainderBound=eta*exp(phaseBound)*(max(abs(t))*scale*2.5)^order*cfg.Ucut^(order-1)/(factorial(order)*(order-1));
assert(remainderBound<1e-9,'Weak-field exponential evaluation requires additional terms');
out=struct('Component',component,'Porosity',f,'Mean',mu,'Thresholds',thresholds,...
 'Probability',survival,'Probability_GP',gp,'Q99',q,'Q99_GP_Probability',qCheck,...
 'Kappa2',k2,'Kappa3',k3,'Kappa4',k4,'TargetKappa2',target2,'TargetKappa4',target4,...
 'Kappa2_RelativeError',k2/target2-1,'Kappa4_RelativeError',k4/target4-1,...
 'NegativeProbabilityMass',negativeMass,'BoundaryProbabilityMass',boundaryMass,...
 'MaxGP_CDF_Difference',max(abs(gp-survival)),'SeriesRemainderExponentBound',remainderBound,...
 'ZeroMeanShiftApplied',firstMoment,'ElapsedSeconds',toc(timer),'ConfigName',cfg.Name,...
 'dx',dx,'N',N,'L',L,'Ucut',cfg.Ucut,'Nradial',cfg.Nradial,'Ntheta',cfg.Ntheta);
assert(max(abs(gp-survival))<2e-6,'Gil-Pelaez vs FFT CDF discrepancy');
assert(abs(qCheck-.01)<2e-6,'Quantile verification failed');
assert(abs(k2/target2-1)<.002 && abs(k4/target4-1)<.002,'Recovered cumulants disagree');
savedCfg=cfg;save(cache,'out','savedCfg');
fprintf('%s: P=%s, K2/K4 rel.err %.3g / %.3g, %.1f s\n',name,mat2str(survival,7),out.Kappa2_RelativeError,out.Kappa4_RelativeError,out.ElapsedSeconds);
end

function [x,w]=composite_gauss(a,b,n,order)
assert(mod(n,order)==0);[z,ww]=gauss_nodes(order);m=n/order;edges=linspace(a,b,m+1);
x=zeros(n,1);w=x;
for k=1:m,ix=(k-1)*order+(1:order);width=edges(k+1)-edges(k);x(ix)=edges(k)+width*z;w(ix)=width*ww;end
end

function [x,w]=gauss_nodes(n)
k=(1:n-1)';beta=k./sqrt(4*k.^2-1);[V,D]=eig(diag(beta,1)+diag(beta,-1));
[z,order]=sort(diag(D));V=V(:,order);x=(z+1)/2;w=(V(1,:).^2)';
end
