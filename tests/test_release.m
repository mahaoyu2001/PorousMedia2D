function test_release()
% Validate corrected moments, independent kernel quadrature, GND, field IO,
% and consistency with the bundled FEM/Campbell reference comparisons.
root=fileparts(fileparts(mfilename('fullpath')));
oldPath=path;cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'output'),'-begin');
[mu,v,s,k]=cal_stress_statistics([0;1;0],[0 .005 .1]);
assert(abs(mu(2,3)-1/.9)<1e-13);
assert(abs(v(2,3)-.154320987654321)<1e-13);
assert(abs(k(2,3)-7.785714285714286)<1e-12);
assert(all(v(:,1)==0) && all(isnan(k(:,1))) && all(isnan(s(:,1))));
assert(all(s(:,2:3)==0,'all'));
rejected=false;
try,cal_stress_statistics([0;1;1],.01);catch,rejected=true;end
assert(rejected,'Unsupported remote shear should be rejected.');
% Independent angular/radial midpoint integration, not polynomial evaluation.
theta=((0:511)+.5)*2*pi/512;u=((0:8191)'+.5)/8192;
for loadCase=[0 1; 1 0; .4 1; -1 1; 1 1]'
    p=sum(loadCase)/2;q=diff(flipud(loadCase))/2;
    B=q*u.*(3*u-2);
    h={-(p+2*q)*u*cos(2*theta)+B*cos(4*theta), ...
       (p-2*q)*u*cos(2*theta)-B*cos(4*theta), ...
       -p*u*sin(2*theta)+B*sin(4*theta)};
    coefficients=cal_stress_coef(loadCase(1),loadCase(2));
    for i=1:3
        numeric2=mean(mean(h{i}.^2,2)./u.^2);
        numeric4=mean(mean(h{i}.^4,2)./u.^2);
        assert(abs(numeric2-coefficients(i,2))<1e-6);
        assert(abs(numeric4-coefficients(i,4))<2e-6);
    end
end
for f=[.005 .05 .2]
    [mu,v,~,k]=cal_stress_statistics([0;1;0],f);
    for i=1:3
        [a,b]=gnd_parameters([],mu(i),v(i),k(i));
        variance=exp(2*log(a)+gammaln(3/b)-gammaln(1/b));
        kurt=exp(gammaln(5/b)+gammaln(1/b)-2*gammaln(3/b));
        assert(abs(variance/v(i)-1)<1e-10 && abs(kurt/k(i)-1)<1e-8);
        quantile=gnd_upper_quantile(mu(i),v(i),k(i),.01);
        assert(abs(gnd_survival(quantile,mu(i),v(i),k(i))-.01)<1e-10);
    end
end
% Match the published 0.5% S22 > 2 result from a separate reference export.
tableDir=fullfile(root,'output','Round2_Tail_Validation','Results','Tables');
t=readtable(fullfile(tableDir,'R5Q4_Exceedance_All.csv'));
r=t(t.Porosity==.005 & t.Component==2 & t.Threshold==2,:);
assert(height(r)==1 && abs(r.Campbell_Full/r.FEM_Mean-1)<.002);
assert(abs(r.GND_Analytical-gnd_survival(2,1/.995,1.25*.005/.995^2,3+(335/448)/(1.25^2*.005)))<1e-12);
% Headerless legacy exports and commented Abaqus exports must agree.
testDir=fullfile(fileparts(mfilename('fullpath')),'Results');
if ~isfolder(testDir),mkdir(testDir);end
file=fullfile(testDir,'field_reader_fixture.txt');
for header={'','# S11 S22 S12 EVOL'}
    fid=fopen(file,'wt');assert(fid>0);
    if ~isempty(header{1}),fprintf(fid,'%s\n',header{1});end
    fprintf(fid,'-1 2 0 0.25\n1 3 0.5 0.75\n');fclose(fid);
    assert(isequal(read_stress_field(file),[-1 2 0 .25;1 3 .5 .75]));
end
delete(file);
% Parse every shipped MATLAB source without executing Abaqus or large fields.
files=dir(fullfile(root,'**','*.m'));
for j=1:numel(files)
    messages=checkcode(fullfile(files(j).folder,files(j).name),'-id');
    if ~isempty(messages)
        bad=contains(string({messages.message}), ...
            ["Parse error","Invalid expression","Unexpected MATLAB expression","END is missing"],'IgnoreCase',true);
        assert(~any(bad),'MATLAB parse error in %s',files(j).name);
    end
end
fprintf('PorousMedia2D v2.0.0 checks passed: moments, kernel quadrature, GND, tail references, IO and syntax.\n');
end
