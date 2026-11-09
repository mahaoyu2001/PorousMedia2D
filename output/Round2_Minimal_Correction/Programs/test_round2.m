function test_round2()
p=round2_paths();
assert(strcmp(which('cal_stress_statistics'),fullfile(p.programDir,'cal_stress_statistics.m')),...
    'MATLAB resolved a different copy of cal_stress_statistics. Use Start_Round2.');
assert(strcmp(which('gnd_parameters'),fullfile(p.programDir,'gnd_parameters.m')),...
    'MATLAB resolved a different copy of gnd_parameters.');
[mu,v,s,k]=cal_stress_statistics([0 1 0],[0 0.005 0.1]);
assert(isequal(size(v),[3 3]) && abs(v(2,3)-0.154320987654321)<1e-13);
assert(abs(k(2,3)-7.785714285714286)<1e-12 && s(2,3)==0);
assert(all(v(:,1)==0) && all(isnan(s(:,1))) && all(isnan(k(:,1))));
assert(abs(mu(2,3)-1/0.9)<1e-13);
loads=[0 1;1 0;1 1;-1 1;0.4 1];
for j=1:size(loads,1)
    [kn,c,c2]=cal_stress_coef(loads(j,1),loads(j,2));
    assert(max(abs(diag(c)-kn(:,2)))<1e-12);
    assert(max(abs(diag(c2)-kn(:,4)))<1e-10*max(1,max(abs(kn(:,4)))));
    assert(all(kn(:,1)==0) && all(kn(:,3)==0));
end
for f=[0.0001 0.005 0.1 0.2]
    [mu,v,~,k]=cal_stress_statistics([0;1;0],f);
    for i=1:3
        [a,b]=gnd_parameters([],mu(i),v(i),k(i));
        assert(a>0 && b>0 && isfinite(a) && isfinite(b));
        reconstructedK=exp(gammaln(5/b)+gammaln(1/b)-2*gammaln(3/b));
        reconstructedV=exp(2*log(a)+gammaln(3/b)-gammaln(1/b));
        assert(abs(reconstructedK/k(i)-1)<1e-8);
        assert(abs(reconstructedV/v(i)-1)<1e-10);
        q=gnd_upper_quantile(mu(i),v(i),k(i),0.01);
        assert(abs(gnd_survival(q,mu(i),v(i),k(i))-0.01)<1e-10);
    end
end
assert(abs(gnd_upper_quantile(0,1,3,0.01)-2.32634787404084)<1e-10);
assert(abs(gnd_survival(2,0,1,3)-0.5*erfc(sqrt(2)))<1e-12);
assert(abs(gnd_survival(-2,0,1,3)+gnd_survival(2,0,1,3)-1)<1e-12);
for kind={'geom','load','vol20','point'}
    t=round2_read_statistics(kind{1});
    assert(all(isfinite(t.Variance) & t.Variance>0));
    assert(all(isfinite(t.Kurtosis) & t.Kurtosis>0));
    assert(all(t.Variance_SD>=0 & t.Kurtosis_SD>=0));
    fprintf('Original %s cache: %d rows.\n',kind{1},height(t));
end
[~,~,cachedV,~,~,~,~,~]=round2_cached_arrays('geom',[0.005 0.1],["S11","S22","S12"]);
assert(isequal(size(cachedV),[3 2]) && all(isfinite(cachedV(:))));
[~,~,pointV,~,~,~,~,~,counts]=round2_cached_point_arrays([0.005 0.1],...
    {'Hard-core','r_min=0'},["S11","S22","S12"]);
assert(isequal(size(pointV),[2 3 2]) && all(isfinite(pointV(:))) && all(counts(:)>0));
files=dir(fullfile(p.programDir,'*.m')); parseErrors=0;
for j=1:numel(files)
    messages=checkcode(fullfile(files(j).folder,files(j).name),'-id');
    if ~isempty(messages)
        bad=contains(string({messages.message}),...
            ["Parse error","Invalid expression","Unexpected MATLAB expression","END is missing"],'IgnoreCase',true);
        parseErrors=parseErrors+nnz(bad);
    end
end
assert(parseErrors==0,'Code analyzer reported MATLAB syntax errors.');
fprintf('Round-2 checks passed: coefficients, moments, GND inversion, source caches, and paths.\n');
end
