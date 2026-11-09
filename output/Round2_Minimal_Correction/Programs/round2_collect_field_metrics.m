function records = round2_collect_field_metrics(force)
% Read each original field once; compute all three stress components together.
% Binning, thresholds, natural-log JS and FEM||Model KL follow root scripts.
% Per-file caches include source size/time and corrected model moments.
if nargin<1, force=false; end
p=round2_paths(); cacheDir=fullfile(p.rootDir,'Results','Cache','Field_Distributions');
if ~isfolder(cacheDir),mkdir(cacheDir);end
models={'2DVoid_Vol05','2DVoid_Vol1','2DVoid_Vol2','2DVoid_Vol3',...
 '2DVoid_Vol5','2DVoid_Vol7_5','2DVoid_Vol10'};
porosity=[0.005 0.01 0.02 0.03 0.05 0.075 0.1]; suffixes={'a','b','c','d','e'};
records=struct([]); count=0; Sigma=[0;1;0];
test_fast_quantile();
for im=1:numel(models)
 [muModel,varModel,~,kurtModel]=cal_stress_statistics(Sigma,porosity(im));
 for ir=1:numel(suffixes)
  count=count+1; fileName=sprintf('%s_%s_S2.txt',models{im},suffixes{ir});
  source=fullfile(p.fieldDir,fileName); info=dir(source);
  assert(~isempty(info),'Missing FEM source %s',source);
  stamp=struct('Version',2,'Bytes',info.bytes,'Datenum',info.datenum,...
   'Moments',[muModel,varModel,kurtModel],'Quantile',0.99,...
   'TailFactor',6,'MinBinCount',250,'EpsProb',1e-12);
  cache=fullfile(cacheDir,[erase(fileName,'.txt') '.mat']);
  if ~force && isfile(cache)
   saved=load(cache,'stamp','r');
   if isequaln(saved.stamp,stamp)
    r=saved.r; records=[records;r(:)]; %#ok<AGROW>
    fprintf('[%d/35] Cached: %s\n',count,fileName);continue
   end
  end
  timer=tic; fprintf('[%d/35] Reading %.1f MB: %s\n',count,info.bytes/1e6,fileName);
  Data=read_stress_field(source);
  assert(size(Data,2)==4 && ~isempty(Data) && all(isfinite(Data(:))),...
   'Expected four finite columns: S11 S22 S12 integration-area.');
  weights=Data(:,4);assert(all(weights>0));weights=weights/sum(weights);
  fprintf('  Read %d integration points in %.1f s\n',size(Data,1),toc(timer));
  r=struct([]);
  for component=1:3
   values=Data(:,component); mu=sum(values.*weights); delta=values-mu;
   variance=sum(delta.^2.*weights); skew=sum(delta.^3.*weights)/variance^1.5;
   kurt=sum(delta.^4.*weights)/variance^2;clear delta
   thresholds=0.5;if component==2,thresholds=[1.5 2 2.5];end
   prob=arrayfun(@(threshold) sum(weights(values>threshold)),thresholds);
   q99=fast_weighted_quantile(values,weights,0.99);
   if count==1 && component==2
    slow=weighted_upper_quantile(values,weights,0.01);
    assert(q99==slow,'Exact weighted-quantile cross-check failed.');
   end
   lo=min([min(values),mu-6*sqrt(variance),muModel(component)-6*sqrt(varModel(component))]);
   hi=max([max(values),mu+6*sqrt(variance),muModel(component)+6*sqrt(varModel(component))]);
   width=max(sqrt(variance)/20,(hi-lo)/250);finiteEdges=lo:width:hi;
   if finiteEdges(end)<hi,finiteEdges(end+1)=hi;end
   edges=[-Inf finiteEdges(2:end-1) Inf];
   [~,~,bins]=histcounts(values,edges);
   empirical=accumarray(bins,weights,[numel(edges)-1,1]).';empirical=empirical/sum(empirical);
   [~,~,~,cg]=gnd_parameters(edges,muModel(component),varModel(component),kurtModel(component));
   [~,~,~,cn]=gnd_parameters(edges,muModel(component),varModel(component));
   [jg,kg]=distance(empirical,diff(cg),1e-12);[jn,kn]=distance(empirical,diff(cn),1e-12);
   a=struct('FileName',fileName,'SourcePath',source,'Porosity',porosity(im),...
    'Realization',suffixes{ir},'Component',component,'NumPoints',size(Data,1),...
    'Mean',mu,'Variance',variance,'Skewness',skew,'Kurtosis',kurt,...
    'Q99',q99,'Thresholds',thresholds,'Exceedance',prob,...
    'JS_GND',jg,'JS_Gaussian',jn,'KL_GND',kg,'KL_Gaussian',kn,...
    'Edges',edges,'FEM_BinProbability',empirical,...
    'AnalyticalMean',muModel(component),'AnalyticalVariance',varModel(component),...
    'AnalyticalKurtosis',kurtModel(component));
   r=[r;a]; %#ok<AGROW>
   fprintf('  S%d: JS %.6g / %.6g; q99 %.9g; %.1f s elapsed\n',...
    component,jg,jn,q99,toc(timer));
  end
  clear Data weights values bins
  save(cache,'r','stamp');records=[records;r]; %#ok<AGROW>
  fprintf('  Cached all components: %.1f s\n',toc(timer));
 end
end
save(fullfile(cacheDir,'Round2_Field_Metrics.mat'),'records','-v7');
end

function q=fast_weighted_quantile(x,w,target)
% Exact inverse weighted empirical CDF. Sort only the bin containing target.
lo=min(x);hi=max(x);if lo==hi,q=lo;return,end
edges=linspace(lo,hi,513);[~,~,bin]=histcounts(x,edges);
mass=accumarray(bin,w,[512 1]);cum=cumsum(mass);ib=find(cum>=target,1);
if isempty(ib),q=hi;return,end
below=sum(mass(1:ib-1));inside=bin==ib;
[selected,idx]=sort(x(inside));ww=w(inside);cdf=below+cumsum(ww(idx));
k=find(cdf>=target,1);if isempty(k),k=numel(selected);end
q=selected(k);
end

function test_fast_quantile()
x=[-1;-1;0;0.3;1;2;5];w=[1;2;1;3;2;1;4];w=w/sum(w);
for target=[0.01 0.1 0.5 0.9 0.99]
 assert(fast_weighted_quantile(x,w,target)==weighted_upper_quantile(x,w,1-target));
end
end

function [js,kl]=distance(p,q,epsilon)
p=max(p(:).',0);p=p/sum(p);q=max(q(:).',0);q=q/sum(q);
assert(all(isfinite(p))&&all(isfinite(q)));
qFloor=max(q,epsilon);qFloor=qFloor/sum(qFloor);kl=relative_entropy(p,qFloor);
m=0.5*(p+q);js=0.5*relative_entropy(p,m)+0.5*relative_entropy(q,m);
assert(isfinite(js)&&isfinite(kl)&&js>=-1e-12&&js<=log(2)+1e-12&&kl>=-1e-12);
end

function value=relative_entropy(p,q)
positive=p>0;value=sum(p(positive).*log(p(positive)./q(positive)));
end
