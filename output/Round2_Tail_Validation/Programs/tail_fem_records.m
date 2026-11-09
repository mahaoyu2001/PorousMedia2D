function records=tail_fem_records(force)
% Read source fields or reuse only byte/time-validated per-source field metrics.
if nargin<1,force=false;end
p=tail_paths();records=struct([]);sourceRows=struct([]);
for im=1:numel(p.models)
 for ir=1:numel(p.realizations)
  file=sprintf('%s_%s_S2.txt',p.models{im},p.realizations{ir});
  source=fullfile(p.fieldDir,file);info=dir(source);assert(~isempty(info),'Missing %s',source);
  ownCache=fullfile(p.cacheDir,['FEM_Source_' erase(file,'.txt') '.mat']);
  cacheCandidates={ownCache,fullfile(p.legacyCacheDir,[erase(file,'.txt') '.mat'])};valid=false;
  if ~force
   for ic=1:numel(cacheCandidates)
    if ~isfile(cacheCandidates{ic}),continue,end
    saved=load(cacheCandidates{ic},'stamp','r');
    valid=saved.stamp.Version==2 && saved.stamp.Bytes==info.bytes && saved.stamp.Datenum==info.datenum;
    if valid
     r=saved.r;
     assert(numel(r)==3 && all([r.Porosity]==p.f(im)) && isequal([r.Component],[1 2 3]));
     break
    end
   end
  end
  if ~valid
   fprintf('Reading source field: %s (%.1f MB)\n',file,info.bytes/1e6);
   data=read_stress_field(source);
   assert(size(data,2)==4 && all(isfinite(data(:))) && all(data(:,4)>0));
   w=data(:,4)/sum(data(:,4));r=struct([]);
   for c=1:3
    x=data(:,c);mu=sum(w.*x);d=x-mu;v=sum(w.*d.^2);k=sum(w.*d.^4)/v^2;sk=sum(w.*d.^3)/v^1.5;
    thresholds=.5;if c==2,thresholds=[1.5 2 2.5];end
    a=struct('FileName',file,'SourcePath',source,'Porosity',p.f(im),'Realization',p.realizations{ir},...
      'Component',c,'NumPoints',numel(x),'Mean',mu,'Variance',v,'Skewness',sk,'Kurtosis',k,...
      'Q99',weighted_q(x,w,.99),'Thresholds',thresholds,'Exceedance',arrayfun(@(z)sum(w(x>z)),thresholds));
    r=[r;a]; %#ok<AGROW>
   end
   clear data w x d
  end
  perSource=struct([]);
  for c=1:3
   a=r(c); % Keep only source statistics, never prior model predictions.
   slim=struct('FileName',file,'Porosity',a.Porosity,'Realization',p.realizations{ir},'Component',c,...
    'NumPoints',a.NumPoints,'Mean',a.Mean,'Variance',a.Variance,'Skewness',a.Skewness,...
    'Kurtosis',a.Kurtosis,'Q99',a.Q99,'Thresholds',a.Thresholds,'Exceedance',a.Exceedance);
   records=[records;slim]; %#ok<AGROW>
   perSource=[perSource;slim]; %#ok<AGROW>
  end
  r=perSource;stamp=struct('Version',2,'Bytes',info.bytes,'Datenum',info.datenum);save(ownCache,'r','stamp');
  method='raw_field_read';if valid,method='validated_source_statistics_cache';end
  sourceRows=[sourceRows;struct('FileName',file,'FullPath',source,'Bytes',info.bytes,...
   'SourceDatenum',info.datenum,'Method',method,'NumPoints',r(1).NumPoints)]; %#ok<AGROW>
  fprintf('[%d/35] %s: %s\n',(im-1)*5+ir,file,method);
 end
end
writetable(struct2table(sourceRows),fullfile(p.tableDir,'R5Q4_Source_Manifest.csv'));
save(fullfile(p.cacheDir,'FEM_Statistics.mat'),'records','sourceRows');
end

function q=weighted_q(x,w,target)
lo=min(x);hi=max(x);edges=linspace(lo,hi,513);[~,~,bin]=histcounts(x,edges);
mass=accumarray(bin,w,[512 1]);cum=cumsum(mass);ib=find(cum>=target,1);
below=sum(mass(1:ib-1));inside=bin==ib;[xx,order]=sort(x(inside));ww=w(inside);cdf=below+cumsum(ww(order));
q=xx(find(cdf>=target,1));
end
