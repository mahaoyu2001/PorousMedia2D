function run_round2_distributions(mode)
% Corrected Eq.21--23 PDF/CDF, JS/KL, exceedance and quantile single panels.
% all: raw FEM/cache + figures; plots: cache + figures; collect: FEM/cache only.
% All styles/empirical graphical objects come from read-only original FIGs.
% Output names omit _S2.txt; empty original templates are not exported.
if nargin<1,mode='plots';end
p=round2_paths();oldPath=path;restore=onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(p.programDir,'-begin');
figDir=fullfile(p.figureDir,'Distributions_and_Tails');if ~isfolder(figDir),mkdir(figDir);end
diary(fullfile(p.logDir,['Distribution_Plots_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.log']));
cleanupDiary=onCleanup(@() diary('off')); %#ok<NASGU>
assert(strcmp(which('cal_stress_statistics'),fullfile(p.programDir,'cal_stress_statistics.m')));
assert(strcmp(which('gnd_parameters'),fullfile(p.programDir,'gnd_parameters.m')));
switch lower(mode)
 case {'all','collect'},records=round2_collect_field_metrics(false);
 case 'plots'
  saved=load(fullfile(p.rootDir,'Results','Cache','Field_Distributions','Round2_Field_Metrics.mat'),'records');records=saved.records;
 otherwise,error('Use all, collect or plots.');
end
assert(numel(records)==105,'Expected 7 porosities x 5 RVEs x 3 components.');
[jsTable,pTable,qTable,parameterTable]=export_numeric_results(records,p);
if strcmpi(mode,'collect'),return,end
report=struct([]); mapping=struct([]);plotData=struct([]);
for f=[0.01 0.05]
 for component=1:3
  for kind={'pdf','cdf'}
   for view={'','L','E'}
    originalName=sprintf('%s_2DVoid_Vol%d_a_S2.txt_stress%d%s',kind{1},100*f,component,view{1});
    name=erase(originalName,'_S2.txt');
    template=fullfile(p.styleDir,[originalName '.fig']);
    fig=openfig(template,'new','invisible');closeFig=onCleanup(@() close(fig));
    ax=one_axes(fig);round2_distribution_axis_override(ax,name);
    freeze_style(fig,ax);before=style_snapshot(fig,ax);
    empirical=empirical_snapshot(ax);
    [mu,v,~,k]=cal_stress_statistics([0;1;0],f);
    lines=findall(ax,'Type','line');gridVisible=false;
    for il=1:numel(lines)
     if ~isempty(model_kind(lines(il).DisplayName))
      gridVisible=gridVisible||any(lines(il).XData>=ax.XLim(1)&lines(il).XData<=ax.XLim(2));
     end
    end
    if ~gridVisible
     fprintf('Skipped empty original template: %s\n',originalName);
     clear closeFig
     continue
    end
    for il=1:numel(lines)
     line=lines(il);whichModel=model_kind(line.DisplayName);if isempty(whichModel),continue,end
     x=line.XData;[a,b,pdf,cdf]=gnd_parameters(x,mu(component),v(component),...
      choose_kurtosis(whichModel,k(component)));
     if strcmp(kind{1},'pdf'),line.YData=pdf;else,line.YData=cdf;end
     gridVisible=gridVisible||any(x>=ax.XLim(1)&x<=ax.XLim(2));
     plotData=[plotData;curve_records(name,whichModel,x,line.YData)]; %#ok<AGROW>
     assert(a>0&&b>0);
    end
    assert(isequaln(empirical,empirical_snapshot(ax)),'FEM curve/histogram changed.');
    entry=finish_figure(fig,ax,before,name,template,figDir,gridVisible,NaN);
    report=[report;entry]; %#ok<AGROW>
    source='hisogram.m';if ~isempty(view{1}),source='hisogram_ExtVal.m';end
    mapping=[mapping;mapping_row(name,source,template,f,component,NaN)]; %#ok<AGROW>
    clear closeFig
   end
  end
 end
end
for component=1:3
 for metric={'JS','KL'}
  name=sprintf('%s_Divergence_vol_stress%d',metric{1},component);
  template=fullfile(p.styleDir,[name '.fig']);
  fig=openfig(template,'new','invisible');closeFig=onCleanup(@() close(fig));
  ax=one_axes(fig);round2_distribution_axis_override(ax,name);freeze_style(fig,ax);before=style_snapshot(fig,ax);
  rows=jsTable(jsTable.Component==component,:);lines=findall(ax,'Type','line');
  for il=1:numel(lines)
   line=lines(il);whichModel=model_kind(line.DisplayName);assert(~isempty(whichModel));
   assert(max(abs(line.XData(:)-rows.Porosity))<1e-12);
   line.YData=rows.([metric{1} '_' whichModel '_Mean']).';
   plotData=[plotData;curve_records(name,whichModel,line.XData,line.YData)]; %#ok<AGROW>
  end
  report=[report;finish_figure(fig,ax,before,name,template,figDir,true,NaN)]; %#ok<AGROW>
  mapping=[mapping;mapping_row(name,'ExtVal_Divergence.m',template,NaN,component,NaN)]; %#ok<AGROW>
  clear closeFig
 end
end
for component=1:3
 thresholds=0.5;if component==2,thresholds=[2 1.5 2.5];end
 for threshold=thresholds
  suffix='';if component==2&&threshold~=2,suffix=sprintf('_SigmaC%d',10*threshold);end
  name=sprintf('P_EaxtVal_vol_stress%d_log_errorband%s',component,suffix);
  template=fullfile(p.styleDir,[name '.fig']);
  fig=openfig(template,'new','invisible');closeFig=onCleanup(@() close(fig));
  ax=one_axes(fig);round2_distribution_axis_override(ax,name);freeze_style(fig,ax);before=style_snapshot(fig,ax);empirical=empirical_snapshot(ax);
  rows=pTable(pTable.Component==component & pTable.Threshold==threshold,:);
  probabilityFloor=1e-4;if component~=2,probabilityFloor=4e-5;end
  femDifference=validate_fem_objects(ax,rows.Porosity,rows.FEM_Mean,rows.FEM_SD,probabilityFloor);
  lines=findall(ax,'Type','line');
  for il=1:numel(lines)
   line=lines(il);whichModel=model_kind(line.DisplayName);if isempty(whichModel),continue,end
   x=line.XData;[mu,v,~,k]=cal_stress_statistics([0;1;0],x);y=zeros(size(x));
   for ix=1:numel(x),y(ix)=gnd_survival(threshold,mu(component,ix),v(component,ix),choose_kurtosis(whichModel,k(component,ix)));end
   line.YData=y;plotData=[plotData;curve_records(name,whichModel,x,y)]; %#ok<AGROW>
  end
  assert(isequaln(empirical,empirical_snapshot(ax)));
  report=[report;finish_figure(fig,ax,before,name,template,figDir,true,femDifference)]; %#ok<AGROW>
  source='ExtVal_Probability_Log_ErrorBand_Sigma1112.m';
  if component==2,source='ExtVal_Probability_Log_ErrorBand.m';if threshold~=2,source='ExtVal_Probability_Log_ErrorBand_SigmaC.m';end,end
  mapping=[mapping;mapping_row(name,source,template,NaN,component,threshold)]; %#ok<AGROW>
  clear closeFig
 end
 name=sprintf('Q_EaxtVal_vol_stress%d_errorband',component);
 template=fullfile(p.styleDir,[name '.fig']);
 fig=openfig(template,'new','invisible');closeFig=onCleanup(@() close(fig));
 ax=one_axes(fig);round2_distribution_axis_override(ax,name);freeze_style(fig,ax);before=style_snapshot(fig,ax);empirical=empirical_snapshot(ax);
 rows=qTable(qTable.Component==component,:);
 femDifference=validate_fem_objects(ax,rows.Porosity,rows.FEM_Mean,rows.FEM_SD,NaN);
 lines=findall(ax,'Type','line');
 for il=1:numel(lines)
  line=lines(il);whichModel=model_kind(line.DisplayName);if isempty(whichModel),continue,end
  x=line.XData;[mu,v,~,k]=cal_stress_statistics([0;1;0],x);y=zeros(size(x));
  for ix=1:numel(x),y(ix)=gnd_upper_quantile(mu(component,ix),v(component,ix),choose_kurtosis(whichModel,k(component,ix)),0.01);end
  line.YData=y;plotData=[plotData;curve_records(name,whichModel,x,y)]; %#ok<AGROW>
 end
 assert(isequaln(empirical,empirical_snapshot(ax)));
 report=[report;finish_figure(fig,ax,before,name,template,figDir,true,femDifference)]; %#ok<AGROW>
 source='ExtVal_Quantile_ErrorBand_Sigma1112.m';if component==2,source='ExtVal_Quantile_ErrorBand.m';end
 mapping=[mapping;mapping_row(name,source,template,NaN,component,NaN)]; %#ok<AGROW>
 clear closeFig
end
assert(numel(report)==44,'Expected 30 PDF/CDF + 6 JS/KL + 8 tail panels; 6 empty templates skipped.');
writetable(struct2table(report),fullfile(p.logDir,'Distribution_Style_Verification.csv'));
issueRows=~cellfun(@isempty,{report.OriginalAxisIssue});
writetable(struct2table(report(issueRows)),fullfile(p.logDir,'Original_Axis_Issues.csv'));
writetable(struct2table(mapping),fullfile(p.tableDir,'Round2_Distribution_Source_Map.csv'));
writetable(struct2table(plotData),fullfile(p.tableDir,'Round2_Distribution_Plot_Curves.csv'));
save(fullfile(p.tableDir,'Round2_Distribution_Results.mat'),'records','jsTable','pTable','qTable','parameterTable','mapping','report');
fprintf('Completed 44 single panels (600 dpi PNG + editable FIG); empty templates skipped. Output: %s\n',figDir);
end

function [jt,pt,qt,params]=export_numeric_results(records,p)
fValues=unique([records.Porosity]);jRows=struct([]);pRows=struct([]);qRows=struct([]);paramRows=struct([]);
raw=rmfield(records,{'Thresholds','Exceedance','Edges','FEM_BinProbability'});
writetable(struct2table(raw),fullfile(p.tableDir,'Round2_Distribution_Per_Realization.csv'));
pRaw=struct([]);
for i=1:numel(records)
 r=records(i);for j=1:numel(r.Thresholds)
  pRaw=[pRaw;struct('Porosity',r.Porosity,'Component',r.Component,'Realization',r.Realization,...
   'Threshold',r.Thresholds(j),'FEM_Probability',r.Exceedance(j),'SourceFile',r.SourcePath)]; %#ok<AGROW>
 end
end
writetable(struct2table(pRaw),fullfile(p.tableDir,'Round2_Exceedance_Per_Realization.csv'));
for component=1:3
 for f=fValues
  r=records([records.Component]==component & abs([records.Porosity]-f)<1e-12);assert(numel(r)==5);
  [mu,v,~,k]=cal_stress_statistics([0;1;0],f);[a,b]=gnd_parameters([],mu(component),v(component),k(component));
  row=struct('Porosity',f,'Component',component);
  for metric={'JS','KL'}
   for model={'GND','Gaussian'}
    field=[metric{1} '_' model{1}];values=[r.(field)];
    row.([field '_Mean'])=mean(values);row.([field '_SD'])=std(values,0);
   end
  end
  jRows=[jRows;row]; %#ok<AGROW>
  fq=[r.Q99];gq=gnd_upper_quantile(mu(component),v(component),k(component),0.01);
  nq=gnd_upper_quantile(mu(component),v(component),3,0.01);
  qRows=[qRows;struct('Porosity',f,'Component',component,'FEM_Mean',mean(fq),'FEM_SD',std(fq,0),...
   'GND',gq,'Gaussian',nq,'GND_Error_Percent',100*abs(gq-mean(fq))/abs(mean(fq)),...
   'Gaussian_Error_Percent',100*abs(nq-mean(fq))/abs(mean(fq)))]; %#ok<AGROW>
  for j=1:numel(r(1).Thresholds)
   threshold=r(1).Thresholds(j);fp=arrayfun(@(s)s.Exceedance(j),r);
   gp=gnd_survival(threshold,mu(component),v(component),k(component));np=gnd_survival(threshold,mu(component),v(component),3);
   pgErr=NaN;pnErr=NaN;if mean(fp)>0,pgErr=100*abs(gp-mean(fp))/mean(fp);pnErr=100*abs(np-mean(fp))/mean(fp);end
   pRows=[pRows;struct('Porosity',f,'Component',component,'Threshold',threshold,...
    'FEM_Mean',mean(fp),'FEM_SD',std(fp,0),'GND',gp,'Gaussian',np,...
    'GND_Error_Percent',pgErr,'Gaussian_Error_Percent',pnErr)]; %#ok<AGROW>
  end
  paramRows=[paramRows;struct('Porosity',f,'Component',component,'AnalyticalMean',mu(component),...
   'AnalyticalVariance',v(component),'AnalyticalKurtosis',k(component),'GND_Alpha',a,'GND_Beta',b,...
   'FEM_Mean_Mean',mean([r.Mean]),'FEM_Variance_Mean',mean([r.Variance]),...
   'FEM_Skewness_Mean',mean([r.Skewness]),'FEM_Kurtosis_Mean',mean([r.Kurtosis]))]; %#ok<AGROW>
 end
end
jt=struct2table(jRows);pt=struct2table(pRows);qt=struct2table(qRows);params=struct2table(paramRows);
cached=readtable(fullfile(p.tableDir,'Round2_Comparison_geom.csv'));
labels={'S11','S22','S12'};
for i=1:height(params)
 row=cached(abs(cached.Porosity-params.Porosity(i))<1e-12 & ...
  strcmp(string(cached.StressComponent),labels{params.Component(i)}),:);
 assert(height(row)==1);
 assert(abs(params.FEM_Variance_Mean(i)-row.Variance)<1e-9,...
  'Raw FEM variance differs from the source used for Fig.4.');
 assert(abs(params.FEM_Kurtosis_Mean(i)-row.Kurtosis)<1e-6,...
  'Raw FEM kurtosis differs from the source used for Fig.4.');
end
writetable(jt,fullfile(p.tableDir,'Round2_JS_KL_Comparison.csv'));
writetable(pt,fullfile(p.tableDir,'Round2_Exceedance_Comparison.csv'));
writetable(qt,fullfile(p.tableDir,'Round2_Q99_Comparison.csv'));
writetable(params,fullfile(p.tableDir,'Round2_Distribution_Parameters.csv'));
end

function ax=one_axes(fig)
ax=findall(fig,'Type','axes');assert(numel(ax)==1,'Expected one original axes.');
end

function freeze_style(fig,ax)
drawnow;
ax.XLim=ax.XLim;ax.YLim=ax.YLim;ax.XTick=ax.XTick;ax.YTick=ax.YTick;
legends=findall(fig,'Type','legend');
for j=1:numel(legends),legends(j).Position=legends(j).Position;legends(j).AutoUpdate='off';end
drawnow;
ax.XLabel.Position=ax.XLabel.Position;ax.YLabel.Position=ax.YLabel.Position;
drawnow;
end

function snapshot=style_snapshot(fig,ax)
snapshot.Figure=get(fig,{'Units','Position','Color','Renderer','PaperUnits','PaperSize','PaperPosition'});
snapshot.Axes=get(ax,{'Units','Position','XLim','YLim','XTick','YTick','XScale','YScale',...
 'XDir','YDir','FontName','FontSize','LineWidth','TickLength','TickDir','XMinorTick','YMinorTick','Box'});
snapshot.Labels={get(ax.XLabel,{'String','FontName','FontSize','FontWeight','Interpreter','Position'}),...
 get(ax.YLabel,{'String','FontName','FontSize','FontWeight','Interpreter','Position'})};
legends=findall(fig,'Type','legend');snapshot.Legends=cell(numel(legends),1);
for j=1:numel(legends),snapshot.Legends{j}=get(legends(j),{'Position','String','FontName','FontSize','Box','ItemTokenSize'});end
types={'line','scatter','bar','patch'};snapshot.Objects=cell(numel(types),1);
for j=1:numel(types)
 objs=findall(ax,'Type',types{j});data=cell(numel(objs),1);
 for k=1:numel(objs)
  switch types{j}
   case 'line',properties={'DisplayName','LineWidth','LineStyle','Color','Marker','MarkerSize','MarkerFaceColor','MarkerEdgeColor'};
   case 'scatter',properties={'DisplayName','LineWidth','SizeData','Marker','MarkerFaceColor','MarkerEdgeColor'};
   case 'bar',properties={'DisplayName','LineWidth','BarWidth','FaceColor','FaceAlpha','EdgeColor','EdgeAlpha'};
   case 'patch',properties={'DisplayName','LineWidth','FaceColor','FaceAlpha','EdgeColor'};
  end
  data{k}=get(objs(k),properties);
 end
 snapshot.Objects{j}=data;
end
end

function snapshot=empirical_snapshot(ax)
snapshot=struct([]);types={'line','scatter','bar','patch'};
for i=1:numel(types)
 objs=findall(ax,'Type',types{i});for j=1:numel(objs)
  obj=objs(j);if ~isempty(model_kind(obj.DisplayName)),continue,end
  snapshot=[snapshot;struct('Type',obj.Type,'Name',obj.DisplayName,'X',obj.XData,'Y',obj.YData)]; %#ok<AGROW>
 end
end
end

function result=finish_figure(fig,ax,before,name,template,figDir,visible,femDifference)
drawnow;after=style_snapshot(fig,ax);
if ~isequaln(before,after)
 fields=fieldnames(before);
 for i=1:numel(fields)
  if ~isequaln(before.(fields{i}),after.(fields{i}))
   fprintf('Style difference %s / %s\n',name,fields{i});disp(before.(fields{i}));disp(after.(fields{i}));
  end
 end
end
assert(isequaln(before,after),'Original style changed: %s',name);
base=fullfile(figDir,name);exportgraphics(fig,[base '.png'],'Resolution',600);savefig(fig,[base '.fig']);
oldPng=[erase(template,'.fig') '.png'];old=imfinfo(oldPng);new=imfinfo([base '.png']);
issue='';parts=regexp(name,'^(pdf|cdf)_2DVoid_Vol(?<volume>\d+)_a_stress(?<component>[123])(?<view>[LE]?)$','names','once');
if ~visible
 issue='Original XLim excludes the model data grid; original empty-plot problem retained';
elseif ~isempty(parts) && isempty(parts.view)
 center=0;if str2double(parts.component)==2,center=1/(1-str2double(parts.volume)/100);end
 if center<ax.XLim(1)||center>ax.XLim(2)
  issue='Original main-plot XLim excludes the distribution center; original range retained';
 end
end
result=struct('Figure',name,'StyleCheck','PASS','FEMObjects','UNCHANGED',...
 'ModelPointsInOriginalXLim',visible,'FEMNumericMaxDifference',femDifference,...
 'OriginalPixelWidth',old.Width,'OriginalPixelHeight',old.Height,...
 'OutputPixelWidth',new.Width,'OutputPixelHeight',new.Height,'OriginalAxisIssue',issue);
fprintf('Original-style panel: %s\n',name);
end

function which=model_kind(name)
name=lower(strtrim(name));which='';
if contains(name,'gen'),which='GND';elseif contains(name,'normal')||contains(name,'gauss'),which='Gaussian';end
end

function kurt=choose_kurtosis(which,kurt)
if strcmp(which,'Gaussian'),kurt=3;end
end

function rows=curve_records(name,model,x,y)
rows=struct('Figure',{},'Model',{},'X',{},'Y',{});
for j=1:numel(x),rows(j,1)=struct('Figure',name,'Model',model,'X',x(j),'Y',y(j));end
end

function row=mapping_row(name,script,template,f,component,threshold)
row=struct('Figure',name,'OriginalScript',script,'OriginalFIG',template,...
 'Porosity',f,'Component',component,'Threshold',threshold,...
 'ModelFormula','Revised manuscript Eqs.21-23; rho*A0 approximately f; analytical GND and Gaussian moments');
end

function difference=validate_fem_objects(ax,x,meanValue,sd,probabilityFloor)
scatter=findall(ax,'Type','scatter');assert(numel(scatter)==1);
assert(max(abs(scatter.XData(:)-x(:)))<1e-12);
difference=max(abs(scatter.YData(:)-meanValue(:)));
assert(difference<1e-8,'Original FEM markers differ from original raw-data recalculation: %.3g',difference);
patches=findall(ax,'Type','patch');assert(numel(patches)==1);
assert(numel(patches.XData)==2*numel(x),'Unexpected original SD band geometry.');
lower=meanValue-sd;upper=meanValue+sd;
if isfinite(probabilityFloor)
 lower=max(lower,probabilityFloor);upper=max(upper,probabilityFloor);
end
expected=[upper(:);flipud(lower(:))];
bandDifference=max(abs(patches.YData(:)-expected));
assert(bandDifference<1e-8,'Original SD band differs from raw FEM: %.3g',bandDifference);
difference=max(difference,bandDifference);
end
