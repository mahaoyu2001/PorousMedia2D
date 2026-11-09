function export_tail_validation()
p=tail_paths();saved=load(fullfile(p.cacheDir,'FEM_Statistics.mat'),'records');records=saved.records;
per=struct([]);qper=struct([]);mom=struct([]);
for k=1:numel(records)
 r=records(k);[ma,va,~,ka]=cal_stress_statistics([0;1;0],r.Porosity);c=r.Component;
 qa=gnd_upper_quantile(ma(c),va(c),ka(c),.01);
 qe=gnd_upper_quantile(r.Mean,r.Variance,r.Kurtosis,.01);
 m=rmfield(r,{'Thresholds','Exceedance'});mom=[mom;m]; %#ok<AGROW>
 qper=[qper;struct('Porosity',r.Porosity,'Component',c,'Realization',r.Realization,...
  'FEM_Q99',r.Q99,'GND_Analytical',qa,'GND_FEMmoments',qe)]; %#ok<AGROW>
 for j=1:numel(r.Thresholds)
  threshold=r.Thresholds(j);
  a=struct('Porosity',r.Porosity,'Component',c,'Realization',r.Realization,'Threshold',threshold,...
   'FEM',r.Exceedance(j),'GND_Analytical',gnd_survival(threshold,ma(c),va(c),ka(c)),...
   'GND_FEMmoments',gnd_survival(threshold,r.Mean,r.Variance,r.Kurtosis),...
   'GND_FEMmean_only',gnd_survival(threshold,r.Mean,va(c),ka(c)),...
   'GND_FEMvariance_only',gnd_survival(threshold,ma(c),r.Variance,ka(c)),...
   'GND_FEMkurtosis_only',gnd_survival(threshold,ma(c),va(c),r.Kurtosis),...
   'Gaussian_Analytical',gnd_survival(threshold,ma(c),va(c)));
  per=[per;a]; %#ok<AGROW>
 end
end
P=struct2table(per);Q=struct2table(qper);M=struct2table(mom);
summary=struct([]);qsummary=struct([]);convergence=struct([]);numerical=struct([]);
configNames={'Base','SpatialFine','GridFine','DomainWide','WeakFine','SpatialVeryFine'};
for f=p.f
 for c=1:3
  outs=cell(1,numel(configNames));
  for k=1:numel(configNames)
   name=sprintf('%s_S%d_f%04d.mat',configNames{k},c,round(f*10000));
   s=load(fullfile(p.cacheDir,name),'out');outs{k}=s.out;
   o=s.out;diag=rmfield(o,{'Thresholds','Probability','Probability_GP'});numerical=[numerical;diag]; %#ok<AGROW>
  end
  final=outs{6};[ma,va,~,ka]=cal_stress_statistics([0;1;0],f);
  empirical=records([records.Porosity]==f & [records.Component]==c);
  for j=1:numel(final.Thresholds)
   threshold=final.Thresholds(j);chosen=P(P.Porosity==f & P.Component==c & P.Threshold==threshold,:);
   assert(height(chosen)==5);pf=mean(chosen.FEM);ga=chosen.GND_Analytical(1);ge=mean(chosen.GND_FEMmoments);pc=final.Probability(j);
   allP=cellfun(@(o)o.Probability(j),outs);
   spatial=abs(allP(6)-allP(3));grid=abs(allP(3)-allP(2));domain=abs(allP(4)-allP(3));weak=abs(allP(5)-allP(3));
   estimate=max([spatial grid domain weak final.MaxGP_CDF_Difference]);
   passed=estimate<=max(1e-8,.005*pc);
   conv=struct('Porosity',f,'Component',c,'Threshold',threshold,'Campbell_Final',pc,...
    'CoarseSpatial_AbsChange',abs(allP(2)-allP(1)),...
    'Spatial_AbsChange',spatial,'Grid_AbsChange',grid,'Domain_AbsChange',domain,'WeakRegion_AbsChange',weak,...
    'GP_vs_FFT_AbsDifference',abs(final.Probability_GP(j)-pc),'EstimatedAbsSensitivity',estimate,...
    'EstimatedRelSensitivity_Percent',100*estimate/pc,'Passed',passed);
   convergence=[convergence;conv]; %#ok<AGROW>
   pooled=gnd_survival(threshold,mean([empirical.Mean]),mean([empirical.Variance]),mean([empirical.Kurtosis]));
   a=struct('Porosity',f,'Component',c,'Threshold',threshold,'FEM_Mean',pf,'FEM_SD',std(chosen.FEM),...
    'FEM_SE',std(chosen.FEM)/sqrt(5),'GND_Analytical',ga,'GND_FEMmoments_Mean',ge,...
    'GND_FEMmoments_SD',std(chosen.GND_FEMmoments),'GND_AveragedMoments',pooled,...
    'Campbell_Full',pc,'Gaussian_Analytical',chosen.Gaussian_Analytical(1),...
    'GND_FEMmean_only',mean(chosen.GND_FEMmean_only),...
    'GND_FEMvariance_only',mean(chosen.GND_FEMvariance_only),...
    'GND_FEMkurtosis_only',mean(chosen.GND_FEMkurtosis_only),...
    'Error_GND_Analytical_Percent',100*abs(ga-pf)/pf,...
    'Error_GND_FEMmoments_Percent',100*abs(ge-pf)/pf,'Error_Campbell_Percent',100*abs(pc-pf)/pf,...
    'GND_MomentInput_Contribution',ga-ge,'GND_Representation_Residual',ge-pf,...
    'Poisson_Closure_Contribution',ga-pc,'Poisson_FEM_Residual',pc-pf,...
    'Total_Signed_Error',ga-pf,'CF_NumericalSensitivity',estimate,'CF_Passed',passed);
   assert(abs(a.Total_Signed_Error-a.GND_MomentInput_Contribution-a.GND_Representation_Residual)<1e-14);
   assert(abs(a.Total_Signed_Error-a.Poisson_Closure_Contribution-a.Poisson_FEM_Residual)<1e-14);
   summary=[summary;a]; %#ok<AGROW>
  end
  chosen=Q(Q.Porosity==f & Q.Component==c,:);qf=mean(chosen.FEM_Q99);qg=chosen.GND_Analytical(1);qe=mean(chosen.GND_FEMmoments);
  qvalues=cellfun(@(o)o.Q99,outs);qSensitivity=max([abs(qvalues(6)-qvalues(3)),abs(qvalues([2 4 5])-qvalues(3))]);
  qsummary=[qsummary;struct('Porosity',f,'Component',c,'FEM_Mean',qf,'FEM_SD',std(chosen.FEM_Q99),...
   'FEM_SE',std(chosen.FEM_Q99)/sqrt(5),'GND_Analytical',qg,'GND_FEMmoments_Mean',qe,...
   'GND_FEMmoments_SD',std(chosen.GND_FEMmoments),'Campbell_Full',final.Q99,...
   'Error_GND_Analytical_Percent',100*abs(qg-qf)/abs(qf),...
   'Error_GND_FEMmoments_Percent',100*abs(qe-qf)/abs(qf),...
   'Error_Campbell_Percent',100*abs(final.Q99-qf)/abs(qf),'CF_NumericalSensitivity',qSensitivity)]; %#ok<AGROW>
 end
end
T=struct2table(summary);QT=struct2table(qsummary);C=struct2table(convergence);D=struct2table(numerical);
main=T(T.Component==2 & T.Threshold==2 & ismember(T.Porosity,[.005 .01 .02 .05 .1]),:);
main=main(:,{'Porosity','FEM_Mean','FEM_SD','GND_Analytical','GND_FEMmoments_Mean','Campbell_Full',...
 'Error_GND_Analytical_Percent','Error_GND_FEMmoments_Percent','Error_Campbell_Percent','CF_NumericalSensitivity'});
writetable(main,fullfile(p.tableDir,'R5Q4_Main_Sigma22_Threshold2.csv'));
writetable(T,fullfile(p.tableDir,'R5Q4_Exceedance_All.csv'));
writetable(P,fullfile(p.tableDir,'R5Q4_Exceedance_Per_Realization.csv'));
writetable(QT,fullfile(p.tableDir,'R5Q4_Q99_All.csv'));
writetable(Q,fullfile(p.tableDir,'R5Q4_Q99_Per_Realization.csv'));
writetable(M,fullfile(p.tableDir,'R5Q4_FEM_Moments.csv'));
writetable(C,fullfile(p.tableDir,'R5Q4_Convergence.csv'));
writetable(D,fullfile(p.tableDir,'R5Q4_Numerical_Diagnostics.csv'));
workbook=fullfile(p.tableDir,'R5Q4_Tail_Validation.xlsx');
for pair={ {'Main_S22',main},{'Exceedance',T},{'Q99',QT},{'Per_Realization',P},{'FEM_Moments',M},{'Convergence',C},{'Numerical_Checks',D} }
 item=pair{1};writetable(item{2},workbook,'Sheet',item{1});
end
manifest=readtable(fullfile(p.tableDir,'R5Q4_Source_Manifest.csv'));writetable(manifest,workbook,'Sheet','Source_Files');
save(fullfile(p.tableDir,'R5Q4_Tail_Validation.mat'),'T','QT','P','Q','M','C','D','main');
write_tex(main,p);
fprintf('\nPrimary R5Q4 table:\n');disp(main);
fprintf('Convergence checks passed %d/%d; maximum relative sensitivity %.6g percent.\n',sum(C.Passed),height(C),max(C.EstimatedRelSensitivity_Percent));
assert(all(C.Passed),'Some probabilities require additional convergence checks; see R5Q4_Convergence.csv.');
end

function write_tex(main,p)
fid=fopen(fullfile(p.tableDir,'R5Q4_Comparison_Table.tex'),'w');cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%% Input fragment, uses booktabs. Stress thresholds are normalized by remote Sigma2; probabilities are dimensionless.\n');
fprintf(fid,'\\begin{tabular}{@{}ccccc@{}}\n\\toprule\n');
fprintf(fid,'$f$ (\\%%) & FEM mean $\\pm$ SD & \\shortstack{GND\\\\analytical moments} & \\shortstack{GND\\\\FEM moments} & \\shortstack{Full\\\\Campbell}\\\\\n\\midrule\n');
for k=1:height(main)
 fprintf(fid,'%.1f & $%s\\pm%s$ & $%s$ & $%s$ & $%s$ \\\\\n',...
 100*main.Porosity(k),texnum(main.FEM_Mean(k),5),texnum(main.FEM_SD(k),3),texnum(main.GND_Analytical(k),5),texnum(main.GND_FEMmoments_Mean(k),5),texnum(main.Campbell_Full(k),5));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
end

function value=texnum(x,digits)
if x==0,value='0';return,end
parts=split(sprintf('%.*e',digits-1,x),'e');
value=[parts{1} '\times10^{' sprintf('%d',str2double(parts{2})) '}'];
end
