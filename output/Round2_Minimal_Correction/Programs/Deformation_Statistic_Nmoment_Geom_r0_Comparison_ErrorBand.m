clear;
r2=round2_paths(); % Isolated round-2 outputs; original FEM inputs are read-only.

% Compare analytical predictions with two FEM microstructure generators:
%   1) the original hard-core model;
%   2) the overlap-allowed model with minimum centre distance r_min = 0.
%
% The existing plotting scripts and field-value files are not modified.
% Missing files are skipped. Vol2 is currently omitted from the two aligned
% lists below; add its model name and porosity back after the r0 files exist.

suffixList={'_a','_b','_c','_d','_e'};
baseModelNameList={'2DVoid_Vol05','2DVoid_Vol1','2DVoid_Vol2',...
    '2DVoid_Vol3','2DVoid_Vol5','2DVoid_Vol7_5','2DVoid_Vol10'};
vol=[0.5 1 2 3 5 7.5 10]/100;

methodFileTags={'','_r0'};
methodNames={'Hard-core FEM','  Poisson FEM'};
methodShortNames={'Hard-core','r_min=0'};

Sigma=[0;1;0];
components=[1 2 3];
componentNames=["S11","S22","S12"];
fieldName='S';
step=2;

% Publication-oriented vertical-axis limits for S22.

varYLim=[0 0.045;...
         0 0.20;...
         0 0.04];
skeYLim=repmat([-1.5 1.5],3,1);
kurYLim=[0 100;...
         0 115;...
         0 90];

scriptDir=fileparts(mfilename('fullpath'));
fieldDir=r2.fieldDir;
picDir=r2.figureDir;
if ~isfolder(picDir)
    mkdir(picDir);
end
addpath(scriptDir);

nMethod=length(methodFileTags);
nRealization=length(suffixList);
nModel=length(baseModelNameList);
nComponent=length(components);

muRVE=nan(nMethod,nRealization,nModel,nComponent);
moment2=nan(nMethod,nRealization,nModel,nComponent);
skewnessRVE=nan(nMethod,nRealization,nModel,nComponent);
kurtosisRVE=nan(nMethod,nRealization,nModel,nComponent);
fileAvailable=false(nMethod,nRealization,nModel);

useCachedFEM=true; % Set false to recalculate from original FieldValue files.
if useCachedFEM
    [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD,sampleCount]=...
        round2_cached_point_arrays(vol,methodShortNames,componentNames);
else
%% Read each large field file once and calculate all three components
fprintf('Reading field files from:\n  %s\n',fieldDir);
for im=1:nMethod
    for nf=1:nModel
        modelName=[baseModelNameList{nf} methodFileTags{im}];
        for ir=1:nRealization
            fileName=[modelName suffixList{ir} '_' fieldName num2str(step) '.txt'];
            filePath=fullfile(fieldDir,fileName);
            if ~isfile(filePath)
                warning('Skipping missing field file: %s',filePath);
                continue
            end

            Data=importdata(filePath);
            if isstruct(Data)
                Data=Data.data;
            end
            if ~isnumeric(Data) || size(Data,2)<4
                error('Unexpected data format in %s.',filePath);
            end

            Evol=Data(:,end);
            totalArea=sum(Evol);
            if ~isfinite(totalArea) || totalArea<=0 || any(~isfinite(Evol))
                error('Invalid element-area column in %s.',filePath);
            end
            weights=Evol/totalArea;
            fileAvailable(im,ir,nf)=true;

            for ia=1:nComponent
                component=components(ia);
                phaseValue=Data(:,component);
                mu=sum(phaseValue.*weights);
                centralValue=phaseValue-mu;
                m2=sum(centralValue.^2.*weights);
                m3=sum(centralValue.^3.*weights);
                m4=sum(centralValue.^4.*weights);

                muRVE(im,ir,nf,ia)=mu;
                moment2(im,ir,nf,ia)=m2;
                if m2>0
                    skewnessRVE(im,ir,nf,ia)=m3/m2^1.5;
                    kurtosisRVE(im,ir,nf,ia)=m4/m2^2;
                end
            end
            fprintf('Loaded: %s\n',fileName);
            clear Data Evol weights phaseValue centralValue
        end
    end
end

%% Means and sample SD across the available realizations
meanMean=nan(nMethod,nComponent,nModel);
meanSD=nan(nMethod,nComponent,nModel);
varMean=nan(nMethod,nComponent,nModel);
varSD=nan(nMethod,nComponent,nModel);
skeMean=nan(nMethod,nComponent,nModel);
skeSD=nan(nMethod,nComponent,nModel);
kurMean=nan(nMethod,nComponent,nModel);
kurSD=nan(nMethod,nComponent,nModel);
sampleCount=zeros(nMethod,nComponent,nModel);

for im=1:nMethod
    for ia=1:nComponent
        [valueMean,valueSD,valueN]=column_mean_sd(...
            squeeze(muRVE(im,:,:,ia)));
        meanMean(im,ia,:)=reshape(valueMean,1,1,[]);
        meanSD(im,ia,:)=reshape(valueSD,1,1,[]);
        sampleCount(im,ia,:)=reshape(valueN,1,1,[]);

        [valueMean,valueSD]=column_mean_sd(squeeze(moment2(im,:,:,ia)));
        varMean(im,ia,:)=reshape(valueMean,1,1,[]);
        varSD(im,ia,:)=reshape(valueSD,1,1,[]);

        [valueMean,valueSD]=column_mean_sd(...
            squeeze(skewnessRVE(im,:,:,ia)));
        skeMean(im,ia,:)=reshape(valueMean,1,1,[]);
        skeSD(im,ia,:)=reshape(valueSD,1,1,[]);

        [valueMean,valueSD]=column_mean_sd(...
            squeeze(kurtosisRVE(im,:,:,ia)));
        kurMean(im,ia,:)=reshape(valueMean,1,1,[]);
        kurSD(im,ia,:)=reshape(valueSD,1,1,[]);
    end
end

end % useCachedFEM

%% Analytical predictions
% cal_stress_statistics returns variance, standardized skewness and
% kurtosis. The analytical curve starts at the smallest simulated porosity.
xAnalytical=min(vol):0.0001:0.1125;
[~,varAnalytical,skeAnalytical,kurAnalytical]=...
    cal_stress_statistics(Sigma,xAnalytical);

[curLineWidth,axLineWidth,~,axTickFontSize,...
    ~,lgdFontSize,axTickFontName,~,...
    lgdFontName,~]=plotParams('curLineWidth',1.25,'markerSize',5.5);

hardCoreColor=[0.6350 0.0780 0.1840];
r0Color=[0.3010 0.7450 0.9330];
methodColors=[hardCoreColor;r0Color];
methodBandColors=[0.86 0.36 0.43;...
                  0.55 0.80 0.94];
methodMarkers={'+','o'};
bandAlpha=0.18;

%% Variance, skewness and kurtosis comparison figures
for ia=1:nComponent
    component=components(ia);

    plot_comparison(vol,squeeze(varMean(:,ia,:)),squeeze(varSD(:,ia,:)),...
        xAnalytical,varAnalytical(component,:),methodNames,methodColors,...
        methodBandColors,methodMarkers,bandAlpha,varYLim(ia,:),...
        fullfile(picDir,['var_vol_stress' num2str(component) '_r0_comparison_errorband']),...
        false,false,80,1.6,ia==1,ia==nComponent,curLineWidth,axLineWidth,...
        axTickFontSize,lgdFontSize,axTickFontName,lgdFontName);

    plot_comparison(vol,squeeze(skeMean(:,ia,:)),squeeze(skeSD(:,ia,:)),...
        xAnalytical,skeAnalytical(component,:),methodNames,methodColors,...
        methodBandColors,methodMarkers,bandAlpha,skeYLim(ia,:),...
        fullfile(picDir,['ske_vol_stress' num2str(component) '_r0_comparison_errorband']),...
        false,false,100,1.8,ia==1,ia==nComponent,curLineWidth,axLineWidth,...
        axTickFontSize,lgdFontSize,axTickFontName,lgdFontName);

    plot_comparison(vol,squeeze(kurMean(:,ia,:)),squeeze(kurSD(:,ia,:)),...
        xAnalytical,kurAnalytical(component,:),methodNames,methodColors,...
        methodBandColors,methodMarkers,bandAlpha,kurYLim(ia,:),...
        fullfile(picDir,['kur_vol_stress' num2str(component) '_r0_comparison_errorband']),...
        true,true,80,1.6,ia==1,ia==nComponent,curLineWidth,axLineWidth,...
        axTickFontSize,lgdFontSize,axTickFontName,lgdFontName);
end

%% Export summary statistics for the reviewer-response analysis
nSummary=nMethod*nComponent*nModel;
Method=strings(nSummary,1);
StressComponent=strings(nSummary,1);
Porosity=zeros(nSummary,1);
NumRealizations=zeros(nSummary,1);
MeanStress=nan(nSummary,1);
MeanStress_SD=nan(nSummary,1);
Variance=nan(nSummary,1);
Variance_SD=nan(nSummary,1);
Skewness=nan(nSummary,1);
Skewness_SD=nan(nSummary,1);
Kurtosis=nan(nSummary,1);
Kurtosis_SD=nan(nSummary,1);

row=0;
for im=1:nMethod
    for ia=1:nComponent
        for nf=1:nModel
            row=row+1;
            Method(row)=methodShortNames{im};
            StressComponent(row)=componentNames{ia};
            Porosity(row)=vol(nf);
            NumRealizations(row)=sampleCount(im,ia,nf);
            MeanStress(row)=meanMean(im,ia,nf);
            MeanStress_SD(row)=meanSD(im,ia,nf);
            Variance(row)=varMean(im,ia,nf);
            Variance_SD(row)=varSD(im,ia,nf);
            Skewness(row)=skeMean(im,ia,nf);
            Skewness_SD(row)=skeSD(im,ia,nf);
            Kurtosis(row)=kurMean(im,ia,nf);
            Kurtosis_SD(row)=kurSD(im,ia,nf);
        end
    end
end

summaryTable=table(Method,StressComponent,Porosity,NumRealizations,...
    MeanStress,MeanStress_SD,Variance,Variance_SD,Skewness,Skewness_SD,...
    Kurtosis,Kurtosis_SD);

nDetail=nMethod*nComponent*nModel*nRealization;
Method=strings(nDetail,1);
StressComponent=strings(nDetail,1);
Porosity=zeros(nDetail,1);
Realization=strings(nDetail,1);
FileAvailable=false(nDetail,1);
MeanStress=nan(nDetail,1);
Variance=nan(nDetail,1);
Skewness=nan(nDetail,1);
Kurtosis=nan(nDetail,1);

row=0;
for im=1:nMethod
    for ia=1:nComponent
        for nf=1:nModel
            for ir=1:nRealization
                row=row+1;
                Method(row)=methodShortNames{im};
                StressComponent(row)=componentNames{ia};
                Porosity(row)=vol(nf);
                Realization(row)=erase(suffixList{ir},'_');
                FileAvailable(row)=fileAvailable(im,ir,nf);
                MeanStress(row)=muRVE(im,ir,nf,ia);
                Variance(row)=moment2(im,ir,nf,ia);
                Skewness(row)=skewnessRVE(im,ir,nf,ia);
                Kurtosis(row)=kurtosisRVE(im,ir,nf,ia);
            end
        end
    end
end
realizationTable=table(Method,StressComponent,Porosity,Realization,...
    FileAvailable,MeanStress,Variance,Skewness,Kurtosis);

excelFile=fullfile(r2.tableDir,...
    'Deformation_Statistics_Geom_r0_Comparison_ErrorBand.xlsx');
if isfile(excelFile)
    delete(excelFile);
end
writetable(summaryTable,excelFile,'Sheet','Summary');
if useCachedFEM
    realizationTable=readtable(fullfile(r2.sourceOutputDir,...
        'Deformation_Statistics_Geom_r0_Comparison_ErrorBand.xlsx'),'Sheet','Realizations');
end
writetable(realizationTable,excelFile,'Sheet','Realizations');

fprintf('\nFinished. Figures were saved in:\n  %s\n',picDir);
fprintf('Statistics were saved in:\n  %s\n',excelFile);


function [valueMean,valueSD,valueN]=column_mean_sd(values)
% Mean and sample SD by column while ignoring unavailable realizations.
valid=isfinite(values);
valueN=sum(valid,1);
valuesForSum=values;
valuesForSum(~valid)=0;
valueMean=sum(valuesForSum,1)./valueN;
valueMean(valueN==0)=nan;

deviation=bsxfun(@minus,values,valueMean);
deviation(~valid)=0;
valueSD=sqrt(sum(deviation.^2,1)./(valueN-1));
valueSD(valueN<2)=nan;
end


function plot_comparison(xData,meanData,sdData,xAnalytical,yAnalytical,...
    methodNames,methodColors,methodBandColors,methodMarkers,bandAlpha,...
    yLimits,picBase,nonnegativeBand,showGaussianLine,scatterSize,...
    scatterLineWidth,showLegend,showXTickLabels,curLineWidth,...
    axLineWidth,axTickFontSize,lgdFontSize,axTickFontName,lgdFontName)

f=figure;
f.Units='centimeters';
ax=gca;
ax.Units='centimeters';
ax.Position=[2,2,5.8,4.25];
hold(ax,'on');

nMethod=size(meanData,1);
hFEM=gobjects(1,nMethod);
for im=1:nMethod
    validBand=isfinite(meanData(im,:)) & isfinite(sdData(im,:));
    lowerBound=meanData(im,validBand)-sdData(im,validBand);
    if nonnegativeBand
        lowerBound=max(0,lowerBound);
    end
    upperBound=meanData(im,validBand)+sdData(im,validBand);
    draw_error_band(xData(validBand),lowerBound,upperBound,...
        methodBandColors(im,:),bandAlpha);
end

hAnalytical=plot(ax,xAnalytical,yAnalytical,'--',...
    'LineWidth',curLineWidth,'Color',[0.45 0.45 0.45],...
    'DisplayName','Analytical');

for im=1:nMethod
    validMean=isfinite(meanData(im,:));
    hFEM(im)=scatter(ax,xData(validMean),meanData(im,validMean),...
        scatterSize,methodMarkers{im},'LineWidth',scatterLineWidth,...
        'MarkerEdgeColor',methodColors(im,:),...
        'DisplayName',methodNames{im});
end

if showGaussianLine
    hGaussian=yline(ax,3,':','LineWidth',1.2,'Color',[0 0 0],...
        'DisplayName',sprintf('    Kurtosis = 3'));
end

set(ax,'FontSize',axTickFontSize,'LineWidth',axLineWidth,...
    'FontName',axTickFontName,'Layer','top','box','on');
xlim(ax,[0 0.1125]);
ylim(ax,yLimits);
ax.XTick=0:0.025:0.125;
if ~showXTickLabels
    ax.XTickLabel=[];
end

if showLegend
    if showGaussianLine
        legendHandles=hGaussian;
        legendLabels={sprintf('   Kurtosis = 3')};
    else
        legendHandles=[hAnalytical,hFEM(1),hFEM(2)];
        legendLabels={'    Analytical',[methodNames{1}],...
            [methodNames{2}]};
    end
    legend(ax,legendHandles,legendLabels,'FontSize',lgdFontSize,...
        'FontName',lgdFontName,'Location','best','Box','off','Interpreter','tex');
else
    legend(ax,'off');
end

exportgraphics(f,[picBase '.png'],'Resolution',600);
saveas(f,[picBase '.fig']);
end


function h=draw_error_band(x,lowerBound,upperBound,bandColor,bandAlpha)
x=x(:).';
lowerBound=lowerBound(:).';
upperBound=upperBound(:).';
if isempty(x)
    h=gobjects(1);
    return
end

% Smooth the discrete SD envelope with a shape-preserving piecewise cubic
% interpolation. This avoids the visibly straight polygonal band while
% reducing the oscillatory overshoot common to unrestricted high-order fits.
if length(x)>=3
    xSmooth=linspace(min(x),max(x),400);
    lowerSmooth=interp1(x,lowerBound,xSmooth,'pchip');
    upperSmooth=interp1(x,upperBound,xSmooth,'pchip');
    % Guard against a locally inverted envelope after interpolation.
    bandLower=min(lowerSmooth,upperSmooth);
    bandUpper=max(lowerSmooth,upperSmooth);
else
    xSmooth=x;
    bandLower=lowerBound;
    bandUpper=upperBound;
end

h=fill([xSmooth fliplr(xSmooth)],[bandUpper fliplr(bandLower)],bandColor,...
    'FaceAlpha',bandAlpha,'EdgeColor','none');
end
