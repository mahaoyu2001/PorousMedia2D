clear;
r2=round2_paths(); % Isolated round-2 outputs; original FEM inputs are read-only.

% Extended porosity study (0.5%-20%) based on
% Deformation_Statistic_Nmoment_Geom_ErrorBand.m.
% The reference script and field files are not modified.

suffixList={'_a','_b','_c','_d','_e'};
modelNameList={'2DVoid_Vol05','2DVoid_Vol1','2DVoid_Vol2',...
    '2DVoid_Vol3','2DVoid_Vol5','2DVoid_Vol7_5','2DVoid_Vol10',...
    '2DVoid_Vol15','2DVoid_Vol20'};
vol=[0.5 1 2 3 5 7.5 10 15 20]/100;

% Data source switch:
%   'field' - calculate from FieldValue files and refresh the Excel cache.
%   'excel' - read cached statistics and draw figures without field files.
dataSource='excel';

Sigma=[0,1,0];
[kn,~,~]=cal_stress_coef(Sigma(1),Sigma(2));

components=[1 2 3];
componentNames=["S11","S22","S12"];
fieldName='S';
step=2;

% Centralized vertical-axis limits. Rows correspond to S11, S22 and S12.
% The variance limits are extended to include porosities up to 20%.
varYLim=[0 0.135;...
         0 0.55;...
         0 0.11];
skeYLim=repmat([-1.5 1.5],3,1);
kurYLim=[0 100;...
         0 115;...
         0 90];

bandColor=[0.86 0.36 0.43];
bandAlpha=0.18;

scriptDir=fileparts(mfilename('fullpath'));
fieldDir=r2.fieldDir;
picDir=r2.figureDir;
excelFile=fullfile(r2.tableDir,...
    'Deformation_Statistics_Geom_Vol20_ErrorBand.xlsx');
if ~isfolder(picDir)
    mkdir(picDir);
end
addpath(scriptDir);

if length(modelNameList)~=length(vol)
    error('modelNameList and vol must have the same length.');
end

nComponent=length(components);
nModel=length(modelNameList);
nRealization=length(suffixList);
dataSource=lower(strtrim(dataSource));

if strcmp(dataSource,'field')
    meanMean=zeros(nComponent,nModel);
    meanSD=zeros(nComponent,nModel);
    varMean=zeros(nComponent,nModel);
    varSD=zeros(nComponent,nModel);
    skeMean=zeros(nComponent,nModel);
    skeSD=zeros(nComponent,nModel);
    kurMean=zeros(nComponent,nModel);
    kurSD=zeros(nComponent,nModel);

    mu_RVE=zeros(nRealization,nModel,nComponent);
    moment2=zeros(nRealization,nModel,nComponent);
    moment3=zeros(nRealization,nModel,nComponent);
    moment4=zeros(nRealization,nModel,nComponent);
    ske_RVE=zeros(nRealization,nModel,nComponent);
    kur_RVE=zeros(nRealization,nModel,nComponent);

    %% Check all inputs before reading the large field files
    for nf=1:nModel
        for is=1:nRealization
            fileName=[modelNameList{nf} suffixList{is} '_' ...
                fieldName num2str(step) '.txt'];
            filePath=fullfile(fieldDir,fileName);
            if ~isfile(filePath)
                error('Missing field file: %s',filePath);
            end
        end
    end

    %% Import every field file once and calculate all stress components
    fprintf('Reading %d field files from:\n  %s\n',...
        nModel*nRealization,fieldDir);
    for nf=1:nModel
        for is=1:nRealization
            fileName=[modelNameList{nf} suffixList{is} '_' ...
                fieldName num2str(step) '.txt'];
            filePath=fullfile(fieldDir,fileName);
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

            for ia=1:nComponent
                component=components(ia);
                phaseValue=Data(:,component);
                mu=sum(phaseValue.*weights);
                centralValue=phaseValue-mu;

                mu_RVE(is,nf,ia)=mu;
                moment2(is,nf,ia)=sum(centralValue.^2.*weights);
                moment3(is,nf,ia)=sum(centralValue.^3.*weights);
                moment4(is,nf,ia)=sum(centralValue.^4.*weights);
                ske_RVE(is,nf,ia)=moment3(is,nf,ia)/...
                    moment2(is,nf,ia)^1.5;
                kur_RVE(is,nf,ia)=moment4(is,nf,ia)/...
                    moment2(is,nf,ia)^2;
            end
            fprintf('Loaded: %s\n',fileName);
            clear Data Evol weights phaseValue centralValue
        end
    end

    for ia=1:nComponent
        meanMean(ia,:)=mean(mu_RVE(:,:,ia),1);
        meanSD(ia,:)=std(mu_RVE(:,:,ia),0,1);
        varMean(ia,:)=mean(moment2(:,:,ia),1);
        varSD(ia,:)=std(moment2(:,:,ia),0,1);
        skeMean(ia,:)=mean(ske_RVE(:,:,ia),1);
        skeSD(ia,:)=std(ske_RVE(:,:,ia),0,1);
        kurMean(ia,:)=mean(kur_RVE(:,:,ia),1);
        kurSD(ia,:)=std(kur_RVE(:,:,ia),0,1);
    end

    %% Store plot statistics and realization-level data in Excel
    Porosity=repmat(vol(:),nComponent,1);
    ModelName=repmat(string(modelNameList(:)),nComponent,1);
    StressComponent=strings(nComponent*nModel,1);
    for ia=1:nComponent
        rowIndex=(ia-1)*nModel+(1:nModel);
        StressComponent(rowIndex)=componentNames(ia);
    end
    NumRealizations=repmat(nRealization,nComponent*nModel,1);
    Mean=reshape(meanMean.',[],1);
    Mean_SD=reshape(meanSD.',[],1);
    Variance=reshape(varMean.',[],1);
    Variance_SD=reshape(varSD.',[],1);
    Skewness=reshape(skeMean.',[],1);
    Skewness_SD=reshape(skeSD.',[],1);
    Kurtosis=reshape(kurMean.',[],1);
    Kurtosis_SD=reshape(kurSD.',[],1);
    summaryTable=table(Porosity,ModelName,StressComponent,NumRealizations,...
        Mean,Mean_SD,Variance,Variance_SD,Skewness,Skewness_SD,...
        Kurtosis,Kurtosis_SD);

    Porosity=repmat(repelem(vol(:),nRealization),nComponent,1);
    ModelName=repmat(repelem(string(modelNameList(:)),nRealization),...
        nComponent,1);
    Realization=repmat(string(suffixList(:)),nModel*nComponent,1);
    StressComponent=repelem(componentNames(:),nRealization*nModel);
    Mean=mu_RVE(:);
    Variance=moment2(:);
    Skewness=ske_RVE(:);
    Kurtosis=kur_RVE(:);
    realizationTable=table(Porosity,ModelName,Realization,...
        StressComponent,Mean,Variance,Skewness,Kurtosis);

    if isfile(excelFile)
        delete(excelFile);
    end
    writetable(summaryTable,excelFile,'Sheet','Summary');
    writetable(realizationTable,excelFile,'Sheet','Realizations');
    fprintf('Statistics cache refreshed:\n  %s\n',excelFile);

elseif strcmp(dataSource,'excel')
    if ~isfile(fullfile(r2.sourceOutputDir,'Deformation_Statistics_Geom_Vol20_ErrorBand.xlsx'))
        error(['Statistics cache not found: %s\nRun once with ' ...
            'dataSource=''''field'''' first.'],excelFile);
    end
    summaryTable=round2_read_statistics('vol20');
    [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD]=...
        read_summary_cache(summaryTable,vol,modelNameList,componentNames);
    fprintf('Statistics loaded from cache:\n  %s\n',excelFile);

else
    error('dataSource must be either ''field'' or ''excel''.');
end

[curLineWidth,axLineWidth,~,axTickFontSize,~,lgdFontSize,...
    axTickFontName,~,lgdFontName,~]=plotParams();

%% Draw variance, skewness and kurtosis for all components
for ia=1:nComponent
    component=components(ia);
    i=component;

    % Moment 2: variance
    f=figure;
    f.Units='centimeters';
    picName=fullfile(picDir,...
        ['var_vol_stress' num2str(i) '_Vol20_errorband']);
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.25];

    m2=varMean(ia,:);
    sd2=varSD(ia,:);
    hBand=draw_error_band(vol,max(0,m2-sd2),m2+sd2,...
        bandColor,bandAlpha);
    hold on

    k2=kn(i,2);
    x=0:0.001:0.225;
    [~,r2Var,~,~]=cal_stress_statistics(Sigma,x);
    y_fit=r2Var(component,:);
    hFit=plot(x,y_fit,LineWidth=curLineWidth,LineStyle="--",...
        Color=[0.0471 0.3059 0.6078 0.6],DisplayName='  Analytical');
    hFEM=scatter(vol,m2,80,'+',LineWidth=1.6,...
        MarkerEdgeColor=[0.6350 0.0780 0.1840],DisplayName='FEM results');

    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,...
        FontName=axTickFontName,Layer='top');
    xlim([0 0.2125]);
    ylim(varYLim(ia,:));
    ax.XTick=0:0.05:0.20;
    if ia<nComponent
        ax.XTickLabel=[];
    end
    if ia==1
        legend([hFit,hFEM,hBand],...
            {'  Analytical','FEM results','   SD band'},...
            FontSize=lgdFontSize,FontName=lgdFontName,...
            Location='best',Box='off');
    end

    exportgraphics(f,[picName '.png'],Resolution=600);
    saveas(f,[picName '.fig']);

    % Moment 3: skewness
    f=figure;
    f.Units='centimeters';
    picName=fullfile(picDir,...
        ['ske_vol_stress' num2str(i) '_Vol20_errorband']);
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.25];

    m3=skeMean(ia,:);
    sd3=skeSD(ia,:);
    hBand=draw_error_band(vol,m3-sd3,m3+sd3,bandColor,bandAlpha);
    hold on

    k3=kn(i,3);
    x=0:0.001:0.225;
    [~,~,r2Skew,~]=cal_stress_statistics(Sigma,x);
    y_fit=r2Skew(component,:);
    hFit=plot(x,y_fit,LineWidth=curLineWidth,LineStyle="--",...
        Color=[0.0471 0.3059 0.6078 0.6],DisplayName='  Analytical');
    hFEM=scatter(vol,m3,100,'+',LineWidth=1.8,...
        MarkerEdgeColor=[0.6350 0.0780 0.1840],DisplayName='FEM results');

    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,...
        FontName=axTickFontName,Layer='top');
    xlim([0 0.2125]);
    ylim(skeYLim(ia,:));
    ax.XTick=0:0.05:0.20;
    ax.YTick=-10:5:10;
    if ia<nComponent
        ax.XTickLabel=[];
    end
    if ia==1
        legend([hFit,hFEM,hBand],...
            {'  Analytical','FEM results','   SD band'},...
            FontSize=lgdFontSize,FontName=lgdFontName,...
            Location='best',Box='off');
    end

    exportgraphics(f,[picName '.png'],Resolution=600);
    saveas(f,[picName '.fig']);

    % Moment 4: kurtosis
    f=figure;
    f.Units='centimeters';
    picName=fullfile(picDir,...
        ['kur_vol_stress' num2str(i) '_Vol20_errorband']);
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.25];

    m4=kurMean(ia,:);
    sd4=kurSD(ia,:);
    draw_error_band(vol,max(0,m4-sd4),m4+sd4,bandColor,bandAlpha);
    hold on

    k4=kn(i,4);
    k2=kn(i,2);
    x=0:0.001:0.5;
    [~,~,~,r2Kurt]=cal_stress_statistics(Sigma,x);
    y_fit=r2Kurt(component,:);
    plot(x,y_fit,LineWidth=curLineWidth,LineStyle="--",...
        Color=[0.0471 0.3059 0.6078 0.6],DisplayName='  Analytical');
    scatter(vol,m4,80,'+',LineWidth=1.6,...
        MarkerEdgeColor=[0.6350 0.0780 0.1840],DisplayName='FEM results');
    yl=yline(3,LineStyle=":",LineWidth=1.2,Color=[0 0 0],...
        DisplayName=sprintf('         Kurtosis = 3\n(Gaussian distribution)'));

    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,...
        FontName=axTickFontName,Layer='top');
    xlim([0 0.2125]);
    ylim(kurYLim(ia,:));
    ax.XTick=0:0.05:0.20;
    if ia<nComponent
        ax.XTickLabel=[];
    end
    if ia==1
        lgd=legend(yl,FontSize=lgdFontSize,FontName=lgdFontName,...
            Location='north',Box='off');
        drawnow;
        lgdText=findall(lgd,'Type','Text');
        set(lgdText,'HorizontalAlignment','center');
    end

    exportgraphics(f,[picName '.png'],Resolution=600);
    saveas(f,[picName '.fig']);
end

fprintf('\nFigures saved in:\n  %s\n',picDir);


function [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD]=...
    read_summary_cache(summaryTable,vol,modelNameList,componentNames)
% Validate cached rows against the active model/component configuration.
requiredVariables={'Porosity','ModelName','StressComponent',...
    'NumRealizations','Mean','Mean_SD','Variance','Variance_SD',...
    'Skewness','Skewness_SD','Kurtosis','Kurtosis_SD'};
missingVariables=setdiff(requiredVariables,summaryTable.Properties.VariableNames);
if ~isempty(missingVariables)
    error('Summary sheet is missing column(s): %s',...
        strjoin(missingVariables,', '));
end

nComponent=length(componentNames);
nModel=length(modelNameList);
meanMean=nan(nComponent,nModel);
meanSD=nan(nComponent,nModel);
varMean=nan(nComponent,nModel);
varSD=nan(nComponent,nModel);
skeMean=nan(nComponent,nModel);
skeSD=nan(nComponent,nModel);
kurMean=nan(nComponent,nModel);
kurSD=nan(nComponent,nModel);

for ia=1:nComponent
    for nf=1:nModel
        rowMask=abs(summaryTable.Porosity-vol(nf))<1e-12 & ...
            string(summaryTable.ModelName)==string(modelNameList{nf}) & ...
            string(summaryTable.StressComponent)==componentNames(ia);
        if nnz(rowMask)~=1
            error(['Cache must contain exactly one matching Summary row for ' ...
                '%s, %s and porosity %g; found %d. Run with ' ...
                'dataSource=''''field'''' to refresh it.'],...
                modelNameList{nf},componentNames(ia),vol(nf),nnz(rowMask));
        end
        meanMean(ia,nf)=summaryTable.Mean(rowMask);
        meanSD(ia,nf)=summaryTable.Mean_SD(rowMask);
        varMean(ia,nf)=summaryTable.Variance(rowMask);
        varSD(ia,nf)=summaryTable.Variance_SD(rowMask);
        skeMean(ia,nf)=summaryTable.Skewness(rowMask);
        skeSD(ia,nf)=summaryTable.Skewness_SD(rowMask);
        kurMean(ia,nf)=summaryTable.Kurtosis(rowMask);
        kurSD(ia,nf)=summaryTable.Kurtosis_SD(rowMask);
    end
end

cachedValues=[meanMean(:);meanSD(:);varMean(:);varSD(:);...
    skeMean(:);skeSD(:);kurMean(:);kurSD(:)];
if any(~isfinite(cachedValues)) || any(meanSD(:)<0) || ...
        any(varSD(:)<0) || any(skeSD(:)<0) || any(kurSD(:)<0)
    error('Summary sheet contains invalid statistic or SD values.');
end
end


function h=draw_error_band(x,lowerBound,upperBound,bandColor,bandAlpha)
x=x(:).';
lowerBound=lowerBound(:).';
upperBound=upperBound(:).';

if length(x)>=3
    xSmooth=linspace(min(x),max(x),400);
    lowerSmooth=interp1(x,lowerBound,xSmooth,'pchip');
    upperSmooth=interp1(x,upperBound,xSmooth,'pchip');
    bandLower=min(lowerSmooth,upperSmooth);
    bandUpper=max(lowerSmooth,upperSmooth);
else
    xSmooth=x;
    bandLower=lowerBound;
    bandUpper=upperBound;
end

h=fill([xSmooth fliplr(xSmooth)],[bandUpper fliplr(bandLower)],...
    bandColor,FaceAlpha=bandAlpha,EdgeColor='none',...
    HandleVisibility='off');
end
