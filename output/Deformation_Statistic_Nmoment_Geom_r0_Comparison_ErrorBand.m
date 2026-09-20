clear;

% Compare analytical predictions with two FEM microstructure generators:
%   (1) hard-core pores, whose centres are separated by at least 2a;
%   (2) Poisson pores with r_min = 0, for which overlap is allowed.
%
% DATA_SOURCE controls how statistics are obtained:
%   'auto'  - use all raw FieldValue files when present, otherwise cache;
%   'field' - recompute from raw files and refresh the Excel cache;
%   'cache' - draw directly from the bundled processed Excel workbook.

DATA_SOURCE = 'auto';

suffixList = {'_a','_b','_c','_d','_e'};
baseModelNameList = {'2DVoid_Vol05','2DVoid_Vol1','2DVoid_Vol2', ...
    '2DVoid_Vol3','2DVoid_Vol5','2DVoid_Vol7_5','2DVoid_Vol10'};
vol = [0.5 1 2 3 5 7.5 10] / 100;

methodFileTags = {'','_r0'};
methodNames = {'Hard-core FEM','Poisson FEM'};
methodShortNames = {'Hard-core','r_min=0'};

Sigma = [0;1;0];
components = [1 2 3];
componentNames = ["S11","S22","S12"];
fieldName = 'S';
outputFrameTag = 2;

% Publication-oriented vertical-axis limits for S11, S22 and S12.
varYLim = [0 0.045; 0 0.20; 0 0.04];
skeYLim = [-10 10; -10 10; -10 10];
kurYLim = [0 100; 0 115; 0 90];

scriptDir = fileparts(mfilename('fullpath'));
fieldDir = fullfile(scriptDir,'FieldValue');
picDir = fullfile(scriptDir,'pic');
cacheFile = fullfile(scriptDir, ...
    'Deformation_Statistics_Geom_r0_Comparison_ErrorBand.xlsx');
if ~isfolder(picDir)
    mkdir(picDir);
end
addpath(scriptDir);

nMethod = numel(methodFileTags);
nRealization = numel(suffixList);
nModel = numel(baseModelNameList);
nComponent = numel(components);

expectedFiles = expected_field_files(fieldDir,baseModelNameList, ...
    methodFileTags,suffixList,fieldName,outputFrameTag);
allRawFilesExist = all(isfile(expectedFiles));
switch lower(DATA_SOURCE)
    case 'auto'
        useCache = ~allRawFilesExist;
    case 'cache'
        useCache = true;
    case 'field'
        useCache = false;
    otherwise
        error('DATA_SOURCE must be ''auto'', ''field'' or ''cache''.')
end

if useCache
    if ~isfile(cacheFile)
        error(['Processed cache not found: %s\nSet DATA_SOURCE=''field'' ' ...
            'after generating the FieldValue files.'],cacheFile);
    end
    fprintf('Reading processed statistics from:\n  %s\n',cacheFile);
    summaryTable = readtable(cacheFile,'Sheet','Summary', ...
        'VariableNamingRule','preserve');
    [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD, ...
        sampleCount] = summary_to_arrays(summaryTable,methodShortNames, ...
        componentNames,vol);
else
    fprintf('Computing statistics from:\n  %s\n',fieldDir);
    [muRVE,moment2,skewnessRVE,kurtosisRVE,fileAvailable] = ...
        read_field_statistics(fieldDir,baseModelNameList,methodFileTags, ...
        suffixList,components,fieldName,outputFrameTag);
    if ~any(fileAvailable,'all')
        error('No readable FieldValue files were found in %s.',fieldDir);
    end
    [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD, ...
        sampleCount] = aggregate_realizations(muRVE,moment2, ...
        skewnessRVE,kurtosisRVE);

    [summaryTable,realizationTable] = build_output_tables( ...
        methodShortNames,componentNames,vol,suffixList,fileAvailable, ...
        muRVE,moment2,skewnessRVE,kurtosisRVE,meanMean,meanSD, ...
        varMean,varSD,skeMean,skeSD,kurMean,kurSD,sampleCount);
    if isfile(cacheFile)
        delete(cacheFile);
    end
    writetable(summaryTable,cacheFile,'Sheet','Summary');
    writetable(realizationTable,cacheFile,'Sheet','Realizations');
    fprintf('Updated processed cache:\n  %s\n',cacheFile);
end

%% Analytical predictions and figure style
xAnalytical = 0.0001:0.0001:0.1125;
[~,varAnalytical,skeAnalytical,kurAnalytical] = ...
    cal_stress_statistics(Sigma,xAnalytical);

[curLineWidth,axLineWidth,~,axTickFontSize,~,lgdFontSize, ...
    axTickFontName,~,lgdFontName,~] = plotParams( ...
    'curLineWidth',1.25,'markerSize',5.5);

methodColors = [0.6350 0.0780 0.1840; ...
                0.3010 0.7450 0.9330];
methodBandColors = [0.86 0.36 0.43; ...
                    0.55 0.80 0.94];
methodMarkers = {'+','o'};
bandAlpha = 0.18;

%% Variance, skewness and kurtosis comparison figures
for ia = 1:nComponent
    component = components(ia);

    plot_comparison(vol,squeeze(varMean(:,ia,:)),squeeze(varSD(:,ia,:)), ...
        xAnalytical,varAnalytical(component,:),methodNames,methodColors, ...
        methodBandColors,methodMarkers,bandAlpha,varYLim(ia,:), ...
        fullfile(picDir,['var_vol_stress' num2str(component) ...
        '_r0_comparison_errorband']),false,false,80,1.6,ia==1, ...
        ia==nComponent,curLineWidth,axLineWidth,axTickFontSize, ...
        lgdFontSize,axTickFontName,lgdFontName);

    plot_comparison(vol,squeeze(skeMean(:,ia,:)),squeeze(skeSD(:,ia,:)), ...
        xAnalytical,skeAnalytical(component,:),methodNames,methodColors, ...
        methodBandColors,methodMarkers,bandAlpha,skeYLim(ia,:), ...
        fullfile(picDir,['ske_vol_stress' num2str(component) ...
        '_r0_comparison_errorband']),false,false,100,1.8,ia==1, ...
        ia==nComponent,curLineWidth,axLineWidth,axTickFontSize, ...
        lgdFontSize,axTickFontName,lgdFontName);

    % As in the paper layout, only the first kurtosis panel has a legend,
    % and that legend contains only the Gaussian reference line.
    plot_comparison(vol,squeeze(kurMean(:,ia,:)),squeeze(kurSD(:,ia,:)), ...
        xAnalytical,kurAnalytical(component,:),methodNames,methodColors, ...
        methodBandColors,methodMarkers,bandAlpha,kurYLim(ia,:), ...
        fullfile(picDir,['kur_vol_stress' num2str(component) ...
        '_r0_comparison_errorband']),true,true,80,1.6,ia==1, ...
        ia==nComponent,curLineWidth,axLineWidth,axTickFontSize, ...
        lgdFontSize,axTickFontName,lgdFontName);
end

fprintf('Finished. Figures were saved in:\n  %s\n',picDir);


function files = expected_field_files(fieldDir,baseNames,methodTags, ...
    suffixes,fieldName,frameTag)
%EXPECTED_FIELD_FILES Build the complete set of raw-data paths.
n = numel(baseNames)*numel(methodTags)*numel(suffixes);
files = strings(n,1);
row = 0;
for im = 1:numel(methodTags)
    for nf = 1:numel(baseNames)
        for ir = 1:numel(suffixes)
            row = row+1;
            files(row) = fullfile(fieldDir,[baseNames{nf} methodTags{im} ...
                suffixes{ir} '_' fieldName num2str(frameTag) '.txt']);
        end
    end
end
end


function [muRVE,moment2,skewnessRVE,kurtosisRVE,fileAvailable] = ...
    read_field_statistics(fieldDir,baseNames,methodTags,suffixes, ...
    components,fieldName,frameTag)
%READ_FIELD_STATISTICS Calculate area-weighted moments for every realization.
nMethod = numel(methodTags);
nRealization = numel(suffixes);
nModel = numel(baseNames);
nComponent = numel(components);
muRVE = nan(nMethod,nRealization,nModel,nComponent);
moment2 = nan(size(muRVE));
skewnessRVE = nan(size(muRVE));
kurtosisRVE = nan(size(muRVE));
fileAvailable = false(nMethod,nRealization,nModel);

for im = 1:nMethod
    for nf = 1:nModel
        modelName = [baseNames{nf} methodTags{im}];
        for ir = 1:nRealization
            fileName = [modelName suffixes{ir} '_' fieldName ...
                num2str(frameTag) '.txt'];
            filePath = fullfile(fieldDir,fileName);
            if ~isfile(filePath)
                warning('Skipping missing field file: %s',filePath);
                continue
            end
            imported = importdata(filePath);
            if isstruct(imported)
                data = imported.data;
            else
                data = imported;
            end
            if ~isnumeric(data) || size(data,2) < 4
                error('Unexpected data format in %s.',filePath);
            end

            elementArea = data(:,end);
            totalArea = sum(elementArea);
            if ~isfinite(totalArea) || totalArea <= 0 || ...
                    any(~isfinite(elementArea))
                error('Invalid EVOL column in %s.',filePath);
            end
            weights = elementArea/totalArea;
            fileAvailable(im,ir,nf) = true;

            for ia = 1:nComponent
                value = data(:,components(ia));
                mu = sum(value.*weights);
                centered = value-mu;
                m2 = sum(centered.^2.*weights);
                m3 = sum(centered.^3.*weights);
                m4 = sum(centered.^4.*weights);
                muRVE(im,ir,nf,ia) = mu;
                moment2(im,ir,nf,ia) = m2;
                if m2 > 0
                    skewnessRVE(im,ir,nf,ia) = m3/m2^1.5;
                    kurtosisRVE(im,ir,nf,ia) = m4/m2^2;
                end
            end
            fprintf('Loaded: %s\n',fileName);
        end
    end
end
end


function [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD, ...
    sampleCount] = aggregate_realizations(muRVE,moment2,skewnessRVE, ...
    kurtosisRVE)
%AGGREGATE_REALIZATIONS Mean and sample SD across available realizations.
[nMethod,~,nModel,nComponent] = size(muRVE);
meanMean = nan(nMethod,nComponent,nModel);
meanSD = nan(size(meanMean));
varMean = nan(size(meanMean));
varSD = nan(size(meanMean));
skeMean = nan(size(meanMean));
skeSD = nan(size(meanMean));
kurMean = nan(size(meanMean));
kurSD = nan(size(meanMean));
sampleCount = zeros(size(meanMean));

for im = 1:nMethod
    for ia = 1:nComponent
        [m,s,n] = column_mean_sd(squeeze(muRVE(im,:,:,ia)));
        meanMean(im,ia,:) = reshape(m,1,1,[]);
        meanSD(im,ia,:) = reshape(s,1,1,[]);
        sampleCount(im,ia,:) = reshape(n,1,1,[]);
        [m,s] = column_mean_sd(squeeze(moment2(im,:,:,ia)));
        varMean(im,ia,:) = reshape(m,1,1,[]);
        varSD(im,ia,:) = reshape(s,1,1,[]);
        [m,s] = column_mean_sd(squeeze(skewnessRVE(im,:,:,ia)));
        skeMean(im,ia,:) = reshape(m,1,1,[]);
        skeSD(im,ia,:) = reshape(s,1,1,[]);
        [m,s] = column_mean_sd(squeeze(kurtosisRVE(im,:,:,ia)));
        kurMean(im,ia,:) = reshape(m,1,1,[]);
        kurSD(im,ia,:) = reshape(s,1,1,[]);
    end
end
end


function [valueMean,valueSD,valueN] = column_mean_sd(values)
%COLUMN_MEAN_SD Mean and sample SD by column, ignoring missing realizations.
valid = isfinite(values);
valueN = sum(valid,1);
values(~valid) = 0;
valueMean = sum(values,1)./valueN;
valueMean(valueN==0) = nan;
deviation = values-valueMean;
deviation(~valid) = 0;
valueSD = sqrt(sum(deviation.^2,1)./(valueN-1));
valueSD(valueN<2) = nan;
end


function [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD, ...
    sampleCount] = summary_to_arrays(t,methodNames,componentNames,porosity)
%SUMMARY_TO_ARRAYS Restore plotting arrays from the processed workbook.
nMethod = numel(methodNames);
nComponent = numel(componentNames);
nModel = numel(porosity);
arraySize = [nMethod nComponent nModel];
meanMean = nan(arraySize); meanSD = nan(arraySize);
varMean = nan(arraySize); varSD = nan(arraySize);
skeMean = nan(arraySize); skeSD = nan(arraySize);
kurMean = nan(arraySize); kurSD = nan(arraySize);
sampleCount = zeros(arraySize);

methods = string(t.Method);
components = string(t.StressComponent);
for im = 1:nMethod
    for ia = 1:nComponent
        for nf = 1:nModel
            row = methods==string(methodNames{im}) & ...
                components==componentNames(ia) & ...
                abs(t.Porosity-porosity(nf)) < 1e-12;
            if nnz(row) ~= 1
                error('Cache row not unique for %s, %s, porosity %.4g.', ...
                    methodNames{im},componentNames(ia),porosity(nf));
            end
            meanMean(im,ia,nf) = t.MeanStress(row);
            meanSD(im,ia,nf) = t.MeanStress_SD(row);
            varMean(im,ia,nf) = t.Variance(row);
            varSD(im,ia,nf) = t.Variance_SD(row);
            skeMean(im,ia,nf) = t.Skewness(row);
            skeSD(im,ia,nf) = t.Skewness_SD(row);
            kurMean(im,ia,nf) = t.Kurtosis(row);
            kurSD(im,ia,nf) = t.Kurtosis_SD(row);
            sampleCount(im,ia,nf) = t.NumRealizations(row);
        end
    end
end
end


function [summaryTable,realizationTable] = build_output_tables( ...
    methodNames,componentNames,porosity,suffixes,fileAvailable,muRVE, ...
    moment2,skewnessRVE,kurtosisRVE,meanMean,meanSD,varMean,varSD, ...
    skeMean,skeSD,kurMean,kurSD,sampleCount)
%BUILD_OUTPUT_TABLES Make compact summary and realization-level tables.
nMethod = numel(methodNames);
nComponent = numel(componentNames);
nModel = numel(porosity);
nRealization = numel(suffixes);

n = nMethod*nComponent*nModel;
Method = strings(n,1); StressComponent = strings(n,1);
Porosity = zeros(n,1); NumRealizations = zeros(n,1);
MeanStress = nan(n,1); MeanStress_SD = nan(n,1);
Variance = nan(n,1); Variance_SD = nan(n,1);
Skewness = nan(n,1); Skewness_SD = nan(n,1);
Kurtosis = nan(n,1); Kurtosis_SD = nan(n,1);
row = 0;
for im = 1:nMethod
    for ia = 1:nComponent
        for nf = 1:nModel
            row = row+1;
            Method(row) = methodNames{im};
            StressComponent(row) = componentNames(ia);
            Porosity(row) = porosity(nf);
            NumRealizations(row) = sampleCount(im,ia,nf);
            MeanStress(row) = meanMean(im,ia,nf);
            MeanStress_SD(row) = meanSD(im,ia,nf);
            Variance(row) = varMean(im,ia,nf);
            Variance_SD(row) = varSD(im,ia,nf);
            Skewness(row) = skeMean(im,ia,nf);
            Skewness_SD(row) = skeSD(im,ia,nf);
            Kurtosis(row) = kurMean(im,ia,nf);
            Kurtosis_SD(row) = kurSD(im,ia,nf);
        end
    end
end
summaryTable = table(Method,StressComponent,Porosity,NumRealizations, ...
    MeanStress,MeanStress_SD,Variance,Variance_SD,Skewness,Skewness_SD, ...
    Kurtosis,Kurtosis_SD);

n = nMethod*nComponent*nModel*nRealization;
Method = strings(n,1); StressComponent = strings(n,1);
Porosity = zeros(n,1); Realization = strings(n,1);
FileAvailable = false(n,1); MeanStress = nan(n,1);
Variance = nan(n,1); Skewness = nan(n,1); Kurtosis = nan(n,1);
row = 0;
for im = 1:nMethod
    for ia = 1:nComponent
        for nf = 1:nModel
            for ir = 1:nRealization
                row = row+1;
                Method(row) = methodNames{im};
                StressComponent(row) = componentNames(ia);
                Porosity(row) = porosity(nf);
                Realization(row) = erase(suffixes{ir},'_');
                FileAvailable(row) = fileAvailable(im,ir,nf);
                MeanStress(row) = muRVE(im,ir,nf,ia);
                Variance(row) = moment2(im,ir,nf,ia);
                Skewness(row) = skewnessRVE(im,ir,nf,ia);
                Kurtosis(row) = kurtosisRVE(im,ir,nf,ia);
            end
        end
    end
end
realizationTable = table(Method,StressComponent,Porosity,Realization, ...
    FileAvailable,MeanStress,Variance,Skewness,Kurtosis);
end


function plot_comparison(xData,meanData,sdData,xAnalytical,yAnalytical, ...
    methodNames,methodColors,methodBandColors,methodMarkers,bandAlpha, ...
    yLimits,picBase,nonnegativeBand,showGaussianLine,scatterSize, ...
    scatterLineWidth,showLegend,showXTickLabels,curLineWidth, ...
    axLineWidth,axTickFontSize,lgdFontSize,axTickFontName,lgdFontName)
%PLOT_COMPARISON Draw analytical line, FEM markers and smooth +/-1 SD bands.
f = figure('Units','centimeters');
ax = axes(f,'Units','centimeters','Position',[2,2,5.8,4.25]);
hold(ax,'on');

nMethod = size(meanData,1);
hFEM = gobjects(1,nMethod);
for im = 1:nMethod
    valid = isfinite(meanData(im,:)) & isfinite(sdData(im,:));
    lower = meanData(im,valid)-sdData(im,valid);
    if nonnegativeBand
        lower = max(0,lower);
    end
    upper = meanData(im,valid)+sdData(im,valid);
    draw_error_band(ax,xData(valid),lower,upper, ...
        methodBandColors(im,:),bandAlpha);
end

hAnalytical = plot(ax,xAnalytical,yAnalytical,'--', ...
    'LineWidth',curLineWidth,'Color',[0.45 0.45 0.45], ...
    'DisplayName','Analytical');
for im = 1:nMethod
    valid = isfinite(meanData(im,:));
    hFEM(im) = scatter(ax,xData(valid),meanData(im,valid), ...
        scatterSize,methodMarkers{im},'LineWidth',scatterLineWidth, ...
        'MarkerEdgeColor',methodColors(im,:), ...
        'DisplayName',methodNames{im});
end
if showGaussianLine
    hGaussian = yline(ax,3,':','LineWidth',1.2,'Color',[0 0 0], ...
        'DisplayName','Kurtosis = 3');
end

set(ax,'FontSize',axTickFontSize,'LineWidth',axLineWidth, ...
    'FontName',axTickFontName,'Layer','top','Box','on');
xlim(ax,[0 0.1125]);
ylim(ax,yLimits);
ax.XTick = 0:0.025:0.125;
if ~showXTickLabels
    ax.XTickLabel = [];
end
if showLegend
    if showGaussianLine
        legend(ax,hGaussian,{'   Kurtosis = 3'}, ...
            'FontSize',lgdFontSize,'FontName',lgdFontName, ...
            'Location','best','Box','off','Interpreter','tex');
    else
        legend(ax,[hAnalytical,hFEM], ...
            {'    Analytical',methodNames{1},methodNames{2}}, ...
            'FontSize',lgdFontSize,'FontName',lgdFontName, ...
            'Location','best','Box','off','Interpreter','tex');
    end
else
    legend(ax,'off');
end

exportgraphics(f,[picBase '.png'],'Resolution',600);
saveas(f,[picBase '.fig']);
end


function h = draw_error_band(ax,x,lower,upper,bandColor,bandAlpha)
%DRAW_ERROR_BAND Draw a shape-preserving cubic interpolation of the envelope.
x = x(:).'; lower = lower(:).'; upper = upper(:).';
if isempty(x)
    h = gobjects(1);
    return
end
if numel(x) >= 3
    xSmooth = linspace(min(x),max(x),400);
    lower = interp1(x,lower,xSmooth,'pchip');
    upper = interp1(x,upper,xSmooth,'pchip');
    x = xSmooth;
end
bandLower = min(lower,upper);
bandUpper = max(lower,upper);
h = fill(ax,[x fliplr(x)],[bandUpper fliplr(bandLower)],bandColor, ...
    'FaceAlpha',bandAlpha,'EdgeColor','none','HandleVisibility','off');
end
