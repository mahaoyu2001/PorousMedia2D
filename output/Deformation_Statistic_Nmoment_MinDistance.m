clear;

% Sensitivity of S22 variance and kurtosis to minimum pore-centre distance
% at 5% porosity. DATA_SOURCE has the same meaning as in the main comparison
% script: 'auto', 'field' or 'cache'.
DATA_SOURCE = 'auto';

suffixList = {'_a','_b','_c','_d','_e'};
modelNameList = {'2DVoid_Vol5_r0','2DVoid_Vol5'};
minDistanceOverRadius = [0 2];  % d_min/a; a is the pore radius.
porosity = 0.05;
component = 2;                  % S22 under uniaxial loading.
fieldName = 'S';
outputFrameTag = 2;

if numel(modelNameList) ~= numel(minDistanceOverRadius)
    error('modelNameList and minDistanceOverRadius must have equal lengths.');
end
[minDistanceOverRadius,order] = sort(minDistanceOverRadius);
modelNameList = modelNameList(order);

scriptDir = fileparts(mfilename('fullpath'));
fieldDir = fullfile(scriptDir,'FieldValue');
picDir = fullfile(scriptDir,'pic');
cacheFile = fullfile(scriptDir,'Minimum_Distance_Statistics.xlsx');
if ~isfolder(picDir)
    mkdir(picDir);
end
addpath(scriptDir);

nCase = numel(modelNameList);
nRealization = numel(suffixList);
rawFiles = strings(nCase*nRealization,1);
row = 0;
for ic = 1:nCase
    for ir = 1:nRealization
        row = row+1;
        rawFiles(row) = fullfile(fieldDir,[modelNameList{ic} ...
            suffixList{ir} '_' fieldName num2str(outputFrameTag) '.txt']);
    end
end
allRawFilesExist = all(isfile(rawFiles));

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
        error('Processed cache not found: %s',cacheFile);
    end
    summaryTable = readtable(cacheFile,'Sheet','Summary', ...
        'VariableNamingRule','preserve');
    [found,location] = ismember(minDistanceOverRadius, ...
        summaryTable.MinimumCenterDistanceOverRadius);
    if ~all(found)
        error('The cache does not contain every requested distance case.');
    end
    varianceMean = summaryTable.Variance(location).';
    varianceSD = summaryTable.Variance_SD(location).';
    kurtosisMean = summaryTable.Kurtosis(location).';
    kurtosisSD = summaryTable.Kurtosis_SD(location).';
    fprintf('Reading processed statistics from:\n  %s\n',cacheFile);
else
    varianceByRealization = nan(nRealization,nCase);
    kurtosisByRealization = nan(nRealization,nCase);
    for ic = 1:nCase
        for ir = 1:nRealization
            filePath = fullfile(fieldDir,[modelNameList{ic} ...
                suffixList{ir} '_' fieldName num2str(outputFrameTag) '.txt']);
            if ~isfile(filePath)
                error('Missing field file: %s',filePath);
            end
            imported = importdata(filePath);
            if isstruct(imported)
                data = imported.data;
            else
                data = imported;
            end
            if ~isnumeric(data) || size(data,2) < component+1
                error('Unexpected data format in %s.',filePath);
            end
            elementArea = data(:,end);
            totalArea = sum(elementArea);
            if ~isfinite(totalArea) || totalArea <= 0
                error('Invalid EVOL column in %s.',filePath);
            end
            weights = elementArea/totalArea;
            value = data(:,component);
            mu = sum(value.*weights);
            m2 = sum((value-mu).^2.*weights);
            m4 = sum((value-mu).^4.*weights);
            varianceByRealization(ir,ic) = m2;
            kurtosisByRealization(ir,ic) = m4/m2^2;
        end
    end

    varianceMean = mean(varianceByRealization,1);
    varianceSD = std(varianceByRealization,0,1);
    kurtosisMean = mean(kurtosisByRealization,1);
    kurtosisSD = std(kurtosisByRealization,0,1);

    MinimumCenterDistanceOverRadius = minDistanceOverRadius(:);
    ModelName = string(modelNameList(:));
    Porosity = repmat(porosity,nCase,1);
    NumRealizations = repmat(nRealization,nCase,1);
    Variance = varianceMean(:); Variance_SD = varianceSD(:);
    Kurtosis = kurtosisMean(:); Kurtosis_SD = kurtosisSD(:);
    summaryTable = table(MinimumCenterDistanceOverRadius,ModelName,Porosity, ...
        NumRealizations,Variance,Variance_SD,Kurtosis,Kurtosis_SD);

    MinimumCenterDistanceOverRadius = repelem(minDistanceOverRadius(:), ...
        nRealization,1);
    ModelName = repelem(string(modelNameList(:)),nRealization,1);
    Realization = repmat(string(suffixList(:)),nCase,1);
    Variance = varianceByRealization(:);
    Kurtosis = kurtosisByRealization(:);
    realizationTable = table(MinimumCenterDistanceOverRadius,ModelName, ...
        Realization,Variance,Kurtosis);

    if isfile(cacheFile)
        delete(cacheFile);
    end
    writetable(summaryTable,cacheFile,'Sheet','Summary');
    writetable(realizationTable,cacheFile,'Sheet','Realizations');
    fprintf('Updated processed cache:\n  %s\n',cacheFile);
end

%% Variance and kurtosis versus normalized minimum centre distance
[curLineWidth,axLineWidth,markerSize,axTickFontSize, ...
    axLabelFontSize,lgdFontSize,axTickFontName,axLabelFontName, ...
    lgdFontName,~] = plotParams('curLineWidth',1.2,'markerSize',3.8);

f = figure('Units','centimeters');
ax = axes(f,'Units','centimeters','Position',[2,2,6.4,4.95]);
set(ax,'FontSize',axTickFontSize,'LineWidth',axLineWidth, ...
    'FontName',axTickFontName,'Box','on');

varianceColor = [0.6350 0.0780 0.1840];
kurtosisColor = [0.0471 0.3059 0.6078];
yyaxis(ax,'left');
hVariance = errorbar(ax,minDistanceOverRadius,varianceMean,varianceSD,'o-', ...
    'LineWidth',curLineWidth,'Color',varianceColor, ...
    'MarkerSize',markerSize,'MarkerFaceColor',varianceColor,'CapSize',8, ...
    'DisplayName','  Var($\sigma_{\rm{22}}$)');
ylabel(ax,'Variance','FontSize',axLabelFontSize, ...
    'FontName',axLabelFontName,'Color','k');

yyaxis(ax,'right');
hKurtosis = errorbar(ax,minDistanceOverRadius,kurtosisMean,kurtosisSD,'s-', ...
    'LineWidth',curLineWidth,'Color',kurtosisColor, ...
    'MarkerSize',markerSize+0.7,'MarkerFaceColor',kurtosisColor, ...
    'CapSize',8,'DisplayName','Kurt($\sigma_{\rm{22}}$)');
ylabel(ax,'Kurtosis','FontSize',axLabelFontSize, ...
    'FontName',axLabelFontName,'Color','k');

xlim(ax,padded_x_limits(minDistanceOverRadius));
ax.XTick = minDistanceOverRadius;
xlabel(ax,'Minimum center distance $d_{\min}/a$','Interpreter','latex', ...
    'FontSize',axLabelFontSize,'FontName',axLabelFontName);
legend(ax,[hVariance,hKurtosis],'FontSize',lgdFontSize, ...
    'Interpreter','latex','FontName',lgdFontName,'Location','best', ...
    'Box','off');
yyaxis(ax,'left'); ax.YColor = varianceColor;
yyaxis(ax,'right'); ax.YColor = kurtosisColor;

picBase = fullfile(picDir,'var_kur_minimum_distance');
exportgraphics(f,[picBase '.png'],'Resolution',600);
saveas(f,[picBase '.fig']);


function limits = padded_x_limits(x)
%PADDED_X_LIMITS Add a small horizontal margin around discrete cases.
x = x(:).';
if numel(x) == 1
    margin = max(0.5,0.1*max(1,abs(x)));
else
    margin = 0.12*(max(x)-min(x));
end
limits = [min(x)-margin max(x)+margin];
end
