function t = round2_read_statistics(kind)
% Read-only access to the exact original statistics workbooks.
p=round2_paths();
switch lower(kind)
    case 'geom'
        name='Deformation_Statistics_Geom_ErrorBand.xlsx'; sheet='Statistics';
    case 'load'
        name='Deformation_Statistics_Load_ErrorBand.xlsx'; sheet='Statistics';
    case 'vol20'
        name='Deformation_Statistics_Geom_Vol20_ErrorBand.xlsx'; sheet='Summary';
    case 'point'
        name='Deformation_Statistics_Geom_r0_Comparison_ErrorBand.xlsx'; sheet='Summary';
    otherwise
        error('Unknown statistics dataset: %s',kind);
end
source=fullfile(p.sourceOutputDir,name);
assert(isfile(source),'Original FEM statistics workbook is missing: %s',source);
t=readtable(source,'Sheet',sheet);
if ismember('MeanStress',t.Properties.VariableNames)
    t=renamevars(t,{'MeanStress','MeanStress_SD'},{'Mean','Mean_SD'});
end
required={'Porosity','StressComponent','Variance','Variance_SD','Skewness','Skewness_SD','Kurtosis','Kurtosis_SD'};
if strcmpi(kind,'load'), required{1}='LoadRatio'; end
assert(all(ismember(required,t.Properties.VariableNames)),...
    'Unexpected source table columns in %s.',source);
t.SourceFile=repmat(string(source),height(t),1);
t.SourceSheet=repmat(string(sheet),height(t),1);
end
