function run_round2(mode)
% Modes: single (default), tables, check, field_distributions.
% cached and legacy_cached are retained aliases for single-panel output.
if nargin<1, mode='single'; end
p=round2_paths(); oldDir=pwd; oldPath=path;
restore=onCleanup(@() restore_environment(oldDir,oldPath)); %#ok<NASGU>
cd(p.programDir); addpath(p.programDir,'-begin');
logName=fullfile(p.logDir,['Round2_' char(mode) '_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.log']);
diary(logName); cleanupDiary=onCleanup(@() diary('off')); %#ok<NASGU>
test_round2();
switch lower(mode)
    case 'check'
        % Validation only; no large field files are read.
    case 'tables'
        export_round2_comparisons();
    case {'single','cached','legacy_cached'}
        export_round2_comparisons(); % Numeric exports only; no grouped figures.
        jobs={'Deformation_Statistic_Nmoment_Geom_ErrorBand.m',...
            'Deformation_Statistic_Nmoment_Geom_load_ErrorBand.m',...
            'Deformation_Statistic_Nmoment_Geom_Vol20_ErrorBand.m',...
            'Deformation_Statistic_Nmoment_Geom_r0_Comparison_ErrorBand.m'};
        for j=1:numel(jobs)
            fprintf('Original-style single panels: %s\n',jobs{j});
            before=findall(groot,'Type','figure');
            run_isolated(jobs{j});
            created=setdiff(findall(groot,'Type','figure'),before);
            close(created);
        end
        round2_restore_original_styles();
    case {'field_distributions','distributions'}
        run_round2_distributions('plots');
    otherwise
        error('Unknown mode. Choose single, tables, check, or field_distributions.');
end
fprintf('Round-2 mode %s completed. Results: %s\n',mode,fullfile(p.rootDir,'Results'));
end

function run_isolated(job)
% Each legacy script contains clear; isolate it from this driver's workspace.
run(job);
end

function restore_environment(oldDir,oldPath)
cd(oldDir); path(oldPath);
end
