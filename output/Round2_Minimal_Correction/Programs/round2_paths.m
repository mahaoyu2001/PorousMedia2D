function p = round2_paths()
% Round-2 outputs are isolated; all FEM inputs remain in the original folder.
p.programDir = fileparts(mfilename('fullpath'));
p.rootDir = fileparts(p.programDir);
p.sourceOutputDir = getenv('POROUSMEDIA_OUTPUT_DIR');
if isempty(p.sourceOutputDir),p.sourceOutputDir=fileparts(p.rootDir);end
p.styleDir=fullfile(fileparts(p.rootDir),'Style_Templates');
p.fieldDir = fullfile(p.sourceOutputDir,'FieldValue');
p.figureDir = fullfile(p.rootDir,'Results','Figures');
p.tableDir = fullfile(p.rootDir,'Results','Tables');
p.logDir = fullfile(p.rootDir,'Results','Logs');
for target = {p.figureDir,p.tableDir,p.logDir}
    if ~isfolder(target{1}), mkdir(target{1}); end
end
% Raw fields are optional when reproducing bundled metrics and tables.
end
