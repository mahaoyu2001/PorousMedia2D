function p=tail_paths()
p.programDir=fileparts(mfilename('fullpath'));
p.rootDir=fileparts(p.programDir);
p.sourceDir=getenv('POROUSMEDIA_OUTPUT_DIR');
if isempty(p.sourceDir),p.sourceDir=fileparts(p.rootDir);end
p.fieldDir=fullfile(p.sourceDir,'FieldValue');
p.legacyCacheDir=fullfile(p.sourceDir,'Round2_Minimal_Correction','Results','Cache','Field_Distributions');
p.tableDir=fullfile(p.rootDir,'Results','Tables');
p.figureDir=fullfile(p.rootDir,'Results','Figures');
p.cacheDir=fullfile(p.rootDir,'Results','Cache');
p.logDir=fullfile(p.rootDir,'Results','Logs');
for s={'tableDir','figureDir','cacheDir','logDir'}
 if ~isfolder(p.(s{1})),mkdir(p.(s{1}));end
end
p.f=[.005 .01 .02 .03 .05 .075 .1];
p.models={'2DVoid_Vol05','2DVoid_Vol1','2DVoid_Vol2','2DVoid_Vol3','2DVoid_Vol5','2DVoid_Vol7_5','2DVoid_Vol10'};
p.realizations={'a','b','c','d','e'};
end
