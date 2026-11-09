function run_tail_validation(stage)
% Main entry point. All outputs are isolated and legacy FEM files are read-only.
if nargin<1,stage='plots';end
oldPath=path;cleanupPath=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fileparts(mfilename('fullpath')),'-begin');
assert(any(strcmp(stage,{'plots','pilot','fem','base','convergence','export','all'})),'Unknown stage.');
p=tail_paths();diary(fullfile(p.logDir,['Run_' stage '.log']));cleanup=onCleanup(@()diary('off'));
fprintf('R5Q4 tail validation; %s\n',datestr(now,31));
cfg=struct('Name','Base','N',2^19,'L',16,'Nradial',1024,'Ntheta',2048,'Ucut',1e-5,'SeriesOrder',40);
if strcmp(stage,'pilot')
 out=campbell_fft(2,.005,cfg);disp(out);return
end
if any(strcmp(stage,{'fem','all'})),tail_fem_records(false);end
if any(strcmp(stage,{'base','all'}))
 for f=p.f,for component=1:3,campbell_fft(component,f,cfg);end,end
end
if any(strcmp(stage,{'convergence','all'}))
 configs={cfg};
 c=cfg;c.Name='SpatialFine';c.Nradial=2048;c.Ntheta=4096;configs{end+1}=c;
 c=cfg;c.Name='GridFine';c.N=2^20;c.Nradial=2048;c.Ntheta=4096;configs{end+1}=c;
 c=cfg;c.Name='DomainWide';c.N=2^21;c.L=32;c.Nradial=2048;c.Ntheta=4096;configs{end+1}=c;
 c=cfg;c.Name='WeakFine';c.N=2^20;c.Nradial=2048;c.Ntheta=4096;c.Ucut=5e-6;configs{end+1}=c;
 c=cfg;c.Name='SpatialVeryFine';c.N=2^20;c.Nradial=4096;c.Ntheta=8192;configs{end+1}=c;
 for j=1:numel(configs)
  for f=p.f,for component=1:3,campbell_fft(component,f,configs{j});end,end
 end
 save(fullfile(p.cacheDir,'Convergence_Configurations.mat'),'configs');
end
if any(strcmp(stage,{'export','all'})),export_tail_validation();end
if any(strcmp(stage,{'plots','export','all'})),plot_tail_validation();end
fprintf('Stage %s complete.\n',stage);
end
