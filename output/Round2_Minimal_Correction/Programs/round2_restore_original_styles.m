function round2_restore_original_styles()
% Clone saved original FIGs; transfer only revised numerical plot data.
% Do not reset Color: RGB getters can discard original line opacity.
p=round2_paths(); sourcePic=p.styleDir;
files=dir(fullfile(p.figureDir,'*.fig'));
names=strings(0,1); status=strings(0,1); dimensions=strings(0,1);
for j=1:numel(files)
    name=files(j).name;
    if isempty(regexp(name,'^(var|ske|kur)_(vol|load)_stress[123].*errorband\.fig$','once'))
        continue
    end
    originalPath=fullfile(sourcePic,name); newPath=fullfile(p.figureDir,name);
    assert(isfile(originalPath),'Original style reference FIG missing: %s',originalPath);
    original=openfig(originalPath,'new','invisible');
    revised=openfig(newPath,'new','invisible');
    cleanup=onCleanup(@() close_pair(original,revised));
    oldAxes=findall(original,'Type','axes'); newAxes=findall(revised,'Type','axes');
    assert(numel(oldAxes)==1 && numel(newAxes)==1,'Expected one axes in %s.',name);
    a=oldAxes(1); b=newAxes(1);
    figurePosition=original.Position; axesPosition=a.Position;
    fontName=a.FontName; fontSize=a.FontSize; lineWidth=a.LineWidth;
    xLimits=a.XLim; yLimits=a.YLim;
    legends=findall(original,'Type','legend');
    for k=1:numel(legends), legends(k).Position=legends(k).Position; end
    for type={'line','scatter','patch'}
        templates=findall(a,'Type',type{1}); objects=findall(b,'Type',type{1});
        assert(numel(templates)==numel(objects),'Object count differs: %s, %s.',type{1},name);
        for k=1:numel(templates)
            assert(strcmp(templates(k).DisplayName,objects(k).DisplayName),...
                'Object ordering differs: %s, %s.',type{1},name);
            templates(k).XData=objects(k).XData;
            templates(k).YData=objects(k).YData;
            % R2023a FIG serialization drops the fourth RGBA component.
            % Reinstate the exact RGBA literal used by the original scripts.
            if strcmp(type{1},'line') && ~contains(name,'r0_comparison') && ...
                    contains(strtrim(templates(k).DisplayName),'Analytical')
                templates(k).Color=[0.0471 0.3059 0.6078 0.6];
            end
        end
    end
    drawnow;
    assert(isequal(original.Position,figurePosition) && isequal(a.Position,axesPosition),...
        'Original figure/axes geometry changed in %s.',name);
    assert(strcmp(a.FontName,fontName) && a.FontSize==fontSize && a.LineWidth==lineWidth,...
        'Original typography changed in %s.',name);
    assert(isequal(a.XLim,xLimits) && isequal(a.YLim,yLimits),'Original limits changed in %s.',name);
    if startsWith(name,'ske_'),ylim(a,[-1.5 1.5]);end
    base=fullfile(p.figureDir,erase(name,'.fig'));
    exportgraphics(original,[base '.png'],'Resolution',600);
    savefig(original,newPath);
    names(end+1,1)=string(name);
    status(end+1,1)="Original FIG cloned; numerical data only replaced";
    originalPng=fullfile(sourcePic,[erase(name,'.fig') '.png']);
    if isfile(originalPng)
        oldInfo=imfinfo(originalPng); newInfo=imfinfo([base '.png']);
        dimensions(end+1,1)=sprintf('Original %dx%d; revised %dx%d',...
            oldInfo.Width,oldInfo.Height,newInfo.Width,newInfo.Height);
    else
        dimensions(end+1,1)="No original PNG; original FIG template retained";
    end
    fprintf('Cloned original FIG style: %s\n',name);
    clear cleanup
end
assert(numel(names)==36,'Expected 36 original-style panels; got %d.',numel(names));
report=table(names,status,dimensions,'VariableNames',{'Figure','StyleCheck','RasterDimensions'});
writetable(report,fullfile(p.logDir,'Original_Style_Verification.csv'));
end

function close_pair(a,b)
if isgraphics(a), close(a); end
if isgraphics(b), close(b); end
end
