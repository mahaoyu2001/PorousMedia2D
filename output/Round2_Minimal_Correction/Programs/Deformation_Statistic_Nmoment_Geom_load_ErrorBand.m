clear;
r2=round2_paths(); % Isolated round-2 outputs; original FEM inputs are read-only.

% Error-band version of Deformation_Statistic_Nmoment_Geom_load.m.
% The original script is not modified. A single run processes all three
% stress components and exports the statistics to one Excel workbook.

suffixList={'_a','_b','_c','_d','_e'};
modelNameList={'2DVoid_Vol05-B-1','2DVoid_Vol05-B-075',...
    '2DVoid_Vol05-B-05','2DVoid_Vol05-B-025','2DVoid_Vol05',...
    '2DVoid_Vol05-B025','2DVoid_Vol05-B05','2DVoid_Vol05-B075',...
    '2DVoid_Vol05-B'};
RL=[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];

Sigma2=1;
vol=0.005;
components=[1 2 3];
componentNames=["S11","S22","S12"];
fieldName='S';
step=2;

% Centralized vertical-axis limits. Rows correspond to S11, S22 and S12.
% Modify these matrices only if the panel ranges need to be adjusted.
varYLim=[-0.0010 0.0150;...
          0.0015 0.0135;...
          0.0010 0.0030];
skeYLim=repmat([-1.5 1.5],3,1);
kurYLim=[25 125;...
         40 125;...
         25 110];

bandColor=[0.86 0.36 0.43];
bandAlpha=0.18;

nComponent=length(components);
nModel=length(modelNameList);
nRealization=length(suffixList);

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

useCachedFEM=true; % Set false to recalculate from the original FieldValue files.
if useCachedFEM
    [meanMean,meanSD,varMean,varSD,skeMean,skeSD,kurMean,kurSD]=...
        round2_cached_arrays('load',RL,componentNames);
else
%% Import each field file once and calculate all three stress components
for nf=1:nModel
    for is=1:nRealization
        suffix=suffixList{is};
        fileName=[modelNameList{nf} suffix '_' fieldName num2str(step) '.txt'];
        Data=importdata(fullfile(r2.fieldDir,fileName));

        Evol=Data(:,end);
        weights=Evol/sum(Evol);
        for ia=1:nComponent
            component=components(ia);
            phaseValue=Data(:,component);
            mu=sum(phaseValue.*weights);

            mu_RVE(is,nf,ia)=mu;
            moment2(is,nf,ia)=sum((phaseValue-mu).^2.*weights);
            moment3(is,nf,ia)=sum((phaseValue-mu).^3.*weights);
            moment4(is,nf,ia)=sum((phaseValue-mu).^4.*weights);
            ske_RVE(is,nf,ia)=moment3(is,nf,ia)./moment2(is,nf,ia).^1.5;
            kur_RVE(is,nf,ia)=moment4(is,nf,ia)./moment2(is,nf,ia).^2;
        end
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

end % useCachedFEM

[curLineWidth, axLineWidth, markerSize, ...
    axTickFontSize, axLabelFontSize, lgdFontSize, ...
    axTickFontName, axLabelFontName, lgdFontName, colorList] = plotParams();
%%
for ia=1:nComponent
    component=components(ia);

    % Analytical curves
    x=-1.2:0.05:1.2;
    y_fit_var=zeros(1,length(x));
    y_fit_ske=zeros(1,length(x));
    y_fit_kur=zeros(1,length(x));
    for ir=1:length(x)
        [~, var, ske, kur] = cal_stress_statistics(...
            [x(ir)*Sigma2;Sigma2;0],vol);
        y_fit_var(ir)=var(component);
        y_fit_ske(ir)=ske(component);
        y_fit_kur(ir)=kur(component);
    end

    % Moment 2: variance
    f=figure;
    f.Units='centimeters';
    i=component;
    picName=[r2.figureDir filesep 'var_load_stress' num2str(i) '_errorband'];
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.25];

    m2=varMean(ia,:);
    sd2=varSD(ia,:);
    hBand=draw_error_band(RL,max(0,m2-sd2),m2+sd2,bandColor,bandAlpha);
    hold on

    hFit=plot(x,y_fit_var,LineWidth=curLineWidth,LineStyle="--",...
        Color=[0.0471 0.3059 0.6078 0.6],DisplayName='  Analytical');
    hFEM=scatter(RL,m2,80,'+',LineWidth=1.5,...
        MarkerEdgeColor=[0.6350 0.0780 0.1840],DisplayName='FEM results');

    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,...
        FontName=axTickFontName,Layer='top');
    xlim([-1.1 1.1]);
    ylim(varYLim(ia,:));
    if ia<nComponent
        ax.XTickLabel=[];
    end
    if ia==1
        legend([hFit,hFEM,hBand],{'  Analytical','FEM results','   SD band'},...
            FontSize=lgdFontSize,FontName=lgdFontName,...
            Location='best',Box='off');
    end

    exportgraphics(gcf,[picName '.png'],Resolution=600);
    saveas(gcf,[picName '.fig']);

    % Moment 3: skewness
    f=figure;
    f.Units='centimeters';
    picName=[r2.figureDir filesep 'ske_load_stress' num2str(i) '_errorband'];
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.25];

    m3=skeMean(ia,:);
    sd3=skeSD(ia,:);
    draw_error_band(RL,m3-sd3,m3+sd3,bandColor,bandAlpha);
    hold on

    plot(x,y_fit_ske,LineWidth=curLineWidth,LineStyle="--",...
        Color=[0.0471 0.3059 0.6078 0.6],DisplayName='  Analytical');
    scatter(RL,m3,100,'+',LineWidth=1.8,...
        MarkerEdgeColor=[0.6350 0.0780 0.1840],DisplayName='FEM results');

    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,...
        FontName=axTickFontName,Layer='top');
    xlim([-1.1 1.1]);
    ylim(skeYLim(ia,:));
    ax.YTick=-10:5:10;
    if ia<nComponent
        ax.XTickLabel=[];
    end
    % xlabel('Volume fraction',FontSize=axLabelFontSize,FontName=axLabelFontName)
    % ylabel('Skewness',FontSize=axLabelFontSize,FontName=axLabelFontName)
    % if ia==1
    %     legend(FontSize=lgdFontSize,FontName=lgdFontName,...
    %         Location='best',Box='off');
    % end

    exportgraphics(gcf,[picName '.png'],Resolution=600);
    saveas(gcf,[picName '.fig']);

    % Moment 4: kurtosis
    f=figure;
    f.Units='centimeters';
    picName=[r2.figureDir filesep 'kur_load_stress' num2str(i) '_errorband'];
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.25];

    m4=kurMean(ia,:);
    sd4=kurSD(ia,:);
    draw_error_band(RL,max(0,m4-sd4),m4+sd4,bandColor,bandAlpha);
    hold on

    plot(x,y_fit_kur,LineWidth=curLineWidth,LineStyle="--",...
        Color=[0.0471 0.3059 0.6078 0.6],DisplayName='  Analytical');
    scatter(RL,m4,80,'+',LineWidth=1.5,...
        MarkerEdgeColor=[0.6350 0.0780 0.1840],DisplayName='FEM results');

    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,...
        FontName=axTickFontName,Layer='top');
    xlim([-1.1 1.1]);
    ylim(kurYLim(ia,:));
    if ia<nComponent
        ax.XTickLabel=[];
    end

    exportgraphics(gcf,[picName '.png'],Resolution=600);
    saveas(gcf,[picName '.fig']);
end

%% Export mean statistics and their realization-to-realization SD to Excel
LoadRatio=repmat(RL(:),nComponent,1);
StressComponent=strings(nComponent*nModel,1);
for ia=1:nComponent
    rowIndex=(ia-1)*nModel+(1:nModel);
    StressComponent(rowIndex)=componentNames(ia);
end

Mean=reshape(meanMean.',[],1);
Mean_SD=reshape(meanSD.',[],1);
Variance=reshape(varMean.',[],1);
Variance_SD=reshape(varSD.',[],1);
Skewness=reshape(skeMean.',[],1);
Skewness_SD=reshape(skeSD.',[],1);
Kurtosis=reshape(kurMean.',[],1);
Kurtosis_SD=reshape(kurSD.',[],1);

statisticsTable=table(LoadRatio,StressComponent,Mean,Mean_SD,...
    Variance,Variance_SD,Skewness,Skewness_SD,Kurtosis,Kurtosis_SD);
excelFile=fullfile(r2.tableDir,'Deformation_Statistics_Load_ErrorBand.xlsx');
if isfile(excelFile)
    delete(excelFile);
end
writetable(statisticsTable,excelFile,Sheet='Statistics');

function h=draw_error_band(x,lowerBound,upperBound,bandColor,bandAlpha)
x=x(:).';
lowerBound=lowerBound(:).';
upperBound=upperBound(:).';

% Smooth the discrete uncertainty envelope with shape-preserving cubic
% interpolation, avoiding oscillations typical of unrestricted splines.
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

h=fill([xSmooth fliplr(xSmooth)],[bandUpper fliplr(bandLower)],bandColor,...
    FaceAlpha=bandAlpha,EdgeColor='none',HandleVisibility='off');
end
