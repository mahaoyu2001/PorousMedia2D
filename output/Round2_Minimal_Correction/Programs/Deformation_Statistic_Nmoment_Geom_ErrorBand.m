clear;
r2=round2_paths(); % Isolated round-2 outputs; original FEM inputs are read-only.

% Error-band version of Deformation_Statistic_Nmoment_Geom.m.
% The original script is not modified. A single run processes all three
% stress components and exports the statistics to one Excel workbook.

suffixList={'_a','_b','_c','_d','_e'};
modelNameList={'2DVoid_Vol05','2DVoid_Vol1','2DVoid_Vol2',...
    '2DVoid_Vol3','2DVoid_Vol5','2DVoid_Vol7_5','2DVoid_Vol10'};
vol=[0.5 1 2 3 5 7.5 10]/100;

Sigma=[0,1,0];
[kn, corr, corr2] = cal_stress_coef(Sigma(1),Sigma(2));
kcova=[corr(1,2),corr(1,3),corr(2,3)];

components=[1 2 3];
componentNames=["S11","S22","S12"];
fieldName='S';
step=2;

% Centralized vertical-axis limits. Rows correspond to S11, S22 and S12.
% Modify these matrices only if the panel ranges need to be adjusted.
varYLim=[0 0.0425;...
         0 0.20;...
         0 0.04];
skeYLim=repmat([-1.5 1.5],3,1);
kurYLim=[0 100;...
         0 115;...
         0 90];

% Pale warm red band associated with the red FEM symbols.
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
        round2_cached_arrays('geom',vol,componentNames);
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

% Mean curves and sample standard deviations across the five RVEs.
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

    % Moment 2: variance
    f=figure;
    f.Units='centimeters';
    i=component;
    Sigma_i=Sigma(i);
    picName=[r2.figureDir filesep 'var_vol_stress' num2str(i) '_errorband'];
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.25];

    m2=varMean(ia,:);
    sd2=varSD(ia,:);
    hBand=draw_error_band(vol,max(0,m2-sd2),m2+sd2,bandColor,bandAlpha);
    hold on

    k2=kn(i,2);
    x=0:0.001:0.125;
    [~,r2Var,~,~]=cal_stress_statistics(Sigma,x);
    y_fit=r2Var(component,:);
    hFit=plot(x,y_fit,LineWidth=curLineWidth,LineStyle="--",...
        Color=[0.0471 0.3059 0.6078 0.6],DisplayName='  Analytical');
    hFEM=scatter(vol,m2,80,'+',LineWidth=1.6,...
        MarkerEdgeColor=[0.6350 0.0780 0.1840],DisplayName='FEM results');

    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,...
        FontName=axTickFontName,Layer='top');
    xlim([0 0.1125]);
    ylim(varYLim(ia,:));
    ax.XTick=0:0.025:0.125;
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
    picName=[r2.figureDir filesep 'ske_vol_stress' num2str(i) '_errorband'];
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.25];

    m3=skeMean(ia,:);
    sd3=skeSD(ia,:);
    hBand=draw_error_band(vol,m3-sd3,m3+sd3,bandColor,bandAlpha);
    hold on

    k3=kn(i,3);
    x=0:0.001:0.125;
    [~,~,r2Skew,~]=cal_stress_statistics(Sigma,x);
    y_fit=r2Skew(component,:);
    hFit=plot(x,y_fit,LineWidth=curLineWidth,LineStyle="--",...
        Color=[0.0471 0.3059 0.6078 0.6],DisplayName='  Analytical');
    hFEM=scatter(vol,m3,100,'+',LineWidth=1.8,...
        MarkerEdgeColor=[0.6350 0.0780 0.1840],DisplayName='FEM results');

    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,...
        FontName=axTickFontName,Layer='top');
    xlim([0 0.1125]);
    ylim(skeYLim(ia,:));
    ax.XTick=0:0.025:0.125;
    ax.YTick=-10:5:10;
    if ia<nComponent
        ax.XTickLabel=[];
    end
    % xlabel('Volume fraction',FontSize=axLabelFontSize,FontName=axLabelFontName)
    % ylabel('Skewness',FontSize=axLabelFontSize,FontName=axLabelFontName)
    if ia==1
        legend([hFit,hFEM,hBand],{'  Analytical','FEM results','   SD band'},FontSize=lgdFontSize,FontName=lgdFontName,...
            Location='best',Box='off');
    end

    exportgraphics(gcf,[picName '.png'],Resolution=600);
    saveas(gcf,[picName '.fig']);

    % Moment 4: kurtosis
    f=figure;
    f.Units='centimeters';
    picName=[r2.figureDir filesep 'kur_vol_stress' num2str(i) '_errorband'];
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
    xlim([0 0.1125]);
    ylim(kurYLim(ia,:));
    ax.XTick=0:0.025:0.125;
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

    exportgraphics(gcf,[picName '.png'],Resolution=600);
    saveas(gcf,[picName '.fig']);
end

%% Export mean statistics and their realization-to-realization SD to Excel
Porosity=repmat(vol(:),nComponent,1);
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

statisticsTable=table(Porosity,StressComponent,Mean,Mean_SD,...
    Variance,Variance_SD,Skewness,Skewness_SD,Kurtosis,Kurtosis_SD);
excelFile=fullfile(r2.tableDir,'Deformation_Statistics_Geom_ErrorBand.xlsx');
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
