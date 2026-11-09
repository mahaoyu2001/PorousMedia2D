clear;
r2=round2_paths(); % Isolated round-2 outputs; original FEM inputs are read-only.

% Error-band version of ExtVal_Quantile.m.
% The original script is not modified. The shaded band represents the FEM
% mean +/- one sample standard deviation across five random RVE realizations.

suffixList={'_a','_b','_c','_d','_e'};
modelNameList={'2DVoid_Vol05','2DVoid_Vol1','2DVoid_Vol2',...
    '2DVoid_Vol3','2DVoid_Vol5','2DVoid_Vol7_5','2DVoid_Vol10'};
vol=[0.5 1 2 3 5 7.5 10]/100;

Sigma=[0;1;0];
[kn, corr, corr2] = cal_stress_coef(Sigma(1),Sigma(2));
kcova=[corr(1,2),corr(1,3),corr(2,3)];

mu_RVE=zeros(length(suffixList),length(modelNameList));
var_RVE=zeros(length(suffixList),length(modelNameList));
ske_RVE=zeros(length(suffixList),length(modelNameList));
kur_RVE=zeros(length(suffixList),length(modelNameList));
moment2=zeros(length(suffixList),length(modelNameList));
moment3=zeros(length(suffixList),length(modelNameList));
moment4=zeros(length(suffixList),length(modelNameList));
q_fem=zeros(length(suffixList),length(modelNameList));

acr=0.01;
components=2;
fieldName='S';
step=2;
bandColor=[0.55 0.55 0.55];
bandAlpha=0.20;

for ia=1:length(components)
    component=components(ia);

    %% Calculate FEM upper quantiles for every RVE realization
    for nf=1:length(modelNameList)
        for is=1:length(suffixList)
            suffix=suffixList{is};
            fileName=[modelNameList{nf} suffix '_' fieldName num2str(step) '.txt'];
            Data=importdata(fullfile(r2.fieldDir,fileName));

            Evol=Data(:,end);
            weights=Evol/sum(Evol);
            phaseValue=Data(:,component);
            mu=sum(phaseValue.*weights);
            var=sum((phaseValue-mu).^2.*weights);

            moment2(is,nf)=var;
            moment3(is,nf)=sum((phaseValue-mu).^3.*weights);
            moment4(is,nf)=sum((phaseValue-mu).^4.*weights);
            mu_RVE(is,nf)=mu;
            var_RVE(is,nf)=var;
            ske_RVE(is,nf)=moment3(is,nf)./moment2(is,nf).^1.5;
            kur_RVE(is,nf)=moment4(is,nf)./moment2(is,nf).^2;
            q_fem(is,nf)=weighted_upper_quantile(phaseValue,weights,acr);
        end
    end

    [curLineWidth, axLineWidth, markerSize, ...
        axTickFontSize, axLabelFontSize, lgdFontSize, ...
        axTickFontName, axLabelFontName, lgdFontName, colorList] = plotParams();

    %% Upper 1% quantile of stress
    f=figure;
    f.Units='centimeters';
    i=component;
    Sigma_i=Sigma(i);
    picName=[r2.figureDir filesep 'Q_EaxtVal_vol_stress' num2str(i) '_errorband'];
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.65];

    extv=mean(q_fem,1);
    extvSD=std(q_fem,0,1);
    hBand=draw_error_band(vol,extv-extvSD,extv+extvSD,...
        bandColor,bandAlpha);
    hold on

    x=0.0001:0.0025:0.125;
    [mu_a, var_a, ske_a, kur_a] = cal_stress_statistics(Sigma,x);
    q_gnd=zeros(1,length(x));
    q_n=zeros(1,length(x));
    for iv=1:length(x)
        mu=mu_a(component,iv);
        var=var_a(component,iv);
        kur=kur_a(component,iv);
        q_gnd(iv)=gnd_upper_quantile(mu,var,kur,acr);
        q_n(iv)=gnd_upper_quantile(mu,var,3,acr);
    end

    hFEM=scatter(vol,extv,60,'o',LineWidth=1.25,...
        MarkerEdgeColor=[1 1 1]*0.25,DisplayName='     FEM');
    hGND=plot(x,q_gnd,LineWidth=curLineWidth,LineStyle="--",...
        Color=[0.0471 0.3059 0.6078],DisplayName='GenNormal');
    hNormal=plot(x,q_n,LineWidth=curLineWidth,LineStyle=":",...
        Color=[0.7804 0.1333 0.1569],DisplayName='   Normal');

    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,...
        FontName=axTickFontName,box='on',Layer='top');
    xlim([0 0.1125]);
    ylim([1 2.55]);
    ax.XTick=0:0.025:0.125;
    % xlabel('Pore volume fraction',FontSize=axLabelFontSize,FontName=axLabelFontName)
    % ylabel('upper 1% quantile of stress',FontSize=axLabelFontSize,FontName=axLabelFontName)
    legend([hFEM,hGND,hNormal,hBand],...
        {'     FEM','GenNormal','  Gaussian','  SD band'},...
        FontSize=lgdFontSize,FontName=lgdFontName,...
        Location='southeast',Box='off');

    exportgraphics(gcf,[picName '.png'],Resolution=600);
    saveas(gcf,[picName '.fig']);
end

function h=draw_error_band(x,lowerBound,upperBound,bandColor,bandAlpha)
x=x(:).';
lowerBound=lowerBound(:).';
upperBound=upperBound(:).';
h=fill([x fliplr(x)],[upperBound fliplr(lowerBound)],bandColor,...
    FaceAlpha=bandAlpha,EdgeColor='none');
end
