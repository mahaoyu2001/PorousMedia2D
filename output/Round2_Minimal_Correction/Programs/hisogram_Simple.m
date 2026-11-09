clear;
r2=round2_paths(); % Isolated round-2 outputs; original FEM inputs are read-only.
%%
fileNameList={'2DVoid_Vol1_a_S2.txt'};
[curLineWidth, axLineWidth, markerSize, ...
axTickFontSize, axLabelFontSize, lgdFontSize, ...
axTickFontName, axLabelFontName, lgdFontName, colorList] = ...
plotParams(curLineWidth=1.65, axLineWidth=0.85, ...
axTickFontSize=10.5, axLabelFontSize=13.5, lgdFontSize=12.5);

comp_list=[2];
% solution for n-moments

Sigma=[0;1;0]; % 远场应力载荷
vf=[1]/100;
[mu_a, var_a, ske_a, kur_a] = cal_stress_statistics(Sigma, vf);

XLimL=[0.5 1.5];
XLimE=[1.25 1.45];

%%
xlabelName={'\it{\sigma}_{\rm{11}}';'\it{\sigma}_{\rm{22}}';'\it{\sigma}_{\rm{12}}'};
for nf = 1:length(fileNameList)
    fileName=fileNameList{nf};
    Data=load(fullfile(r2.fieldDir,fileName));
    for component=comp_list
        phaseValue=Data(:,component);
        weights=Data(:,end)/(sum(Data(:,end)));
        mu0=sum(phaseValue.*weights);
        moment2=sum((phaseValue-mu0).^2.*weights);
        moment3=sum((phaseValue-mu0).^3.*weights);
        moment4=sum((phaseValue-mu0).^4.*weights);
        var0=moment2;
        ske0=moment3./moment2.^1.5;
        kur0=moment4./moment2.^2;
        std0=sqrt(var0);

        mu=mu_a(component);
        var=var_a(component);
        kur=kur_a(component);

        xLim=[mu0-5*std0, mu0+5*std0];
        binWidth=std0/5;

        % weighted histogram & cdf
        [hist_data, cdf_data] = weighted_hist_cdf(phaseValue, weights, ...
            BinWidth=binWidth,XLim=xLim);
        % semi-analytical solution for pdf & cdf
        x=xLim(1):0.5*binWidth:xLim(2);
        [~, ~, pdf_gnd,cdf_gnd] = gnd_parameters(x, mu, var, kur);
        [~, ~, pdf_n,cdf_n] = gnd_parameters(x, mu0, var0);
        %% histogram/pdf
        picName=[r2.figureDir filesep 'pdf_' fileName '_stress' num2str(component) '_simple'];
        f=figure;
        ax=gca;
        ax.Units='centimeters';
        ax.Position=[2,2,10.8,4.2];
        bar(hist_data.binCenters, hist_data.pdf, ...
            FaceColor=[0.0471 0.3059 0.6078], FaceAlpha=0.4,...
            EdgeColor=[0 0 0], EdgeAlpha=0.5, LineWidth=0.5, ...
            DisplayName='      FEM');
        hold on;
        plot(x, pdf_n,lineStyle=':',LineWidth=1.8,Color=[0 0 0], ...
            DisplayName='Gaussian fit');

        lgd=legend(FontSize=lgdFontSize,FontName=lgdFontName,Location='northwest',Box='off');
        set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,FontName=axTickFontName);
        xlim(XLimL);
        ylim([0 15]);
        ax.XTick=0:0.2:2;
        % xlabel(['Stress ' xlabelName{component}],FontSize=axLabelFontSize,FontName=axLabelFontName);
        % ylabel('Probobility density',FontSize=axLabelFontSize,FontName=axLabelFontName);
        exportgraphics(gcf, [picName '.png'],Resolution=600); saveas(gcf, [picName '.fig']);

        %%
        picName=[r2.figureDir filesep 'pdf_' fileName '_stress' num2str(component) '_simpleE'];
        f=figure;
        ax=gca;
        ax.Units='centimeters';
        ax.Position=[2,2,3,2];
        bar(hist_data.binCenters, hist_data.pdf, ...
            FaceColor=[0.0471 0.3059 0.6078], FaceAlpha=0.4,...
            EdgeColor=[0 0 0], EdgeAlpha=0.5, LineWidth=0.5, ...
            DisplayName='      FEM');
        hold on;
        plot(x, pdf_n,lineStyle=':',LineWidth=1.5,Color=[0 0 0], ...
            DisplayName='Gaussian fit');

        % lgd=legend(FontSize=lgdFontSize,FontName=lgdFontName,Location='northwest',Box='off');
        set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,FontName=axTickFontName);
        xlim(XLimE);
        ylim([0 0.1]);
        ax.XTick=0:0.1:2;
        % xlabel(['Stress ' xlabelName{component}],FontSize=axLabelFontSize,FontName=axLabelFontName);
        % ylabel('Probobility density',FontSize=axLabelFontSize,FontName=axLabelFontName);
        exportgraphics(gcf, [picName '.png'],Resolution=600); saveas(gcf, [picName '.fig']);
    end
end

%%
