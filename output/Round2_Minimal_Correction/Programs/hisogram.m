clear;
r2=round2_paths(); % Isolated round-2 outputs; original FEM inputs are read-only.
%%
fileNameList={'2DVoid_Vol1_a_S2.txt'};
[curLineWidth, axLineWidth, markerSize, ...
axTickFontSize, axLabelFontSize, lgdFontSize, ...
axTickFontName, axLabelFontName, lgdFontName, colorList] = ...
plotParams(curLineWidth=1.5, axLineWidth=0.75, ...
axTickFontSize=10.5, axLabelFontSize=13.5, lgdFontSize=12);

comp_list=[1];
% solution for n-moments

Sigma=[0;1;0]; % 远场应力载荷
vf=[1]/100;
[mu_a, var_a, ske_a, kur_a] = cal_stress_statistics(Sigma, vf);

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
        std=sqrt(var);

        xLim=[mu0-5*std0, mu0+5*std0];
        binWidth=std0/20;

        % weighted histogram & cdf
        [hist_data, cdf_data] = weighted_hist_cdf(phaseValue, weights, ...
            BinWidth=binWidth,XLim=xLim);
        % semi-analytical solution for pdf & cdf
        x=xLim(1):0.5*binWidth:xLim(2);
        [~, ~, pdf_gnd,cdf_gnd] = gnd_parameters(x, mu, var, kur);
        [~, ~, pdf_n,cdf_n] = gnd_parameters(x, mu, var);
        %% example: critical value
        scr=mu+3*std;
        [~, ~, ~,cdf_cr_gnd] = gnd_parameters(scr, mu, var, kur);
        [~, ~, ~,cdf_cr_n] = gnd_parameters(scr, mu, var);
        acr_gnd=1-cdf_cr_gnd;acr_n=1-cdf_cr_n;
        acr_fem=sum(weights(phaseValue>scr));
        %% stress_invariants
        % [sigma1, sigma2, sigma_vm, ~, ~, ~] = stress_invariants(Data(:,1:3));
        % [hist_data_s1, cdf_data_s1] = weighted_hist_cdf(sigma1, weights, ...
        % BinWidth=binWidth,XLim=xLim);
        % [hist_data_sm, cdf_data_sm] = weighted_hist_cdf(sigma_vm, weights, ...
        % BinWidth=binWidth,XLim=xLim);
        % figure;
        % plot(cdf_data.x, cdf_data.y,lineStyle='--', ...
        %     Color=[1 1 1]*0.25, LineWidth=curLineWidth, DisplayName='s11');
        % hold on;
        % plot(cdf_data_s1.x, cdf_data_s1.y,lineStyle='-',LineWidth=curLineWidth,Color=[0.0471 0.3059 0.6078], ...
        %     DisplayName=' s1');
        % hold on
        % plot(cdf_data_sm.x, cdf_data_sm.y,lineStyle='-',LineWidth=curLineWidth,Color=[0.7804 0.1333 0.1569], ...
        %     DisplayName=' sm');
        % lgd=legend(FontSize=lgdFontSize,FontName=lgdFontName,Location='best',Box='off');
        % ax=gca;
        % set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,FontName=axTickFontName);
        %% histogram/pdf
        picName=[r2.figureDir filesep 'pdf_' fileName '_stress' num2str(component)];
        f=figure;
        ax=gca;
        ax.Units='centimeters';
        ax.Position=[2,2,6.75,5.4];
        bar(hist_data.binCenters, hist_data.pdf, ...
            FaceColor=[0.5 0.5 0.5], FaceAlpha=0.4,...
            EdgeColor=[0 0 0], EdgeAlpha=0.2, LineWidth=0.1, ...
            DisplayName='    FEM');
        hold on;
        plot(x, pdf_gnd,lineStyle='-',LineWidth=curLineWidth,Color=[0.0471 0.3059 0.6078], ...
            DisplayName='GenNorm');
        hold on
        plot(x, pdf_n,lineStyle='-',LineWidth=curLineWidth,Color=[0.7804 0.1333 0.1569], ...
            DisplayName='  Normal');

        lgd=legend(FontSize=lgdFontSize,FontName=lgdFontName,Location='northeast',Box='off');
        set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,FontName=axTickFontName);
        xlim([0.61 1.41]);
        ylim([0 20]);
        % ax.XTick=[-0.2 0 0.2];
        xlabel(['Stress ' xlabelName{component}],FontSize=axLabelFontSize,FontName=axLabelFontName);
        ylabel('Probobility density',FontSize=axLabelFontSize,FontName=axLabelFontName);

        exportgraphics(gcf, [picName '.png'],Resolution=600); saveas(gcf, [picName '.fig']);
        %% cdf
        picName=[r2.figureDir filesep 'cdf_' fileName '_stress' num2str(component)];
        f=figure;
        ax=gca;
        ax.Units='centimeters';
        ax.Position=[2,2,6.75,5.4];
        plot(cdf_data.x, cdf_data.y, lineStyle='--',...
            Color=[1 1 1]*0.25, LineWidth=curLineWidth, DisplayName='    FEM');
        hold on;
        plot(x, cdf_gnd,lineStyle='-',LineWidth=curLineWidth,Color=[0.0471 0.3059 0.6078], ...
            DisplayName='GenNorm');
        hold on
        plot(x, cdf_n,lineStyle='-',LineWidth=curLineWidth,Color=[0.7804 0.1333 0.1569], ...
            DisplayName='  Normal');

        ax=gca;
        lgd=legend(FontSize=lgdFontSize,FontName=lgdFontName,Location='southeast',Box='off');
        set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,FontName=axTickFontName);
        xlim([0.73 1.29]);
        % ylim([0 80]);
        % ax.XTick=[-0.2 0 0.2];
        xlabel(['Stress ' xlabelName{component}],FontSize=axLabelFontSize,FontName=axLabelFontName)
        ylabel('Cumulative distribution',FontSize=axLabelFontSize,FontName=axLabelFontName)

        exportgraphics(gcf, [picName '.png'],Resolution=600);
        saveas(gcf, [picName '.fig']);
    end
end

%%
