function export_round2_comparisons()
% Refresh analytic predictions, using the unchanged original FEM caches.
p=round2_paths();
book=fullfile(p.tableDir,'Round2_Analytical_FEM_Comparison.xlsx');
kinds={'geom','load','vol20','point'};
sheetNames={'Geom','Load','Vol20','Point_Process'};
for ik=1:numel(kinds)
    kind=kinds{ik}; t=round2_read_statistics(kind); n=height(t);
    if strcmp(kind,'load')
        t.Porosity=repmat(0.005,n,1);
    else
        t.LoadRatio=zeros(n,1);
    end
    t.AnalyticalMean=nan(n,1); t.AnalyticalVariance=nan(n,1);
    t.AnalyticalSkewness=nan(n,1); t.AnalyticalKurtosis=nan(n,1);
    t.GND_Alpha=nan(n,1); t.GND_Beta=nan(n,1);
    for ir=1:n
        c=find(["S11","S22","S12"]==string(t.StressComponent(ir)));
        assert(isscalar(c),'Unrecognized stress component in cached row %d.',ir);
        [mu,v,s,k]=cal_stress_statistics([t.LoadRatio(ir);1;0],t.Porosity(ir));
        t.AnalyticalMean(ir)=mu(c); t.AnalyticalVariance(ir)=v(c);
        t.AnalyticalSkewness(ir)=s(c); t.AnalyticalKurtosis(ir)=k(c);
        [t.GND_Alpha(ir),t.GND_Beta(ir)]=gnd_parameters([],mu(c),v(c),k(c));
    end
    t.Variance_Error_Percent=relative_error(t.AnalyticalVariance,t.Variance);
    t.Kurtosis_Error_Percent=relative_error(t.AnalyticalKurtosis,t.Kurtosis);
    t.Skewness_AbsoluteError=abs(t.AnalyticalSkewness-t.Skewness);
    writetable(t,book,'Sheet',sheetNames{ik});
    writetable(t,fullfile(p.tableDir,['Round2_Comparison_' kind '.csv']));
    fprintf('Updated %s: %d original FEM rows retained.\n',kind,n);
end
metadata=table(["Model";"Porosity identification";"Background stress";...
    "Relative error reference";"FEM uncertainty";"Vol20 source";"Created"],...
    ["Poisson single-pore superposition; minimal statistical correction";...
    "rho*A0 approximately f; no silent replacement by -log(1-f)";...
    "Sigma/(1-f), retained mean-field approximation";...
    "abs(Analytical-FEM)/abs(FEM)*100";...
    "Original realization-to-realization sample SD";...
    "Original Vol20 workbook, Summary sheet";string(datetime('now'))],...
    'VariableNames',{'Item','Convention'});
writetable(metadata,book,'Sheet','Model_Info');
export_model_curves(p);
fprintf('Comparison workbook: %s\n',book);
end

function e=relative_error(a,b)
e=nan(size(a)); valid=isfinite(a)&isfinite(b)&b~=0;
e(valid)=100*abs(a(valid)-b(valid))./abs(b(valid));
end

function plot_comparison(t,kind,base)
fig=figure('Visible','off','Color','w','Position',[100 100 1000 900]);
layout=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','loose');
names=["S11","S22","S12"]; labels={'\sigma_{11}','\sigma_{22}','\sigma_{12}'};
isLoad=strcmp(kind,'load'); key='Porosity'; if isLoad, key='LoadRatio'; end
if isLoad
    x=linspace(min(t.LoadRatio),max(t.LoadRatio),241);
    va=zeros(3,numel(x)); ka=va;
    for j=1:numel(x), [~,va(:,j),~,ka(:,j)]=cal_stress_statistics([x(j);1;0],0.005); end
else
    x=linspace(min(t.Porosity),max(t.Porosity),401);
    [~,va,~,ka]=cal_stress_statistics([0;1;0],x);
end
methods="FEM"; if ismember('Method',t.Properties.VariableNames), methods=unique(string(t.Method),'stable'); end
colors=[0.75 0.12 0.16;0.05 0.45 0.72];
for i=1:3
    for statistic=1:2
        ax=nexttile(layout); hold(ax,'on');
        field='Variance'; ya=va(i,:); if statistic==2, field='Kurtosis'; ya=ka(i,:); end
        for im=1:numel(methods)
            mask=string(t.StressComponent)==names(i);
            if ismember('Method',t.Properties.VariableNames), mask=mask & string(t.Method)==methods(im); end
            rows=t(mask,:); [xx,idx]=sort(rows.(key));
            yy=rows.(field)(idx); sd=rows.([field '_SD'])(idx);
            xx=xx(:).'; yy=yy(:).'; sd=sd(:).';
            if ~isLoad, xx=100*xx; end
            color=colors(min(im,size(colors,1)),:);
            fill(ax,[xx fliplr(xx)],[max(0,yy-sd) fliplr(yy+sd)],color,...
                'FaceAlpha',0.14,'EdgeColor','none','HandleVisibility','off');
            plot(ax,xx,yy,'o','Color',color,'MarkerFaceColor',color,'MarkerSize',4,...
                'DisplayName',char(methods(im)));
        end
        xa=x; if ~isLoad, xa=100*x; end
        plot(ax,xa,ya,'--','Color',[0.1 0.25 0.5],'LineWidth',1.5,'DisplayName','Analytical');
        ylabel(ax,[field ' (' labels{i} ')'],'Interpreter','tex');
        if isLoad, xlabel(ax,'Biaxial stress ratio'); else, xlabel(ax,'Porosity (%)'); end
        set(ax,'FontName','Times New Roman','FontSize',11,'Box','on');
        if i==1 && statistic==1, legend(ax,'Location','best','Box','off'); end
    end
end
exportgraphics(fig,[base '.png'],'Resolution',300);
exportgraphics(fig,[base '.jpg'],'Resolution',300);
savefig(fig,[base '.fig']); close(fig);
end

function export_model_curves(p)
vol=[0.005 0.01 0.02 0.03 0.05 0.075 0.10 0.15 0.20];
curveF=linspace(min(vol),max(vol),401);
[mu,v,s,k]=cal_stress_statistics([0;1;0],curveF);
curves=table();
componentNames=["S11","S22","S12"];
for i=1:3
    part=table(repmat(componentNames(i),numel(curveF),1),...
        curveF(:),mu(i,:).',v(i,:).',s(i,:).',k(i,:).',...
        'VariableNames',{'StressComponent','Porosity','Mean','Variance','Skewness','Kurtosis'});
    curves=[curves;part]; %#ok<AGROW>
end
writetable(curves,fullfile(p.tableDir,'Round2_Analytical_Curves.csv'));
rows=cell(0,10); labels=["S11","S22","S12"];
for f=vol
    [mu,v,~,k]=cal_stress_statistics([0;1;0],f);
    for i=1:3
        thresholds=0.5;
        if i==2, thresholds=[1.5 2.0 2.5]; end
        [a,b]=gnd_parameters([],mu(i),v(i),k(i));
        for threshold=thresholds
            rows(end+1,:)={labels(i),f,threshold,a,b,...
                gnd_survival(threshold,mu(i),v(i),k(i)),...
                gnd_survival(threshold,mu(i),v(i),3),...
                gnd_upper_quantile(mu(i),v(i),k(i),0.01),...
                gnd_upper_quantile(mu(i),v(i),3,0.01),"Analytical only; FEM tail probabilities not inferred from moments"}; %#ok<AGROW>
        end
    end
end
tails=cell2table(rows,'VariableNames',{'StressComponent','Porosity','Threshold',...
    'GND_Alpha','GND_Beta','P_GND','P_Gaussian','Q99_GND','Q99_Gaussian','DataScope'});
writetable(tails,fullfile(p.tableDir,'Round2_Model_Tail_Predictions.csv'));
end
