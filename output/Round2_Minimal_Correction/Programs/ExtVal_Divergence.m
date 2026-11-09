clear;
r2=round2_paths(); % Isolated round-2 outputs; original FEM inputs are read-only.
baseDir=fileparts(mfilename('fullpath'));
if isempty(baseDir)
    baseDir=pwd;
end
fieldDir=r2.fieldDir;
picDir=r2.figureDir;
suffixList={'_a','_b','_c','_d','_e'};
modelNameList={'2DVoid_Vol05','2DVoid_Vol1','2DVoid_Vol2',...
    '2DVoid_Vol3','2DVoid_Vol5','2DVoid_Vol7_5','2DVoid_Vol10'};
vol=[0.5 1 2 3 5 7.5 10]/100;

Sigma=[0;1;0];
components=[1 2 3];
fieldName='S';
step=2;
tailFactor=6;
minBinCount=250;
epsProb=1e-12;
klDirection='FEM||Model';
lgdItemTokenSize=[20 8];
stressLabelList={'\it\sigma_{\rm{11}}\rm','\it\sigma_{\rm{22}}\rm','\it\sigma_{\rm{12}}\rm'};

js_gnd=zeros(length(suffixList),length(modelNameList));
js_n=zeros(length(suffixList),length(modelNameList));
kl_gnd=zeros(length(suffixList),length(modelNameList));
kl_n=zeros(length(suffixList),length(modelNameList));

for ia=1:length(components)
    component=components(ia);
    [mu_a, var_a, ske_a, kur_a] = cal_stress_statistics(Sigma, vol);

    for nf = 1:length(modelNameList)
        mu=mu_a(component,nf);
        var=var_a(component,nf);
        kur=kur_a(component,nf);
        stdModel=sqrt(var);

        for is=1:length(suffixList)
            suffix=suffixList{is};
            fileName=[modelNameList{nf} suffix '_' fieldName num2str(step) '.txt'];
            filePath=fullfile(fieldDir, fileName);
            if ~isfile(filePath)
                error('Cannot find field-value file: %s', filePath);
            end
            Data=importdata(filePath);

            Evol=Data(:,end);
            weights=Evol/sum(Evol);
            phaseValue=Data(:,component);
            mu0=sum(phaseValue.*weights);
            var0=sum((phaseValue-mu0).^2.*weights);
            std0=sqrt(var0);

            xLower=min([min(phaseValue), mu0-tailFactor*std0, mu-tailFactor*stdModel]);
            xUpper=max([max(phaseValue), mu0+tailFactor*std0, mu+tailFactor*stdModel]);
            binWidth=max(std0/20, (xUpper-xLower)/minBinCount);
            finiteEdges=xLower:binWidth:xUpper;
            if finiteEdges(end)<xUpper
                finiteEdges(end+1)=xUpper;
            end
            edges=[-Inf finiteEdges(2:end-1) Inf];

            p_fem=weighted_bin_prob(phaseValue, weights, edges);

            [~, ~, ~, cdf_gnd_edges] = gnd_parameters(edges, mu, var, kur);
            [~, ~, ~, cdf_n_edges] = gnd_parameters(edges, mu, var);
            q_gnd=diff(cdf_gnd_edges(:))';
            q_n=diff(cdf_n_edges(:))';

            [js_gnd(is,nf), kl_gnd(is,nf)] = distribution_distance(p_fem, q_gnd, epsProb);
            [js_n(is,nf), kl_n(is,nf)] = distribution_distance(p_fem, q_n, epsProb);
        end
    end

    js_gnd_mean=mean(js_gnd,1);
    js_n_mean=mean(js_n,1);
    kl_gnd_mean=mean(kl_gnd,1);
    kl_n_mean=mean(kl_n,1);
    js_gnd_std=std(js_gnd,0,1);
    js_n_std=std(js_n,0,1);
    kl_gnd_std=std(kl_gnd,0,1);
    kl_n_std=std(kl_n,0,1);

    results=table(vol(:), js_gnd_mean(:), js_n_mean(:), ...
        kl_gnd_mean(:), kl_n_mean(:), js_gnd_std(:), js_n_std(:), ...
        kl_gnd_std(:), kl_n_std(:), ...
        'VariableNames', {'VolumeFraction','JS_GenNormal','JS_Normal', ...
        'KL_GenNormal','KL_Normal','JS_GenNormal_Std','JS_Normal_Std', ...
        'KL_GenNormal_Std','KL_Normal_Std'});
    writetable(results, fullfile(r2.tableDir, ['Divergence_JS_KL_stress' num2str(component) '.csv']));
    save(fullfile(r2.tableDir, ['Divergence_JS_KL_stress' num2str(component) '.mat']), ...
        'vol','suffixList','modelNameList','component','Sigma', ...
        'js_gnd','js_n','kl_gnd','kl_n','klDirection','results');

    [curLineWidth, axLineWidth, markerSize, ...
    axTickFontSize, axLabelFontSize, lgdFontSize, ...
    axTickFontName, axLabelFontName, lgdFontName, colorList] = plotParams();

    if ~exist(picDir,'dir')
        mkdir(picDir);
    end
%%
    picName=fullfile(picDir, ['JS_Divergence_vol_stress' num2str(component)]);
    f=figure;
    f.Units='centimeters';
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.65];
    plot(vol,js_gnd_mean,'--o',LineWidth=curLineWidth,MarkerSize=5, ...
        MarkerFaceColor=[1 1 1],Color=[0.0471 0.3059 0.6078], ...
        DisplayName='GenNormal');
    hold on
    plot(vol,js_n_mean,':s',LineWidth=curLineWidth,MarkerSize=6, ...
        MarkerFaceColor=[1 1 1],Color=[0.7804 0.1333 0.1569], ...
        DisplayName='  Gaussian');
    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,FontName=axTickFontName,box='on');
    xlim([0 0.1125]);
    ax.XTick=0:0.025:0.125;
    ylim([-0.005 0.175]);
    % ax.YTick=0:0.01:0.03;
    xlabel('Pore volume fraction',FontSize=axLabelFontSize,FontName=axLabelFontName)
    ylabel(['JS divergence (' stressLabelList{component} ')'], ...
        FontSize=axLabelFontSize,FontName=axLabelFontName,Interpreter='tex')
    lgd=legend(FontSize=lgdFontSize,FontName=lgdFontName,Location='best', ...
        Box='off');
    lgd.ItemTokenSize=lgdItemTokenSize;
    exportgraphics(gcf, [picName '.png'],Resolution=600);
    saveas(gcf, [picName '.fig']);
%%
    picName=fullfile(picDir, ['KL_Divergence_vol_stress' num2str(component)]);
    f=figure;
    f.Units='centimeters';
    ax=gca;
    ax.Units='centimeters';
    ax.Position=[2,2,5.8,4.65];
    plot(vol,kl_gnd_mean,'--o',LineWidth=curLineWidth,MarkerSize=5, ...
        MarkerFaceColor=[1 1 1],Color=[0.0471 0.3059 0.6078], ...
        DisplayName='GenNormal');
    hold on
    plot(vol,kl_n_mean,':s',LineWidth=curLineWidth,MarkerSize=6, ...
        MarkerFaceColor=[1 1 1],Color=[0.7804 0.1333 0.1569], ...
        DisplayName='  Gaussian');
    set(ax,FontSize=axTickFontSize,LineWidth=axLineWidth,FontName=axTickFontName,box='on');
    xlim([0 0.1125]);
    ax.XTick=0:0.025:0.125;
    ylim([-0.02 0.85]);
    % ax.YTick=0:0.01:0.03;
    xlabel('Pore volume fraction',FontSize=axLabelFontSize,FontName=axLabelFontName)
    ylabel(['KL divergence (' stressLabelList{component} ')'], ...
        FontSize=axLabelFontSize,FontName=axLabelFontName,Interpreter='tex')
    lgd=legend(FontSize=lgdFontSize,FontName=lgdFontName,Location='best', ...
        Box='off');
    lgd.ItemTokenSize=lgdItemTokenSize;
    exportgraphics(gcf, [picName '.png'],Resolution=600);
    saveas(gcf, [picName '.fig']);
end

function p=weighted_bin_prob(data, weights, edges)
    data=data(:);
    weights=weights(:);
    weights=weights/sum(weights);
    p=zeros(1,length(edges)-1);
    for ib=1:length(p)
        if ib==length(p)
            inBin=(data>=edges(ib)) & (data<=edges(ib+1));
        else
            inBin=(data>=edges(ib)) & (data<edges(ib+1));
        end
        p(ib)=sum(weights(inBin));
    end
    p=p/sum(p);
end

function [js, kl]=distribution_distance(p, q, epsProb)
    p=clean_prob(p);
    q=clean_prob(q);
    qKL=max(q,epsProb);
    qKL=qKL/sum(qKL);
    kl=kl_divergence(p,qKL);

    m=0.5*(p+q);
    js=0.5*kl_divergence(p,m)+0.5*kl_divergence(q,m);
end

function p=clean_prob(p)
    p=p(:)';
    p(~isfinite(p))=0;
    p=max(p,0);
    p=p/sum(p);
end

function d=kl_divergence(p,q)
    idx=(p>0);
    d=sum(p(idx).*log(p(idx)./q(idx)));
end

function set_legend_line_length(lgd, lineLength, tokenHeight)
    try
        lgd.ItemTokenSize=[lineLength, tokenHeight];
    catch
        try
            set(lgd,'ItemTokenSize',[lineLength, tokenHeight]);
        catch
        end
    end
end
