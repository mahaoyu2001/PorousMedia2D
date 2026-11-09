function plot_round2_relative_errors()
% Existing full-precision comparisons only; no FEM solve or field-file read.
p=round2_paths();t=readtable(fullfile(p.tableDir,'Round2_Comparison_vol20.csv'),'TextType','string');
out=fullfile(p.figureDir,'Relative_Errors');if ~isfolder(out),mkdir(out);end
[lw,aw,~,tickSize,labelSize,legendSize,tickFont,labelFont,legendFont]=plotParams();
components={'S11','S22','S12'};labels={'\sigma_{11}','\sigma_{22}','\sigma_{12}'};
colors=[0.0471 0.3059 0.6078;0.7804 0.1333 0.1569;0.9765 0.5608 0.2039];
styles={'-o','--s',':^'};
for metric={'Variance','Kurtosis'}
 name=metric{1};fig=figure('Visible','off','Units','centimeters');
 ax=axes(fig,'Units','centimeters','Position',[2 2 5.8 4.65]);hold(ax,'on');
 for i=1:3
  rows=sortrows(t(t.StressComponent==components{i},:),'Porosity');assert(height(rows)==9);
  value=100*abs(rows.(['Analytical' name])-rows.(name))./abs(rows.(name));
  assert(max(abs(value-rows.([name '_Error_Percent'])))<1e-9);
  plot(ax,rows.Porosity,value,styles{i},'LineWidth',lw,'MarkerSize',5,...
   'MarkerFaceColor','white','Color',colors(i,:),'DisplayName',labels{i});
 end
 set(ax,'FontName',tickFont,'FontSize',tickSize,'LineWidth',aw,'Box','on',...
  'XLim',[0 0.21],'XTick',0:0.05:0.2,'YLim',[0 40],'YTick',0:10:40);
 if strcmp(name,'Kurtosis'),set(ax,'YLim',[0 60],'YTick',0:10:60);end
 xlabel(ax,'Porosity','FontName',labelFont,'FontSize',labelSize);
 ylabel(ax,'Relative error (%)','FontName',labelFont,'FontSize',labelSize);
 legend(ax,'Location','northwest','Box','off','FontName',legendFont,'FontSize',legendSize);
 base=fullfile(out,['Relative_Error_' name '_vol']);
 exportgraphics(fig,[base '.png'],'Resolution',600);savefig(fig,[base '.fig']);close(fig);
 fprintf('Exported %s (PNG + FIG)\n',base);
end
% Reproducible standard errors for both porosity and loading comparisons.
for kind={'vol20','load'}
 source=fullfile(p.tableDir,['Round2_Comparison_' kind{1} '.csv']);
 data=readtable(source,'TextType','string');
 if strcmp(kind{1},'vol20'),key='Porosity';suffix='Porosity';else,key='LoadRatio';suffix='Loading';end
 cols={key,'StressComponent','Mean','Mean_SD','Variance','Variance_SD',...
  'Skewness','Skewness_SD','Kurtosis','Kurtosis_SD'};u=data(:,cols);
 u.NumRealizations=repmat(5,height(u),1);
 for stat={'Mean','Variance','Skewness','Kurtosis'}
  u.([stat{1} '_SE'])=u.([stat{1} '_SD'])/sqrt(5);
 end
 writetable(u,fullfile(p.tableDir,['Round2_Mean_Uncertainty_' suffix '.csv']));
end
end
