function plot_tail_validation()
p=tail_paths();T=readtable(fullfile(p.tableDir,'R5Q4_Exceedance_All.csv'));
Q=readtable(fullfile(p.tableDir,'R5Q4_Q99_All.csv'));
for c=1:3
 thresholds=unique(T.Threshold(T.Component==c));
 for threshold=thresholds'
  r=T(T.Component==c & T.Threshold==threshold,:);
  tag=sprintf('S%d_Threshold_%g',c,threshold);tag=strrep(tag,'.','p');
  fig=figure('Visible','off','Color','w','Units','centimeters','Position',[5 5 10.3 8.6]);
  ax=axes(fig,'Units','centimeters','Position',[1.7 1.4 7.9 6.7]);hold(ax,'on');
  lo=max(r.FEM_Mean-r.FEM_SD,1e-10);hi=r.FEM_Mean+r.FEM_SD;
  hband=fill(ax,[r.Porosity;flipud(r.Porosity)],[hi;flipud(lo)],[.55 .55 .55],'FaceAlpha',.2,'EdgeColor','none');
  h1=plot(ax,r.Porosity,r.FEM_Mean,'o','Color',[.25 .25 .25],'MarkerSize',5,'LineWidth',1.25);
  h2=plot(ax,r.Porosity,r.GND_Analytical,'--','Color',[.0471 .3059 .6078],'LineWidth',1.5);
  h3=plot(ax,r.Porosity,r.GND_FEMmoments_Mean,'-.','Color',[.78 .133 .157],'LineWidth',1.5);
  h4=plot(ax,r.Porosity,r.Campbell_Full,'-','Color',[.12 .56 .26],'LineWidth',1.5);
  set(ax,'YScale','log','FontName','Times New Roman','FontSize',10,'LineWidth',.75,'Box','on','XTick',0:.025:.1);
  ylim(ax,[.8*min([lo;r.GND_Analytical;r.GND_FEMmoments_Mean;r.Campbell_Full]),1.25*max([hi;r.GND_Analytical;r.GND_FEMmoments_Mean;r.Campbell_Full])]);
  xlim(ax,[0 .1125]);xlabel(ax,'Porosity','FontSize',11);ylabel(ax,'Exceedance probability','FontSize',11);
  legend(ax,[h1 h2 h3 h4 hband],{'FEM','GND: analytical moments','GND: FEM moments','Full Campbell','SD band'},'Location','best','Box','off','FontName','Times New Roman','FontSize',8);
  exportgraphics(fig,fullfile(p.figureDir,['Exceedance_' tag '.png']),'Resolution',600);savefig(fig,fullfile(p.figureDir,['Exceedance_' tag '.fig']));close(fig);
 end
 r=Q(Q.Component==c,:);
 fig=figure('Visible','off','Color','w','Units','centimeters','Position',[5 5 10.3 8.6]);
 ax=axes(fig,'Units','centimeters','Position',[1.7 1.4 7.9 6.7]);hold(ax,'on');
 errorbar(ax,r.Porosity,r.FEM_Mean,r.FEM_SD,'o','Color',[.25 .25 .25],'LineWidth',1.25,'MarkerSize',5);
 plot(ax,r.Porosity,r.GND_Analytical,'--','Color',[.0471 .3059 .6078],'LineWidth',1.5);
 plot(ax,r.Porosity,r.GND_FEMmoments_Mean,'-.','Color',[.78 .133 .157],'LineWidth',1.5);
 plot(ax,r.Porosity,r.Campbell_Full,'-','Color',[.12 .56 .26],'LineWidth',1.5);
 set(ax,'FontName','Times New Roman','FontSize',10,'LineWidth',.75,'Box','on','XTick',0:.025:.1);xlim(ax,[0 .1125]);
 xlabel(ax,'Porosity','FontSize',11);ylabel(ax,'0.99 quantile','FontSize',11);
 legend(ax,{'FEM: mean +/- SD','GND: analytical moments','GND: FEM moments','Full Campbell'},'Location','best','Box','off','FontName','Times New Roman','FontSize',8);
 exportgraphics(fig,fullfile(p.figureDir,sprintf('Q99_S%d.png',c)),'Resolution',600);savefig(fig,fullfile(p.figureDir,sprintf('Q99_S%d.fig',c)));close(fig);
end
end
