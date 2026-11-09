function round2_distribution_axis_override(ax,name)
% User axis override; horizontal range matches cdf_2DVoid_Vol1_a_stress2L.
if strcmp(name,'pdf_2DVoid_Vol1_a_stress2L')
 set(ax,'XLim',[0.6 1.42],'YLim',[0 20],'XTick',[0.6 0.7999999999999999 1 1.2 1.4],'YTick',0:5:20);
elseif strcmp(name,'Q_EaxtVal_vol_stress1_errorband')
 set(ax,'YLim',[0 0.55]);
end
end
