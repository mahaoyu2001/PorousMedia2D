%% 参数配置函数（严格匹配输出参数顺序）
function [curLineWidth, axLineWidth, markerSize, ...
    axTickFontSize, axLabelFontSize, lgdFontSize, ...
    axTickFontName, axLabelFontName, lgdFontName, colorList] = plotParams(varargin)
    % 默认参数值（按输出顺序排列）
    defaults = struct(...
        'curLineWidth', 1.5, ...        % 当前线宽
        'axLineWidth', 0.75, ...            % 坐标轴线宽
        'markerSize', 1, ...             % 标记大小
        'axTickFontSize', 9.5, ...      % 坐标轴刻度字体大小
        'axLabelFontSize', 11, ...       % 坐标轴标签字体大小
        'lgdFontSize', 11, ...           % 图例字体大小
        'axTickFontName', 'Times New Roman', ...      % 坐标轴刻度字体名称
        'axLabelFontName', 'Times New Roman', ...     % 坐标轴标签字体名称
        'lgdFontName', 'Times New Roman', ...         % 图例字体名称
        'colorGroup', 'default' ...
    );

    colorSchemes = containers.Map(...
    {'default', 'matlab', 'multi'}, ...  % 颜色组名称
    {...
        { ... % default
          [0.7804 0.1333 0.1569];
          [0.0471 0.3059 0.6078];
          [0.9765 0.5608 0.2039];}; ...

        { ... % matlab default
        [0 0.4470 0.7410];    % 蓝色 (#0072BD)
        [0.8500 0.3250 0.0980]; % 橙色 (#D95319)
        [0.9290 0.6940 0.1250]; % 黄色 (#EDB120)
        [0.4940 0.1840 0.5560]; % 紫色 (#7E2F8E)
        [0.4660 0.6740 0.1880]; % 树绿 (#77AC30)
        [0.3010 0.7450 0.9330]; % 浅天蓝 (#4DBEEE)
        [0.6350 0.0780 0.1840]; % 苝褐 (#A2142F)
        };...

        { ... % multi
          [0.7804 0.1333 0.1569];  % 1. 深红色
          [0.0471 0.3059 0.6078];  % 2. 深蓝色
          [0.1176 0.5647 0.2627];  % 3. 深绿色
          [0.8902 0.4471 0.1333];  % 4. 橙色
          [0.5804 0.1569 0.5569];  % 5. 紫红色
          [0.1725 0.6275 0.7922];  % 6. 天蓝色
          [0.8392 0.6784 0.1529];  % 7. 金色
          [0.5490 0.3373 0.2941];  % 8. 棕色
          [0.0902 0.7451 0.8118];   % 10. 青蓝色
          [0.4980 0.4980 0.4980];  % 9. 中性灰
        };
    });

    % 解析输入参数
    p = inputParser;
    addParameter(p, 'curLineWidth', defaults.curLineWidth, @(x) isscalar(x)>0);
    addParameter(p, 'axLineWidth', defaults.axLineWidth, @(x) isscalar(x)>0);
    addParameter(p, 'markerSize', defaults.markerSize, @(x) isscalar(x)>0);
    addParameter(p, 'axTickFontSize', defaults.axTickFontSize, @(x) isscalar(x)>0);
    addParameter(p, 'axLabelFontSize', defaults.axLabelFontSize, @(x) isscalar(x)>0);
    addParameter(p, 'lgdFontSize', defaults.lgdFontSize, @(x) isscalar(x)>0);
    addParameter(p, 'axTickFontName', defaults.axTickFontName, @(x) isscalar(x)>0);
    addParameter(p, 'axLabelFontName', defaults.axLabelFontName, @(x) isscalar(x)>0);
    addParameter(p, 'lgdFontName', defaults.lgdFontName, @(x) isscalar(x)>0);
    addParameter(p, 'colorGroup', defaults.colorGroup, @ischar);  % 新增参数
    parse(p, varargin{:});

    % 按输出顺序赋值
    curLineWidth = p.Results.curLineWidth;
    axLineWidth = p.Results.axLineWidth;
    markerSize = p.Results.markerSize;
    axTickFontSize = p.Results.axTickFontSize;
    axLabelFontSize = p.Results.axLabelFontSize;
    lgdFontSize = p.Results.lgdFontSize;
    axTickFontName = p.Results.axTickFontName;
    axLabelFontName = p.Results.axLabelFontName;
    lgdFontName = p.Results.lgdFontName;
    colorList = colorSchemes(p.Results.colorGroup);
end
