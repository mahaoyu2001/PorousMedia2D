function [hist_data, cdf_data] = weighted_hist_cdf(data, weights, varargin)
% WEIGHTED_HIST_CDF 计算带权重的直方图和累积分布函数（CDF）
%
% 输入:
%   data     - 观测数据向量 (N×1)
%   weights  - 对应权重向量 (N×1)，自动归一化
%   'BinWidth' - 分箱宽度 (默认: 0.005)
%   'XLim'     - x轴范围 [xmin, xmax] (默认: [min(data), max(data)])
%
% 输出:
%   hist_data - 结构体: binEdges, binCenters, pdf, weightedCounts
%   cdf_data  - 结构体:
%                 .x_grid    → 等间距网格 (binEdges)
%                 .y_grid    → 对应 CDF 值（右连续）
%                 .x_step    → 阶梯函数 x（原始，可选）
%                 .y_step    → 阶梯函数 y（原始，可选）

    % 参数解析
    p = inputParser;
    addRequired(p, 'data', @isnumeric);
    addRequired(p, 'weights', @isnumeric);
    addParameter(p, 'BinWidth', 0.005, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'XLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    parse(p, data, weights, varargin{:});

    data = data(:);
    weights = weights(:);

    % 归一化权重
    weights = weights / sum(weights);

    % 设置 x 范围
    if isempty(p.Results.XLim)
        xlim = [min(data), max(data)];
    else
        xlim = p.Results.XLim;
    end

    % 构建 bin edges（确保覆盖整个 xlim）
    binWidth = p.Results.BinWidth;
    binEdges = xlim(1):binWidth:xlim(2);

    % 如果最后一个 bin 不覆盖 xlim(2)，手动添加
    if binEdges(end) < xlim(2)
        binEdges(end+1) = xlim(2);
    end

    nBins = length(binEdges) - 1;
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

    %% 计算加权直方图
    weightedCounts = zeros(nBins, 1);
    for i = 1:nBins
        if i == nBins
            % 最后一个 bin 包含右边界（闭区间）
            inBin = (data >= binEdges(i)) & (data <= binEdges(i+1));
        else
            inBin = (data >= binEdges(i)) & (data < binEdges(i+1));
        end
        weightedCounts(i) = sum(weights(inBin));
    end

    % PDF（概率密度）
    pdf_weights = weightedCounts / binWidth;

    %% 计算 CDF 在 binEdges 上的值（右连续）
    % 方法：对每个 binEdge，CDF = P(X <= edge)
    cdf_at_edges = zeros(size(binEdges));
    for i = 1:length(binEdges)
        cdf_at_edges(i) = sum(weights(data <= binEdges(i)));
    end

    %% （可选）保留原始阶梯 CDF 用于精确绘图
    cdf_x_step = [];
    cdf_y_step = [];
    % [sortedData, sortIdx] = sort(data);
    % sortedWeights = weights(sortIdx);
    % cumulativeWeights = cumsum(sortedWeights);
    %
    % cdf_x_step = [xlim(1); sortedData; xlim(2)];
    % cdf_y_step = [0; cumulativeWeights; 1];

    %% 输出结构
    hist_data = struct(...
        'binEdges',    binEdges, ...
        'binCenters',  binCenters, ...
        'pdf',         pdf_weights, ...
        'weightedCounts', weightedCounts);

    cdf_data = struct(...
        'x',      binEdges, ...          % 等间距网格（推荐用于分析）
        'y',      cdf_at_edges, ...      % CDF 在 binEdges 上的值（右连续）
        'x_step',      cdf_x_step, ...        % 原始阶梯 x（用于 plot(cdf_x_step, cdf_y_step, 'Steps','post')）
        'y_step',      cdf_y_step);
end
