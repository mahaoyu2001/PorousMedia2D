function q = weighted_upper_quantile(x, w, alpha)
%WEIGHTED_UPPER_QUANTILE  计算加权离散数据的上alpha分位数 q_{1-alpha}
%
% 输入:
%   x     - 离散数据向量 (Nx1 或 1xN)
%   w     - 权重向量，长度与 x 相同，可不归一化，允许 0
%   alpha - 尾部概率 (0<alpha<1)，例如 0.01 / 1e-3
%
% 输出:
%   q          - 上alpha分位数，使得 P(X > q) = alpha (按加权经验分布定义)
%   x_sorted   - 排序后的样本（可选）
%   cdf_sorted - 对应加权经验CDF（可选）
%
% 说明:
%   这里采用“左连续/最小满足”定义：
%     q = inf{x: F(x) >= 1-alpha}
%   对离散分布，这等价于取达到累计概率阈值的第一个样本点。
%
% 示例:
%   q = weighted_upper_quantile(x, w, 1e-3);

    % ---- 输入检查 ----
    if nargin < 3
        error('Need inputs: x, w, alpha.');
    end
    x = x(:);
    w = w(:);
    if numel(x) ~= numel(w)
        error('x and w must have the same length.');
    end
    if ~(isscalar(alpha) && alpha > 0 && alpha < 1)
        error('alpha must be a scalar in (0,1).');
    end

    % 去除 NaN/Inf 和非正权重
    valid = isfinite(x) & isfinite(w) & (w > 0);
    x = x(valid);
    w = w(valid);

    if isempty(x)
        error('No valid samples after filtering (need finite x and positive weights).');
    end

    % ---- 权重归一化 ----
    W = sum(w);
    if W <= 0
        error('Sum of weights must be positive.');
    end
    w = w / W;

    % ---- 排序并计算加权CDF ----
    [x_sorted, idx] = sort(x, 'ascend');
    w_sorted = w(idx);
    cdf_sorted = cumsum(w_sorted);

    % ---- 上alpha分位数: 找到最小x使得CDF >= 1-alpha ----
    p = 1 - alpha;
    k = find(cdf_sorted >= p, 1, 'first');
    if isempty(k)
        % 理论上不会发生（cdf_sorted(end)=1），除非数值问题
        q = x_sorted(end);
    else
        q = x_sorted(k);
    end
end
