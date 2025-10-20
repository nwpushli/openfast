function [Hase, G, sys, omega, inputIdx, outputIdx] = compute_hase_transfer(A, B, C, D, varargin)
%COMPUTE_HASE_TRANSFER 组装 OpenFAST 线性模型并计算 H_ase(j*omega)。
%   [HASE, G, SYS, OMEGA, INPUTIDX, OUTPUTIDX] = COMPUTE_HASE_TRANSFER(A, B, C, D)
%   会根据 OpenFAST 线性化得到的状态空间矩阵生成系统模型，并返回频率响应矩阵。
%   通过名称-数值对参数可以筛选指定的输入/输出通道，并设置频率点。
%
%   必选输入参数：
%     A,B,C,D     - OpenFAST 线性化输出的连续时间状态、输入、输出与旁路矩阵。
%
%   名称-数值对参数：
%     'Inputs'            : 要保留的输入通道索引或名称（默认保留全部）。
%     'Outputs'           : 要保留的输出通道索引或名称（默认保留全部）。
%     'InputDescriptions' : 输入通道的描述字符串（如 FASTLinearizationFile.udescr()）。
%                           当通过名称选择输入通道时需要提供。
%     'OutputDescriptions': 输出通道的描述字符串（如 FASTLinearizationFile.ydescr()）。
%                           当通过名称选择输出通道时需要提供。
%     'Frequencies'       : 指定的角频率向量（单位 rad/s），用于计算 H_ase。
%     'FrequencySpan'     : 长度为 2 的向量 [wMin wMax]（单位 rad/s），在未指定
%                           'Frequencies' 时用于自动生成对数频率网格，范围默认为
%                           根据 A 的特征值估计得到。
%     'NumFrequencyPoints': 自动生成频率网格时的点数（默认 200）。
%
%   输出参数：
%     Hase      - ny×nu×nw 的复数频率响应数组。
%     G         - 控制系统工具箱的传递函数对象。
%     sys       - 选定通道构成的状态空间对象。
%     omega     - 计算时使用的角频率列向量（单位 rad/s）。
%     inputIdx  - 保留下来的输入通道索引。
%     outputIdx - 保留下来的输出通道索引。
%
%   使用示例：
%     file = FASTLinearizationFile('MyTurbine.1.lin');
%     [Hase, ~, sys, omega] = compute_hase_transfer(file.A, file.B, file.C, file.D, ...
%         'Inputs', {'HWindSpeed'}, 'Outputs', {'TwrBsFys'}, ...
%         'InputDescriptions', file.udescr(), 'OutputDescriptions', file.ydescr(), ...
%         'FrequencySpan', [0.01 10]);
%     bode(sys, omega);
%
%   需要 MATLAB 控制系统工具箱。
%
%   相关函数：FASTLINEARIZATIONFILE，SS，TF，FREQRESP，BODE，LSIM。

arguments
    A double
    B double
    C double
    D double
end

arguments (Repeating)
    varargin
end

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Inputs', [], @(x) isnumeric(x) || isstring(x) || ischar(x) || iscellstr(x));
addParameter(p, 'Outputs', [], @(x) isnumeric(x) || isstring(x) || ischar(x) || iscellstr(x));
addParameter(p, 'InputDescriptions', {}, @(x) iscellstr(x) || isstring(x));
addParameter(p, 'OutputDescriptions', {}, @(x) iscellstr(x) || isstring(x));
addParameter(p, 'Frequencies', [], @(x) isnumeric(x));
addParameter(p, 'FrequencySpan', [], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'NumFrequencyPoints', 200, @(x) isnumeric(x) && isscalar(x) && x > 1);
parse(p, varargin{:});
opts = p.Results;

% 基本维度检查
[nx, nAcols] = size(A);
if nAcols ~= nx
    error('矩阵 A 必须为方阵。');
end

[nBrows, nu] = size(B);
[nCrows, nCcols] = size(C);
[nDrows, nDcols] = size(D);

if nBrows ~= nx || nCcols ~= nx || nDrows ~= nCrows || nDcols ~= nu
    error('输入的矩阵维度不一致。');
end

% 解析输入/输出通道选择
inputIdx = resolve_selection(opts.Inputs, nu, opts.InputDescriptions, 'input');
if isempty(inputIdx)
    inputIdx = 1:nu;
end

outputIdx = resolve_selection(opts.Outputs, nCrows, opts.OutputDescriptions, 'output');
if isempty(outputIdx)
    outputIdx = 1:nCrows;
end

% 提取选定通道的子模型
Bsel = B(:, inputIdx);
Csel = C(outputIdx, :);
Dsel = D(outputIdx, inputIdx);

% 构建状态空间对象与传递函数对象
sys = ss(A, Bsel, Csel, Dsel);
G = tf(sys);

% 确定需要计算的频率点
omega = prepare_frequency_grid(A, opts.Frequencies, opts.FrequencySpan, opts.NumFrequencyPoints);

% 调用控制系统工具箱计算频率响应
Hase = freqresp(sys, omega.');

% 将角频率转换为列向量
omega = omega(:);

end

%--------------------------------------------------------------------------
function idx = resolve_selection(selection, nAvailable, descriptors, kind)
    if isnumeric(selection)
        idx = selection(:).';
        if any(idx < 1) || any(idx > nAvailable)
            error('请求的 %s 索引超出了可用通道数量（%d）。', kind, nAvailable);
        end
        idx = unique(idx, 'stable');
        return;
    end

    if isempty(selection)
        idx = [];
        return;
    end

    descriptors = cellstr(descriptors);
    if isempty(descriptors)
        error('为 %s 选择提供了通道名称，但未给出通道描述。', kind);
    end

    if ischar(selection) || isstring(selection)
        selection = cellstr(selection);
    else
        selection = cellstr(string(selection));
    end

    idx = zeros(1, numel(selection));
    for ii = 1:numel(selection)
        match = find(strcmpi(strtrim(selection{ii}), strtrim(descriptors)), 1);
        if isempty(match)
            error('未能找到名为 "%s" 的 %s 通道。', selection{ii}, kind);
        end
        idx(ii) = match;
    end
    idx = unique(idx, 'stable');
end

%--------------------------------------------------------------------------
function omega = prepare_frequency_grid(A, explicit, span, nPoints)
    if ~isempty(explicit)
        omega = explicit(:);
        return;
    end

    if isempty(span)
        eigVals = eig(A);
        eigVals = eigVals(~isnan(eigVals) & ~isinf(eigVals));
        wn = abs(eigVals);
        wn(wn == 0) = [];
        if isempty(wn)
            span = [1e-3, 1e2];
        else
            wmin = max(min(wn)/10, 1e-3);
            wmax = max(max(wn)*10, wmin*10);
            span = [wmin, wmax];
        end
    else
        span = sort(abs(span(:).'));
        if span(1) <= 0
            span(1) = max(span(2)/1e3, 1e-4);
        end
        if span(2) <= span(1)
            span(2) = span(1) * 10;
        end
    end

    omega = logspace(log10(span(1)), log10(span(2)), nPoints).';
end
