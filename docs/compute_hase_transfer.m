function [Hase, G, sys, omega, inputIdx, outputIdx] = compute_hase_transfer(A, B, C, D, varargin)
%COMPUTE_HASE_TRANSFER Assemble OpenFAST linear models and evaluate H_ase(j*omega).
%   [HASE, G, SYS, OMEGA, INPUTIDX, OUTPUTIDX] = COMPUTE_HASE_TRANSFER(A, B, C, D)
%   builds the state-space system defined by the linearized OpenFAST matrices and returns
%   its frequency-response matrix. Optional name-value arguments let you pick specific
%   input/output channels and define the frequency grid.
%
%   Required inputs:
%     A,B,C,D     - Continuous-time state, input, output, and feedthrough matrices taken
%                   from an OpenFAST linearization.
%
%   Name-value arguments:
%     'Inputs'            : Indices or names of input channels to retain (default: all).
%     'Outputs'           : Indices or names of output channels to retain (default: all).
%     'InputDescriptions' : Cell array of character vectors describing each input (as
%                           provided by FASTLinearizationFile.udescr()). Needed when
%                           selecting inputs by name.
%     'OutputDescriptions': Cell array of character vectors describing each output (from
%                           FASTLinearizationFile.ydescr()). Needed when selecting outputs
%                           by name.
%     'Frequencies'       : Vector of frequencies (rad/s) at which to evaluate H_ase.
%     'FrequencySpan'     : Two-element vector [wMin wMax] (rad/s). Used to generate a
%                           logarithmic grid when 'Frequencies' is omitted. Defaults to a
%                           span inferred from the eigenvalues of A.
%     'NumFrequencyPoints': Number of points in the generated frequency grid (default 200).
%
%   Outputs:
%     Hase      - ny-by-nu-by-nw array containing the complex frequency response values.
%     G         - Transfer-function model (Control System Toolbox TF object).
%     sys       - State-space model (SS object) constructed from the selected channels.
%     omega     - Column vector of frequencies (rad/s) used for the evaluation.
%     inputIdx  - Numeric indices of the retained input channels.
%     outputIdx - Numeric indices of the retained output channels.
%
%   Example:
%     file = FASTLinearizationFile('MyTurbine.1.lin');
%     [Hase, ~, sys, omega] = compute_hase_transfer(file.A, file.B, file.C, file.D, ...
%         'Inputs', {'HWindSpeed'}, 'Outputs', {'TwrBsFys'}, ...
%         'InputDescriptions', file.udescr(), 'OutputDescriptions', file.ydescr(), ...
%         'FrequencySpan', [0.01 10]);
%     bode(sys, omega);
%
%   Requires MATLAB Control System Toolbox.
%
%   See also FASTLINEARIZATIONFILE, SS, TF, FREQRESP, BODE, LSIM.

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

% Basic dimension checks
[nx, nAcols] = size(A);
if nAcols ~= nx
    error('Matrix A must be square.');
end

[nBrows, nu] = size(B);
[nCrows, nCcols] = size(C);
[nDrows, nDcols] = size(D);

if nBrows ~= nx || nCcols ~= nx || nDrows ~= nCrows || nDcols ~= nu
    error('Inconsistent matrix dimensions.');
end

% Resolve channel selections
inputIdx = resolve_selection(opts.Inputs, nu, opts.InputDescriptions, 'input');
if isempty(inputIdx)
    inputIdx = 1:nu;
end

outputIdx = resolve_selection(opts.Outputs, nCrows, opts.OutputDescriptions, 'output');
if isempty(outputIdx)
    outputIdx = 1:nCrows;
end

% Extract the sub-model
Bsel = B(:, inputIdx);
Csel = C(outputIdx, :);
Dsel = D(outputIdx, inputIdx);

% Assemble state-space and transfer-function objects
sys = ss(A, Bsel, Csel, Dsel);
G = tf(sys);

% Determine the evaluation frequencies
omega = prepare_frequency_grid(A, opts.Frequencies, opts.FrequencySpan, opts.NumFrequencyPoints);

% Evaluate frequency response (Control System Toolbox)
Hase = freqresp(sys, omega.');

% Return omega as a column vector
omega = omega(:);

end

%--------------------------------------------------------------------------
function idx = resolve_selection(selection, nAvailable, descriptors, kind)
    if isnumeric(selection)
        idx = selection(:).';
        if any(idx < 1) || any(idx > nAvailable)
            error('Requested %s index exceeds available channels (%d).', kind, nAvailable);
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
        error('Channel names were provided for the %s selection, but no descriptors were supplied.', kind);
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
            error('Unable to find %s channel named "%s".', kind, selection{ii});
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
