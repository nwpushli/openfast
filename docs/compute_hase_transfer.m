function [Hase, G, sys] = compute_hase_transfer(A, B, C, D, inputChannels, outputChannels, omega)
%COMPUTE_HASE_TRANSFER Compute frequency-domain transfer function from ABCD matrices.
%   [HASE, G, SYS] = COMPUTE_HASE_TRANSFER(A, B, C, D, INPUTCHANNELS, OUTPUTCHANNELS, OMEGA)
%   builds the state-space model defined by the linearized OpenFAST matrices and returns
%   the multi-input multi-output transfer function values between the selected channels.
%
%   INPUTS:
%     A, B, C, D         - Linearized system matrices from an OpenFAST .lin file.
%     inputChannels     - Vector of indices selecting which input channels to keep. These
%                         indices correspond to the order returned by FASTLinearizationFile.udescr().
%     outputChannels    - Vector of indices selecting which outputs to observe. These indices
%                         correspond to FASTLinearizationFile.ydescr().
%     omega             - Column vector of frequencies (rad/s) at which to evaluate the
%                         frequency response. If empty or omitted, a logarithmic grid spanning
%                         the system dynamics is generated automatically.
%
%   OUTPUTS:
%     Hase - 3-D array of complex frequency-response values with dimensions
%            [numOutputs x numInputs x numFrequencies]. For SISO selections, use squeeze(Hase)
%            to obtain a row vector over omega.
%     G    - Transfer function model (TF object) of the selected channels.
%     SYS  - State-space model (SS object) of the selected channels.
%
%   EXAMPLE:
%     % Load matrices from a .mat file exported from FASTLinearizationFile
%     load('5MW_linear_model.mat', 'A', 'B', 'C', 'D', 'u_desc', 'y_desc');
%     windInput = find(strcmp(u_desc, 'HWindSpeed'));   % wind-speed input index
%     towerBaseFy = find(strcmp(y_desc, 'TwrBsFys'));   % tower-base shear force output index
%     omega = logspace(-2, 2, 200);                     % rad/s grid
%     [Hase, G] = compute_hase_transfer(A, B, C, D, windInput, towerBaseFy, omega);
%     bode(G, omega);                                   % plot magnitude/phase
%
%   This helper follows the workflow described in the fatigue-load estimation paper: it
%   extracts the subset of inputs/outputs, forms the transfer function H_ase(j*omega), and
%   evaluates the frequency response on the requested grid. Additional blocks such as the
%   spectral shaping of wind or wave inputs can be applied by multiplying the resulting
%   frequency response with the corresponding PSDs.
%
%   Requires MATLAB Control System Toolbox.
%
%   See also SS, TF, FREQRESP, BODE.

arguments
    A double
    B double
    C double
    D double
    inputChannels (1,:) {mustBeInteger, mustBePositive}
    outputChannels (1,:) {mustBeInteger, mustBePositive}
    omega double = []
end

% Validate dimensions
nx = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

if size(A,2) ~= nx
    error('Matrix A must be square.');
end
if size(B,1) ~= nx || size(D,1) ~= ny || size(D,2) ~= nu || size(C,2) ~= nx
    error('Inconsistent dimensions among A, B, C, D matrices.');
end

% Ensure channel indices are within bounds
if any(inputChannels > nu)
    error('inputChannels index exceeds the number of available inputs (%d).', nu);
end
if any(outputChannels > ny)
    error('outputChannels index exceeds the number of available outputs (%d).', ny);
end

% Extract the selected sub-system
Bsel = B(:, inputChannels);
Csel = C(outputChannels, :);
Dsel = D(outputChannels, inputChannels);

% Construct state-space model and equivalent transfer function
sys = ss(A, Bsel, Csel, Dsel);
G = tf(sys);

% Generate a default frequency grid if needed
if isempty(omega)
    % Use eigenvalues of A to determine an appropriate span
    eigVals = eig(A);
    wn = abs(eigVals(real(eigVals) < 0));
    if isempty(wn)
        wn = abs(eigVals);
    end
    wn(wn == 0) = [];
    if isempty(wn)
        omega = logspace(-2, 2, 200);
    else
        wmin = max(min(wn)/10, 1e-4);
        wmax = max(wn)*10;
        if wmax <= wmin
            wmax = wmin*1e3;
        end
        omega = logspace(log10(wmin), log10(wmax), 200);
    end
else
    omega = omega(:)';
end

% Evaluate the frequency response Hase(j*omega)
Hase = freqresp(sys, omega);

end
