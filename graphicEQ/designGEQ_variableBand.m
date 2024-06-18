function [sos, targetF] = designGEQ_variableBand ( targetG, fx_band, R, fs, fftLen, method)
%designGEQ - Graphic EQ with 10 bands
% All About Audio Equalization: Solutions and Frontiers
% by V. Välimäki and J. Reiss, 
% in Applied Sciences, vol. 6, no. 5, p. 129, May 2016.
%
% Syntax:  [sos, targetF] = designGEQ ( targetG )
%
% Inputs:
%    targetG - target magnitude response in dB of size [length(fx_bnd)+1,1]
%    fx_band - crossover frequencies of the bands. Note the lowest and highest bands are shelf filters
%    fs - sampling rate
%    R - relates to bandwidth Q by Q = sqrt(R) / (R-1). 10-band default was R = 2.7;
%    fftLen = 
%    
%    method = 'constrained' or 'unconstrained'
% Outputs:
%    sos - Filter coefficients as SOS of size [10,6]
%    targetF - Band Center frequencies of size [1,10]
%
% Example: 
%    [sos, targetF] = designGEQ ( 1:10 )
%
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Dr.-Ing. Sebastian Jiro Schlecht, 
% Aalto University, Finland
% email address: sebastian.schlecht@aalto.fi
% Website: sebastianjiroschlecht.com
% 23. October 2020; Last revision: 23. October 2020


% Initialization
% Setup variables such as sampling frequency, center frequencies and
% control frequencies, and command gains
if nargin < 5 || isempty(method)
    method = 'constrained';
end
if nargin < 5 || isempty(R)
    R = 2.7;
end
if nargin < 4 || isempty(fftLen)
    fftLen = 2^16;
end
if nargin < 3 || isempty(fs)
    fs = 48000;
end
if nargin < 2 || isempty(fx_band)
    fx_band = [62.5, 250, 1000, 4000, 16000];
end

if (length(targetG) ~= length(fx_band)+1), error("targetG must have 1 element more than fx_band"); end

% center freqs of peak eq filters
nPeakFilts = length(fx_band)-1;
fc_band = zeros(1, nPeakFilts);
for i = 2:length(fx_band)
    fc_band(i-1) = sqrt(fx_band(i-1) * fx_band(i));
end
% centerFrequencies = [ 63, 125, 250, 500, 1000, 2000, 4000, 8000]; % Hz (default)
centerFrequencies = fc_band; % Hz

% ShelvingCrossover = [46 11360]; % Hz
ShelvingCrossover = fx_band([1 end]); % Hz

numFreq = length(centerFrequencies) + length(ShelvingCrossover);
shelvingOmega = hertz2rad(ShelvingCrossover, fs);
centerOmega = hertz2rad(centerFrequencies, fs);
% R = 2.7; % 10-band default

% control frequencies are spaced logarithmically
numControl = 100;
controlFrequencies = round(logspace(log10(1), log10(fs/2.1), numControl+1));

% target magnitude response via command gains
targetF = [1, centerFrequencies fs];
% targetG = 10 * [1; -1; 1; -1; 1; -1; 1; -1; 1; 1]; % dB
targetInterp = interp1(targetF, targetG, controlFrequencies)';

% desgin prototype of the biquad sections
prototypeGain = 10; % dB
prototypeGainArray = prototypeGain * ones(numFreq+1,1);
prototypeSOS = graphicEQ(centerOmega, shelvingOmega, R, prototypeGainArray);
[G,prototypeH,prototypeW] = probeSOS (prototypeSOS, controlFrequencies, fftLen, fs);
G = G / prototypeGain; % interaction matrix: dB vs control frequencies

switch method
    case 'constrained'
        % compute optimal parametric EQ gains
        % Either you can use a unconstrained linear solver or introduce gain bounds
        % at [-20dB,+20dB] with acceptable deviation from the self-similarity
        % property. The plot shows the deviation between design curve and actual
        % curve.
        upperBound = [Inf, 2 * prototypeGain * ones(1,numFreq)];
        lowerBound = -upperBound;

        opts = optimset('Display','off');
        optG = lsqlin(G, targetInterp, [],[],[],[], lowerBound, upperBound, [], opts);
    case 'unconstrained'
        optG = G\targetInterp; % unconstrained solution
    otherwise
        error("method must be 'constrained' or 'unconstrained'.")
end
% optG

sos = graphicEQ( centerOmega, shelvingOmega, R, optG );

% disp('targetG size'); disp(size(targetG))
% disp('G size'); disp(size(G))
% disp('targetInterp size'); disp(size(targetInterp))

% optG2 = pinv(G)*targetInterp;
% disp('optG2 size'); disp(size(optG2))
% optG2
% optG2 - optG 

% optG3 = lsqlin( G./targetInterp, ones(size( targetInterp)), [],[],[],[], lowerBound, upperBound, [], opts);
% optG3 = optG3 .* [targetG(1); targetG']
% disp('optG size'); disp(size(optG3))
% optG3