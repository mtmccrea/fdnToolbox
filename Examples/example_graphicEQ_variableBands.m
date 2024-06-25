% Example for graphic EQ filter design, with variable number of bands,
%
% All About Audio Equalization: Solutions and Frontiers
% by V. Välimäki and J. Reiss, 
% in Applied Sciences, vol. 6, no. 5, p. 129, May 2016.
%
%
% *Reproduce in Code*
% (c) Sebastian Jiro Schlecht:  Monday, 7. January 2019
% 

% close all; clear; clc;

fftLen = 2^16;
fs = 48000;

nBand = 5;

% Target gains and crossover frequencies of the bands. 
% Note the lowest and highest bands are shelf filters
switch nBand
    case 10
        targetG = [-1 1 -1 1 -1 1 -1 1 -1 1].';
        f_shelf = [42 16000];
        R = 2.7;
    case 6
        targetG = [0 1 -1 1 -1 1].';
        f_shelf = [42 16000];
        R = 3;
    case 5
        targetG = [-1 1 -1 1 -1].';
        f_shelf = [42 6000];
        R = 3.5;
    case 3
        targetG = [-1 1 -1].';
        f_shelf = [42 16000];
        R = 5;
end

% fx = [63, 125, 250, 500, 1000, 2000, 4000, 8000, 11360];
fx = round(logspace(log10(f_shelf(1)), log10(f_shelf(end)), nBand-1));

targetG = 10 * targetG; % scale dynamic range of gains
% targetG = targetG * -1; % invert staggered gains

% constrained uses a lsq solver so is more accurate, but has limited
% applicability in the real world
method = 'unconstrained'; % constrained unconstrained

% center freqs, just for reference (calculated in designGEQ_variableBand)
fc_bnd = zeros(1, length(fx)-1);
for i = 2:length(fx)
    fc_bnd(i-1) = sqrt(fx(i-1) * fx(i));
end
fc_bnd

[optimalSOS, targetF] = designGEQ_variableBand(targetG, fx, R, fs, fftLen, method);

[hOpt,wOpt] = freqz(optimalSOS,fftLen,fs);

% plot
figure; hold on; grid on; 
plot(targetF, targetG);
plot(wOpt,mag2db(abs(hOpt)))
set(gca, 'xScale', 'log')
ylim([-12 12])
xlim([1 fs/1.9])
title('Approximation Magnitude Response')
legend('Target', 'Actual EQ', 'Location','SouthEast');

%% Test: Graphic EQ design
assert( 1 == 1); % script runs without error  
