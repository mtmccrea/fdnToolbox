% Comparing multi-shelf GEQ to symmetric peaking GEQ (Liski, 2019)
% Modified based on Figure8_3Case_1st2ndPeak of the MultiShelfGEQ library
%
% Tantep Sinjanakhom, 20 September 2023
% Sebastian J. Schlecht, Friday, 17 November 2023
clear; clc; %close all;

fs = 44100;

designShelfFilter = @(gain, fc,fs,k) designHigherOrderShelvingFilter(fc , fs, gain, k);

gainInterpType = 'spline'; % makima (>2 pts) spline (>4 pts) linear (any num pts) cubic (>3) pchip (>4)
numfilts = 4;
osFactor = 2; % oversample factor for interaction frequencies (integer aligns to all control freqs, fractional aligns to endpoints at least)
numInterx = osFactor*(numfilts-1) + 1; % for equal distribution within control freq endpoints
lowFreq = 31.25;
highFreq = 16e3; %fs/2-1;

% band center freqs, which neighbor shelf break freqs (these are SGE peak filter centers)
% frqControl = [31.25/4, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1]';
frqControl = logspace(log10(lowFreq), log10(highFreq), numfilts)';

% Design prototype shelving filters

% shelf transition/"break" frequencies are the geometric mean of the adjacent band centers (control freqs)
frqShelf = geometricMeanPoints([0 ; frqControl]);

shelfFilterOrder = 2;
prototypeGain = 1; % 1 dB is recommended
numShelfFilters = numel(frqShelf);

% "interaction" frequencies, at which gains are optimized.
% They're the control freqs interleaved with the shelf freqs (midpoints between control freqs)
% Note: lowest and highest interaction freqs comes from control freq list
% frqInterx(1:2:numShelfFilters*2-1) = frqControl;
% frqInterx(2:2:numShelfFilters*2-1) = frqShelf(2:end); % discards lowest
frqInterx = round(logspace(log10(lowFreq), log10(highFreq), numInterx))';

% shelf filter coefficients
for j = 1 : numShelfFilters
    [B(j,:),A(j,:)] = designShelfFilter(db2mag(prototypeGain), frqShelf(j), fs, shelfFilterOrder);
end

% Evaluate the frequency response of the shelf filters at the interaction
% freqs, normalize by the prototype gain.
[P_interx,w] = freqzVec(B,A,frqInterx,fs);
% Normalize the prototype shelf filters' response by the prototype gain.
gainInterx_proto = mag2db(abs(P_interx)) ./ prototypeGain;

% Inspect filter interaction freq response curves
figure; semilogx(frqInterx, gainInterx_proto); legend

% Design target filter response

% Frequencies at which target gains are defined
% frqTarget = [0, 31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1];
frqTarget = [0, frqControl'];

figure();

% Compute different cases
tl = tiledlayout(3, 1);
tiles = cell(3, 1);
gainCase = 1;


% The spec used to define target gains
% Different from command gains of the generated filters, these can be
% arbitrary resolution, but will be resampled to the
targetFreqSpec = [0, 31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1];
switch gainCase
    case 1 % like Two Stage reverb example
        gainTargetSpec = [0 -1 -3 -10 -16 -18 -17 -12 -13 -15 -17 -20];
    case 2 % high dynamic range
        gainTargetSpec = linspace(0, -60, 12);
    case 3 % ZigZag
        gainTargetSpec = 5 * [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];
end
gainTarget = interp1(targetFreqSpec, gainTargetSpec, frqTarget, 'makima','extrap').';

% Interpolate the requested target gains (at target frequencies)
% to the control freqs (for SGE) and interaction matrix freqs
% gainInterx_target = interp1(frqTarget, gainTarget, frqInterx, gainInterpType,'extrap').';
% Increase accuracy by setting interaction matrix target gains by
% interpolating in _linear_ space
% Note the interpolation used here assumes frqInterx maps linearly between the endpoints of frqTarget
gainInterx_target = interp1( ...
    linspace(0,1,numel(frqControl)), gainTarget(2:end), ...
    linspace(0,1,numel(frqInterx)), ...
    'linear','extrap').';
% % Note for simple 2 x interaction freq oversampling, these gains can just
% % be the control gains, interleaved with the mean of neighboring gains
% gainInterx_target = zeros(numInterx,1);
% gainInterx_target([1 end]) = gainTarget([2 end]);
% for i = 1:numfilts-1
%     gainInterx_target(2*i) = mean(gainTarget([i i+1]+1)); % +1 because gainTarget(1) is not used
%     gainInterx_target(2*i+1) = gainTarget(i+1 +1);
% end


% Solve for filter gains by:
% constrained least squares
maxGain = 100;
ub = maxGain * ones(numShelfFilters,1);
lb = -ub;
ub(1) = Inf; % broadband gain
lb(1) = -Inf;
gains = lsqlin(gainInterx_proto, gainInterx_target,[],[],[],[],lb,ub);
% or unconstrained least squares
% gains = gainProtoInterx\frqTargetInterx_interp; % see designGEQ.m


% 1st order shelving
clear B A;
for j = 1 : length(frqShelf)
    [B(j,:),A(j,:)] = designShelfFilter(db2mag(gains(j)), frqShelf(j),fs,1);
end
[F1,w] = freqzVec(B,A,fs,fs);
FFShelf1 = mag2db(abs(prod(F1,2)));

% 2nd order shelving
clear B A;
for j = 1 : length(frqShelf)
    [B(j,:),A(j,:)] = designShelfFilter(db2mag(gains(j)), frqShelf(j),fs,2);
end
[F2,w] = freqzVec(B,A,fs,fs);
FFShelf2 = mag2db(abs(prod(F2,2))); % mag response of the combined filters


% target: interpolated from target freq/gains across the whole frequency range
% targetPCHIP = makima(frqTarget, gainTarget, w);
% targetPCHIP = interp1(frqTarget, gainTarget, w, gainInterpType, 'extrap');
% for error calculation, extend the first/last target gain to DC/Nyquist to
% avoid misleading error values at end of spectrum
targetPCHIP = interp1( ...
    [0; frqInterx; fs/2], [gainInterx_target(1); gainInterx_target; gainInterx_target(end)], ...
    w, gainInterpType, 'extrap'); % todo: this should probably also be interpolated in log space

% PLOT


% plot magnitude of individual and aggregate shelves (second order)
tiles{1} = nexttile; hold on; box on;

plot(w, mag2db(abs(F2)))
plot(w, mag2db(abs(prod(F2,2))), LineWidth=2)
xline(frqControl, color=[0.98 0.1 0.10], alpha=0.1, LineWidth=1);
title(gca, "second order shelves")
ylabel('Magnitude [dB]')


% plot magnitude response with control and interaction points
tiles{2} = nexttile; hold on; box on;

plot(w, FFShelf1,    ':', Color=[0.3 0.65 0.35], DisplayName='first-order cascade');
plot(w, FFShelf2,    '-', Color=[1 0.45 0.0],    DisplayName='second-order cascade');
plot(w, targetPCHIP, '-', Color=0.8*[1 1 1],     DisplayName='target response', LineWidth=1);
% gain control markers
plot(frqControl, gainTarget(2:end), 'o', Color=[0.98 0.1 0.10], MarkerSize=5, LineWidth=1.5, DisplayName='command gain');
plot(frqInterx,  gainInterx_target, 'x', Color=0.6*[1 1 1],     MarkerSize=5, LineWidth=1.2, DisplayName='command gain (interx)');
xline(frqControl, color=[0.98 0.1 0.10], alpha=0.1, LineWidth=1, HandleVisibility="off");

set(gca,'YLim', [min(gainTarget) max(gainTarget)] + 5*[-1 1]);
ylabel('Magnitude [dB]')
legend


% plot error
tiles{3} = nexttile; hold on; box on;

errShelf1 = targetPCHIP - FFShelf1;
errShelf2 = targetPCHIP - FFShelf2;

maxError = max([ ...
    max(abs(errShelf1)), ...
    max(abs(errShelf2))  ]);
maxError = ceil(maxError * 4) / 4;

% freq indices over which the filter response and error was evalutated
% (at fs number of points)
controlResponseIdc = round(frqControl * 2)+1;

% first-order shelf error
plot(w,abs(errShelf1),':',Color=[0.3 0.65 0.35])
plot(w(controlResponseIdc), abs(errShelf1(controlResponseIdc)), 'o', color=[0.3 0.65 0.35], MarkerSize=4, LineWidth=1);

% second-order shelf error
plot(w,abs(errShelf2),'-',Color=[1.00 0.45 0.00])
plot(w(controlResponseIdc), abs(errShelf2(controlResponseIdc)), 'o', color=[0.98 0.1 0.10], MarkerSize=5, LineWidth=1.5);
xline(frqControl, color=[0.98 0.1 0.10], alpha=0.1, LineWidth=1);

xlabel('Frequency [Hz]')
ylabel('Absolute error [dB]')
set(gca,'YLim', [-0.1 4]); % db error


title(tl, sprintf("Number of - filters: %d, interaction freqs: %d", numShelfFilters, length(frqInterx)))

for i = 1:numel(tiles)
    tile = tiles{i};
    set(tile,'XScale','log')
    set(tile,'XTick',[30 100 300 1000 3000 10000 2e4]);
    set(tile,'XTicklabels',{'30','100','300','1k','3k','10k','20k'});
    set(tile,'XLim', [3 sqrt(fs/2 * fs)]);
end

% [F2,w] = freqzVec(B,A,fs,fs);
% FFShelf2 = mag2db(abs(prod(F2,2))); % mag response of the combined filters

%% Print Figures
set(gcf,'Units', 'inches', 'Position', [0 0 7 4.5]);
exportgraphics(gcf,'./Figures/Figure8_3_Cases_1st_2nd_peak.pdf')
