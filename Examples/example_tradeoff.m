% Example for FDN Tradeoff
%
% FDN design typically needs to balance modal and echo density with
% computational complexity:
% - the longer the delays, the more modes, but less echo density
% - the more delays, the higher modal and echo density, but more expensive
%
% Schlecht, S. (2020). FDNTB: The Feedback Delay Network Toolbox,
% Proceedings of the 23rd International Conference on Digital Audio Effects
% (DAFx-20)
%
% Sebastian J. Schlecht, Monday, 06 March 2023
clear; clc; close all;

rng(2)
fs = 48000;
irLen = 2*fs;
gainPerSample = db2mag(RT602slope(2,fs));


% We compare 3x3 settings
Ns = [4 8 16]; % FDN size: small, medium, large
Ds = round((rand(1,16) + 0.5) .* [300; 1000; 3000]); % delays: short, medium, long

% {'orthogonal','Hadamard','circulant','Householder','parallel','series', ...
%       'diagonalConjugated','tinyRotation','allpassInFDN','nestedAllpass', ...
%           'SchroederReverberator','AndersonMatrix'};
matrixType = 'Hadamard';
matrixType = 'allpassInFDN';
matrixType = 'nestedAllpass';
matrixType = 'AndersonMatrix';
matrixType = 'circulant';
matrixType = 'orthogonal';

COMPLETE = 1; % to generate IRs of completed random orthogonal matrices


% generate impulse responses
for itN = 1:size(Ns,2)
    for itD = 1:size(Ds,1)
        % Define FDN
        N = Ns(itN);
        delays = Ds(itD,1:N);
        numInput = 1;
        numOutput = 1;
        g = diag(gainPerSample.^delays);
        A = fdnMatrixGallery(N,matrixType) * g;
        B = ones(N,numInput);
        C = ones(numOutput,N);
        D = zeros(numOutput,numInput);

        rir(:,itN,itD) = dss2impz(irLen,delays,A,B,C,D); 
        
        dir = sprintf('renders/FDNs/%s/', matrixType);
        
        if ~exist(dir,'dir'), mkdir(dir); end
        audiowrite( ...
            sprintf('%stradeoff_N%d_D%d.wav',dir,N,round(mean(delays))), ...
            0.5*rir(:,itN,itD),fs);
    end
    % sprintf('/renders/FDNs/%s/tradeoff_N%d_D%d.wav',matrixType,N,round(mean(delays))), ...
end

if COMPLETE && strcmp(matrixType, 'orthogonal')
    % generate impulse responses
    for itN = 1:size(Ns,2)
        for itD = 1:size(Ds,1)
            % Define FDN
            N = Ns(itN);
            delays = Ds(itD,1:N);
            numInput = 1;
            numOutput = 1;
            g = diag(gainPerSample.^delays);
            % A = fdnMatrixGallery(N,matrixType) * g;
            % B = ones(N,numInput);
            % C = ones(numOutput,N);
            % D = zeros(numOutput,numInput);

            numIO = numInput;
            V = randomOrthogonal(N + numIO);
            A = V(1:N,1:N);% * g; % g: include delay proportional gain
            
            % A = randomOrthogonal(N);

            [B,C,D,PP] = completeOrthogonal(A, numIO, 'verbose', true);
            % VV = [A,b;c d];

            A = A*g; % include delay-proportional gains

            % [B,C,D,V] = completeOrthogonal(A, numInput);
            % A = V(1:N,1:N);
            % C=C';
            % B=B';
            % C = C * N; 
            D = 0;
            C = C * db2mag(8);

            rir(:,itN,itD) = dss2impz(irLen,delays,A,B,C,D);

            dir = sprintf('renders/FDNs/%s/', [matrixType '_COMPLETE']);

            if ~exist(dir,'dir'), mkdir(dir); end
            audiowrite( ...
                sprintf('%stradeoff_N%d_D%d.wav',dir,N,round(mean(delays))), ...
                0.5*rir(:,itN,itD),fs);

        end
    end
    % sprintf('/renders/FDNs/%s/tradeoff_N%d_D%d.wav',matrixType,N,round(mean(delays))), ...
end




%% Test: Script completed
assert( true ) 

