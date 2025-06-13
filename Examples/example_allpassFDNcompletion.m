% Example for allpass FDN completion problem. For a given feedback matrix
% A, the goal is to construct b,c, and d, such that the FDN is uniallpass.
%
% see "Allpass Feedback Delay Networks", Sebastian J. Schlecht, submitted
% to IEEE TRANSACTIONS ON SIGNAL PROCESSING.
%
% Sebastian J. Schlecht, Tuesday, 9. June 2020
%
clear; clc; close all;
rng(2);


%% Test: complete orthogonal
N = 3;
numIO = 2;
V = randomOrthogonal(N + numIO);
A = V(1:N,1:N);

[b,c,d,PP] = completeOrthogonal(A, numIO, 'verbose', true);
VV = [A,b;c d];

disp('is Orthogonal?')
VV*VV' % is orthogonal if it's an identity matrix
assert(isAlmostZero(VV*VV' - eye(size(VV)),'tol',10^-6))

%% Test-mtm mod: complete orthogonal
clear rir

N = 4;
numIO = 1;
V = randomOrthogonal(N + numIO);
Aortho = fdnMatrixGallery(N, 'orthogonal'); % randomOrthogonal(N);
Acirculant = fdnMatrixGallery(N, 'circulant');

% rng(2)
fs = 48000;
T60 = 2.0;
irLen = T60*fs;
gainPerSample = db2mag(RT602slope(T60,fs));

Ds = round((rand(1,16) + 0.5) .* [500; 1000; 3000]); % delays: short, medium, long

for str = ["ortho", "circulant"]
    for itD = 1:size(Ds,1)
        delays = Ds(itD,1:N);
        g = diag(gainPerSample.^delays);

        % % Complete A B C D
        % See section V.A of the reference paper.
        % A is with unilossless matrix U and diagonal matrix Γ ...
        % ... In (54), U can be a unilossless triangular matrix, i.e., with
        % ... a diagonal of ones [6]. In Section VI, we revisit this structure
        % ... for series allpasses. In the following, we focus on the more
        % ... intricate case of orthogonal U .

        % TODO:
        % Need to try some other matrix than random orthogonal?

        % A = V(1:N,1:N);
        % Ag = A * g; % gives incorrect decay times. Seems can't be random orthogonal with delay-proportional gains
        % [b,c,d,PP] = completeOrthogonal(Ag, numIO, 'verbose', true);
        % B = b;
        % C = c;
        % D = 0; %d;

        % Non-complete A B C D (this produces the correct t60)
        numInput = 1;
        numOutput = N;
        switch str
            case "ortho"
                Ag = Aortho * g;
            case "circulant"
                Ag = Acirculant * g;
            otherwise
                error("Unsupported matrix case")
        end
        B = ones(N,numInput);
        % C = ones(numOutput,N);
        C = diag(ones(N,1));
        D = zeros(numOutput,numInput);

        % rir(:,itD) = dss2impz(irLen,delays,Ag,B,C,D);
        rir(:,:,itD) = dss2impz(irLen,delays,Ag,B,C,D);

        dir = sprintf('renders/FDNs/%s/', str);

        if ~exist(dir,'dir'), mkdir(dir); end
        audiowrite( ...
            sprintf('%stradeoff_N%d_D%d.wav',dir,N,round(mean(delays))), ...
            0.05*rir(:,:,itD),...0.05*rir(:,itD), ... %0.5*rir(:,itD), ...
            fs);
    end
end

%%
fn = '/Users/michaelmccrea/src/toolkits/fdnToolbox/Examples/renders/FDNs/ortho_complete/tradeoff_N16_';
[sig, fs] = audioread([fn 'D269.wav']);
% [sig, fs] = audioread([fn 'D897.wav']);
[sig, fs] = audioread([fn 'D2692.wav']);

edc = getEDC(sig(:,2))
t60 = getRT60FromEDC(edc, fs)

%% Test: complete diagonally similar to orthogonal
N = 3;
numIO = 1;
X1 = blkdiag( diag(rand(N,1)), 1);
V = randomOrthogonal(N+numIO);
XVX = X1 \ V * X1;
XAX = XVX(1:N,1:N);

% [b,c,d,X,V] = completeAllpassFDN(XAX, 'verbose', true);
[b,c,d,X,V] = completeAllpassFDN(XAX, 'verbose', false);

V*V'
assert(isAlmostZero(V*V' - eye(size(V)),'tol',10^-6))

%% Test: complete series allpass
N = 4;
g = rescale(rand(N,1),0.5,0.99);

[A, b, c, d] = seriesAllpass(g);

[b,c,d,X,V] = completeAllpassFDN(A, 'verbose', true);

V*V'
assert(isAlmostZero(V*V' - eye(size(V)),'tol',10^-6))

%% Test: nested allpass - TODO sometimes fails due to poor low rank solution
N = 3;
g = rescale(rand(N,1),0.5,0.99);
[A, b, c, d] = nestedAllpass(g);

[b,c,d,X,V] = completeAllpassFDN(A, 'verbose', true);

V*V'
assert(isAlmostZero(V*V' - eye(size(V)),'tol',10^-6))

%% Test: homogeneous allpass
delays = [32 19 13];
g = 0.99;
G = diag( g.^delays )
X = -diag([0.4, 0.6, .85])
[A, b, c, d] = homogeneousAllpassFDN(G, X, 'verbose', true);
[isA, X] = isUniallpass(A, b, c, d);

[b,c,d,X,V] = completeAllpassFDN(A, 'verbose', true);

V*V'
assert(isAlmostZero(V*V' - eye(size(V)),'tol',10^-6))

