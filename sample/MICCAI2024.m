clear, close all
cd('/v/raid10/users/sadams/2023/2023-07-12_SimulatedMPME_NN');
addpath(genpath('./'))
addpath(genpath('/v/raid10/users/sadams/Resources/Code/BlochSim/'))

datasetDir = './SimDatasets/';
scratchDir = '/working/sadams/';
simlabel = 'mpmesim_FA=1-45';

datestr = char(datetime('today', 'Format', 'yyyyMMdd'));

%% Parameters
disp("Setting parameters...")
Nsamp = 1024; % Number of samples to simulate in each iteration
Niter = 256; % Number of iterations to simulate
% Output has shape (Npwy, Necho, Nsamp, Niter)

FArng = [5, 15];
TRrng = [20, 20];

B1maprng = [.5, 1.5]; % scale factor
B0maprng = [-3, 3]; % uT

% Breast fat: T1 367+/-8 ms, T2 53.0+/-1.5 ms
% Breast gland: T1 1440+/-90 ms, T2 54+/-9 ms
T1rng = [200, 1800]; % T1 range
T2rng = [25, 250]; % T2 range [25, min(T1, 250)]
% Range of T2', log scale (overweights short T2')
T2prng = log10([1, 10000]);

% Simulation parameters
Nx = 128;
Nf = 128;
nRep = 250;

% Constants
B0 = 3; % [T] Siemens spec is 2.89==123e6/(267.52218744e6/(2*pi))
gamma = 267.52218744e6; % [rad/(s*T)]

%% Define pulse sequence
% For all definitions, times are in [us]
disp('Initializing sequence...')
seq = get_mpme_seq(10, 20);
seq.plot();

%% Iterate with random configurations

if isempty(gcp('nocreate'))
    parpool('local')
end

% Setup dataset file
datasetName = [datestr, '_', simlabel, '.mat'];
fprintf('Creating dataset %s\n', datasetName)

if exist(fullfile(datasetDir, datasetName), 'file')
    warning('File already exists. Attempting to resume simulation.')
    dfile = matfile(fullfile(datasetDir, datasetName), 'Writable', true);
else
    dfile = matfile(fullfile(datasetDir, datasetName), 'Writable', true);
    dfile.dimLabels = {'pwy', 'echo', 'Nsamp', 'Niter'};
    dfile.T1(Nsamp, Niter) = single(0);
    dfile.T2(Nsamp, Niter) = single(0);
    dfile.T2s(Nsamp, Niter) = single(0);
    dfile.FA(1, Niter) = single(0);
    dfile.TR(1, Niter) = single(0);
    dfile.B0map(Nsamp, Niter) = single(0);
    dfile.B1map(Nsamp, Niter) = single(0);
    dfile.pwy = [1, 0, -1, -2];
    dfile.seq = get_mpme_seq(max(FArng), max(TRrng));
    dfile.Mxy(4, 3, Nsamp, Niter) = complex(0);
    dfile.Mz(4, 3, Nsamp, Niter) = 0;
    dfile.simTime(1, Niter) = 0;
    % Helper properties for partial computation
    dfile.iter_saved = 0;
end

% Simulation inputs have shape 
% (:, :, Nsamp, Nf, Nx)
disp('Starting simulations...')
pb = ProgressBar('Generating dataset', Niter, 'iterBefore', false);

dq = parallel.pool.DataQueue;
afterEach(dq, @(varargin) pb.iter());

% Number of iterations to compute in parallel before saving them to the
% dataset and clearing the scratch directory
maxIterAtOnce = 64;
iterRemain = Niter - dfile.iter_saved;

while iterRemain > 0
    cIter = min(iterRemain, maxIterAtOnce);
    fprintf('%d simulations remain. Starting cycle of %d simulations.\n', iterRemain, cIter)
    parfor iter = 1:cIter
        % Repetition time [ms]
        TR = randi(round(TRrng), 1);
        FA = randi(round(FArng), 1);
        seq = get_mpme_seq(FA, TR);
        Tmax = nRep*TR;

        %%% Randomly set parameters, vary along dim3
        T1 = rand(1, 1, Nsamp) * (max(T1rng) - min(T1rng)) + min(T1rng); % [ms]
        T2 = rand(1, 1, Nsamp) .* (min(T1, max(T2rng)) - min(T2rng)) + min(T2rng); % [ms]
        T2p = 10.^(rand(1, 1, Nsamp) * (max(T2prng) - min(T2prng)) + min(T2prng)); % [ms]
        B1map = rand(1, 1, Nsamp) * (max(B1maprng) - min(B1maprng)) + min(B1maprng); % [scale factor]
        B0map = rand(1, 1, Nsamp) * (max(B0maprng) - min(B0maprng)) + min(B0maprng); % [uT]

        %%% Set up T2*
        T2s = 1./(1./T2 + 1./T2p);

        coeff = Tmax./T2s - log(1e-3); % Replica error limiting coefficient
        df = min(1/(2*Tmax), 1./(coeff.*T2s))/1e-3;

        % Create spectral line with specified number of samples
        % Because we are assuming symmetric spectrum, we simulate only
        % positive off-resonance values
        f = permute(0:(Nf-1), [1, 4, 3, 2]); % delta varies along dim4
        f = df.*f; % Evenly spaced, starting at zero
        % Convert f to delta in ppm
        delta = 2*pi*f./(gamma*B0)*1e6;

        %%% Set up configuration states

        % Configure phase coherence pathways
        [tmp, ~] = seq.gradient_moment('cumulative'); % in T/m
        Gx_moment = tmp(1, end);

        xmax = 2*pi./( gamma*(1 + delta*1e-6)*Gx_moment*TR*1e-3 ); % [m]
        % Establish position based on number of samples
        dx = xmax/Nx; % Sample spacing
        x = permute((-Nx/2:Nx/2-1), [1, 5, 3, 4, 2]); % x varies along dim5
        x = dx.*x;
        pos = [x; zeros(size(x)); zeros(size(x))];

        %%% Prepare and run simulation
        saveName = sprintf("%s_%s_iter=%d.mat", simlabel, datestr, iter);

        % Make map dimensions match
        T1tmp = T1.*ones(size(delta)).*ones(size(x));
        T2tmp = T2.*ones(size(T1tmp));
        B1tmp = B1map.*ones(size(T1tmp));
        B0tmp = B0map.*ones(size(T1tmp));
        delta = delta.*ones(size(T1tmp));
        pos = pos.*ones(size(T1tmp));

        simulate_acquisition(T1tmp, T2tmp, pos, seq, 'delta', delta, ...
            'B1map', B1tmp, 'B0map', B0tmp, 'numRepetitions', nRep, ...
            'savePath', fullfile(scratchDir, saveName), 'pos_units', 'm', ...
            'showProgressBar', false);

        % Load simulation data and combine dimensions to give T2* and configuration state signals
        sfile = matfile(fullfile(scratchDir, saveName), 'Writable', true);
        sfile.T1 = T1; sfile.T2 = T2; sfile.FA = FA; sfile.T2p = T2p;
        sfile.B1map = B1map; sfile.B0map = B0map; sfile.T2s = T2s;
        sfile.TR = TR;

        % Get Lorentzian weighting
        sz = size(f); sz(4) = Nf*2 - 1;
        f_full = zeros(sz);
        f_full(:,:,:,1:Nf) = -flip(f,4);
        f_full(:,:,:,Nf:end) = f;
        % Single-sided distribution
        weight = spectral_line_shape(f, 'T2prime', T2p, 'spectAxis', 4);
        % 'wings' are underweighted for single-sided distribution
        weight(:,:,:,2:end) = 2*weight(:,:,:,2:end);
        % Double-sided symmetric distribution
        weight_full = spectral_line_shape(f_full, 'T2prime', T2p, 'spectAxis', 4);

        % Get longitudinal and transverse components
        [Mxy, Mz] = split_magnetization(sfile.magnetization);

        % Fill full data from the simulated half of the spectrum
        % (T, 1, T2p, 2*Nf-1, Nx)
        % Allocate space needed
        sz = size(Mxy); sz(4) = 2*Nf - 1;
        Mxy_full = zeros(sz);
        % Get on-resonance magnetization
        Mxy_ctr = Mxy(:,:,:,1,:);
        % Get phasor of rotation angle for that magnetization
        phrot = exp(1i*angle(Mxy_ctr));
        % Fill out unsampled half of spectrum
        Mxy_full(:,:,:,1:Nf,:) = phrot.*conj(flip(Mxy, 4).*conj(phrot));
        % Fill out sampled half of spectrum
        Mxy_full(:,:,:,Nf:end,:) = Mxy;

        % Add up spectrum
        Mxy = sum(sum(Mxy_full.*weight_full, 4), 5)./Nx;
        Mz = sum(sum(Mz.*weight, 4), 5)./Nx;

        % Reshape to give pathways and echoes
        Mxy = reshape(permute(Mxy, [3, 1, 2]), [], 4, 3);
        Mz = reshape(permute(Mz, [3, 1, 2]), [], 4, 3);

        sfile.Mxy = Mxy;
        sfile.Mz = Mz;

        send(dq, iter)
    end

    for iter = 1:cIter
        saveName = sprintf("%s_%s_iter=%d.mat", simlabel, datestr, iter);
        sfile = matfile(fullfile(scratchDir, saveName), 'Writable', false);

        cidx = dfile.iter_saved + 1;

        %%% Compile dataset
        % dimlabels {'pwy', 'echo', 'Nsamp', 'Niter'};
        dfile.Mxy(:,:,:,cidx) = permute(sfile.Mxy, [2, 3, 1]);
        dfile.Mz(:,:,:,cidx) = permute(sfile.Mz, [2, 3, 1]);
        dfile.T1(:,cidx) = single(squeeze(sfile.T1));
        dfile.T2(:,cidx) = single(squeeze(sfile.T2));
        dfile.T2s(:,cidx) = single(squeeze(sfile.T2s));
        dfile.FA(1,cidx) = single(sfile.FA);
        dfile.TR(1,cidx) = single(sfile.TR);
        dfile.B0map(:,cidx) = single(squeeze(sfile.B0map));
        dfile.B1map(:,cidx) = single(squeeze(sfile.B1map));
        dfile.simTime(1,cidx) = sfile.telapsed;

        % Delete this iteration
        clear sfile
        delete(fullfile(scratchDir, saveName));

        % Increment iteration counter
        dfile.iter_saved = cidx;
    end

    % Update number of iterations remaining
    iterRemain = iterRemain - cIter;
end
pb.close();

totalTime = sum(dfile.simTime);
avgTime = mean(dfile.simTime);

% Convert file to v7; remove helper variables
load(fullfile(datasetDir, datasetName));
save(fullfile(datasetDir, datasetName), 'avgTime', 'dimLabels', 'FA', 'Mxy', 'Mz', 'pwy', 'seq', 'T1', 'T2', 'T2s', 'B0map', 'B1map', 'totalTime', 'TR', '-v7')

function seq = get_mpme_seq(FA, TR)
seq = PulseSequence(PulseSequence.rf_specs(0, {'hard'}, FA, 320, 0), ...
    PulseSequence.gradient_specs( ... % Gx
    [   320;  880;   4760; 5680;   9560; 10480;  14360], ...
    -[-34.56; 6.95; -33.64; 6.95; -33.64;  6.95; -34.58], ...
    [   280;  250;    170;  250;    170;   250;    280], ...
    [   280; 3630;    750; 3630;    750;  3630;    280], ...
    [   280;  250;    170;  250;    170;   250;    280]), ...
    [], [], ... % Gy, Gz
    PulseSequence.adc_specs([1130; 5930; 10730], [4; 4; 4], [3379.2; 3379.2; 3379.2]), ...
    TR*1e3);

% Set readout times based on gradient timing (when particular pathways
% become coherent)
[ts, rogs] = seq.pathway_readout_times(-2:1, 1);
seq.ADC.sampleTime = ts - rogs(:, 1);
end
