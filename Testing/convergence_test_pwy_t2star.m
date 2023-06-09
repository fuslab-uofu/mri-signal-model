clear, close all

[file, pth] = uiputfile('*.mat', 'Select a location to save test files', 'mpme_pwy+t2s.mat');
savePath = fullfile(pth, file(1:end-4));

% Testing parameters
Nxs = [64, 32];
Nfs = [32, 16, 64];
T2p = 200; % T2' [ms]
fmax = 25; % Maximum off resonance [Hz]
nRep = 50; % 100

% Constants
B0 = 3; % [T]
gamma = 267.52218744e6; % [rad/(s*T)]

% Parameters vary along dim3
T1 = permute(linspace(200, 1800, 25), [1, 3, 2]);
T2 = permute(linspace(50, 200, 25), [1, 3, 2]);

seq = PulseSequence(PulseSequence.rf_specs(0, {'hard'}, 15, 320, 0), ...
    PulseSequence.gradient_specs( ... % Gx
    [   320;  880;   4760; 5680;   9560; 10480;  14360], ...
    -[-34.56; 6.95; -33.64; 6.95; -33.64;  6.95; -34.58], ...
    [   280;  250;    170;  250;    170;   250;    280], ...
    [   280; 3630;    750; 3630;    750;  3630;    280], ...
    [   280;  250;    170;  250;    170;   250;    280]), ...
    [], [], ... % Gy, Gz
    PulseSequence.adc_specs([1130; 5930; 10730], [4; 4; 4], [3380; 3380; 3380]), ...
    20000);

% Set readout times based on gradient timing (when particular pathways
% become coherent)
[ts, rogs] = seq.pathway_readout_times(-2:1, 1);
seq.ADC.sampleTime = ts - rogs(:, 1);

seq.plot()

% Configure phase coherence pathways
[tmp, ~] = seq.gradient_moment('cumulative'); % in T/m
Gx_moment = tmp(1, end);

for fidx = 1:length(Nfs)
    Nf = Nfs(fidx);
    % Create spectral line with specified number of samples 
    f = permute(linspace(-fmax, fmax, Nf), [1, 3, 4, 2]); % delta varies along dim4
    weight = spectral_line_shape(f, 'T2prime', T2p);
    % Convert f to delta in ppm
    delta = 2*pi*f./(gamma*B0)*1e6;

    for xidx = 1:length(Nxs)
        Nx = Nxs(xidx);
        % Determine how far out we have to go to get 1 cycle of dephasing over 1 TR
        % xmax = 2*pi/( gamma*Gx_moment )
        % [m] = [rad]/([rad/(s*T)] * [T/m])
        TR = seq.blockDuration*1e-6; % Repetition time [s]
        xmax = 2*pi./( gamma*(1 + delta*1e-6)*Gx_moment*TR ); % [m]
        % Establish position based on number of samples
        dx = xmax/Nx; % Sample spacing
        x = permute(dx.*(-Nx/2:Nx/2-1), [1, 5, 3, 4, 2]); % x varies along dim5
        pos = [x; zeros(size(x)); zeros(size(x))];

        % Create temporary objects that match the required dimensions
        % Simulation object shape (3, 1, Nparam, Nf, Nx)
        T1tmp = T1.*ones(size(delta)).*ones(size(x));
        T2tmp = T2.*ones(size(T1tmp));
        deltatmp = delta.*ones(size(T1tmp));
        pos = pos.*ones(size(T2tmp));

        tmppath = char(sprintf("%s_Nf=%d_Nx=%d.mat", savePath, Nf, Nx));
        simulate_acquisition(T1tmp, T2tmp, pos, seq, 'numRepetitions', nRep, 'pos_units', 'm', 'delta', deltatmp, 'savePath', tmppath)
        
        % Load simulation data and combine dimensions to give T2* and
        % phase coherence pathway signals
        file = load(tmppath);
        [meas, Mz] = split_magnetization(file.magnetization);
        meas = sum(sum(meas.*weight, 4), 5)./Nx;
        Mz = sum(sum(Mz.*weight, 4), 5)./Nx;
        save(tmppath, 'meas', 'Mz', 'Nx', 'Nf', 'weight', '-append');
        
        % Create overivew figure
        fig = figure();
        xax =  squeeze(file.T1map(1,1,:,1,1)./file.T2map(1,1,:,1,1));
        plot(xax, 4*abs(squeeze(meas(4, 1, :)))); hold on
        plot(xax, 2*abs(squeeze(meas(3, 1, :)))); hold on
        plot(xax,   abs(squeeze(meas(2, 1, :)))); hold on
        plot(xax, 2*abs(squeeze(meas(1, 1, :)))); hold off
        legend('4\times|F_{-2}|', '2\times|F_{-1}|', '|F_{0}|', '2\times|F_{+1}|')
        grid on
        saveas(fig, [tmppath(1:end-4), '.png']);
    end
end