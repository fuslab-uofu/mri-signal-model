%% Convergence test - coherence pathways with symmetric splitting
% Simulates MPME sequence with 4 pathways using the symmetric splitting
% bloch solver
%
% To simulate phase coherence pathways, spins are simulated with a range of
% offsets along the unbalanced gradient direction (x). Once simulated, all
% the spins along that axis for a particular voxel are added up to give the
% observed magnetization.
%
% T1=200 ms, T2=T2*=100 ms, FA \in [1, 45], Nx \in [16, 128]
%
% 2023-06-07 Samuel Adams-Tew

[file, pth] = uiputfile('*.mat', 'Select a location to save test files', 'mpme_fixedT1.mat');
savePath = fullfile(pth, file(1:end-4));

Nxs = [16, 32, 64, 128];

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

FAs = logspace(0,log10(45),35);

%% Make range of x values for gradient-induced pathway splitting

[tmp, ~] = seq.gradient_moment('cumulative'); % in T/m
Gx_moment = tmp(1, end);

% Determine how far out we have to go to get 1 cycle of dephasing over 1 TR
% xmax = 2*pi/( gamma*Gx_moment )
% [m] = [rad]/([rad/(s*T)] * [T/m])
TR = seq.blockDuration*1e-6; % Repetition time [s]
gamma = 267.52218744e6; % [rad/(s*T)]
xmax = 2*pi/( gamma*Gx_moment*TR ); % [m]

for xidx = 1:length(Nxs)
    % Establish position based on number of samples
    Nx = Nxs(xidx); % Number of samples along x
    dx = xmax/Nx; % Sample spacing
    x = dx*(-Nx/2:Nx/2-1);
    y = zeros(size(x));
    z = zeros(size(x));
    
    T1map = 200;
    T2map = 100;
    [T1map, T2map] = meshgrid(T1map, T2map);
    
    pos = ones([3, 1, size(T2map), Nx]);
    pos = permute([x; y; z], [1, 3:5, 2]).*pos;
    
    % Make maps match
    T1map = permute(repmat(T1map, 1, 1, Nx), [4, 5, 1:3]);
    T2map = permute(repmat(T2map, 1, 1, Nx), [4, 5, 1:3]);
    
    pb = ProgressBar(sprintf('Nx=%d MPME flip angle', Nx), length(FAs)); 
    for idx = 1:length(FAs)
        pb.iter();
        seq.RF.FA(1) = FAs(idx); % Update the flip angle
        tmppath = sprintf("%s_FA=%.1f_Nx=%d.mat", savePath, FAs(idx), Nx);
        
        simulate_acquisition(T1map, T2map, pos, seq, 'numRepetitions', 100, 'savePath', tmppath, 'pos_units', 'm');
    end
    pb.close();
    
    %% Compile into single file
    
    file = load(sprintf("%s_FA=%.1f_Nx=%d.mat", savePath, FAs(idx), Nx));
    sz = size(file.magnetization);
    
    im = complex(zeros(sz(5), length(FAs), 4, 3)); % (N_T1, N_FA, N_pwy, N_echo)
    T1 = squeeze(file.T1map(1, 1, 1, :, 1));
    T2 = squeeze(file.T2map(1, 1, 1, :, 1));
    
    for idx = 1:length(FAs)
        tmppath = sprintf("%s_FA=%.1f_Nx=%d.mat", savePath, FAs(idx), Nx);
        
        file = load(tmppath);

        % Add up the transverse magnetization along the unbalanced
        % dimension to get multi-pathway signals
        [meas, ~] = split_magnetization(file.magnetization);
        meas = sum(meas, 5)./Nx;
        meas = permute(meas, [4, 1, 2, 3]);
        
        im(:, idx, :, :) = reshape(meas, [sz(5), 1, 4, 3]);
    end
    
    save(sprintf("%s_Nx=%d_compiled.mat", savePath, Nx), 'T1', 'T2', 'FAs', 'im');

    % Clean up temporary files
    for idx = 1:length(FAs)
        delete(sprintf("%s_FA=%.1f_Nx=%d.mat", savePath, FAs(idx), Nx));
    end
end

%% Show convergence test results

% Plot most accurate one first, show first echo
Nxs = flip(Nxs); ech = 1;
leg = {};
for xidx = 1:length(Nxs)
    Nx = Nxs(xidx);
    file = load(sprintf("%s_Nx=%d_compiled.mat", savePath, Nx));

    for pwyIdx = 1:4
        figure(pwyIdx + 1)
        plot(file.FAs, abs(file.im(1, :, pwyIdx, ech))); hold on
    end
    leg = [leg, sprintf("Nx=%d", Nx)];
end

figure(2); title('|F_{+1}|'); xlim([1, 45]);
figure(3); title('|F_{0}|'); xlim([1, 45]);
figure(4); title('|F_{-1}|'); xlim([1, 45]);
figure(5); title('|F_{-2}|'); xlim([1, 45]);

for pwyIdx = 1:4
    figure(pwyIdx + 1);
    xlabel('FA (\circ)'); ylabel('Signal')
    grid on; hold off
    legend(leg, 'Location','southeast')
    saveas(gcf, sprintf(sprintf("%s_pwyIdx=%d.png", savePath, pwyIdx)))
end