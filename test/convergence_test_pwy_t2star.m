close all

[file, pth] = uiputfile('*.mat', 'Select a location to save test files', 'mpme_pwy+t2s.mat');
savePath = fullfile(pth, file(1:end-4));

% Testing parameters
Nxs = 128;
Nfs = [511, 512, 1023, 1024];
T2p = 200; % T2' [ms]

TR = 20; % Repetition time [ms]

T1 = 200; % [ms]
T2s = 100; % [ms] Observed transverse relaxation time

% 10*T2/TR repetitions to reach steady state
% Add one, readout after steady state
nRep = 200; % ceil(10*T2/TR) + 1;
Tmax = nRep*TR;
T2 = 1/(1/T2s - 1/T2p);

coeff = Tmax/T2s - log(1e-3); % Replica error limiting coefficient
df = min(1/(2*Tmax), 1/(coeff*T2s))/1e-3;

% Constants
B0 = 3; % [T]
gamma = 267.52218744e6; % [rad/(s*T)]

seq = PulseSequence(PulseSequence.rf_specs(0, {'hard'}, 15, 320, 0), ...
    PulseSequence.gradient_specs( ... % Gx
    [   320;  880;   4760; 5680;   9560; 10480;  14360], ...
    -[-34.56; 6.95; -33.64; 6.95; -33.64;  6.95; -34.58], ...
    [   280;  250;    170;  250;    170;   250;    280], ...
    [   280; 3630;    750; 3630;    750;  3630;    280], ...
    [   280;  250;    170;  250;    170;   250;    280]), ...
    [], [], ... % Gy, Gz
    PulseSequence.adc_specs([1130; 5930; 10730], [4; 4; 4], [3380; 3380; 3380]), ...
    TR*1e3);

% Set readout times based on gradient timing (when particular pathways
% become coherent)
[ts, rogs] = seq.pathway_readout_times(-2:1, 1);
seq.ADC.sampleTime = ts - rogs(:, 1);

seq.plot()

FAs = logspace(0,log10(45),20);%35);

% Configure phase coherence pathways
[tmp, ~] = seq.gradient_moment('cumulative'); % in T/m
Gx_moment = tmp(1, end);

for fidx = 1:length(Nfs)
    Nf = Nfs(fidx);
    % Create spectral line with specified number of samples 
    f = df*(-(Nf-1)/2:(Nf-1)/2); % Evenly spaced, centered about zero
    f = permute(f, [1, 3, 4, 2]); % delta varies along dim4
    weight = spectral_line_shape(f, 'T2prime', T2p, 'spectAxis', 4);
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
        
        pb = ProgressBar(sprintf('Nf=%d, Nx=%d MPME flip angle', Nf, Nx), length(FAs));
        for idx = 1:length(FAs)
            pb.iter();
            seq.RF.FA(1) = FAs(idx); % Update the flip angle
            tmppath = sprintf("%s_FA=%.1f_Nf=%d_Nx=%d.mat", savePath, FAs(idx), Nf, Nx);
            
            simulate_acquisition(T1tmp, T2tmp, pos, seq, 'delta', deltatmp, 'numRepetitions', nRep, 'savePath', tmppath, 'pos_units', 'm');
            
            % Load simulation data and combine dimensions to give T2* and
            % phase coherence pathway signals
            file = load(tmppath);
            [meas, Mz] = split_magnetization(file.magnetization);
            meas = sum(sum(meas.*weight, 4), 5)./Nx;
            Mz = sum(sum(Mz.*weight, 4), 5)./Nx;
            save(tmppath, 'meas', 'Mz', 'Nx', 'Nf', 'weight', '-append');
        end
        pb.close();
        
        % Compile into single file
        file = load(sprintf("%s_FA=%.1f_Nf=%d_Nx=%d.mat", savePath, FAs(idx), Nf, Nx));
        sz = size(file.meas);
        
        im = complex(zeros(sz(2), length(FAs), 4, 3));

        for idx = 1:length(FAs)
            file = load(sprintf("%s_FA=%.1f_Nf=%d_Nx=%d.mat", savePath, FAs(idx), Nf, Nx));
            meas = permute(file.meas, [2, 3, 1]);
            im(:, idx, :, :) = reshape(meas, [sz(2), 1, 4, 3]);
        end
        
        save(sprintf("%s_Nf=%d_Nx=%d_compiled.mat", savePath, Nf, Nx), 'im', 'FAs', 'Nf', 'Nx', 'T1', 'T2', 'T2s', 'T2p');
        
        % Clean up temporary files
        for idx = 1:length(FAs)
            delete(sprintf("%s_FA=%.1f_Nf=%d_Nx=%d.mat", savePath, FAs(idx), Nf, Nx));
        end
    end
end

%%
Nfs = flip([127, 128, 255, 256, 511, 512, 1023, 1024]);
ech = 1; leg = {};
for fidx = 1:length(Nfs)
    Nf = Nfs(fidx);
    for xidx = 1:length(Nxs)
        Nx = Nxs(xidx);
        file = load(sprintf("%s_Nf=%d_Nx=%d_compiled.mat", savePath, Nf, Nx));
        
        for pwyIdx = 1:4
            figure(pwyIdx + 1)
            plot(file.FAs, abs(file.im(1, :, pwyIdx, ech))); hold on
        end
        leg = [leg, sprintf("Nf=%d, Nx=%d", Nf, Nx)];
    end
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