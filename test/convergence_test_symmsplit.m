%% symmetric splitting convergence test

clear, close all

%%% Temporally varying
Gz = 30; % mT/m
dz = 5; % mm
FA = 90; % deg
T = 2; % ms
z = -5:0.02:5;

dts = [5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 1e-4];
meass = NaN(length(dts), length(z));
Mzs = NaN(length(dts), length(z));

for iter = 1:length(dts)
    dt = dts(iter);

    % RF field, mT
    [B1, t] = b1_sliceselect(Gz, dz, FA, T, dt);%, 'pulse', 'hard', 'filter', 'none');
    % Gradient, mT/m
    grad = gradient_trap(30, 0.150, T-0.150, 0.150, dt)';

    % Add rephasing gradient
    tmp = gradient_trap(-30, 0.150, 0.95, 0.15, dt)';
    B1 = [B1, zeros(size(tmp))];
    grad = [grad, tmp];
    t = [t, ( T + (0:(length(tmp)-1))*dt )];

    grad = [zeros(2, length(grad)); grad];

    if iter == length(dts)
        % Show sequence
        figure(); plot_sequence(t, B1, grad);
        xlabel('Time (ms)')
    end

    %% Construct phantom and assign properties

    % Position, mm
    pos = permute([zeros(size(z))
        zeros(size(z))
        z], [1, 3, 2]);
    % Chemical shift, ppm
    delta = 0;
    % Single-compartment T1 relaxation time, ms
    T1 = 600*ones(size(z));
    % Single-compartment T2 relaxation time, ms
    T2 = 100*ones(size(z));
    % B0 drift / variation
    dB0 = zeros(size(z));
    % Equilibrium magnetization
    M0 = ones(size(z));
    
    %% Simulate
    [Mfinal, ~] = bloch_symmetric_splitting(dt*1e-3, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0);

    [meas, Mz] = split_magnetization(Mfinal);

    figure(); plot_magnetization(z, meas, Mz)
    xlabel('Z-position (mm)')
    sgtitle(sprintf('dt=%g', dt*1e-3))

    meass(iter, :) = meas(end, :);
    Mzs(iter, :) = Mz(end, :);
end

rmse_meas = abs(sqrt(mean((meass.' - meas(end, :).').^2)))
rmse_Mz = sqrt(mean((Mzs.' - Mz(end, :).').^2))

figure()
plot(z, flipud(abs(meass)))
legend('100 ns', '1 us', '2 us', '5 us', '10 us', '20 us', '50 us')
grid on
xlabel('Position (mm)')
ylabel('Transverse magnitude')
