clear, close all
addpath(genpath('./'))
%% Construct pulse sequence
Gz = 30; % mT/m
dz = 5; % mm
FA = 90; % deg
T = 2; % ms
dt = 4e-3; % ms
samplePeriod = 4e-3; % ms

% RF field, mT 
[B1, t] = b1_sliceselect(Gz, dz, FA, T, dt);%, 'pulse', 'hard', 'filter', 'none');
% Gradient, mT/m
grad = gradient_trap(30, 0.150, T-0.150, 0.150, dt)';

% Add rephasing gradient
tmp = gradient_trap(-30, 0.150, 0.95, 0.15, dt)';
B1 = [B1, zeros(size(tmp))];
grad = [grad, tmp];
t = [t, ( T + (0:(length(tmp)-1))*dt )];

% Make other gradient directions zero magnitude
grad = [zeros(2, length(grad)); grad];

figure()
plot_sequence(t, B1, grad);

%% Construct phantom and assign properties

% Position, mm
z = -5:0.1:5;
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
[Mfinal, Msample] = bloch_symmetric_splitting(dt*1e-3, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0, samplePeriod=samplePeriod*1e-3);

[meas, Mz] = split_magnetization(Mfinal);

figure()
plot_magnetization(z, meas, Mz);
xlabel('Z-position (mm)')

figure()
tmp = ( 0:(size(Msample, 2)-1) )*samplePeriod;
plot(tmp, Msample(1, :, round(length(z)/2)), tmp, Msample(2, :, round(length(z)/2)), tmp, Msample(3, :, round(length(z)/2)))
legend('M_x', 'M_y', 'M_z'); grid on
xlabel('Time (ms)'); ylabel('Magnetization (fracton of M_0)')
