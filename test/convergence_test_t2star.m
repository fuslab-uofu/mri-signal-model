%% Convergence test - T2*
% Simulates spin echo sequence with T2* using the symmetric splitting bloch
% solver
%
% To simulate T2*, create a set of isochromats with varying levels of off
% resonance and add them together with Lorentzian weighting.
%
% Note that I did test whether an even or odd number of isochromats was
% required (set Nf to include both even and odd values), and observed no
% real difference in performance between the two.
%
% 2023-06-08 Samuel Adams-Tew
clear, close all

% Set T2* testing parameters
dt = 10e-3; % [ms]
% List containing how many frequency offsets to use in each test. This
% determines how finely sampled your frequency space is. Larger values are
% more accurate, and generally, larger number of offsets means the
% simulation is accurate for a longer time following RF pulse.
Nf = 2.^(3:8); 
T2p = 50; % [ms]
% Maximum off-resonance for any isochromat. The larger fmax, the closer the
% decay is to exponenetial
fmax = 40; % [Hz] 

% Configure spin echo testing sequence
TE = 30;
pulseDur = 1; % ms
% 90 deg
[B1, t] = b1_hardpulse(90, pulseDur, dt);
% Wait TE/2
tmp = t(end) + (0:dt:(TE/2 - pulseDur));
B1 = [B1, zeros(size(tmp))];
t = [t, tmp];
% 180 deg pulse
[B1tmp, tmp] = b1_hardpulse(180, pulseDur, dt);
B1 = [B1, B1tmp];
t = [t, tmp + t(end) + dt];
% Wait TE
tmp = t(end) + (0:dt:(TE - pulseDur/2));
B1 = [B1, zeros(size(tmp))];
t = [t, tmp];
% No gradients
grad = zeros(3, length(B1));

% Set constants
gamma = 267.52218744e6; % [rad/(s*T)]
B0 = 3; % [T]

% Set sample parameters
T2 = permute(20:10:200, [1, 3, 2]);
T1 = 200*ones(size(T2));
pos = zeros(3, 1, length(T2));

leg = {};
for iter = 1:length(Nf)
    % Create spectral line with specified number of samples 
    f = permute(linspace(-fmax, fmax, Nf(iter)), [1, 3, 4, 2]);
    weight = spectral_line_shape(f, 'T2prime', T2p);
    % Convert f to delta in ppm
    delta = 2*pi*f./(gamma*B0)*1e6;

    % Create parameter maps for this simulation
    postmp = pos.*ones(size(delta));
    T1tmp = T1.*ones(size(delta));
    T2tmp = T2.*ones(size(delta));
    delta = delta.*ones(size(T2));

    [Mfinal, Msignal] = bloch_symmetric_splitting(dt*1e-3, B1*1e-3, grad*1e-3, postmp*1e-3, T1tmp*1e-3, T2tmp*1e-3, delta=delta, samplePeriod=10e-6);
    % Add up isochromats with appropriate weighting to get T2* decay
    Mfinal = sum(weight.*Mfinal, 4); Msignal = sum(weight.*Msignal, 4);
    [meastmp, Mztmp] = split_magnetization([Msignal, Mfinal]);
    
    t = 10e-3*((1:length(meastmp)) - 1);
    % Time evolution
    figure(1)
    plot(t, abs(meastmp(:, T2==100))); hold on

    % Signal vs T2
    figure(2)
    plot(squeeze(T2), abs(meastmp(end, :))); hold on

    leg = [leg sprintf('Nf=%d',Nf(iter))];
end

figure(1);
hold off; grid on
legend(leg)
ylabel('Signal magnitude')
xlabel('Time (ms)')
xlim([0, max(t)]);

figure(2)
hold off; grid on
legend(leg, 'Location', 'northwest')
ylabel('Signal magnitude')
xlabel('T2 (ms)')
