classdef PulseSequenceTest < matlab.unittest.TestCase

    methods(TestClassSetup)
        % Shared setup for the entire test class
        function setup(testCase)
            addpath(genpath('./'))
        end
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods
        function sliceSelectTest(testCase)
            %% Setup object and sequence for slice-selective excitation test
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            T1 = 600*T1; T2 = 100*T2;
            [~, B1, grad, ~] = seq.event_waveforms(1, 2, dt);

            % Symmetric splitting simulation
            [Msim, ~] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0, B1map=B1map, Meq=M0);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('Slice select test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, Mop) + b;

            % Compare
            testCase.verifyLessThan(mean( (Msim(:) - Mop(:)).^2 ), 1e-6);
        end

        function ssADCTest(testCase)
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            T1 = 600*T1; T2 = 100*T2;
            seq.ADC = PulseSequence.adc_specs(2000+150, 8, 800);
            [~, B1, grad, sampleComb] = seq.event_waveforms(1, seq.numEvents, dt);

            % Symmetric splitting simulation
            [Msim, Msamp] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0, B1map=B1map, Meq=M0, sampleComb=sampleComb);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('Slice select test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            Msampop = zeros(size(Msamp));
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, Mop) + b;
            [A, b] = seq.get_event_operator(3, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            for s = 1:8
                Msampop(:,s,:) = pagemtimes(A{s}, Mop) + b{s};
            end
            Mop = pagemtimes(A{end}, Mop) + b{end};
            [A, b] = seq.get_event_operator(4, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, Mop) + b;

            % Compare
            testCase.verifyLessThan(mean( (Msamp(:) - Msampop(:)).^2 ), 1e-6);
            testCase.verifyLessThan(mean( (Msim(:) - Mop(:)).^2 ), 1e-6);
        end

        function ssDeltaTest(testCase)
            %% Setup object and sequence for slice-selective excitation test
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            T1 = 600*T1; T2 = 100*T2;
            [~, B1, grad, ~] = seq.event_waveforms(1, 2, dt);
            delta = reshape(linspace(-50, 50, length(delta)), size(delta));

            % Symmetric splitting simulation
            [Msim, ~] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0, B1map=B1map, Meq=M0);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('SS delta test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, Mop) + b;

            [measop, Mzop] = split_magnetization(Mop);

            % Compare
            testCase.verifyLessThan(mean( ((Mz(:) - Mzop(:)).^2 + (abs(meas(:)) - abs(measop(:))).^2)/2 ), 1e-6);
        end

        function ssT1Test(testCase)
            %% Setup object and sequence for slice-selective excitation test
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            [~, B1, grad, ~] = seq.event_waveforms(1, 2, dt);
            T1 = reshape(linspace(1, 10, length(T1)), size(T1)); T2 = 100*T2;

            % Symmetric splitting simulation
            [Msim, ~] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0, B1map=B1map, Meq=M0);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('SS T1 test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, Mop) + b;

            % Compare
            testCase.verifyLessThan(mean( (Msim(:) - Mop(:)).^2 ), 1e-6);
        end

        function ssT2Test(testCase)
            %% Setup object and sequence for slice-selective excitation test
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            [~, B1, grad, ~] = seq.event_waveforms(1, 2, dt);
            T1 = 600*T1; T2 = reshape(linspace(1, 10, length(T2)), size(T2));

            % Symmetric splitting simulation
            [Msim, ~] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0, B1map=B1map, Meq=M0);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('SS T2 test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, Mop) + b;

            % Compare
            testCase.verifyLessThan(mean( (Msim(:) - Mop(:)).^2 ), 1e-6);
        end

        function ssB0mapTest(testCase)
            %% Setup object and sequence for slice-selective excitation test
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            T1 = 600*T1; T2 = 100*T2;
            [~, B1, grad, ~] = seq.event_waveforms(1, 2, dt);
            dB0 = reshape(linspace(0, 100, length(dB0)), size(dB0));

            % Symmetric splitting simulation
            [Msim, ~] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0*1e-6, B1map=B1map, Meq=M0);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('SS B0 test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0*1e-6, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0*1e-6, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, Mop) + b;

            [measop, Mzop] = split_magnetization(Mop);

            % Compare
            testCase.verifyLessThan(mean( ((Mz(:) - Mzop(:)).^2 + (abs(meas(:)) - abs(measop(:))).^2)/2 ), 1e-6);
        end

        function ssB1mapTest(testCase)
            %% Setup object and sequence for slice-selective excitation test
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            T1 = 600*T1; T2 = 100*T2;
            [~, B1, grad, ~] = seq.event_waveforms(1, 2, dt);
            B1map = B1map.*reshape(linspace(0, 2, length(B1map)), size(B1map));

            % Symmetric splitting simulation
            [Msim, ~] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0*1e-6, B1map=B1map, Meq=M0);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('SS B1 map test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, Mop) + b;

            % Compare
            testCase.verifyLessThan(mean( (Msim(:) - Mop(:)).^2 ), 1e-6);
        end

        function ssM0Test(testCase)
            %% Setup object and sequence for slice-selective excitation test
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            T1 = 600*T1; T2 = 100*T2;
            [~, B1, grad, ~] = seq.event_waveforms(1, 2, dt);
            M0 = reshape(linspace(0.5, 1.5, length(M0)), size(M0));

            % Symmetric splitting simulation
            [Msim, ~] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0*1e-6, B1map=B1map, Meq=M0);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('SS M0 test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map);
            Mop = pagemtimes(A, Mop) + b;
            % Compare
            testCase.verifyLessThan(mean( (Msim(:) - Mop(:)).^2 ), 1e-6);
        end

        function ssB0Test(testCase)
            %% Setup object and sequence for slice-selective excitation test
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            T1 = 600*T1; T2 = 100*T2;
            [~, B1, grad, ~] = seq.event_waveforms(1, 2, dt);
            delta = delta+5;

            % Symmetric splitting simulation
            [Msim, ~] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0*1e-6, B1map=B1map, Meq=M0, B0=2.9);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('SS B0 test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map, 'B0', 2.9);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map, 'B0', 2.9);
            Mop = pagemtimes(A, Mop) + b;

            [measop, Mzop] = split_magnetization(Mop);

            % Compare
            testCase.verifyLessThan(mean( ((Mz(:) - Mzop(:)).^2 + (abs(meas(:)) - abs(measop(:))).^2)/2 ), 1e-6);
        end

        function ssGammaTest(testCase)
            %% Setup object and sequence for slice-selective excitation test
            [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup();
            T1 = 600*T1; T2 = 100*T2;
            [~, B1, grad, ~] = seq.event_waveforms(1, 2, dt);

            % Symmetric splitting simulation
            [Msim, ~] = bloch_symmetric_splitting(dt*1e-6, B1*1e-3, grad*1e-3, pos*1e-3, T1*1e-3, T2*1e-3, delta=delta, B0map=dB0*1e-6, B1map=B1map, Meq=M0, gamma=200e6);

            figure();
            [meas, Mz] = split_magnetization(Msim);
            plot_magnetization(pos(3,:),meas,Mz);
            sgtitle('SS gamma test simulation'); xlabel('Z-position (mm)')

            % Compute and apply operator
            [A, b] = seq.get_event_operator(1, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map, 'gamma', 200e6);
            Mop = pagemtimes(A, [zeros(size(T1)); zeros(size(T1)); M0]) + b;
            [A, b] = seq.get_event_operator(2, 'pos', pos*1e-3, 'T1', T1*1e-3, 'T2', T2*1e-3, 'B0map', dB0, 'delta', delta, 'B1map', B1map, 'gamma', 200e6);
            Mop = pagemtimes(A, Mop) + b;

            % Compare
            testCase.verifyLessThan(mean( (Msim(:) - Mop(:)).^2 ), 1e-6);
        end
    end

end

function [dt, seq, pos, delta, T1, T2, dB0, B1map, M0] = sssetup(Gz, dz, FA, T, dt)
arguments
    Gz = 30; % mT/m
    dz = 5; % mm
    FA = 90; % deg
    T = 2000; % us
    dt = 10; % us
end
%% Construct pulse sequence

% Make pulse sequence
seq = PulseSequence(PulseSequence.rf_specs(0, {'sinc'}, FA, T, [], dz), ...
    [], [], ...
    PulseSequence.gradient_specs([0; T], ...
    [Gz; -Gz], ...
    [150; 150], ...
    [T-150; 950], ...
    [150; 150]), []);

% Construct phantom and assign properties

% Position, mm
z = permute(-5:0.1:5, [1, 3, 2]);
pos = [zeros(size(z)); zeros(size(z)); z];
% Chemical shift, ppm
delta = zeros(size(z));
% Single-compartment T1 relaxation time, ms
T1 = ones(size(z));
% Single-compartment T2 relaxation time, ms
T2 = ones(size(z));
% B0 drift / variation
dB0 = zeros(size(z));
% B1 map
B1map = ones(size(z));
% Equilibrium magnetization
M0 = ones(size(z));

end