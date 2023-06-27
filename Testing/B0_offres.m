%% Unit Tests that check that B0 and gamma inputs work as expected for the function bloch_symmetric_splitting
% confirms that changes in B0 or gamma have the correct impact on
% simulations, also confirms that B0map creates the correct off resonance.
%% 2023-06-26 Addison Powell

classdef B0_offres < matlab.unittest.TestCase
    properties 
        reference
    end

    methods (Test)
        % function for testing different magnetic field strengths
        function diff_B0(testCase)
            addpath(genpath('./'))

            T = 2; % s
            z = -5:5;
            sz = size(z);
            dt = 1e-3; % ms
            pos = permute([zeros(size(z))
                    zeros(size(z))
                    z], [1, 3, 2]);
            grad = zeros(3, ceil(T/dt));
            T1 = inf*ones(sz);
            T2 = inf*ones(sz);
            [B1e, ~] = b1_hardpulse(90, T, dt, t_units='s', B1_units='T');

            for B0=5:20
                expected = [1 0 0; 0 cos(-pi/2) -sin(-pi/2); 0 sin(-pi/2) cos(-pi/2)]*[0; 0; 1]; %analytical solution
    
                [result, ~] = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2, B0=B0, delta=5); %computed solution
                testCase.verifyLessThan(abs(result - expected), 1e-10);
            end
        end

        %function for testisng different gyromagnetic ratios
        function diff_gamma(testCase)
            addpath(genpath('./'))

            T = 2;
            z = -5:5;
            sz = size(z);
            dt = 1e-3; % ms
            pos = permute([zeros(size(z))
                    zeros(size(z))
                    z], [1, 3, 2]);
            grad = zeros(3, ceil(T/dt));
            T1 = inf*ones(sz);
            T2 = inf*ones(sz);
            [B1e, ~] = b1_hardpulse(90, T, dt, t_units='s', B1_units='T');

            for gamma = 200:10:300
                expected =0;
                result = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2, gamma=gamma);

                testCase.verifyLessThan(abs(result - expected), 1e-10);
            end
        end

        % function for testing a magnetic field with varying strengths along the z-axis
        function B0_map(testCase)
            addpath(genpath('./'))

            T = 2;
            z = -5:5;
            sz = size(z);
            dt = 1e-3; % ms
            pos = permute([zeros(size(z))
                    zeros(size(z))
                    z], [1, 3, 2]);
            grad = zeros(3, ceil(T/dt));
            T1 = inf*ones(sz);
            T2 = inf*ones(sz);
            [B1e, ~] = b1_hardpulse(90, T, dt, t_units='s', B1_units='T');

            B0map = 2.5:.1:3.5;
            expected = 0;
            result = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2, B0map=B0map);

            testCase.verifyLessThan(abs(result - expected), 1e-10);

        end
    end
end