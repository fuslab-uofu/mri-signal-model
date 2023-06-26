%% Unit Tests that check that off resonance inputs work as expected for the function bloch_symmetric_splitting
% 
%% 2023-06-26 Addison Powell

classdef off_resonance_test < matlab.unittest.TestCase
    properties 
        reference
    end

    methods (Test)
        function diff_delta(testCase)
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
            Min = [0; 0; 1];
            [B1e, ~] = b1_hardpulse(90, T, dt, t_units='s', B1_units='T');

            for ppm=1:20
                expected = (1+ppm*3e-6)*[1 0 0; 0 cos(-pi/2) -sin(-pi/2); 0 sin(-pi/2) cos(-pi/2)]*[0; 0; 1]; %analytical solution
    
                [result, ~] = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2, delta=ppm, Minit=[0; 0; 1]); %computed solution
                testCase.verifyLessThan(abs(result - expected), 1e-10);
            end
            
            
        end
    end
end