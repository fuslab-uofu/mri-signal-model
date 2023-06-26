%% Unit Test that compares simulated T2 decay to an analytical solution for the function bloch_symmetric_splitting
% test bloch_sym on 11 different T2 values with no B1 or gradient field
%% 2023-06-26 Addison Powell

classdef T2_decay < matlab.unittest.TestCase
    properties
    end

    methods (Test)
        function final_state_err(testCase)
            addpath(genpath('./'));

            T = 2; % s
            z = -5:5;
            sz = size(z);
            dt = 1e-3; % ms
            pos = permute([zeros(size(z))
                    zeros(size(z))
                    z], [1, 3, 2]);
            grad = zeros(3, ceil(T/dt));
            B1 = zeros(1, ceil(T/dt));
            T1 = 600*ones(sz);
            T2 = 20:10:120*ones(sz);
            
            [result, ~] = bloch_symmetric_splitting(dt, B1, grad, pos, T1, T2, Minit = [0; 1; 0]); %computed solution
            expected = permute([zeros(1, 11); exp(-T*1./T2); ones(1, 11)-exp(-T./T1)], [1, 3, 2]); %analytical solution
            
            figure()
            plot(T2, permute(result(2, 1, :), [1, 3, 2]))
            hold on
            plot(T2, permute(expected(2, 1, :), [1, 3, 2]))
            hold off
            grid on
            xlabel('T2 time (ms)')
            ylabel('transverse magnitude')
            testCase.verifyLessThan(abs(result(2, 1, :) - expected(2, 1, :)), 1e-10)
        end
    end
end