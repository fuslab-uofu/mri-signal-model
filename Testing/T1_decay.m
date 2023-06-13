%% Unit Test that compares simulated T1 decay to an analytical solution
%test bloch_sym on 11 different T1 values with no B1 or gradient field
classdef T1_decay < matlab.unittest.TestCase
    properties
    end

    methods (Test)
        function final_state_err(testCase)
            clear, close all
            addpath(genpath('./'))

            T = 2; % s
            z = -5:5;
            sz = size(z);
            dt = 1e-3; % ms
            pos = permute([zeros(size(z))
                    zeros(size(z))
                    z], [1, 3, 2]);
            grad = zeros(3, ceil(T/dt));
            B1 = zeros(1, ceil(T/dt));
            T1 = 100:100:1100*ones(sz);
            T2 = 50*ones(sz);
            
            [result, ~] = bloch_symmetric_splitting(dt, B1, grad, pos, T1, T2, Minit = [0; 1; 0]); %computed solution
            expected = permute([zeros(1, 11); exp(-T*1./T2); ones(1, 11)-exp(-T./T1)], [1, 3, 2]); %analytical solution
            
            figure()
            plot(T1, permute(result(3, 1, :), [1, 3, 2]))
            hold on
            plot(T1, permute(expected(3, 1, :), [1, 3, 2]))
            hold off
            grid on
            xlabel('T1 time (ms)')
            ylabel('horizontal magnitude')  
            testCase.verifyEqual(result(3, 1, :), expected(3, 1, :))
        end
    end
end