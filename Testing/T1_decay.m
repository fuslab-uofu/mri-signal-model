%% Unit Test that compares simulated T1 decay to an analytical solution
%test bloch_sym on 11 different T1 values with no B1 or gradient field
classdef T1_decay < matlab.unittest.TestCase
    properties 
        reference
    end

    methods (Test)
        function final_state_relax(testCase)
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
            Min = [0; 1; 0];
            expected = permute([zeros(1, 11); exp(-T*1./T2); ones(1, 11)-exp(-T*1./T1)], [1, 3, 2]); %analytical solution

            [result, ~] = bloch_symmetric_splitting(dt, B1, grad, pos, T1, T2, Minit = Min); %computed solution
            
            figure()
            plot(T1, permute(result(3, 1, :), [1, 3, 2]))
            hold on
            plot(T1, permute(expected(3, 1, :), [1, 3, 2]))
            hold off
            grid on
            xlabel('T1 time (ms)')
            ylabel('logitudinal magnitude')  
            testCase.verifyLessThan((abs(result(3, 1, :) - expected(3, 1, :))), 1e-10)
        end
        
        function final_state_rec(testCase)
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
            Min = [0; 0; -1];
            expected = permute([zeros(1, 11); exp(-T*1./T2); ones(1, 11)-2*exp(-T./T1)], [1, 3, 2]); %analytical solution

            [result, ~] = bloch_symmetric_splitting(dt, B1, grad, pos, T1, T2, Minit = Min); %computed solution
            
            figure()
            plot(T1, permute(result(3, 1, :), [1, 3, 2]))
            hold on
            plot(T1, permute(expected(3, 1, :), [1, 3, 2]))
            hold off
            grid on
            xlabel('T1 time (ms)')
            ylabel('logitudinal magnitude')  
            testCase.verifyLessThan((abs(result(3, 1, :) - expected(3, 1, :))), 1e-10)
                
        end
    end
end