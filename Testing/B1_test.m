%% Unit test that compares B1 and B1map to analytical solutions
% test #1: test that the flip angle matches the area underneath the B1
% curve
% test #2: test that a spatially varying non uniform B1 map functions as expected
classdef B1_test < matlab.unittest.TestCase
    properties
    end

    methods (Test)
        %tests different flip angles 1-90
        function diff_FA(testCase)
            addpath(genpath('./'));
            
            Gz = 30;
            dz = 5;
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

            for FA=1:90
                [B1e, ~] = b1_sliceselect(Gz, dz, FA, T, dt);
                
                [result, ~] = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2); %computed solution
                expected = permute([zeros(1, 11); exp(-T*1./T2); ones(1, 11)-exp(-T./T1)], [1, 3, 2]); %analytical solution
                testCase.verifyLessThan(abs(result(3, 1, :) - expected(3, 1, :)), 1e-10);
            end
            
            figure()
            plot(T1, permute(result(3, 1, :), [1, 3, 2]))
            hold on
            plot(T1, permute(expected(3, 1, :), [1, 3, 2]))
            hold off
            grid on
            xlabel('T1 time (ms)')
            ylabel('horizontal magnitude')  
            
        end
        
        % non uniform B1map test case
        function non_uniform_B1map(testCase)
            addpath(genpath('./'));
            
            Gz = 30;
            dz = 5;
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
            [B1e, ~] = b1_sliceselect(Gz, dz, 90, T, dt);
            B1map = .9001:.0001:1.1;
            
            [result, ~] = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2, B1map=B1map); %computed solution
            expected = permute([zeros(1, 11); exp(-T*1./T2); ones(1, 11)-exp(-T./T1)], [1, 3, 2]); %analytical solution
            testCase.verifyLessThan(abs(result(3, 1, :) - expected(3, 1, :)), 1e-10);
            
            
        end
    end
end