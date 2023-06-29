%% Unit Tests that check that off gradients and position cause changes in resonance frequency as expected for the function bloch_symmetric_splitting
% all other variables are held constant
%% 2023-06-26 Addison Powell

classdef gradient_test < matlab.unittest.TestCase
    properties 
        reference
    end

    methods (Test)
        function grad_test(testCase)
            addpath(genpath('./'))

            T = 2; % s
            z = -5:5;
            sz = size(z);
            dt = 1e-3; % ms
            pos = permute([zeros(size(z))
                    zeros(size(z))
                    z], [1, 3, 2]);
            T1 = inf*ones(sz);
            T2 = inf*ones(sz);
            grad = gradient_trap(30, 0.150, T-0.150, 0.150, dt)';
            tmp = gradient_trap(-30, 0.150, 0.95, 0.15, dt)';
            grad = [grad, tmp];
            B1 = zeros(1, length(grad));
            grad = [zeros(2, length(grad)); grad];
            zerograd = zeros(size(grad));
            Mloop = [0;1;0];
            

            for iter=length(B1)
                tmp = sum(grad(3, iter).*pos);
                Mloop = [cos(-pi/2) -sin(-pi/2) 0; sin(-pi/2) cos(-pi/2) 0; 0 0 1]*Mloop; %not rotating
            end
           
            expected = [1 0 0; 0 cos(-pi/2) -sin(-pi/2); 0 sin(-pi/2) cos(-pi/2)]*[0; 1; 0]; %analytical solution

            [result, ~] = bloch_symmetric_splitting(dt, B1, grad, pos, T1, T2, Minit=[0;1;0]); %computed solution
            testCase.verifyLessThan(abs(result - expected), 1e-10);

            
            
        end
    end
end