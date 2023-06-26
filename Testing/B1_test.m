%% Unit test that compares B1 and B1map to analytical solutions for the function bloch_symmetric_splitting
% test #1: test that the flip angle matches the area underneath the B1 curve
% test #2: test that a spatially varying non uniform B1 map functions as expected
% test #3: test both #1 and #2 work at the same time
% test #4: test that different phase angle work as expected
%% 2023-06-26 Addison Powell

classdef B1_test < matlab.unittest.TestCase
    properties
    end

    methods (Test)
        %tests different flip angles 1-90
        function diff_FA(testCase)
            addpath(genpath('./'));
            
            T = 2; % s
            z = 0;
            sz = size(z);
            dt = 1e-3; % ms
            pos = permute([zeros(size(z))
                    zeros(size(z))
                    z], [1, 3, 2]);
            grad = zeros(3, ceil(T/dt));
            T1 = inf*ones(sz);
            T2 = inf*ones(sz);
            Minit = [0;0;1];

            for FA=1:90
                [B1e, ~] = b1_hardpulse(FA, T, dt, t_units='s', B1_units='T'); %calculated by hand that B1 is correct for 10deg 2s 1ms
                [result, ~] = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2, Minit=Minit); %computed solution
                expected = [1 0 0; 0 cos(-FA*pi/180) -sin(-FA*pi/180); 0 sin(-FA*pi/180) cos(-FA*pi/180)]*Minit; %analytical solution rotates clockwise
                testCase.verifyLessThan(abs(result - expected), 1e-10);
            end    
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
            [B1e, ~] = b1_hardpulse(30, T, dt, t_units='s', B1_units='T');
            B1map = .9:.02:1.1;
            expected=[];
            Minit = [0;0;1];
            
            [result, ~] = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2, Minit=Minit, B1map=B1map); %computed solution
            for scale=B1map
                expected = [expected, [1 0 0; 0 cos(-scale*30*pi/180) -sin(-scale*30*pi/180); 0 sin(-scale*30*pi/180) cos(-scale*30*pi/180)]*[0;0;1]]; %analytical solution
            end
            testCase.verifyLessThan(abs(permute(result, [1 3 2]) - expected), 1e-10);
            
            
        end
            
        %test different flip angles and a nonuniform B1map
        function FA_and_B1map(testCase)
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
            
            B1map = .9:.02:1.1;
            Minit = [0;0;1];

            for FA=1:90
                [B1e, ~] = b1_hardpulse(FA, T, dt, t_units='s', B1_units='T');
                expected=[];
                for scale=B1map
                    expected = [expected, [1 0 0; 0 cos(-scale*FA*pi/180) -sin(-scale*FA*pi/180); 0 sin(-scale*FA*pi/180) cos(-scale*FA*pi/180)]*[0;0;1]]; %analytical solution
                end
                [result, ~] = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2, Minit=Minit, B1map=B1map); %computed solution
                testCase.verifyLessThan(abs(permute(result, [1 3 2]) - expected), 1e-10);
            end
        end
            
            %test that phase angle works correctly
        function phase_angle_imag(testCase)
            addpath(genpath('./'));
        
            T = 2; % s
            z = 0;1
            sz = size(z);
            dt = 1e-3; % ms
            pos = permute([zeros(size(z))
                    zeros(size(z))
                    z], [1, 3, 2]);
            grad = zeros(3, ceil(T/dt));
            T1 = inf*ones(sz);
            T2 = inf*ones(sz);
            Minit = [0;0;1];
            FA = pi/2;

            for PA=1:90
                [B1e, ~] = b1_hardpulse(90+PA*1i, T, dt, t_units='s', B1_units='T'); %calculated by hand that B1 is correct for 10deg 2s 1ms
                [result, ~] = bloch_symmetric_splitting(dt, B1e, grad, pos, T1, T2, Minit=Minit); %computed solution
                expected = [1 0 0; 0 cos(-(FA+PA*1i)) -sin(-(FA+PA*1i)); 0 sin(-(FA+PA*1i)) cos(-(FA+PA*1i))]*Minit; %analytical solution rotates clockwise
                testCase.verifyLessThan(abs(result - expected), 1e-10);
            end
        end
    end
end