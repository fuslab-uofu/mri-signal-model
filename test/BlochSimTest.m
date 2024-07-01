classdef BlochSimTest < matlab.unittest.TestCase
    %% Unit tests for BlochSim package
    %
    % 2023-05-22 Samuel Adams-Tew

    properties
        reference
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function setup(testCase)
            addpath(genpath('./'))
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        function sliceSelectFATest(testCase)
            FA = 90;
            gamma = 267.5153151e6;
            dt = 4e-3;
            [B1e, ~] = b1_sliceselect(5, 10, 90, 6, dt, gamma);

            % Adjust units
            B1e = B1e*1e-3;
            dt = dt*1e-3;

            % Compute flip angle
            FA_result = gamma*sum(B1e*dt)*180/pi;

            % Give test result
            testCase.verifyTrue(sum(abs(FA(:) - FA_result(:))) < 1e-12);
        end

        function relaxationTest(testCase)
            t = permute([1, 10, 100, 1000], [1, 3, 2]);
            T1 = permute([200, 400, 600, 800], [1, 3, 2]);
            T2 = permute([50, 60, 70, 80], [1, 3, 2]);
            M0 = permute([1, 3, 0.5, 1], [1, 3, 2]);

            [Rel, ~] = relaxation(t, T1, T2, M0);

            Rel_true = permute([exp(-t./T2); exp(-t./T2); exp(-t./T1)], [1, 3, 2]);

            testCase.verifyTrue(sum(abs(Rel(:) - Rel_true(:))) < 1e-12);
        end

        function recoveryTest(testCase)
            t = permute([1, 10, 100, 1000], [1, 3, 2]);
            T1 = permute([200, 400, 600, 800], [1, 3, 2]);
            T2 = permute([50, 60, 70, 80], [1, 3, 2]);
            M0 = permute([1, 3, 0.5, 1], [1, 3, 2]);

            [~, Rec] = relaxation(t, T1, T2, M0);

            Rec_true = [zeros(size(t)); zeros(size(t)); M0.*(1 - exp(-t./T1))];

            testCase.verifyTrue(sum(abs(Rec(:) - Rec_true(:))) < 1e-12);
        end

        function directionTest(testCase)
            % Rotate known vectors around known axes and angles,
            % compare results
            ax = permute([0, 0, 1, 0, 1, 0
                1, 0, 0, 0, 0, 1
                0, 1, 0, 1, 0, 0], [1, 3, 2]);
            ang = permute([-pi/2, -pi/2, -pi/2, -pi/2, -pi/2, -pi/2], [1, 3, 2]);
            R = axis_angle_rotation_matrix(ax, ang);

            in = permute([1, 1, 0, 0, 0, 0
                0, 0, 1, 1, 0, 0
                0, 0, 0, 0, 1, 1], [1, 3, 2]);
            truth = permute([0,  0,  0, 1, 0, -1
                0, -1,  0, 0, 1,  0
                1,  0, -1, 0, 0,  0], [1, 3, 2]);

            result = pagemtimes(R, in);
            testCase.verifyTrue(sum(abs(result(:) - truth(:))) < 1e-12);
        end
    end

end