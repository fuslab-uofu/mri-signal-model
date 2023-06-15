function simulate_acquisition(T1map, T2map, pos, seq, options, units)
%% simulate_acquisition(T1map, T2map, pos, seq, options, units)
%
% ~ Input ~
% NOTE: [X] denotes up to 3 dimensions, enabling the manipulation of
% collections that have up to 3 independent variables, such as full 3D.
% * T1map (1,1,[X]): Longitudinal relaxation time for each point
% * T2map (1,1,[X]): Transverse relaxation time for each point
% * pos (3,1,[X]): Spatial coordinates for each point
% * seq (PulseSequence): PulseSequence object that specifies the RF,
% gradients, and sampling.
%
% ~ Options ~
% * numRepetitions (scalar): Number of times to consecutively simulate the
% application of seq, which may be required for achieving steady-state
% condition for some sequences. Default is 1.
% * saveEvery (scalar): How many repetitions should be completed before
% saving the result. Also saves every saveEvery subsequent iterations.
% Default is -1, which saves only the last repetition.
% * savePath (string): File path where the output should be saved. See
% output section below for more information on the saved file.
% * dt_max (scalar): Maximum time step permitted for any part of the
% simulation. If the provided seq requires that a smaller timestep be used
% for a particular portion of the simulation, that timestep will be used.
% If dt_max==Inf, the time specified by seq is always used. Default is Inf.
% * B0map (1,1,[X]): Spatially varying main field inhomogeneities, applied
% as an additive factor. Default is 0 for all points.
% * B1map (1,1,[X]) [scale factor]: Spatially varying B1 field, applied by
% multiplication with the B1 field at each time. Default is 1 for all
% points.
% * delta (1,1,[X]) [ppm]: Frequency offset / chemical shift. Default is 0
% for all points.
% * B0 (scalar) [T]: Main magnetic field strength. Default is 3.
% * gamma (scalar) [rad/(s*T)]: Gyromagnetic ratio for on-resonance
% excitation. Default is that of water protons, 267.5153151e6.
%
% ~ Units ~
% Use the units options to specify what units are being used for each input
% variable. Includes T1_units, T2_units, pos_units, dt_units, B0map_units.
% See CONVERT_UNITS for options.
%
% ~ Output ~
% This function has no returned output. Simulation results are saved to the
% file specified by savePath. The saved file contains copies of the input
% parameters with an additional variable M
% * M (3,S,R,[X]): Magnetization vector at each sample time specified by
% seq (samples along dim=2, with S=seq.numSamples). Each repetition of the
% sequence is along dim=3 (R=numRepetitions).
%
%% 2023-05-22 Samuel Adams-Tew
arguments
    T1map (1,1,:,:,:,:,:)
    T2map (1,1,:,:,:,:,:)
    pos (3,1,:,:,:,:,:)
    seq PulseSequence
    options.numRepetitions (1, 1) = 1;
    options.saveEvery (1, 1) = -1;
    options.savePath = './mrsim.mat';
    options.dt_max (1, 1) double {mustBePositive} = Inf
    options.B1map (1,1,:,:,:,:,:) = 1;
    options.B0map (1,1,:,:,:,:,:) = 0;
    options.delta (1,1,:,:,:,:,:) = 0;
    options.B0 (1,1) = 3;
    options.gamma (1,1) = 267.5153151e6;
    options.showProgressBar (1,1) logical = false;
    units.T1_units = 'ms'
    units.T2_units = 'ms'
    units.pos_units = 'mm'
    units.dt_units = 'us'
    units.B0map_units = 'uT'
end

%% Parse and validate inputs

% Extent of the object being simulated 
fieldSize = size(pos, 3:7);
lastDim = 2 + find(fieldSize > 1, 1, 'last'); % Find the last non-singleton dimension
if isempty(lastDim); lastDim = 2; end

if ~isequal(size(T1map), size(T2map))
    error('T1map and T2map must have same shape')
end
if ~isequal(size(T1map, 3:7), size(pos, 3:7))
    error('T1map and position must have same field size')
end
% Check that the frequency offset is the same size as spatial inputs
delta = options.delta;
if length(delta) == 1
    delta = delta*ones(size(T1map)); % apply same offset to all points
elseif ~isequal(size(delta), size(T1map))
    error('delta must have same shape as T1map')
end
% Check that B0 inhomogeneity is the same size as spatial inputs
B0map = convert_units(options.B0map, units.B0map_units, 'T');
if length(B0map) == 1
    B0map = B0map*ones(size(T1map)); % apply same offset to all points
elseif ~isequal(size(B0map), size(T1map))
    error('delta must have same shape as T1map')
end
% Check that B1 inhomogeneity is the same size as spatial inputs
B1map = options.B1map;
if length(B1map) == 1
    B1map = B1map*ones(size(T1map)); % apply same offset to all points
elseif ~isequal(size(B1map), size(T1map))
    error('delta must have same shape as T1map')
end

% Constants
gamma = options.gamma;
B0 = options.B0;

% Check save parameters
numRepetitions = options.numRepetitions;
saveEvery = options.saveEvery;
if saveEvery == -1
    saveEvery = numRepetitions;
elseif saveEvery < 1
    error('Provided value saveEvery=%d is not allowed. Must be at least 1.', saveEvery)
elseif saveEvery > numRepetitions
    error('Provided value saveEvery=%d is not allowed. Must be less than numRepetitions=%d.', saveEvery, numRepetitions)
end
% Initialize save file
if ~isempty(options.savePath)
    save(options.savePath, 'T1map', 'T2map', 'pos', 'seq', 'options', 'units', '-v7.3')
    file = matfile(options.savePath, 'Writable', true);
    % Initialize the last element in the array to force file system to
    % prepare to accomodate the full data size
    file.magnetization(3, seq.numSamples, floor(numRepetitions/saveEvery), ...
        fieldSize(1), fieldSize(2), fieldSize(3), fieldSize(4), fieldSize(5)) = double(0);
end

% Convert units
T1map = convert_units(T1map, units.T1_units, 's');
T2map = convert_units(T2map, units.T2_units, 's');
pos = convert_units(pos, units.pos_units, 'm');
dt_max = convert_units(options.dt_max, units.dt_units, 's');

% Progress bar option
showpb = options.showProgressBar && exist('ProgressBar', 'class');

clear options units
% Start timing
tstart = tic();

%% Prepare for multiple iterations, if needed

% Compute sequence event operations to accelerate repetitions with no
% sampling
sz = size(T1map); sz(1:2) = 3;
eventOp = cell(1, seq.numEvents);
eventAdd = cell(1, seq.numEvents);
if numRepetitions > 1 && saveEvery ~= 1
    % If there are repetitions without sampling, we can accelerate by
    % precomputing an operator for the full sequence
    seqOp = eye(3).*ones(sz);
    seqAdd = [0; 0; 0].*zeros(size(T1map));
end

%% Compute operators for events in the sequence

if showpb; pb = ProgressBar('Computing operators', seq.numEvents); end

% Iterate over all events in the sequence
for eventNum = 1:seq.numEvents
    if showpb; pb.iter(); end

    % Get waveforms for this event
    [dt, B1, G, sampleComb] = seq.get_event(eventNum, dt_max, 's');

    % Convert units
    B1 = convert_units(B1, 'mT', 'T');
    G = convert_units(G, 'mT/m', 'T/m');

    [eventOp{eventNum}, eventAdd{eventNum}] = event_operator(dt, B1, G, 'pos', pos, 'T1', T1map, 'T2', T2map, 'delta', delta, 'B0map', B0map, 'B1map', B1map, 'sampleComb', sampleComb, 'gamma', gamma, 'B0', B0, 'showProgressBar', showpb);
    % If multiple TRs, speed up by computing operator for full TR
    if numRepetitions > 1 && saveEvery ~= 1
        if any(sampleComb)
            % Event is sampled and has intermediate operators. Use the
            % operator for the final state.
            seqOp = pagemtimes(eventOp{eventNum}{end}, seqOp);
            seqAdd = pagemtimes(eventOp{eventNum}{end}, seqAdd) + eventAdd{eventNum}{end};
        else
            % Event has no intermediate operators
            seqOp = pagemtimes(eventOp{eventNum}, seqOp);
            seqAdd = pagemtimes(eventOp{eventNum}, seqAdd) + eventAdd{eventNum};
        end
    end
end
if showpb; pb.close(); end

%% Iterate over all repetitions

if showpb; pb = ProgressBar('Simulating repetition', numRepetitions); end

% Create initial magnetization
Mloop = [zeros([2, 1, fieldSize]); ones([1, 1, fieldSize])];

for repNum = 1:numRepetitions
    % Update progress bar
    if showpb; pb.iter(repNum); end

    if mod(repNum, saveEvery) ~= 0
        % This iteration is not sampled. Perform as single operation.
        Mloop = pagemtimes(seqOp, Mloop) + seqAdd;
    else
        % This iteration will be sampled. Perform each event separately
        samples = zeros([3, seq.numSamples, fieldSize]);
        
        repSampIdx = 0; % Index counting samples saved for this repetition
        roEvents = seq.readOutEvents;
        % Iterate over all events in the pulse sequence
        for eventNum = 1:seq.numEvents
            if any(roEvents == eventNum)
                % This event is sampled. Use individual operators for each
                % sample.
                nEventSamples = length(eventOp{eventNum}) - 1;
                for eventSampIdx = 1:nEventSamples
                    samples(:, repSampIdx + eventSampIdx,:,:,:) = ...
                        pagemtimes(eventOp{eventNum}{eventSampIdx}, Mloop) + eventAdd{eventNum}{eventSampIdx};
                end

                % Update count of samples for this index
                repSampIdx = repSampIdx + nEventSamples;

                % Get new Mloop using operator for final state
                Mloop = pagemtimes(eventOp{eventNum}{end}, Mloop) + eventAdd{eventNum}{end};
            else
                % This event is not sampled. Perform as single operation.
                Mloop = pagemtimes(eventOp{eventNum}, Mloop) + eventAdd{eventNum};
            end
        end

        % Save out sampled data for this repetition
        if ~isempty(samples)
            switch lastDim
                % Writing to a matfile requires that the shape of the RHS
                % exactly matches that of the LHS, which requires the permute
                % statement included below
                case 2
                    file.magnetization(:, :, floor(repNum/saveEvery)) = permute(samples, [1, 2, lastDim+1, 3:lastDim]);
                case 3
                    file.magnetization(:, :, floor(repNum/saveEvery), :) = permute(samples, [1, 2, lastDim+1, 3:lastDim]);
                case 4
                    file.magnetization(:, :, floor(repNum/saveEvery), :, :) = permute(samples, [1, 2, lastDim+1, 3:lastDim]);
                case 5
                    file.magnetization(:, :, floor(repNum/saveEvery), :, :, :) = permute(samples, [1, 2, lastDim+1, 3:lastDim]);
                case 6
                    file.magnetization(:, :, floor(repNum/saveEvery), :, :, :, :) = permute(samples, [1, 2, lastDim+1, 3:lastDim]);
                case 7
                    file.magnetization(:, :, floor(repNum/saveEvery), :, :, :, :, :) = permute(samples, [1, 2, lastDim+1, 3:lastDim]);
                case 8
                    file.magnetization(:, :, floor(repNum/saveEvery), :, :, :, :, :, :) = permute(samples, [1, 2, lastDim+1, 3:lastDim]);
                case 9
                    file.magnetization(:, :, floor(repNum/saveEvery), :, :, :, :, :, :, :) = permute(samples, [1, 2, lastDim+1, 3:lastDim]);
            end
        end
    end
end

% Close progress bar
if showpb; pb.close(); end
file.telapsed = toc(tstart); % Record how long the simulation took

end
