function [Mfinal, Msample] = bloch_symmetric_splitting(dt, B1, grad, x, T1, T2, options)
%% [Mfinal, Msample] = bloch_symmetric_splitting(dt, B1, grad, x, T1, T2, options)
% Simulate the Bloch equations using symmetric operator splitting.
%
% Units are whole SI units [m, s, T], NOT fractional units [mm, ms, us mT].
% Unless otherwise noted, all inputs are expected to be real numbers.
%
% ~ Input ~
% * dt (scalar): Time step between each sample of the input waveforms
%
% Temporally-varying inputs
% * B1 (1,T) (complex): Demodulated B1 (envelope function for on-resonance
% excitation). Real and imaginary parts correspond to the x- and y-axes of
% the rotating frame, respectively.
% * grad (3,T): Gradient field waveforms along x-, y-, and z-axes
%
% Spatially-varying inputs
% NOTE: [X] denotes up to 4 dimensions, enabling
% the manipulation of collections that have up to 4 independent variables,
% such as full 3D plus varying degrees of frequency shift.
% * x (3,1,[X]): Spatial coordinates for each point
% * T1 (1,1,[X]): Longitudinal relaxation time for each point
% * T2 (1,1,[X]): Transverse relaxation time for each point
%
% ~ Options ~
% * gamma (scalar) [rad/(s*T)]: Gyromagnetic ratio for on-resonance
% excitation. Default is that of water protons, 267.5153151e6.
% * B0 (scalar) [T]: Main magnetic field strength. Default is 3.
% * samplePeriod (scalar): Alternative to sampleComb, specifies how
% frequently states should be added to Msample. Saves simulation state
% after steps s:s:end, where s = round(samplePeriod/dt). Default gives
% Msample = [].
%
% Temporally-varying optional inputs
% * sampleComb (1,T) (boolean): Array that indicates which time points
% should be included in Msample, with 'true' at all times where the state
% should be included and 'false' everywhere else. Default gives Msample =
% [].
% * B0drift (1,T): Temporally varying changes in the main magnetic field,
% applied as an additive factor. Default is 0 for all time.
%
% Spatially-varying optional inputs
% * B0map (1,1,[X]): Spatially varying main field inhomogeneities, applied
% as an additive factor. Default is 0 for all points.
% * B1map (1,1,[X]) [scale factor]: Spatially varying B1 field, applied by
% multiplication with the B1 field at each time. Default is 1 for all
% points.
% * delta (1,1,[X]) [ppm]: Frequency offset / chemical shift. Default is 0
% for all points.
% * Meq (1,1,[X]): Equilibrium magnetization along Z-direction (usually
% denoted by M0, but named differently to avoid confusion with Minit).
% Default is 1 for all points.
% * Minit (3,1,[X]): Initial magnetization vectors, if different from the
% equilibrium magnetization. Default is [0, 0, Meq]' for all points.
%
% ~ Output ~
% * Mfinal (3,1,[X]): The last simulation state for each point in space
% * Msample (3,S,[X]): Intermediate simulation states for each point in
% space. Empty unless sampleComb or samplePeriod are defined.
%
%% 2023-05-15 Samuel Adams-Tew
arguments
    dt (1, 1) {mustBeNumeric, mustBeReal, mustBePositive}
    B1 (1, :) {mustBeNumeric}
    grad (3, :) {mustBeNumeric, mustBeReal}
    x (3, 1, :, :, :, :) {mustBeNumeric, mustBeReal}
    T1 (1, 1, :, :, :, :) {mustBeNumeric, mustBeReal}
    T2 (1, 1, :, :, :, :) {mustBeNumeric, mustBeReal}
    options.gamma (1, 1) {mustBeNumeric, mustBeReal} = 267.5153151e6;
    options.B0 (1, 1) {mustBeNumeric, mustBeReal} = 3;
    options.samplePeriod (1, 1) {mustBePositive, mustBeReal}
    options.sampleComb (1, :) logical = [];
    options.B0drift (1, 1, :, :, :, :) {mustBeNumeric, mustBeReal} = 0;
    options.B0map (1, 1, :, :, :, :) {mustBeNumeric, mustBeReal} = 0;
    options.B1map (1, 1, :, :, :, :) {mustBeNumeric, mustBeReal} = 1;
    options.delta (1, 1, :, :, :, :) {mustBeNumeric, mustBeReal} = 0;
    options.Meq (1, 1, :, :, :, :) {mustBeNumeric, mustBeReal} = 1;
    options.Minit (3, 1, :, :, :, :) {mustBeNumeric, mustBeReal};
end

%% Parse and validate inputs

% Check that the input waveforms have the same number of samples
if length(B1) ~= size(grad, 2) % Validate duration
    error('Temporally varying inputs (B1, grad) must have same number of columns.')
end
% Temporal sampling of results
if isempty(options.sampleComb)
    if isfield(options, 'samplePeriod')
        % Compute how frequently to save states
        s = round(options.samplePeriod/dt);
        % Identify which states to save
        saveIter = s:s:length(B1);
    else
        % Do not save intermediate states
        saveIter = [];
    end
else
    % Identify which states are specified to save
    saveIter = find(options.sampleComb);
end
% Extent of spatial / independent variables
fieldSize = size(x, 3:6);
% Check that spatially varying inputs have same number of entries
if any(fieldSize ~= [size(T1, 3:6); size(T2, 3:6)], "all")
    error('Spatially varying inputs (x, T1, T2) must have same number of independent variables.')
end
% Check that the frequency offset is the same size as spatial inputs
delta = options.delta*1e-6; % convert 1 ppm = 1e-6
if length(delta) == 1
    delta = delta.*ones(size(T1));
elseif any(fieldSize ~= size(delta, 3:6))
    error('delta must have same number of independent variables as x.')
end
% Check that B0 inhomogeneity is the same size as spatial inputs
B0map = options.B0map;
if length(B0map) == 1
    B0map = B0map.*ones(size(T1)); % Apply same weight to all points
elseif any(fieldSize ~= size(B0map, 3:6)) % Validate number of points
    error('dB0 must have same number of independent variables as x.')
end
% Check that B0 drift is the same size as temporal inputs
B0drift = options.B0drift;
if length(B0drift) == 1
    B0drift = B0drift.*ones(1, size(B1, 2));
elseif size(B1, 2) ~= size(B0drift, 2)
    error('B0drift must have same number of columns as B1.')
end
% Check that B1 inhomogeneity is same size as spatial inputs
B1map = options.B1map;
if length(B1map) == 1
    B1map = B1map.*ones(size(T1)); % Apply same weight to all points
elseif any(fieldSize ~= size(B1map, 3:6)) % Validate number of points
    error('dB1 must have same number of independent variables as x.')
end
% Check that equilibrium magnetization is same size as spatial inputs
Meq = options.Meq;
if length(Meq) == 1
    Meq = Meq.*ones(size(T1)); % Apply same weight to all points
elseif any(fieldSize ~= size(Meq, 3:6)) % Validate number of points
    error('delta must have same number of independent variables as x.')
end
% Constants
gamma = options.gamma;
B0 = options.B0;

%% Initialize magnetization state and save variables

if isfield(options, 'Minit')
    if size(options.Minit, 3) == 1
        % Apply single value to all points
        Mloop = options.Minit.*ones([3, 1, fieldSize]); 
    elseif all(fieldSize == size(options.Minit, 3:6))
        % Use provided initial value
        Mloop = options.Minit;
    else
        error('Minit must have same number of independent variables as x.')
    end
else
    % Initialize with equilibrium magnetization
    Mloop = [zeros([2, 1, fieldSize]); Meq.*ones([1, 1, fieldSize])];
end

clear options

% Try to allocate memory for intermediate results.
% If this fails, fall back to saving only the final state
try
    Msample = zeros([3, length(saveIter), fieldSize]);
catch ME
    % Handle error and warn, if appropriate
    switch ME.identifier
        case 'MATLAB:array:SizeLimitExceeded'
            warning('Cannot allocate array larger than maximum array size preference. bloch_symmetric_splitting will return only Mfinal.')
        case 'MATLAB:pmaxsize'
            warning('Cannot allocate array larger than maximum possible variable size. bloch_symmetric_splitting will return only Mfinal.')
        otherwise
            rethrow(ME)
    end
    % Do not sample intermediate states
    saveIter = [];
    Msample = [];
end
saveNum = 1; % Track how many values have been saved for efficient indexing into saveIter

%% Simulation loop
for iter = 1:length(B1)
    % Compute axis of rotation
    axX = (1 + delta).*real(B1(iter)*B1map);
    axY = (1 + delta).*imag(B1(iter)*B1map);
    axZ = ( delta*B0 + (1 + delta).*(B0map + B0drift(iter) + sum(grad(:, iter).*x)) );
    % Compute half of rotation angle; note clockwise rotation
    ang = -gamma*sqrt(axX.^2 + axY.^2 + axZ.^2) * dt/2;

    % Compute rotation matrix
    Rot = axis_angle_rotation_matrix([axX; axY; axZ], ang);
    % Compute relaxation terms
    [Rel, Rec] = relaxation(dt, T1, T2, Meq);

    % Update magnetization state
    Mloop = pagemtimes(Rot, Rel.*pagemtimes(Rot, Mloop) + Rec);

    % If the current iteration is a sampled time, save the magnetization
    % state
    if ~isempty(saveIter) && iter == saveIter(saveNum)
        Msample(:, saveNum, :, :, :, :) = Mloop;
        % Increment state counter, but not beyond the number of elements in
        % saveIter
        saveNum = min(length(saveIter), saveNum + 1);
    end
end

% Return the final simulatin state
Mfinal = Mloop;

end