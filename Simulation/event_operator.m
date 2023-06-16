function [A, b] = event_operator(dt, B1, grad, options)
%% [A, b] = event_operator(dt, B1, grad, options)
% Computes the operator for an event step in MR sequence simulation. If the
% event includes sampling, an operator is created for each sample time
% point as well as the final simulation state.
%
% ~ Input ~
% * dt (scalar) [s]: Time step between each element of the input waveforms.
% * B1 (1,T) [T]: Complex valued B1 field values for each time point.
% * grad (3,T) [T/m]: Gradient field waveforms along x-, y-, and z-axes.
%
% ~ Options ~
% * gamma (scalar) [rad/(s*T)]: Gyromagnetic ratio for on-resonance
% excitation. Default is that of water protons, 267.5153151e6.
% * B0 (scalar) [T]: Main magnetic field strength. Default is 3.
% * sampleComb (1,T) (logical): Array that indicates which time points
% should be sampled. If any(sampleComb), [A, b] are cell arrays with
% operators that give the magnetization state at each sample time as well
% as the final state. Otherwise, [A, b] are tensor fields for a single
% operation that applies the full event to get the final state.
% * showProgressBar (logical): Whether a progress bar indicating progress
% in computing operator should be displayed. Default is false.
%
% Spatially-varying optional inputs. 
% NOTE: [X] denotes up to 4 dimensions, enabling the manipulation of
% collections that have up to 4 independent variables, such as full 3D plus
% varying degrees of frequency shift.
% * x (3,1,[X]) [m]: Spatial coordinates for each point
% * T1 (1,1,[X]) [s]: Longitudinal relaxation time for each point
% * T2 (1,1,[X]) [s]: Transverse relaxation time for each point
% * B0map (1,1,[X]): Spatially varying main field inhomogeneities, applied
% as an additive factor. Default is 0 for all points.
% * B1map (1,1,[X]) [scale factor]: Spatially varying B1 field, applied by
% multiplication with the B1 field at each time. Default is 1 for all
% points.
% * delta (1,1,[X]) [ppm]: Frequency offset / chemical shift. Default is 0
% for all points.
%
% ~ Output ~
% * A (3,3,[X]): Tensor field containing the multiplication part of a
% linear operation applying the excitation, precession, and relaxation for
% the given sequence waveforms and object properties.
% * b (3,1,[X]): Vector field containing the bias for relaxation effects
% during the event.
%
%% 2023-06-12 Samuel Adams-Tew
arguments
    dt (1,1) {mustBeNumeric, mustBeReal, mustBePositive}
    B1 (1,:) {mustBeNumeric}
    grad (3,:) {mustBeNumeric, mustBeReal}
    options.T1 (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = Inf;
    options.T2 (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = Inf;
    options.pos (3,1,:,:,:,:) {mustBeNumeric, mustBeReal} = [0; 0; 1];
    options.B1map (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = 1;
    options.B0map (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = 0;
    options.delta (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = 0;
    options.sampleComb (1,:) logical = [];
    options.Meq (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = 1;
    options.B0 (1,1) {mustBeNumeric, mustBeReal} = 3;
    options.gamma (1,1) {mustBeNumeric, mustBeReal} = 267.5153151e6;
    options.showProgressBar (1,1) logical = false;
end
%% Parse inputs
% Check that the input waveforms have the same number of samples
if length(B1) ~= size(grad, 2) % Validate duration
    error('Temporally varying inputs (B1, grad) must have same number of columns.')
end
% Check whether this is a sampled event and validate sampleComb
if any(options.sampleComb)
    if length(B1) ~= length(options.sampleComb)
        error('Temporally varying inputs (B1, sampleComb) must have same number of columns.')
    end
    sampleIter = find(options.sampleComb);
else
    sampleIter = [];
end
% Read in spatial variables
x = options.pos; T1 = options.T1; T2 = options.T2;
% Extent of spatial / independent variables
fieldSize = size(x, 3:6);
% Check that spatially varying inputs have same number of entries
if any(fieldSize ~= [size(T1, 3:6); size(T2, 3:6)], "all")
    error('Spatially varying inputs (pos, T1, T2) must have same number of independent variables.')
end
% Check that the frequency offset is the same size as spatial inputs
delta = options.delta*1e-6; % convert 1 ppm = 1e-6
if length(delta) == 1
    delta = delta.*ones(size(T1));
elseif any(fieldSize ~= size(delta, 3:6))
    error('delta must have same number of independent variables as pos.')
end
% Check that B0 inhomogeneity is the same size as spatial inputs
B0map = options.B0map;
if length(B0map) == 1
    B0map = B0map.*ones(size(T1)); % Apply same weight to all points
elseif any(fieldSize ~= size(B0map, 3:6)) % Validate number of points
    error('dB0 must have same number of independent variables as pos.')
end
% Check that B1 inhomogeneity is same size as spatial inputs
B1map = options.B1map;
if length(B1map) == 1
    B1map = B1map.*ones(size(T1)); % Apply same weight to all points
elseif any(fieldSize ~= size(B1map, 3:6)) % Validate number of points
    error('dB1 must have same number of independent variables as pos.')
end
% Check that equilibrium magnetization is same size as spatial inputs
Meq = options.Meq;
if length(Meq) == 1
    Meq = Meq.*ones(size(T1)); % Apply same weight to all points
elseif any(fieldSize ~= size(Meq, 3:6)) % Validate number of points
    error('delta must have same number of independent variables as pos.')
end
% Constants
gamma = options.gamma;
B0 = options.B0;

% Progress bar option
showpb = options.showProgressBar && exist('ProgressBar', 'class');

clear options

%% Initialize variables 

if ~isempty(sampleIter)
    % If this is sampled, return an operator that gives the state at each
    % sampled time and an operator for the final state as a cell array.
    A = cell(1, length(sampleIter) + 1);
    b = cell(1, length(sampleIter) + 1);
end

% Initialize with the identity element for the respective operations
Aloop = eye(3).*ones(size(T1)); % Matrix multiplication
bLoop = [0; 0; 0].*ones(size(T1)); % Vector addition

%% Iterate over event steps, compute operator field

nIter = length(B1);
if showpb; pb = ProgressBar('Event operator iteration', nIter); end
saveNum = 1; % Index into sampleIter, (saveNum - 1) is how many values have been saved
% Loop over each step of the event
for iter = 1:nIter
    if showpb; pb.iter(); end
    % Compute axis of rotation
    axX = (1 + delta).*real(B1(iter)*B1map);
    axY = (1 + delta).*imag(B1(iter)*B1map);
    axZ = ( delta*B0 + (1 + delta).*(B0map + sum(grad(:, iter).*x)) );
    % Compute half of rotation angle; note clockwise rotation
    ang = -gamma*sqrt(axX.^2 + axY.^2 + axZ.^2) * dt/2;

    % Compute rotation matrix
    Rot = axis_angle_rotation_matrix([axX; axY; axZ], ang);
    % Compute relaxation terms
    [Rel, Rec] = relaxation(dt, T1, T2, Meq);

    % Compute matrix and offset components of operator
    Aloop = pagemtimes(Rot, Rel.*pagemtimes(Rot, Aloop));
    bLoop = pagemtimes(Rot, Rel.*pagemtimes(Rot, bLoop)) + Rec;

    % If this is a sampled event, save operator for sampled states
    if ~isempty(sampleIter) && iter == sampleIter(saveNum)
        A{saveNum} = Aloop; b{saveNum} = bLoop; % Store operators in output
        saveNum = min(saveNum + 1, length(sampleIter)); % Update save count
    end
end
if showpb; pb.close(); end

if isempty(sampleIter)
    % Return operator for final state
    A = Aloop; b = bLoop;
else
    % Save operator for final state to output arrays
    A{end} = Aloop; b{end} = bLoop;
end

end