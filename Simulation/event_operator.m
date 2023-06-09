function eventOp = event_operator(dt, B1, grad, options)
arguments
    dt (1, 1) {mustBeNumeric, mustBeReal, mustBePositive}
    B1 (1, :) {mustBeNumeric}
    grad (3, :) {mustBeNumeric, mustBeReal}
    options.T1 (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = Inf;
    options.T2 (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = Inf;
    options.pos (3,1,:,:,:,:) {mustBeNumeric, mustBeReal} = [0; 0; 1];
    options.B1map (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = 1;
    options.B0map (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = 0;
    options.delta (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = 0;
    options.Meq (1,1,:,:,:,:) {mustBeNumeric, mustBeReal} = 0;
    options.B0 (1,1) {mustBeNumeric, mustBeReal} = 3;
    options.gamma (1,1) {mustBeNumeric, mustBeReal} = 267.5153151e6;
end

% Check that the input waveforms have the same number of samples
if length(B1) ~= size(grad, 2) % Validate duration
    error('Temporally varying inputs (B1, grad) must have same number of columns.')
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

clear options

% Initialize operators
eventOp = eye(3).*ones(size(T1));

% Loop over each step of the event
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

    eventOp = pagemtimes(Rot, pagemtimes(Rel.*Rot, eventOp) + Rec);
end

end