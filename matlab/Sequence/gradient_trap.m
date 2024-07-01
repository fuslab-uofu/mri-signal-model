function [g, t] = gradient_trap(ampl, rut, dura, rdt, dt)
%% [g, t] = gradient_trap(ampl, rut, dura, rdt, dt)
% Creates trapezoidal gradient waveform following Siemens IDEA convention
% for ramp up time, duration, and ramp down time.
%
% ~ Input ~
% * ampl: Gradient amplitude
% * rut: Ramp up time
% * dura: Duration at peak
% * rdt: ramp down time
% * dt: time step between each sample.
% NOTE: All time parameters should have the same units!
%
% ~ Output ~
% * g: gradient waveform
% * t: time base
%
%% 2023-05-10 Samuel Adams

% Compute number of points
N = round( (dura + rdt)/dt );
Nru = ceil(rut/dt);
Nrd = ceil(rdt/dt);

% Create gradient box with correct amplitude and duration
g = ampl*ones(N, 1);

% Compute ramp functions
ru = linspace(0, ampl, Nru + 1); ru = ru(1:end-1);
rd = linspace(ampl, 0, Nru + 1); rd = rd(1:end-1);

% Apply ramping by replacing the corresponding portions of the waveform
% with ramps
g(1:Nru) = ru; 
g((end - Nrd + 1):end) = rd;

% Generate time base
t = transpose( (0:N-1)*dt );

end