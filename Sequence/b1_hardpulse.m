function [B1, t] = b1_hardpulse(FA, dura, dt, options)
%% [B1, t] = b1_hardpulse(FA, dura, dt, options)
% Creates a rectangular 'hard' excitation pulse
%
% ~ Input ~
% * FA: Desired flip angle
% * dura: Desired pulse duration
% * dt: Desired time step
%
% ~ Options ~
% * gamma: Gyromagnetic ratio used for computing flip angle. Default is
% that of water protons, 267.5153151e6.
% * FA_units: Units for input flip angle. Default is 'deg'.
% * t_units: Units for input dt and output t. Default is 'us'.
% * B1_units: Units for output B1 amplitude. Default is 'mT'.
%
% ~ Output ~
% * B1: RF waveform
% * t: time base
%
%% 2023-05-11 Samuel Adams
arguments
    FA
    dura
    dt
    options.gamma = 267.5153151e6
    options.FA_units = 'deg'
    options.t_units = 'ms'
    options.B1_units = 'mT'
end

switch options.FA_units
    case 'deg'
        FA = FA*pi/180;
    case 'rad' % do nothing
    otherwise
        error('Unrecognized units for FA: %s', options.FA_units)
end
switch options.t_units
    case 'ms'
        dura = dura*1e-3;
        dt = dt*1e-3;
    case 'us'
        dura = dura*1e-6;
        dt = dt*1e-6;
    case 'ns'
        dura = dura*1e-9;
        dt = dt*1e-9;
    case 's' % do nothing
    otherwise
        error('Unrecognized units for t: %s', options.t_units)
end

% Create time vector
t = 0:dt:dura; t = t(1:end-1);
% Create rect envelope of correct size
B1 = ones(size(t));
% Compute amplitude for desired flip angle
A = FA/(options.gamma*sum(B1.*dt));
% Apply amplitude change
B1 = A*B1; 

% Adjust output units
switch options.B1_units
    case 'mT'
        B1 = B1*1e3;
    case 'T' % do nothing
    otherwise
        error('Unrecognized units for B1: %s', options.B1_units)
end
switch options.t_units
    case 'ms'
        t = t*1e3;
    case 'us'
        t = t*1e6;
    case 'ns'
        t = t*1e9;
    case 's' % do nothing
    otherwise
        error('Unrecognized units for t: %s', options.t_units)
end


end