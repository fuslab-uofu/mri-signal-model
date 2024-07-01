function [B1e, t] = b1_sliceselect(Gz, dz, FA, T, dt, gamma, options)
arguments
    Gz
    dz
    FA
    T
    dt = 4e-3
    gamma = 267.5153151e6
    options.pulse = 'sinc'
    options.filter = 'hann'
    options.Gz_units = 'mT/m'
    options.dz_units = 'mm'
    options.FA_units = 'deg'
    options.t_units = 'ms'
    options.B1e_units = 'mT'
end
% Adjust input units
switch options.Gz_units
    case 'mT/m'
        Gz = Gz*1e-3;
    case 'G/cm'
        Gz = Gz*1e-2;
    case 'T/m' % do nothing
    otherwise
        error('Unrecognized units for Gz: %s', options.Gz_units)
end
switch options.dz_units
    case 'mm'
        dz = dz*1e-3;
    case 'cm'
        dz = dz*1e-2;
    case 'm' % do nothing
    otherwise
        error('Unrecognized units for dz: %s', options.dz_units)
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
        T = T*1e-3;
        dt = dt*1e-3;
    case 'us'
        T = T*1e-6;
        dt = dt*1e-6;
    case 's' % do nothing
    otherwise
        error('Unrecognized units for t: %s', options.t_units)
end

df = gamma*Gz*dz/(2*pi);
t = (0:dt:T); % Has N+1 points, truncate later

switch options.pulse
    case 'sinc'
        B1e = sinc(df*(t - T/2));
    case 'hard'
        B1e = ones(size(t));
    otherwise
        error('Unrecognized pulse type: %s', options.pulse)
end

switch options.filter
    case 'hann'
        B1e = B1e.*( hann(length(B1e))' );
    case 'none' % do nothing
    otherwise
        error('Unrecognized filter type: %s', options.filter)
end

% Truncate extra time point
t = t(1:end-1);
B1e = B1e(1:end-1);

% Compute amplitude for desired flip angle
A = FA/(gamma*sum(B1e.*dt));
B1e = A.*B1e; % Apply amplitude change

% Adjust output units
switch options.B1e_units
    case 'mT'
        B1e = B1e.*1e3;
    case 'T' % do nothing
    otherwise
        error('Unrecognized units for B1e: %s', options.B1e_units)
end
switch options.t_units
    case 'ms'
        t = t*1e3;
    case 'us'
        t = t*1e6;
    case 's' % do nothing
    otherwise
        error('Unrecognized units for t: %s', options.t_units)
end

end