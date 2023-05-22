function val = convert_units(val, inUnits, outUnits)

%%% Setup input units
switch inUnits
    case 's' % Time: convert inUnits -> seconds
        quantity = 'time';
    case 'ms'
        val = val*1e-3; % 1 ms = 1e-3 s
        quantity = 'time';
    case 'us'
        val = val*1e-6; % 1 us = 1e-6 s
        quantity = 'time';
    case 'ns'
        val = val*1e-9; % 1 ns = 1e-9 s
        quantity = 'time';
    case 'ps'
        val = val*1e-12; % 1 ps = 1e-12 s
        quantity = 'time';
        
    case 'm' % Length: convert inUnits -> meters
        quantity = 'length';
    case 'cm'
        val = val*1e-2; % 1 cm = 1e-2 m
        quantity = 'length';
    case 'mm'
        val = val*1e-3; % 1 mm = 1e-3 m
        quantity = 'length';
    case 'um'
        val = val*1e-6; % 1 um = 1e-6 m
        quantity = 'length';
        
    case 'T' % Magnetic field strength: convert inUnits -> tesla
        quantity = 'magnetic field';
    case 'mT'
        val = val*1e-3; % 1 mT = 1e=3 T
        quantity = 'magnetic field';
    case 'uT'
        val = val*1e-6; % 1 uT = 1e-6 T
        quantity = 'magnetic field';
    case 'nT'
        val = val*1e-9; % 1 nT = 1e-9 T
        quantity = 'magnetic field';
    case 'pT'
        val = val*1e-12; % 1 nT = 1e-12 T
        quantity = 'magnetic field';
    case 'G'
        val = val*1e-4; % 1 G = 1e-4 T
        quantity = 'magnetic field';

    case 'T/m' % Gradient strength: convert inUnits -> T/m
        quantity = 'gradient strength';
    case 'mT/m'
        val = val*1e-3; % 1 mT/m = 1e-3 T/m
        quantity = 'gradient strength';
    case 'G/cm'
        val = val*1e-2; % 1 G/cm = 1e-2 T/m
        quantity = 'gradient strength';

    otherwise
        error('Unrecognized input units %s', inUnits)
end

% Process out units and make sure it is compatible with input units
switch quantity
    case 'time'
        switch outUnits
            case 's' % do nothing, already in seconds
            case 'ms' 
                val = val*1e3; % 1 s = 1e3 ms
            case 'us'
                val = val*1e6; % 1 s = 1e6 ms
            case 'ns'
                val = val*1e9; % 1 s = 1e9 ms
            case 'ps'
                val = val*1e12; % 1 s = 1e12 ms
            otherwise
                error('Units incompatible with %s (%s) requested: %s', inUnits, quantity, outUnits)
        end
    case 'length'
        switch outUnits
            case 'm' % do nothing, already in meters
            case 'cm'
                val = val*1e2; % 1 m = 1e2 cm
            case 'mm'
                val = val*1e3; % 1 m = 1e3 mm
            case 'um'
                val = val*1e6; % 1m = 1e6 um
            otherwise
                error('Units incompatible with %s (%s) requested: %s', inUnits, quantity, outUnits)
        end
    case 'magnetic field'
        switch outUnits
            case 'T' % do nothing, already in tesla
            case 'mT'
                val = val*1e3; % 1 T = 1e3 mT
            case 'uT'
                val = val*1e6; % 1 T = 1e6 uT
            case 'nT'
                val = val*1e9; % 1 T = 1e9 nT
            case 'pT'
                val = val*1e12; % 1 T = 1e12 pT
            case 'G'
                val = val*1e4; % 1 T = 1e4 G
            otherwise
                error('Units incompatible with %s (%s) requested: %s', inUnits, quantity, outUnits)
        end
    case 'gradient strength'
        switch outUnits
            case 'T/m' % do nothing, already in T/m
            case 'mT/m'
                val = val*1e3; % 1 T/m = 1e3 mT/m
            case 'G/cm'
                val = val*1e2; % 1 T/m = 1e2 G/cm
            otherwise
                error('Units incompatible with %s (%s) requested: %s', inUnits, quantity, outUnits)
        end
    otherwise
        error('Unknown unit type %s', quantity)
end

end