function weight = spectral_line_shape(f, distribution, options)
%% weight = spectral_line_shape(f, distribution, options)
% Creates NMR spectral line shapes with theoretical properties
%
% T2* is conventionally modeled as, 
%   1/T2* == 1/T2 + 1/T2',
% where T2' is a decay constant describing field inhomogeneities and
% off-resonance effects on a sample.
%
% FWHM for Lorentzian with T2' [in s]: FWHM = 1/(pi*T2')
%
% ~ Input ~
% * f [Hz]: Vector of evenly spaced frequency values to compute the line
% shape weight for.
% * distribution (string): Spectral line shape to use. Default is
% 'Lorentzian'
% 
% ~ Options ~
% * normalize: Whether to normalize the output weight to have unit area.
% Default is true.
% * T2prime: Directly defines the field inhomogeneity coefficient T2'.
% Default is 200 [ms].
% * T2, T2star: If both these parameters are defined, T2prime is computed
% from 1/T2prime = 1/T2star - 1/T2. 
% * T2_units: Units for T2prime, T2, and T2star. Default 'ms'.
%
% ~ Output ~
% * weight: Vector with same size as f, with the relative weight of each
% off-resonance frequency in the NMR signal.
%
%% 2023-06-08 Samuel Adams-Tew
arguments
    f
    distribution = 'Lorentzian'
    options.T2
    options.T2star
    options.T2prime = 200
    options.T2_units = 'ms'
    options.normalize = true
end

% Determine what axis f varies along, sets the axis for spectroscopic line
% shape / distribution
spectAxis = find(size(f) == length(f));

if isfield(options, 'T2star') && isfield(options, 'T2')
    % Compute T2prime from given T2 and T2star values
    options.T2prime = 1/(1/T2star - 1/T2);
end

% Compute the line shape and check that frequency range is large enough to
% accurately model this distribution
switch distribution
    case 'Lorentzian'
        T2p = convert_units(options.T2prime, options.T2_units, 's');
        weight = 2*T2p ./ ( (2*pi*T2p.*f).^2 + 1 );

        if any(max(f)*T2p < 2)
            warning('Maximum frequency max(f) is not large enough to accurately model T2prime=%g ms', options.T2prime)
        end
    otherwise
        error('Unknown distribution for creating off-resonance ensemble')
end

if options.normalize
    % Normalize to have unit area
    weight = weight./sum(weight, spectAxis);
end

end