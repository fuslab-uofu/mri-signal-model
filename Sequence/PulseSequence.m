classdef PulseSequence
    %% Object representing an MR pulse sequence
    %
    % Note PulseSequence properties and functions expect times to be given
    % in microseconds [us] unless otherwise specified.
    %
    % ~ Properties ~
    % RF, Gx, Gy, Gz, ADC, blockDuration
    % Dependent properties: startTimes, endTimes, eventTimes, numEvents,
    % numSamples
    %
    % ~ Constructor ~
    % self = PulseSequence(RF, Gx, Gy, Gz, ADC, blockDuration)
    %
    % ~ Instance methods ~
    % [dt, B1, G, sampleComb] = get_event(self, eventNum, dt_max, dtUnits)
    %
    % ~ Static methods ~
    % grad = gradient_specs(startTime, amplitude, rampUpTime, duration, rampDownTime)
    % rf = rf_specs(startTime, type, FA, duration, phase)
    % adc = adc_specs(startTime, numSamples, duration)
    % 
    % Use 'help PulseSequence/<property or method name>' for property and
    % method documentation.
    %
    %% 2023-05-12 Samuel Adams-Tew

    properties
        % RF - Specifications for all RF excitation events as a table with
        % six columns:
        % * startTime [us]: Start time relative to the beginning of this
        % event block
        % * type (string): Type of pulse applied ('hard')
        % * FA [deg]: Flip angle of the pulse
        % * duration [us]: How long the pulse lasts
        % * phase [deg]: Phase angle on transverse plane 0=x-axis,
        % 90=y-axis
        % * endTime [us]: End time relative to the beginning of this event
        % block.
        RF table

        % Gx - Specifications for all Gx events as a table with six
        % columns, following Siemens IDEA convention for duration, ramp up
        % time, and ramp down time of trapezoidal gradient:
        % * startTime [us]: Start time relative to the beginning of this
        % event block
        % * amplitude [mT/m]: Max amplitude of gradient along x.
        % * rampUpTime [us]: Time to reach amplitude from zero.
        % * duration [us]: Time from start of the ramp up to start of ramp
        % down.
        % * rampDownTime [us]: Time to ramp down from amplitude to zero.
        % * endTime [us]: End time relative to the beginning of this event
        % block.
        Gx table

        % Gy - Specifications for all Gy events as a table with six
        % columns, following Siemens IDEA convention for duration, ramp up
        % time, and ramp down time of trapezoidal gradient:
        % * startTime [us]: Start time relative to the beginning of this
        % event block
        % * amplitude [mT/m]: Max amplitude of gradient along y.
        % * rampUpTime [us]: Time to reach amplitude from zero.
        % * duration [us]: Time from start of the ramp up to start of ramp
        % down.
        % * rampDownTime [us]: Time to ramp down from amplitude to zero.
        % * endTime [us]: End time relative to the beginning of this event
        % block.
        Gy table

        % Gz - Specifications for all Gz events as a table with six
        % columns, following Siemens IDEA convention for duration, ramp up
        % time, and ramp down time of trapezoidal gradient:
        % * startTime [us]: Start time relative to the beginning of this
        % event block
        % * amplitude [mT/m]: Max amplitude of gradient along z.
        % * rampUpTime [us]: Time to reach amplitude from zero.
        % * duration [us]: Time from start of the ramp up to start of ramp
        % down.
        % * rampDownTime [us]: Time to ramp down from amplitude to zero.
        % * endTime [us]: End time relative to the beginning of this event
        % block.
        Gz table

        % ADC - Specifications for all ADC (readout) events as a table with
        % four columns:
        % * startTime [us]: Beginning of measurement readout window
        % realtive to the beginning of this event block.
        % * numSamples: Number of samples to collect during the measurement
        % window. Samples are evenly spaced by numSamples/duration and
        % centered within the readout window.
        % * duration [us]: Length of the readout window.
        % * endTime [us]: End time relative to beginning of this event
        % block
        ADC table

        % blockDuration - Scalar value defining the total duration of this
        % event block. Determines amount of time to apply free precession
        % to the system after the last RF or gradient event before
        % simulating another event block.
        blockDuration (1, 1) double {mustBeNonnegative}
    end
    properties (Dependent)
        % startTimes - List of all unique start times for events in this
        % sequence
        startTimes

        % endTimes - List of all unique end times for events in this
        % sequence
        endTimes

        % eventTimes - List of all unique event times in this sequence,
        % including the end of the block
        eventTimes

        % numEvents - Number of simulation events required to simulate the
        % full event block specified by the sequence.
        numEvents

        % numSamples - The number of ADC readout samples collected during a
        % single repetition of the sequence.
        numSamples
    end

    methods
        %% Constructor
        function self = PulseSequence(RF, Gx, Gy, Gz, ADC, blockDuration)
            %% self = PulseSequence(RF, Gx, Gy, Gz, ADC)
            % Default constructor for PulseSequence object, setting each
            % property to the corresponding input value
            %
            % See also GRADIENT_SPECS, RF_SPECS, ADC_SPECS
            %
            %% 2023-05-12 Samuel Adams-Tew

            % Set each property of this instance to the provided input
            % value. If no value is provided, use an empty specification
            % table.
            self.RF = RF;
            if isempty(Gx)
                self.Gx = PulseSequence.gradient_specs([], [], [], [], []);
            else
                self.Gx = Gx;
            end
            if isempty(Gy)
                self.Gy = PulseSequence.gradient_specs([], [], [], [], []);
            else
                self.Gy = Gy;
            end
            if isempty(Gz)
                self.Gz = PulseSequence.gradient_specs([], [], [], [], []);
            else
                self.Gz = Gz;
            end
            self.ADC = ADC;
            self.blockDuration = blockDuration;
        end

        %% Instance methods
        function [dt, B1, G, sampleComb] = get_event(self, eventNum, dt_max, dtUnits)
            %% [dt, B1, G, sampleComb] = get_event(self, eventNum, dt_max, dtUnits)
            % Generates waveforms and simulation parameters for the
            % requested event
            %
            % ~ Input ~
            % * eventNum: Which event to get simulation data for. Must be
            % at least 1 and at most self.numEvents.
            % * dt_max (optional): Maximum timestep permitted during the
            % simulation. The actual timestep used to generate waveforms
            % and returned to be passed to solvers is determined by the
            % types of fields present during the step and the duration of
            % the current event, but will never be larger than dt_max.
            % Default is Inf.
            % * dtUnits (optional): Units for the provided dt_max and
            % returned dt. Default is microseconds, 'us'.
            %
            % ~ Output ~
            % * dt (scalar): Timestep used to generate waveforms for B1 and
            % G.
            % * B1 (1,T): Demodulated B1 (envelope function for
            % on-resonance excitation). Real and imaginary parts correspond
            % to the x- and y-axes of the rotating frame, respectively.
            % * grad (3,T): Gradient field waveforms along x-, y-, and
            % z-axes.
            % * sampleComb (1,T): If ADC is on for this event, times at
            % which to record samples. Samples are evenly spaced by
            % duration/numSamples and centered in the readout window.
            % Default is false(1, T), which saves instructs the solver not
            % to save out any intermediate magnetization states.
            %
            %% 2023-05-12 Samuel Adams-Tew
            arguments
                self
                eventNum
                dt_max = Inf
                dtUnits = 'us'
            end
            if ~isequal(dtUnits, 'us'); dt_max = convert_units(dt_max, dtUnits, 'us'); end

            % This event takes place during the right-open interval
            %   [startTime, endTime)
            startTime = self.eventTimes(eventNum);
            endTime = self.eventTimes(eventNum + 1);

            %%% Initialize components that are part of this event

            % Find events that:
            % start at or after startTime and before endTime
            % - OR -
            % start before startTime and end after startTime
            eventRF = find(...
                (self.RF.startTime >= startTime & self.RF.startTime < endTime) ...
                | (self.RF.startTime < startTime & self.RF.endTime > startTime));
            eventGx = find(...
                (self.Gx.startTime >= startTime & self.Gx.startTime < endTime) ...
                | (self.Gx.startTime < startTime & self.Gx.endTime > startTime));
            eventGy = find(...
                (self.Gy.startTime >= startTime & self.Gy.startTime < endTime) ...
                | (self.Gy.startTime < startTime & self.Gy.endTime > startTime));
            eventGz = find(...
                (self.Gz.startTime >= startTime & self.Gz.startTime < endTime) ...
                | (self.Gz.startTime < startTime & self.Gz.endTime > startTime));
            eventADC = find(...
                (self.ADC.startTime >= startTime & self.ADC.startTime < endTime) ...
                | (self.ADC.startTime < startTime & self.ADC.endTime > startTime));

            % Choose timestep
            dur = endTime - startTime; % [us]
            if isempty(eventRF) && isempty(eventGx) && isempty(eventGy) && isempty(eventGz) && isempty(eventADC)
                dt_max = min(dur, dt_max); % If no applied fields, this event can be a single step
            else
                dt_max = min(10, dt_max); % Apply maximum 10 us steps
            end

            % Enable ADC if appropriate
            if isempty(eventADC)
                % Compute dt that divides evenly
                dt = dur/ceil(dur/dt_max);
                % Create time base
                t = startTime:dt:endTime; % [us] Closed interval
                t = t(1:end-1); % [startTime, endTime)
                % No ADC events
                sampleComb = false(size(t));
            else
                % Get ADC specs
                spec = self.ADC(eventADC, :);
                samplePeriod = dur ./ spec.numSamples; % [us]
                % Determine when samples should be collected
                ADCtimes = startTime + samplePeriod*( (1:(dur/samplePeriod)) - 0.5);

                % Compute an initial sampling rate
                dt = dur/ceil(2*dur/samplePeriod); % [us]
                % Require that dt <= dt_max
                if dt > dt_max
                    dt = dt/ceil(dt/dt_max);
                end

                % Create time base
                t = startTime:dt:endTime; % [us] Closed interval
                t = t(1:end-1); % [startTime, endTime)

                % Find times t closest to sample times 
                % (necessary due to rounding)
                [~, idx] = min(abs(t' - ADCtimes));
                % Create sampling comb
                sampleComb = false(size(t));
                sampleComb(idx) = true;
            end

            % Initialize magnetization components
            if isempty(eventRF)
                B1 = zeros(size(t)); % No RF event
            else
                % Specs for this RF event
                % Start time, type, FA [deg], duration, end time
                spec = self.RF(eventRF, :);
                % Generate waveform of specified type
                if isequal(spec.type{1}, 'hard')
                    [B1, t_B1] = b1_hardpulse(spec.FA, spec.duration, dt, 't_units', 'us');
                    t_B1 = t_B1 + spec.startTime;
                    % Apply specified phase
                    B1 = B1 * exp(1i*180/pi * spec.phase);
                end
                % Choose only interval covered by this event
                B1 = B1(t_B1 >= startTime & t_B1 < endTime);
            end
            if isempty(eventGx)
                gradX = zeros(size(t)); % No Gx event
            else
                spec = self.Gx(eventGx, :); % Specs for this event
                % Generate waveform
                [gradX, t_Gx] = gradient_trap(spec.amplitude, spec.rampUpTime, spec.duration, spec.rampDownTime, dt);
                t_Gx = t_Gx + spec.startTime; % Shift time base to specified start
                % Choose only interval covered by this event
                gradX = transpose(gradX(t_Gx >= startTime & t_Gx < endTime));
            end
            if isempty(eventGy)
                gradY = zeros(size(t)); % No Gy event
            else
                spec = self.Gy(eventGy, :); % Specs for this event
                % Generate waveform
                [gradY, t_Gy] = gradient_trap(spec.amplitude, spec.rampUpTime, spec.duration, spec.rampDownTime, dt);
                t_Gy = t_Gy + spec.startTime; % Shift time base to specified start
                % Choose only interval covered by this event
                gradY = transpose(gradY(t_Gy >= startTime & t_Gy < endTime));
            end
            if isempty(eventGz)
                gradZ = zeros(size(t)); % No Gz event
            else
                spec = self.Gz(eventGz, :); % Specs for this event
                % Generate waveform
                [gradZ, t_Gz] = gradient_trap(spec.amplitude, spec.rampUpTime, spec.duration, spec.rampDownTime, dt);
                t_Gz = t_Gz + spec.startTime; % Shift time base to specified start
                % Choose only interval covered by this event
                gradZ = transpose(gradZ(t_Gz >= startTime & t_Gz < endTime));
            end
            G = [gradX; gradY; gradZ];

            % Convert time base to desired units
            if ~isequal(dtUnits, 'us'); dt = convert_units(dt, 'us', dtUnits); end
        end

        function plot(self, startEvent, endEvent)
            %% plot(self)
            % Creates a plot of the pulse sequence
            % 
            %% 2023-05-22 Samuel Adams-Tew

            if ~exist('startEvent', 'var')
                startEvent = 1;
            end
            if ~exist('endEvent', 'var')
                endEvent = self.numEvents;
            end

            B1 = [];
            gradX = [];
            gradY = [];
            gradZ = [];
            sampleComb = [];
            t = 0;

            for eventNum = startEvent:endEvent
                % Get waveforms for this event
                [dt, b1tmp, gtmp, samptmp] = self.get_event(eventNum);
                t = [t, self.eventTimes(eventNum) + dt.*(1:length(b1tmp))];
                B1 = [B1, b1tmp];
                gradX = [gradX, gtmp(1, :)];
                gradY = [gradY, gtmp(2, :)];
                gradZ = [gradZ, gtmp(3, :)];
                sampleComb = [sampleComb, samptmp];
            end

            B1 = [B1, 0];
            gradX = [gradX, 0];
            gradY = [gradY, 0];
            gradZ = [gradZ, 0];
            sampleComb = [sampleComb, 0];

            figure(); plot_sequence(t*1e-3, B1, [gradX; gradY; gradZ], sampleComb);
            xlabel('Time (ms)')
        end
    end

    %% Dependent property getters
    methods
        function startTimes = get.startTimes(self)
            startTimes = unique([self.RF.startTime; self.Gx.startTime; ...
                self.Gy.startTime; self.Gz.startTime; self.ADC.startTime]);
        end

        function endTimes = get.endTimes(self)
            endTimes = unique([self.RF.endTime; self.Gx.endTime; ...
                self.Gy.endTime; self.Gz.endTime; self.ADC.endTime]);
        end

        function eventTimes = get.eventTimes(self)
            eventTimes = unique([self.startTimes; self.endTimes; self.blockDuration]);
        end

        function numEvents = get.numEvents(self)
            numEvents = length(self.eventTimes) - 1;
        end

        function numSamples = get.numSamples(self)
            numSamples = sum(self.ADC.numSamples);
        end
    end

    methods (Static)
        %% Static methods
        function grad = gradient_specs(startTime, amplitude, rampUpTime, duration, rampDownTime)
            %% grad = gradient_specs(startTime, amplitude, rampUpTime, duration, rampDownTime)
            % Creates a table with values specifying a gradient waveform
            % following Siemens IDEA convention.
            %
            % ~ Input ~
            % * startTime [us]: Start time relative to the beginning of
            % this event block
            % * amplitude [mT/m]: Max amplitude of gradient along x.
            % * rampUpTime [us]: Time to reach amplitude from zero.
            % * duration [us]: Time from start of the ramp up to start of
            % ramp down.
            % * rampDownTime [us]: Time to ramp down from amplitude to
            % zero.
            %
            % ~ Output ~
            % * grad (table): Table with all specifications required to
            % create a gradient for a PulseSequence object.
            %
            %% 2023-05-12 Samuel Adams-Tew
            arguments
                startTime (:, 1)
                amplitude (:, 1)
                rampUpTime (:, 1)
                duration (:, 1)
                rampDownTime (:, 1)
            end
            endTime = startTime + duration + rampDownTime;
            grad = table(startTime, amplitude, rampUpTime, duration, rampDownTime, endTime);
        end

        function rf = rf_specs(startTime, type, FA, duration, phase)
            %% rf = rf_specs(startTime, type, FA, duration, phase)
            % Creates a table with values specifying RF excitations.
            %
            % ~ Input ~
            % * startTime [us]: Start time relative to the beginning of
            % this event block
            % * type (string): Type of pulse applied ('hard')
            % * FA [deg]: Flip angle of the pulse
            % * duration [us]: How long the pulse lasts
            % * phase [deg]: Phase angle on transverse plane 0=x-axis,
            % 90=y-axis
            %
            % ~ Output ~
            % * rf (table): Table with all specifications required to
            % create RF events for a PulseSequence object.
            %
            %% 2023-05-12 Samuel Adams-Tew
            arguments
                startTime (:, 1)
                type (:, 1)
                FA (:, 1)
                duration (:, 1)
                phase (:, 1)
            end

            endTime = startTime + duration;
            rf = table(startTime, type, FA, duration, phase, endTime);
        end

        function adc = adc_specs(startTime, numSamples, duration)
            %% adc = adc_specs(startTime, numSamples, duration)
            % Creates a table with values specifying ADC readout events.
            %
            % ~ Input ~
            % * startTime [us]: Beginning of measurement readout window
            % realtive to the beginning of this event block.
            % * numSamples: Number of samples to collect during the
            % measurement window. Samples are evenly spaced by
            % numSamples/duration and centered within the readout window.
            % * duration [us]: Length of the readout window.
            %
            % ~ Output ~
            % * adc (table): Table with all specifications required to
            % create ADC readout events for a PulseSequence object.
            %
            %% 2023-05-12 Samuel Adams-Tew
            arguments
                startTime (:, 1)
                numSamples (:, 1)
                duration (:, 1)
            end

            endTime = startTime + round(duration);
            adc = table(startTime, numSamples, duration, endTime);
        end
    end
end