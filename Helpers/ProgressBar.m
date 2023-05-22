classdef ProgressBar < handle
    %% ProgressBar
    % Creates an object for easily displaying a progress bar for iterative
    % processes. Tracks the time elapsed and time remaining, as well as
    % iterations completed.
    %
    % ~ Properties ~
    % message, nIter, isAppendFractionEnabled, isTimeRemainingEnabled, 
    % currIter
    % Dependent properties: telapsed, tremain
    % Private properties: w, tstart
    %
    % ~ Constructor ~
    % obj = ProgressBar(msg, nIter, options)
    %
    % ~ Instance Methods ~
    % start(obj)
    % iter(obj, currIter)
    % close(obj)
    % update_waitbar(obj) (Access = private) 
    %
    % Use 'help ProgressBar/<property name or method name>' for property
    % and method documentation
    %
    %% 2023-05-18 Samuel Adams-Tew

    properties
        message % Message describing the process being completed
        % Boolean flag indicating whether to append the fraction of 
        % iterations completed to the message
        isAppendFractionEnabled
        % Boolean flag indicating whether to display the estimated time 
        % remaining in the message
        isTimeRemainingEnabled
        currIter % Current iteration
        nIter % Total number of iterations
    end
    properties (Access = private)
        w = [] % waitbar handle
        tstart = [] % time process was started
    end
    properties (Dependent)
        telapsed % in seconds
        tremain % in seconds
    end

    methods
        %% Getters
        function telapsed = get.telapsed(obj)
            % Get the elapsed time in seconds
            telapsed = toc(obj.tstart);
        end
        function tremain = get.tremain(obj)
            % Get the estimated time remaining in seconds
            tremain = toc(obj.tstart)/(obj.currIter - 1)*(obj.nIter - obj.currIter + 1);
        end
        %% Setters
        function set.currIter(obj, value)
            % Set the current iteration
            if value < 0
                % Error if value is negative
                error('Value of currIter must be non-negative')
            else
                obj.currIter = value;
                % If the current iteration is greater than 0, update the waitbar
                if value > 0
                    obj.update_waitbar();
                end
            end
        end
        %% Constructor
        function obj = ProgressBar(msg, nIter, options)
            %% obj = ProgressBar(msg, nIter, options)
            % Creates a progress bar object for tracking and predicting the
            % time remaining for iterative processes
            %
            % ~ Inputs ~ 
            % - msg: Message describing the process being completed
            % - nIter: Total number of iterations
            % 
            % ~ Options ~
            % - startNow: starts the timer and displays the waitbar
            % immediately after initializing the object. Default is true.
            % - currIter: Which iteration to start counting on. Default is
            % 0.
            % - isAppendFractionEnabled: Boolean flag indicating whether to 
            % append the fraction of iterations completed to the message.
            % Default is true;
            % - isTimeRemainingEnabled: Boolean flag indicating whether to 
            % display the estimated time remaining in the message. Default
            % is true.
            %
            % 2023-05-18 Samuel Adams-Tew
            arguments
                msg
                nIter
                options.startNow = true
                options.currIter = 0
                options.isAppendFractionEnabled = true
                options.isTimeRemainingEnabled = true
            end

            % Assign values to properties
            obj.message = msg;
            obj.nIter = nIter;
            obj.currIter = options.currIter;
            obj.isAppendFractionEnabled = options.isAppendFractionEnabled;
            obj.isTimeRemainingEnabled = options.isTimeRemainingEnabled;

            % If the startNow option is set to true, start the progress bar
            if options.startNow
                obj.start();
            end
        end
        %%
        function start(obj)
            % Start timing the process and display the progress bar
            obj.tstart = tic;
            obj.w = waitbar(0, [obj.message, ' initializing']);
        end
        function iter(obj, currIter)
            % Increment the current iteration OR update to iteration
            % currIter
            if nargin > 1
                % Set the current iteration to the value passed in
                obj.currIter = currIter;
            else
                % Increment the current iteration by 1
                obj.currIter = obj.currIter + 1;
            end

            if obj.currIter > obj.nIter
                % If the current iteration is greater than the total number 
                % of iterations, delete the progress bar
                delete(obj);
            end
        end
        function close(obj)
            % Close the progress bar
            close(obj.w);
        end
    end
    methods (Access = private)
        function update_waitbar(obj)
            % Helper function that updates the waitbar message

            % If the isAppendFractionEnabled flag is set to true, append 
            % the fraction of iterations completed to the message
            if obj.isAppendFractionEnabled
                msg = [obj.message sprintf(' %d/%d', obj.currIter, obj.nIter)];
            else
                msg = obj.message;
            end

            % If the isTimeRemainingEnabled flag is set to true and the 
            % estimated time remaining is less than infinity, append the 
            % estimated time remaining to the message
            if obj.isTimeRemainingEnabled && obj.tremain < inf
                msg = {msg, sprintf('Est. %.1f min remaining', obj.tremain/60)};
            end

            % Update the waitbar with the new message
            waitbar(obj.currIter/obj.nIter, obj.w, msg);
        end
    end
end