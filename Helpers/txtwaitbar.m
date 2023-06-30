classdef txtwaitbar

    properties
        message
        barwidth
        progress
    end
    properties (Access = private)
        currentWidth
    end
    properties (Dependent)
        completeElements
    end

    methods
        function waitbar(x, self, msg)
            self.message = msg;
            self.progress = x;
        end

        function close(self)
            fprintf(repelem('\b', self.currentWidth))
        end
    end

    methods (Access = private)
        function str = make_message(self)

        end
        function str = make_bar(self)
            completecount = 

            % Fill in completed elements
            str = repelem(sprintf('\x2588'), )];
            % Blank spaces for pending elements
            str = [str, repelem(' ', self.barwidth - newCompleteCount)];
            % Closing bracket
            str = [str, ']'];
            % Return carriage
            str = [str, '\n'];
        end
    end
end