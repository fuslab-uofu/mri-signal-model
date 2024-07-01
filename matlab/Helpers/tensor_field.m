function value = tensor_field(value, tensorSize, fieldSize)
%% tf = tensor_field(value, tensorSize, fieldSize)
% Creates higher-order arrays for convenient maniuplation of vectors and
% matrices. Structured so they can be maniuplated similar to a single
% vector or matrix while applying operations to the full collection.
% Simplifies writing vectorized code for large collections of tensors, such
% as a field of spatially-varying vector values
%
% 2023-05-22 Samuel Adams-Tew

switch nargin
    case 0 % No inputs provided
        error('Not enough input arguments!')
    case 1 % Infer tensorSize, fieldSize from value
        sz = size(value);
        if length(sz) == 2
            % Assume that tensors are stored as columns and
            % the field size is determined by the number of
            % column tensors
            tensorSize = [sz(1), 1];
            fieldSize = sz(2);
        else
            % Assume 2nd order tensors given by the first
            % two dimensions and field size is determined
            % by higher dimensions
            tensorSize = sz(1:2);
            fieldSize = sz(3:end);
        end
        % Reshape the input to match the specified TensorField,
        % if needed
        if ~(isequal(size(value), [tensorSize, fieldSize]) ...
                || isequal(size(value), [tensorSize, fieldSize, 1]))
            value = reshape(value, [tensorSize, fieldSize]);
        end
    case 2 % Infer fieldSize from value, tensorSize
        if length(tensorSize) == 1
            % If tensorSize is a scalar, append 1
            tensorSize = [tensorSize, 1];
        end
        if ndims(value) > length(tensorSize)
            % If the value has more dimensions than
            % specified by tensorSize, make the field size
            % the remaining dimensions
            fieldSize = size(value);
            fieldSize(tensorSize ~= 1) = [];
        else
            % Find what field size to reshape to
            fieldSize = numel(value)/prod(tensorSize);
        end

        % Reshape the input to match the specified TensorField,
        % if needed
        if ~(isequal(size(value), [tensorSize, fieldSize]) ...
                || isequal(size(value), [tensorSize, fieldSize, 1]))
            value = reshape(value, [tensorSize, fieldSize]);
        end
    case 3 % All arguments specified
        if length(tensorSize) == 1
            % If tensorSize is a scalar, append 1
            tensorSize = [tensorSize, 1];
        end
        % If shape of value matches those specified by
        % tensorSize and fieldSize, skip other processing in
        % order to permit fast creation of new TensorFields
        % during operations
        if ~(isequal(size(value), [tensorSize, fieldSize]) ...
                || isequal(size(value), [tensorSize, fieldSize, 1]))

            if ndims(value) == length(tensorSize) ...
                    && isequal(size(value), tensorSize)
                % value is a single tensor
                % Create copies to fill the field
                value = value.*ones([tensorSize, fieldSize]);
            end

            % Reshape the input to match the specified
            % TensorField, if needed
            if ~(isequal(size(value), [tensorSize, fieldSize]) ...
                    || isequal(size(value), [tensorSize, fieldSize, 1]))
                value = reshape(value, [tensorSize, fieldSize]);
            end
        end
    otherwise
        error('Too many input arguments!')
end

end