classdef BlockLDUMatrix
    properties
        Data
    end
    methods
        % CONSTRUCTOR
        function obj = BlockLDUMatrix(data,itterations)
            % Memory allocation for cell array header.
            if nargin > 1
                obj.Data = cell(itterations,1);
            else
                obj.Data = cell(1,1);
            end
            % Insert data into the first element of the cell array.
            if nargin == 1
                obj.Data{1} = data;
            end
        end
    end
end