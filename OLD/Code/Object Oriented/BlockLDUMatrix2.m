classdef BlockLDUMatrix2
    methods
        % CONSTRUCTOR
        function obj = BlockLDUMatrix(data,itterations)
            % Memory allocation for cell array header.
            if nargin > 1
                obj = cell(itterations,1);
            else
                obj = cell(1,1);
            end
            % Insert data into the first element of the cell array.
            if nargin == 1
                obj{1,1} = data;
            end
        end
    end
end