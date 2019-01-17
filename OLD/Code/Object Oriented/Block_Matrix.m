classdef Block_Matrix
    
    properties
    
        Data
        Part
        Sym
        Size
        Type
        
    end
    
    methods
        
        % CONSTRUCTOR
        function obj = Block_Matrix( partition, type, size, varargin )
            
            % Input parsing.
            p = inputParser;
            chk = @(x,y) any(validatestring(x,y));
                        
            % Specify parameters, their default value and check if they are
            % correct.
            p.addParameter('sym','false',@(x) chk(x,{'false','true'}));
            p.addParameter('alloc','false',@(x) chk(x,{'false','true'}));
            p.addRequired('partition',@isnumeric);
            p.addRequired('type',@(x) chk(x,{'L','D','U'}));
            p.addRequired('size',@isnumeric);
            p.parse(partition,type,size,varargin{:});
                        
            % Store information.
            obj.Part = p.Results.partition;
            obj.Sym  = p.Results.sym;
            obj.Size = p.Results.size;
            obj.Type = p.Results.type;
                        
            % Calculate memory allocation.
            l = length(obj.Part) - double( ...
                (obj.Size(1) == obj.Size(2) && obj.Type ~= 'D') ||  ...
                (obj.Size(1) > obj.Size(2)  && obj.Type == 'U') ||  ...
                (obj.Size(1) < obj.Size(2)  && obj.Type == 'L'));
            
            % Check for additional memory allocation to append data.
            if strcmp(p.Results.alloc,'false')
                l = [1, l];
            end
            
            % Initialize Data as a cell array.
            obj.Data = cell(l);
            
        end
        
    end
    
end