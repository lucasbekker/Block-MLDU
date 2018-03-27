function [ output_test ] = test_block_MLDU( MLDU_function, test_routine, inarg_MLDU, inarg_test )

output_test_str = '[';
for k = 1:nargout(MLDU_function)
    output_test_str = sprintf('%s\x2Coutput_test_%i',output_test_str,k);    
end
output_test_str = [output_test_str,']'];

eval(sprintf('%s\x3D%s%s',output_test_str,test_routine,)

end

function [ output_string ] = cell2inarg_local( input_cell )

% Determine dimensions of input_cell and perform check.
size_cell = size(input_cell);
if min(size_cell) ~= 1
    error('input cell should be one dimensional.')
end

% Loop over the elements and concatenate them.
output_string = string(input_cell{1});
for k = 2:max(size_cell)
    output_string = sprintf('%s\x2C%s',output_string,string(input_cell{k}));
end

% Format the output string to include input format brackets.
output_string = sprintf('(%s)',output_string);

end

