function [ ] = Add_field_nested( struct_name, field_name, field_value_name )

% This function can only be used in a nested situation, because the struct 
% itself is not part of the scope of the function.

% Check if the input arguments are sane.
% if ~(isstring(struct_name) && isstring(field_name) && isstring(field_value_name)) 
%     error('struct_name, field_name and field_value_name should be strings')
% end

% Check if the field is already present in the struct.
field_check = eval(sprintf('isfield(%s,"%s");',struct_name,field_name));
if field_check
    error('field is already present')
end

% Perform the addition of the field.
eval(sprintf('%s.%s = %s;',struct_name,field_name,field_value_name));

end