function spydeep(S,scale)

% Function like spy, but with the ability to visualize the value of the
% non-zero elements of the matrix. "scale" is an optional input argument
% that sets the partition of the data. 'lin' corresponds with a linear
% partition and 'log' with a logarithmic one. The default behaviour tries
% to automatically set the right type of partition.

% Check for empty matrix.
if isempty(S), error('empty matrix'); end

% Convert to COO.
[m,n]   = size(S);
[I,J,V] = find(S);

% Adjust complex matrix.
if ~isreal(V), V = abs(V); end

% Sort for increased performance.
[W,P]   = sort(V);
I       = I(P);
J       = J(P);

% Automatic partition selection.
if nargin < 2
    if (abs(abs(W(end)) - abs(W(1))) < 100)
        scale = 'lin';
    else
        scale = 'log';
    end
end

% Create partitioning.
if strcmp(scale,'lin')
    W_part = linspace(W(1),W(end),9);
elseif strcmp(scale,'log')
    W_part = logspace(W(1),W(end),9);
else
    error('Wrong scale type argument')
end

% Open a new figure.
figure
axis([0 (m + 1) 0 (n + 1)]);
title('Spy plot with depth of matrix S');
xlabel(['nz = ' int2str(nnz(S))]);
hold on

% Select the color range.
cmp = colormap(jet);
color = cmp((1:9:64),:);

% Initialize variables to start the loop.
t      = false(size(W,1),1);
l      = [false(1);true(8,1)];

% Loop over the selected colors.
for i = 1:8
    
    % Find the relevant non-zero elements to color.
    tt = [ false(nnz(t),1); (W(~t) <= W_part(i + 1)) ];
    t = t | tt;
    
    % Plot the relevant non-zero elements with parameters.
    plot(J(tt),I(tt),'marker','.','markersize',14,...
         'linestyle','none','color',color(i,:));  
    
    if sum(tt) == 0, l(i) = false; end
     
end

% Add a legend for reference.
legend(string(W_part(l)))

hold off

end