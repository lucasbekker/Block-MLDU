function [ L, D, U ] = MLDU_2( A, block_sizes )

[m,n] = size(A);
[i,j,v] = find(A);
nz = nnz(A);

S = struct('i',vertcat(i,zeros(nz,1)),...
           'j',vertcat(j,zeros(nz,1)),...
           'v',vertcat(v,zeros(nz,1)),...
           'm',m,'n',n,'nz',nz);

ind_L = false(2*nz,1); ind_D = ind_L; ind_U = ind_L;

shift = [0,block_sizes(1)];

for k = 1:length(block_sizes)
    
    [ind_L_up,ind_D_up,ind_U_up] = Block_Split_local(S,shift(2), ind_L, ind_D, ind_U);
    
    D_inv = inv(sparse((S.i(ind_D_up) - shift(1)),(S.j(ind_D_up) - shift(1)),...
                S.v(ind_D_up),block_sizes(k),block_sizes(k)));
    
    [i,j,v] = find(sparse(S.i(ind_L_up),(S.j(ind_L_up) - shift(1)),S.v(ind_L_up),m,block_sizes(k))*...
                   -D_inv*...
                   sparse((S.i(ind_U_up) - shift(1)),S.j(ind_U_up),S.v(ind_U_up),block_sizes(k),n));
    
    Block_Join_nested;
    Merge_Index_L_nested;
    Merge_Index_D_nested;
    Merge_Index_U_nested;
    
    shift(1) = shift(2);
    shift(2) = shift(2) + block_sizes(k);
               
end

L = sparse(S.i(ind_L),S.j(ind_L),S.v(ind_L),S.m,S.n);
D = sparse(S.i(ind_D),S.j(ind_D),S.v(ind_D),S.m,S.n);
U = sparse(S.i(ind_U),S.j(ind_U),S.v(ind_U),S.m,S.n);
       
    function [ ] = Block_Join_nested( )
    
    l = length(S.i);
    ll = length(i);
    
    if l < (ll + S.nz)
        S.i = vertcat(S.i,zeros(l,1));
        S.j = vertcat(S.j,zeros(l,1));
        S.v = vertcat(S.v,zeros(l,1));
    end
    
    index = S.nz + (1:ll);
    
    S.i(index) = i; S.j(index) = j; S.v(index) = v;
    S.nz = S.nz + ll;
            
    end

    function [ ] = Merge_Index_L_nested( )

    l = length(ind_L);
    ll = length(ind_L_up);

    if l == ll
       ind_L = ind_L | ind_L_up;
    else
        ind_L = vertcat(ind_L,false((ll - l),1)) | ind_L_up;
    end

    end

    function [ ] = Merge_Index_D_nested( )

    l = length(ind_D);
    ll = length(ind_D_up);

    if l == ll
       ind_D = ind_D | ind_D_up;
    else
        ind_D = vertcat(ind_D,false((ll - l),1)) | ind_D_up;
    end

    end

    function [ ] = Merge_Index_U_nested( )

    l = length(ind_U);
    ll = length(ind_U_up);

    if l == ll
       ind_U = ind_U | ind_U_up;
    else
        ind_U = vertcat(ind_U,false((ll - l),1)) | ind_U_up;
    end

    end

end

function [ ind_L, ind_D, ind_U ] = Block_Split_local( S, a, ind_L, ind_D, ind_U )

R = [(ind_L | ind_D | ind_U); true(length(ind_L) - S.nz,1)];
L = (S.j <= a);
U = (S.i <= a);
ind_D = ~R & L & U;
ind_L = L & ~(R | ind_D);
ind_U = U & ~(R | ind_D);

end