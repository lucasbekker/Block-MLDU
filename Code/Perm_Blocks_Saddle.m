function [ P ] = Perm_Blocks_Saddle( n , m )
%--------------------------------------------------------------------------
%MICRO_BLOCKS: Returns a permutation vector P of length n+m such that the
%           transformed saddle point matrix Xp can be partitioned into the 
%           following block structures:
%                     | 2x2 blocks | 2x1 blocks |
%           Xp(P,P) = |------------|------------|
%                     | 1x2 blocks | 1x1 blocks |
%--------------------------------------------------------------------------

P = [];

for i = 1:m
    P = [ P, [i, n+i] ];
end

P = [ P, m+1:n ];

end