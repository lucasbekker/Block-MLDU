% % Assume D represents a diagonal matrix (missing entries are value 1)
% % Assume that B is an infinite matrix (only possible non-zeros in n x m part)
% %
% % Cases: B stored crs, output crs X = D*B
% %        B stored ccs, output ccs X = B*D
% % One implemenations does 2 different multiplications ...
%
% % A possible check would be nD == |ccsB| - 1. However, such a check is not really required:
% % if a row/col of B has to be scaled which is not stored in B
% % then no action needs to be taken (since the non-stored rows are zero rows).
function [nB, mB, ccsB, iB, vB] = ccs_scale_sparse(nB, mB, ccsB, iB, vB, nD, iD, vD)
for k=1:length(iD)
 l = iD(k);
 if l < length(ccsB) && ccsB(l+1) > ccsB(l)
  c = vD(k);
  vB(ccsB(l):ccsB(l+1)-1) = vB(ccsB(l):ccsB(l+1)-1) * c;
 end;
end;
