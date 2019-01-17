% % read output written in ccs format (written by ccs2ccs())
%
function[n, m, ccs, i, v] = ccs_read(fname, format)
if (nargin < 2) || isempty(format), format = 0; end;
fid = fopen(fname, 'r');
if format == 0
 n = fscanf(fid,'%d', 1);
 m = fscanf(fid,'%d', 1);
 nccs = fscanf(fid,'%d', 1);
 nnz = fscanf(fid,'%d', 1);
 ccs=zeros(nccs,1);
 i = zeros(nnz,1);
 v = zeros(nnz,1);
 ccs = fscanf(fid,'%d',nccs);
 i   = fscanf(fid,'%d',nnz);
 v   = fscanf(fid,'%e',nnz);
else
 n = fscanf(fid,'%d', 1);
 m = fscanf(fid,'%d', 1);
 s = fscanf(fid,'%d', 1);
 ccs=[];
 i = [];
 v = [];
 for k=1:s
  ccs = [ccs; fscanf(fid,'%d',1)];
  i   = [i;   fscanf(fid,'%d',ccs(end))];
  v   = [v;   fscanf(fid,'%e',ccs(end))];
 end;
 ccs = cumsum([1;ccs]);
end
fclose(fid);
