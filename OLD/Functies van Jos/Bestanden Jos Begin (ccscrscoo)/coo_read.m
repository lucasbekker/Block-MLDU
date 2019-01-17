function[n, m, i, j, v] = coo_read(fname)
fid = fopen(fname, 'r');
n = fscanf(fid,'%d', 1);
m = fscanf(fid,'%d', 1);
ne = fscanf(fid,'%d', 1);
i = [];
j = [];
v = [];
for k=1:ne
if mod(k,10000) == 0, fprintf('Entry %d!\n',k); end
 i   = [i; fscanf(fid,'%d',1)];
 j   = [j; fscanf(fid,'%d',1)];
 v   = [v; fscanf(fid,'%e',1)];
end;
fclose(fid);
