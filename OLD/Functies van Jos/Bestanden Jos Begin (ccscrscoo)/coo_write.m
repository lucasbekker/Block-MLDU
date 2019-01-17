% % copyright Joseph M.L. Maubach, 2014 -- reproduction prohibited
% % For REAL-VALUED matrices -- same as matrix market format
% % coo_write outputs the matrix in the coo format:
% % (lines with *... are not output):
%
% % [n, m, length(i);
% %  *for k=1:length(i)
% %   i(k), j(k), v(k);
% %  *end
% % ]
%
function [] = coo_write(n, m, i, j, v, fname)
if (nargin < 6) || (isempty(fname)), fid = 1; else fid = fopen(fname,'w'); end;
fprintf(fid,'%d %d %d\n', n, m, length(i));
for k=1:length(i)
 fprintf(fid,'%d %d ', i(k), j(k)); if abs(v(k) - round(v(k))) == 0, fprintf(fid,'%d. ',round(v(k))), else fprintf(fid,'%22.16e ',v(k)); end; fprintf(fid,'\n');
end
if fid ~=1 fclose(fid); end;
