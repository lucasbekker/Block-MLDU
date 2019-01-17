function [] = ccs_write(n, m, ccs, i, v, fname, format, echo)
if (nargin < 6) || (isempty(fname)), fid = 1; else fid = fopen(fname,'w'); end;
if (nargin < 7) || isempty(format), format = 0; end;
if (nargin < 8) || isempty(echo), echo = 0; end;
if format == 0
 if echo > 0, fprintf('variables: n, m, #stored_cols+1, nnz:\n'); end;
 fprintf(fid,'%d %d %d %d\n', n, m, length(ccs), ccs(end)-1);
 if echo > 0, fprintf('ccs:\n'); end;
 for k=1:length(ccs), fprintf(fid,'%d ',ccs(k)); end; fprintf(fid,'\n');
 if echo > 0, fprintf('i:\n'); end;
 for k=1:length(i), fprintf(fid,'%d ',i(k)); end; fprintf(fid,'\n');
 if echo > 0, fprintf('v:\n'); end;
 for k=1:length(v), if abs(v(k) - round(v(k))) == 0, fprintf(fid,'%d. ',round(v(k))); else fprintf(fid,'%22.16e ',v(k)); end; end; fprintf(fid,'\n');
else
% % harwell software library routine hsl_ma77 (example hsl_ma77ds.data) uses this as input format.
%
% % instead of full ccs output (of ccsOccs)
%
% % [m length(i);
% %  ccs';
% %  i';
% %  v']
%
% % one outputs (lines with *... are not output):
%
% % [m;
% %  *for k=1:m
% %   ccs(k+1)-ccs(k);			"amount of rowindices i_{k_l} for column k";
% %   i(ccs(k):ccs(k+1)-1)';		"rowindices i_{k_l}for column k"
% %   v(ccs(k):ccs(k+1)-1)';		"value A_(i_{k_l},k)"
% %  *end
% %  ]
%
 fprintf(fid,'%d %d %d\n',n,m,length(ccs));
 for k=1:length(ccs)-1
  fprintf(fid,'%d\n',ccs(k+1)-ccs(k));
  for l=ccs(k):ccs(k+1)-1, fprintf(fid,'%d ',i(l)); end; fprintf(fid,'\n');
  for l=ccs(k):ccs(k+1)-1, if abs(v(l) - round(v(l))) == 0, fprintf(fid,'%d. ',round(v(l))); else fprintf(fid,'%22.16e ',v(l)); end; end; fprintf(fid,'\n');
 end
end
if fid ~=1 fclose(fid); end;
