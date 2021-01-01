% mexfunction [data count_read bytes_read] = MEM_fread_fast(fid, count, precision)
%
% A mexed version of MEM_fread that may be used when count is finite (not
% to end of file), skip is zero, endianness is default (set when file was
% opened), and file has loaded TO END OF FILE (ie, no calls to
% MEM_freopen).  If these conditions are not met, all return outputs are
% empty arrays.  Also note, three outputs ARE necessary.  ALSO:  filepos
% DOES NOT ADVANCE when using this function -- that's why bytes_read is
% returned, so you can do it yourself if you like (call
% MEM_fseek(fid,bytes_read,'cof')
%
% CED, 08/05/06

disp('MEM_fread_fast is a mex function -- please compile it (no other libraries needed, just mex MEM_fread_fast.c');