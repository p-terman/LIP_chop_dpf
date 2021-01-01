% mex function GZRead
%
% /* This c code is intended for use as a mex file to open, navigate, read, and
%  * close .gz files.  This is done by passing an fid array back and forth between
%  * this code and Matlab every time this function is called.  These elements are stored
%  * as uint32 (unsigned long), and are as follows:
%  * gzfid[0] - pointer to the FILE object
%  * gzfid[1] - pointer to the BZFILE object
%  * gzfid[2] - latest gzerror code (should be converted to int32 to read,
%  *            codes for different values are listed in gzlib.h
%  * gzfid[3] - current position in file (first byte is zero)
%  * gzfid[4] - max buffer to read when scanning through file -- used for setting
%  *            position only, does not limit the size of buffer read out in 'read'
%  *            mode.
%  *
%  * The code should be used as follows:
%  * To open file:
%  *    gzfid = GZRead('open', filename, uint32(max_buffer_size))
%  * To set position in file:
%  *    gzfid = GZRead('setpos', gzfid, uint32(new_position), [suppress eof warning])
%  *                       (suppress eof warning is an optional third
%  *                        argument, default is false, set to true if you
%  *                        don't want a warning if the position you set is
%  *                        past the end of the file)
%  * To read compressed data from current position:
%  *    [gzfid data] = GZRead('read', gzfid, uint32(num_elements_to_read), 'datatype', swapbytes)
%  *     (note, swapbytes needs to be boolean, meaning you should enter a
%  *      boolean variable, or type out true or false -- typing 0 or 1 is
%  *      interpreted as a double when passed to the mex file)
%  * To close file
%  *    gzfid = GZRead('close', gzfid)
%  *
%  * When reading data, the data is returned as a column-vector of desired datatype.
%  *
%  * CED, 06/14/06
%  * Release Version 1.01
%  * Version 1.0 -> 1.01, fixed "help", no recompile needed
%  * */


