% mex function TypeCastAndSwap
%
% /* This c code is intended for use as a mex file to perform typecasting and
%  * byte swapping -- it may be replaced by the typecast and swapbytes functions
%  * in matlab 7.1 (R14SP3) and later.
%  *
%  * Usage is:
%  * data = TypeCastAndSwap(inputdata, format, byteswap_input, byteswap_output)
%  * inputdata:  Matlab vector (numeric or char)
%  * format:     desired type of output data ([u]int(8|16|32|64), double, or single
%  *                                             also understands c-style type naming)
%  *                     if the input format is blank or not understood, the format is not
%  *                     changed.
%  * byteswap_input:   boolean, true to reverse byte order based on the input element size
%  * byteswap_output:  boolean, true to reverse byte order based on the output element size
%  * data:       output Matlab vector (same singleton dimension as input -- if input is scalar,
%  *                                   output will be a column vector)
%  *
%  * The function creates a matlab array of the desired type and copies the data 
%  * to the new array -- if the input is the wrong length to give an integer number
%  * of output elements, the data is padded with zeros and a warning is printed.
%  * After filling the new array, the bytes in each element are reversed according to the
%  * byteswap arguments.  If both byteswap input and output are true, the input byteswapping
%  * is done first (it's left as an excercise to the reader to show that this doesn't matter).
%  * If the code is instructed to byteswap single-byte data, a warning is issued.
%  *
%  * CED, 06/14/06
%  * Release Version 1.0
%  * */
% 
