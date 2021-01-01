% function machine_format = MEM_SwitchEndianness(MEM_fid)
% 
% Quick function to flip the endianness used for a file already loaded into
% memory -- uesd by ReadNtupleFile
%
% 5/13/08, CED

function machine_format = MEM_SwitchEndianness(MEM_fid)

%% declare global variables
global gbl_MEM_FID_INFO

%% set defaults
machine_format = '';

%% switch swapbytes
gbl_MEM_FID_INFO(MEM_fid).swapbytes = ~gbl_MEM_FID_INFO(MEM_fid).swapbytes;

%% output new machineformat, if requested
if nargout > 0
    endiantest = hex2dec('01020304');
    bytewise_endiantest = double(TypeCastAndSwap(uint32(endiantest),'uint8',false,false));
    currently_big_endian = (bytewise_endiantest(1)==1);

    if xor(currently_big_endian, gbl_MEM_FID_INFO(MEM_fid).swapbytes)
        machine_format = 'l';
    else
        machine_format = 'b';
    end
end
