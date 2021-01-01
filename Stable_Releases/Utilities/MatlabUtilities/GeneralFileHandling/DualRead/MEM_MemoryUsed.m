function memory_used = MEM_MemoryUsed()

%% set defaults
memory_used = 0;

%% set global variables
global gbl_MEM_FID_INFO

%% count memory reserved
for f=1:length(gbl_MEM_FID_INFO)
    if ~isempty(gbl_MEM_FID_INFO(f).maxblocksize)
        memory_used = memory_used + gbl_MEM_FID_INFO(f).maxblocksize;
    end;
end;
