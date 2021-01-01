echo "**************** make distclean *****************"
make distclean
echo "****************    rootcint    *****************"
rootcint eventdict.cxx -c LUX_rq1_channel.h LUX_rq1_pulse.h LUX_rq1_event.h LUX_rq1_header.h
echo "****************      MAKE      *****************"
make