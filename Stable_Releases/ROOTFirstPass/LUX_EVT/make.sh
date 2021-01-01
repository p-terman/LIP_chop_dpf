echo "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
echo "ooooo ooooo ooooo         make distclean                ooooo ooooo ooooo"
echo "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
make distclean
echo "ooooo ooooo ooooo            rootcint                   ooooo ooooo ooooo"
rootcint eventdict.cxx -c Lux_EVT_Event.h Lux_EVT_Channel.h Lux_EVT_Pulse.h Lux_EVT_Header.h
echo "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
echo "ooooo ooooo ooooo              MAKE                     ooooo ooooo ooooo"
echo "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo"
make