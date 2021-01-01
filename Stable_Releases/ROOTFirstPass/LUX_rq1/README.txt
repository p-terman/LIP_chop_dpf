ROOT rq1 reader

The idea is that the reader will read the rq1 files and store them in a ROOT tree.  The tree has three classes:

1) The header contains information common to all events such as the name of the dataset and the number of events in the file.
2) The event class stores values which are unique for each event (an event being a dump of the memory) such as the mean value of the event, the timestamp, the baseline and so on.
3) The pulse class is contained within the event class, so that the pulse information (such as head, tail and so forth) can be obtained for each pulse in each event individually.

The program is set up so that the full chain of execution can be handled by one script, in this case called rq1_build_script.  This is called with one argument, being the name of the .rq1 or .rq1.mrg file like:
rq1_build_script lux01_20091018T2344/lux01_20091019T2344_f0000000001-f0000002008.rq1.mrg

This executes the steps:

1) compile the program crazy_class_builder, which will actually write the classes based on the information contained in the header of the .rq1.mrg file.
2) execute crazy_class_builder to write the classes as described above
3) call the make file, which will compile the reader and the example program, and link them with the class files.
4) call the reader, which will create the root file and fill the tree.  The root file will have the same name as the .rq1.mrg file, with the extension .root replacing that type.

The tree is now filled, so the example file rq1_plot_ex can be called, using the root file name as the argument.

Note that if you want to make your way through the file by hand using the ROOT interface (which the author strongly encourages) you have to manually load the library for the classes, so at the command prompt enter this:
.L libTEvent.so

The intended usage in the final product is that the xml reader will inform the script whether it needs to rebuild the classes and the headers by deciding if this file was built with the same verison of the first pass as the currently built headers.  This is a bit unresolved right now, but will hopefully get sorted out soon.  I'll add it to the list of current/future bugs.
FIXED! The reader now handles the xml headers as of 2010-10-18 - PHP

CURRENT AND FUTURE BUGS

This is the stuff I know about right now.  I'll update this file with dates when/if I get the situation resolved.

1) Right now there is no endianness checking.  This is because the endian swap was not necessary on my computer and it was way easier to code up the trees without worrying about it.  I will implement this almost immediately, but for now it doesn't work.
FIXED!  Still has to be tested on big-endian system, but I have to revive my old iBook to make that happen.
KC - 16/04/10

2) If there are more than 26 variables in the rq1 data, you're in trouble.  There are 26 right now, and more will probably be added.  I will come up with some trivial solution for this shortly, I hope.  
FIXED! Largely fixed by moving from the a, b, c variable names to a1, a2, a3 etc.  Possible there are still some of the old type lurking in the script that will rear their ugly heads in the future.  Error to look for is that a variable ends up named { in the LUX_rq1_pulse.cxx which will cause a crash. PHP-2010/10/16

3) In the really weird situation that you are trying to manually read a file (ie in the CINT interface) and it contains a different version of the variables than the last file that you built, you are screwed.  Loading libTEvent.so will give you the wrong header information, and you would have to go back and build it again from the .rq1 file.  There's no easy way to fix this as far as I can tell, but I'm open to suggestions.


4) There is no versioning in this program, ie it rebuilds the classes every time even if you just rebuilt them using a file with exactly the same structure.  This is unnecessary, and should get sorted out soon, but we don't have an xml reader yet so that will have to come later.
FIXED!  The file crazy_class_builder now checks the analysis version that the classes were built with, and if it doesn't need to rebuild them, it doesn't.  The analysis version is stored in the header LUX_rq1_event.h.  This cannot be changed, which is very important. - KJC 04/05/2010

Ken Clark, 15/04/2010
ken.clark@case.edu

Some upkeep of this code has fallen to Patrick Phelps 2010/10/16
patrick.phelps@case.edu
