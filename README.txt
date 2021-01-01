LUXcode Notes - 2010-10-18 - P. Phelps (patrick.phelps@case.edu)
----------------------------------------------------------------------

As of 2010/10/18 svn organization has limited top level directories

Trunk - This is the main shared branch of code.  Things in trunk are considered in beta, usable but not gaurenteed.  This is the main code sharing space for late level debugging and development. Files in Trunk are expected to be decently stable, for instance to run most of the time and to be out of low level development and on to bug testing across platforms etc.

Tags - READ ONLY! DO NOT EDIT TAGS!!
Tags directory is for snapshots of the code at a given time.  For instance we might take a snapshot of the code to be used to publish a paper and save it always and forever.  If you are caught editing a TAG you get one warning, after that your SVN access will be revoked until the next analysis meeting where it will be discussed.  I'm being serous!

Scratch - Scratch is a place to mess around, a place to test and develope alpha code.  If you decide you'd like to modify a script from the trunk and save all kinds of versions of it as you work scratch is the place to do that.  Please start by making a user directork under scratch using the svn mkdir command.  Then you can use this space for work.  When you get your file up to the standards of Trunk you can move it to there.

Stable_Releases - This directory is meant to keep the most up-to-date stable releases of major chunks of code.  This code isn't being developed on other then to release new stable versions so it is a good place for users just looking to use code to do something and move on.  


