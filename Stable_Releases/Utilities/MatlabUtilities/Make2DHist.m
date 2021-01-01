%	function [xbincenters,ybincenters,counts]= ...
%		make2Dhist(xdata,ydata,xbinedges,ybinedges)
%	
%	Creates 2D histograms.  Works for non-uniform bins.
%
%		x(y)data			= data to bin
%		x(y)binedges	= vectors of N+1 bin edges.  
%								First values is low edge of first bin,
%								last value is high edge of last bin.
%	
%		counts 			= 2D histogram of data.
%					 			y given by column (first index),
%								x given by row (second index).
%		x(y)bincenters	= vectors of N centers of bins.
%
%		bins are lower edge inclusive, upper edge exclusive.
%
%	1/29/97	T.S.	

%	3/97 mod bin arguments.
%	9/97	allow for non-uniform bins.

function [xbincenters,ybincenters,counts]=make2Dhist(xdata,ydata,xbinedges,ybinedges)

%
% 	Check input
%

if nargin~=4
	disp('*** Error in make2Dhist: bad number of input arguements ***');
	return;
end;
if nargout~=3
	disp('*** Error in make2Dhist: bad number of output arguements ***');
	return;
end

if (length(xdata) ~= length(ydata))
	disp('*** Error in make2D hist: data vectors not the same size ***');
	return;
end
if length(xbinedges)<2|length(ybinedges)<2
	disp('*** Error in make2Dhist: bad input ranges ***');
	return;
end


%
%	Set up binning
%

xlow=xbinedges(1);
ylow=ybinedges(1);
xhigh=xbinedges(length(xbinedges));
yhigh=ybinedges(length(ybinedges));

numxbins=length(xbinedges)-1;
numybins=length(ybinedges)-1;

xbinwidth=diff(xbinedges);
ybinwidth=diff(ybinedges);

xbincenters=xbinedges(1:length(xbinedges)-1);
xbincenters=(xbincenters+0.5*xbinwidth);

ybincenters=ybinedges(1:length(ybinedges)-1);
ybincenters=(ybincenters+0.5*ybinwidth);

%
% 	build the histogram
%

counts=zeros(numxbins,numybins);

cut=xdata>=xlow&xdata<xhigh&ydata>=ylow&ydata<yhigh;
xdata=xdata(cut);
ydata=ydata(cut);

for ne=1:length(xdata)
	xbin=max(find(xdata(ne)>=xbinedges));
	ybin=max(find(ydata(ne)>=ybinedges));
	counts(xbin,ybin)=counts(xbin,ybin)+1;
end;


