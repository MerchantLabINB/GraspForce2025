function h = myerrorbar(x,y,varargin)
%myerrorbar
%  h = myerrorbar(x,y,e)
%  h = myerrorbar(x,[lower upper])
%  e should be a single colum vector or a two colum matrix [lowerLim upperLim]
%
%  vhdlf 2005

linestyle  = getArgumentValue('linestyle','-',varargin{:},'warningoff');
color      = getArgumentValue('color',[.5 .5 .5],varargin{:},'warningoff');
interval   = getArgumentValue('interval',0,varargin{:},'warningoff');
alphaValue = getArgumentValue('alpha',.08,varargin{:},'warningoff');
logScale   = getArgumentValue('logScale', 0,varargin{:},'warningoff');
errorOnX   = getArgumentValue('erroronx', 0,varargin{:},'warningoff');

x=x(:); % Make sure it is a column vector.
y=y(:); % Make sure it is a column vector.
xticks = [x'; x'];
if isvector(varargin{1})
   e=varargin{1};
   e=e(:);
   yticks = [y'+e' ; y'-e'];
   if errorOnX
      yticks = [y' ; y'];
      xticks = [x'+e'; x'-e'];
   end
else
   e = varargin{1};
   if size(e,1)< size(e,2), e=flipud(rot90(e));end % Make sure it is a two column matrix
   yticks = [e(:,2)' ; e(:,1)'];
end

if logScale
   xticks = log(xticks);
   yticks = log(yticks);
end

if interval
   x = [xticks(2,:) fliplr(xticks(1,:)) ];
   y = [yticks(2,:) fliplr(yticks(1,:)) ];
   handle = patch(x,y,1,'facecolor',color,'edgecolor','none','facealpha',alphaValue);
else
   handle = line(xticks,yticks,'linestyle',linestyle,'color',color,'linewidth',3.5);
end

if nargout
   h = handle;
end
