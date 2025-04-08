function fig = togglefig(name)
% Finds and activates, or creates, figure with user-specified name.
%
% If no name is provided, creates figure named "untitledn" (where n is
% incremented to result in a unique name).
%
% SYNTAX:
% togglefig('My Figure');
%    If figure named 'My Figure' exists, it will be activated (brought to
%    the front and shown). Otherwise, it will be created.
%
% h = togglefig('My Figure');
%    Also returns the handle to the specified or created figure.
%
% togglefig;
%    Creates and activates new figure named untitled1, untitled2, ...
%    Note: You can subsequently activate these figures with, for instance,
%          togglefig('untitled1').
% 
% OTHER EXAMPLES:
% NOTE: This example requires the Image Processing Toolbox
% im = imread('cameraman.tif');
% for ii = 1:10
%    thresh = ii/20;
%    togglefig('Threshold');
%    imshow(im2bw(im,thresh));
%    title(sprintf('Threshold = %0.2f',thresh));
%    pause(1)
% end
% 
% Motivation:
%   I've found this to be exceptionally useful in algorith-development
%   mode, particularly when iterating on cells in the cell-mode editor. (I
%   use this function in almost every mfile I write these days.)
%
% Coded by Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 01/04/08

% Modifications:
% 01/10/08
% To suppress handle generation if nargout == 0 (per John D'Errico's
% suggestion).

if ~nargin
   callf= dbstack; % Added by victor.
    if length(callf)>1
      name = callf(2).file;
    else
       tmp = get(findall(0,'type','figure'),'name');
       if isempty(tmp)
           tmp = 0;
       elseif ischar(tmp) %one match
           tmp = 1;
       else %should be a cell array of matches
           tmp = sum(cell2mat(regexp(tmp,'untitled')));
       end
       name = ['untitled',num2str(tmp+1)];
   end
end

fig = findall(0,'type','figure','name',name);
if isempty(fig)
   % Modified by victor. ------------------------
   pause(.01)   % So that the next call to clock is different from the previous one.
   FigNum = round(sum(1000*clock)); % Create figure number.
   fig = figure(FigNum);
   set(fig,'numbertitle','off','name',name);%shg;
%   axdrag
%   GTscroll(fig)
   % ----------------------------------------------
else
    figure(fig);
end
% If no output was requested, none should be generated
if ~nargout
    clear fig
end