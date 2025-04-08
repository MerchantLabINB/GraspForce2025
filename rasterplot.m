function [xticks, yticks] = rasterplot(spiketimes,varargin)
%
% Generates a raster plot of the spike times.
%
% [xticks, yticks] = rasterplot(SpikeTimes)
% [xticks, yticks] = rasterplot(SpikeTimes,'QuickPlot','yes','TickMarkLength',0.8)
% [xticks, yticks] = rasterplot(SpikeTimes,'xlim',[-1 2])
% [xticks, yticks] = rasterplot(SpikeTimes,'displace',.5) Use 'displace' to move the spikes along the x-axis
%
% Example:
% SpikeTimes{1}   = [22 34 58 77];          
% SpikeTimes{2}   = [33 44 55 45 78 98 99];
% [xticks,yticks] = rasterplot(SpikeTimes)
% plot(xticks,yticks)
%
% vdl aug2007 - jun2011 -may2020

% Default parameter values:
TickMarkLength = 0.9;        % Proportion of row height occupied by the spike tick marks.
QuickPlot      = 'no';       % 'yes' will plot dots instead of lines, which is much faster.
xlim           = [-inf inf]; % Plots all spikes by default.
displace       = 0;          % Move the spikes along the x-axis. 
fig            = [];         % Figure handle.
spkColor       = [0 0 0]; % Gray is the default spike color.
maxNumSpikes   = 5e5;     % Maximum number of total spikes that will be displayed

% To change the parameter values:
for field_num = 1:2:length(varargin)
   if strcmpi(varargin{field_num},'QuickPlot')
      QuickPlot = varargin{field_num+1};   % Modifies the default value;
   elseif strcmpi(varargin{field_num},'TickMarkLength')
      TickMarkLength = varargin{field_num+1};   % Modifies the default value;
   elseif strcmpi(varargin{field_num},'xlim')
      xlim = varargin{field_num+1};   % Modifies the default value;
   elseif strcmpi(varargin{field_num},'displace')
      displace = varargin{field_num+1};   % Modifies the default value;
   elseif strcmpi(varargin{field_num},'maxNumSpikes')
      maxNumSpikes = varargin{field_num+1};   % Modifies the default value;
   elseif strcmpi(varargin{field_num},'fig')
      fig = varargin{field_num+1};   % Modifies the default value;
   elseif strcmpi(varargin{field_num},'color')
      spkColor = varargin{field_num+1};   % Modifies the default value;
   end
end

% Platform specific parameters. See the "Memory Allocation" help topic.
% We're assuming spike times are stored as double.
%headerSize = 60; % Bytes.
%doubleSize =  8; % Bytes.
% Get the info about the spiketimes cell array.
%info        = whos('spiketimes')
% Calculate the number of spikes stored in the whole array.
%numOfSpikes = (info.bytes-headerSize*length(spiketimes))/doubleSize

if ~iscell(spiketimes) % In case spiketimes is not a cell array.
   spiketimes = {spiketimes};
end

% Obtain the total number of spikes
numOfSpikes = sum(cellfun(@length,spiketimes));

% If number of spikes is larger than maxNumSpikes decimate the rasterplot.
% All spike trains are decimated so that the relative number of spikes
% across them remains unchanged. 
if numOfSpikes>maxNumSpikes
    
    decimateRatio = numOfSpikes/maxNumSpikes;
    
    for n = 1:length(spiketimes)
        N = length(spiketimes{n});
%        keep = round(1:decimateRatio:N);
        keep = round(decimateRatio:decimateRatio:N);
        
        spiketimes{n} = spiketimes{n}(keep);
    end
    
    numOfSpikesDecimated = sum(cellfun(@length,spiketimes));
    disp(['rasterplot.m: Decimating number of spikes in rasterplot. Displaying ' ...
              '1 of ~ ' num2str(round(decimateRatio*10)/10) ' spikes (' ...
              num2str(numOfSpikesDecimated) '/' num2str(numOfSpikes) ').'])
    numOfSpikes = numOfSpikesDecimated;
end

% Allocate memory for the output matrices/vectors.
if strcmpi(QuickPlot,'yes') % For plotting dots instead of lines.
   xticks = zeros(1,numOfSpikes);
   yticks = zeros(1,numOfSpikes);
else
   xticks = zeros(2,numOfSpikes);
   yticks = zeros(2,numOfSpikes);
end

dk = 1; % counter
for trial = 1 : length(spiketimes)

%   spiketimes{k}   = ktrial.spikes(ktrial.spikes>=xlim(1) & ktrial.spikes<=xlim(2)) + displ;
    spkt     = spiketimes{trial}(:);    % Make the spike times a colum vector.
    
    % Select the spikes within xlim and displace them
    spkt     = spkt( spkt>=xlim(1) & spkt<=xlim(2) ) + displace;

    rowNumb  = trial-TickMarkLength/2;  % The spike tickmark beggins here.
    numSpk   = length(spkt);            % Get the number of spikes in the trial.

   % x and y position of each tickmark
   if strcmpi(QuickPlot,'yes') % For plotting dots instead of lines (plotting will be faster).
      xticks ( :, dk:dk + numSpk-1 ) = spkt';
      yticks ( :, dk:dk + numSpk-1 ) = trial*ones(1,numSpk);
   else % For plotting lines (plotting will be slower).
      xticks ( :, dk:dk + numSpk-1 ) = [spkt'; spkt'];
      yticks ( :, dk:dk + numSpk-1 ) = [rowNumb*ones(1,numSpk) ; ...
              rowNumb*ones(1,numSpk)+TickMarkLength];
   end
    % Update the counter.
    dk  = dk  + numSpk;
end

% For plotting
if nargout == 0
   if isempty(fig)
      togglefig('rasterplot'), clf,drawnow
   else
      figure(fig)
   end

   if strcmpi(QuickPlot,'yes') 
      linestyle = 'none'; marker    = '.';
   else
      linestyle = '-'; marker    = 'none';
   end

   % Plot and configure axes.
   line(xticks,yticks,'LineStyle',linestyle,'marker',marker, ...
       'markersize', 10,'linewidth',0.3,'color',spkColor)
   set(gca,'tickdir','out'); box off
end

