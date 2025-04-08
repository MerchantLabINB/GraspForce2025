function fRates = firingrate(SpikeTimes, TimeSamples, varargin)
%firingrate() Calculates the firing rate at the provided TimeSamples.
%
%   FRATES = firingrate(SPIKETIMES, TIMESAMPLES, OPTIONS)
%   Examples:
%   fRates = firingrate(SpikeTimes, TimeSamples, 'FilterType', 'exponential', 'TimeConstant', 0.05);
%   fRates = firingrate(SpikeTimes, TimeSamples, 'FilterType', 'boxcar'     , 'TimeConstant', 0.05);
%   fRates = firingrate(SpikeTimes, TimeSamples, 'Attrit',      [0 10]);
% 
%   SpikeTimes can be either a vector of spike times or a cell array containing one vector per cell.
%   Time unit: seconds.
%   fRates: spikes/second
%   Options: 
%   'FilterType'   ,  'boxcar' or 'exponential'
%   'TimeConstant' ,  [.05]        If 'Filtertype';'exponential' the TimeConstant sets the
%                                  decay rate of the exponential (see equation below).
%                                  If 'Filtertype';'boxcar' the TimeConstant sets the total length of the window.
%   'Attrit'       ,  [tIni tEnd]  Either only one pair, or the number of rows must equal the number of spike vectors. 
%                                  Attrit will put NANs at every TimeSample falling   outside the tIni tEnd bounds.
%   'Normalize'    ,  [TimeConstant Center] Always uses a boxcar
%   vhdlf

if ~iscell(SpikeTimes) % In case SpikeTimes is not a cell array.
   SpikeTimes = {SpikeTimes};
end

% getArgumentValue('ARGUMENT',DEFAULTVALUE, VARARGIN);
FilterType   = getArgumentValue('FilterType'   ,'exponential', varargin{:});
TimeConstant = getArgumentValue('TimeConstant' ,0.05         , varargin{:}); % Time constant in seconds.
attrit       = getArgumentValue('attrit',[TimeSamples(1) TimeSamples(end)], varargin{:},'warningoff'); % Attrition times. 
normWidthCent= getArgumentValue('Normalize',[], varargin{:},'warningoff'); % Normalization with and center. 

if size(attrit,1)>1 && size(attrit,1)~=length(SpikeTimes) % Error checking
   error('The number of rows in attrit must equal the number of SpikeTimes vectors.')
end

ntrials = length(SpikeTimes);                 % Get the number of trials.
fRates  = nan(ntrials,length(TimeSamples));   % Initialize fRates matrix (a row for each trial).

% Loop through each trial.
for k = 1:ntrials
   if size(attrit,1)>1 % If only one pair or attrition times were provided use that pair for all the trials. 
      subtimesamp = (TimeSamples<= attrit(k,2) & TimeSamples>= attrit(k,1));
   else                % Otherwise use a different pair for each trial. 
      subtimesamp = (TimeSamples<= attrit(1,2) & TimeSamples>= attrit(1,1));
   end
   timeSamples = TimeSamples(subtimesamp);
   spkT  = SpikeTimes{k};                    % Get the trial spikes.
   spkT  = spkT(:);                          % Make it a column vector.
   spkT  = spkT*ones(1,size(timeSamples,2));
   timesamples = ones(1,size(spkT,1))'*timeSamples;

   switch lower(FilterType)
      case 'boxcar' % Counts spikes between -TimeConstant/2 and +TimeConstant2 of each time sample (i.e. it extends into the past and future). 
         fRate  = (sum(spkT>=(timesamples-TimeConstant/2) & spkT<(timesamples+TimeConstant/2),1)/TimeConstant);
      case 'exponential' % Only counts spikes into the past. 
         SelectSpikesUpTo = TimeConstant * 7.5; % Time window within which spike will contribute to the firing rate of each time sample.
         SpikesBeforeTimeSample=(spkT>(timesamples-SelectSpikesUpTo) & spkT<=(timesamples));
         DistanceToTimeSample = spkT - timesamples;
         DistanceToTimeSample = abs(SpikesBeforeTimeSample.*DistanceToTimeSample);
         fRate = 1./TimeConstant .* exp(-DistanceToTimeSample/TimeConstant);
         fRate(logical(SpikesBeforeTimeSample~=1))=0; % Remove those distaces that we are not considering.
         fRate = sum(fRate,1);
      otherwise
         error(['FilterType "' FilterType '" not supported'])
   end
   fRates(k,subtimesamp) = fRate;  % This is the quick solution to many trials. Look for better ones.   
end

% If normalization is required
if ~isempty(normWidthCent)
  FRforNorm = firingrate(SpikeTimes, normWidthCent(2), 'FilterType','boxcar','TimeConstant',normWidthCent(1));
%  fRates = (fRates-(mean(FRforNorm)) ) / std(FRforNorm);
  fRates = fRates / mean(FRforNorm);  % Mikes normalization
end

if nargout==0
   %togglefig(mfilename)
   plot(TimeSamples,fRates,'.-')
   xlabel('time (s)')
   ylabel('firing rate (spikes/s)')
end