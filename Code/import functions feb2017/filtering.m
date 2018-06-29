
function [xfilt, yfilt, thetafilt, selected] =...
    filtering (oneflydata, velFilt, freqFilt, alignFilt, borderFilt)
    %stimulusfilt, timefilt, timelinefilt, modefilt]...

% [outputs...] = filtering (oneflydata, velFilt, freqFilt, alignFilt, borderFilt)
%
%   oneflydata -> Structure array containing behavioral data from one fly.
%
%   (see more sintax below)
%
% This function does all the preprocessing to a just imported set of
% behavioral data from flies, returning X, Y and THETA (and other) matrices
% if they are requested already filtered. The different processes to be
% made are selected by setting them as '1' in the input arguments,
% following the order in the sintax and the following guide:
%
%   velFilt -> removes trials with less total movement than 25 mm
%   freqfilt -> filters X and Y coordinates and Theta under 2.5 Hz
%   alignFilt -> warps data to actual stimulus experienced
%   borderFilt -> removes 3 mm around the arena, to exclude trajectories along the borders

    

% Sets the initial values as the raw data
xfilt = oneflydata.x; yfilt = oneflydata.y; thetafilt = oneflydata.theta;
%stimulusfilt = oneflydata.stimulus; timefilt = oneflydata.time;
%timelinefilt = oneflydata.timeline; modefilt = oneflydata.mode;



% (1) FILTERS X AND Y COORDINATES, AND THETA ORIENTATION VECTORS
% -------------------------------------------------------------------------

if freqFilt == 1

% Filters coordinates X and Y
[b, a] = butter(2,2.5/25); % Designs a band-pass filter
xfilt = filtfilt(b,a,xfilt); % Filters
yfilt = filtfilt(b,a,yfilt); 

% Filters Theta orientations
[b, a] = butter(2,2.5/25); % Low-pass filter at 2.5 Hz is designed
thetafilt = filtfilt(b,a,(180/pi) * unwrap(thetafilt * (pi/180))); % Filters

clear a b

end % ----------------------------------





% (2) WARPS DATA TO ALIGN TO TRUE STIMULUS ONSET
% -------------------------------------------------------------------------

if alignFilt == 1
    
% [xfilt,yfilt,thetafilt] = retrowarp_old (oneflydata,xfilt,yfilt,thetafilt);

end % --------------------------------------







% (3) SELECTS TRIALS BASED ON MINIMUM MOVEMENT (25 mm)
% -------------------------------------------------------------------------

if velFilt == 1
    
% Calculates distance moved
del = hypot(diff(yfilt,1,1),diff(xfilt,1,1)); % distance moved 
ix = [];
for i = 1:size(xfilt,2)
%     if nanmean(v(:,i)) < 1 % excludes if average velocity is below 1 mm/s
%         ix = [ix, i];
%     end
    if nansum(del(:,i)) < 25 % excludes if total distance moved is below 25 mm
        ix = [ix, i];
    end
end
% Saves the information of which trials are good and selected
selected = setdiff(1:i,ix);
% Eliminates the non-moving trials
xfilt(:,ix) = []; yfilt(:,ix) = []; thetafilt(:,ix) = [];
%stimulusfilt(:,ix) = []; timefilt(:,ix) = []; timelinefilt(ix) = []; modefilt(ix) = [];
% If there are no data left, it ends the run
if isempty(xfilt)
    xfilt = []; yfilt = []; thetafilt = []; selected = [];
    %stimulusfilt = []; timefilt = []; timelinefilt = []; modefilt = [];
    return
end


else % If no velocity filtering is required, a 'selected' output is generated
    selected = 1:size(xfilt,2);

end % ------------------------------------------------




% ***********************   ACTIVATE THIS FILTER  ************************
% activate = 1; % set to 1 to turn on, 0 to turn off
% % ***********************************************************************
% % Extra custom filter to select the trials in which flies start the ON period in a
% % particular part of the arenas AND the ones in which they don't reach the
% % upwind end.
% 
% if activate == 1
% % Define the limits of the region where flies must be at stimulus onset
% reg = [0, 140/5]; % Downwind fifth of the arenas
% stims = oneflydata.stimulus(:,selected); % Takes the stimuli selected so far
% n = 1; nxfilt = [];
% for trial = 1:size(xfilt,2)
%     stind = find(stims(:,trial) ~= 0); % Finds stimulus timing
%     
%     
%     if yfilt(stind(1),trial) > reg(1) && yfilt(stind(1),trial) < reg(2) % if the fly is within the specified part of the arena...
% %         if isempty(find(yfilt(stind(1):stind(end),trial) > 135)) % and if the fly didn't reach the upwind end...
%             nxfilt(:,n) = xfilt(:,trial);
%             nyfilt(:,n) = yfilt(:,trial);
%             nthetafilt(:,n) = thetafilt(:,trial);
%             nselected(n) = selected(trial);
%             n = n+1;
%             
% %         end
%     end
%     
% end
% if isempty(nxfilt)
%     xfilt = []; yfilt = []; thetafilt = []; selected = [];
%     %stimulusfilt = []; timefilt = []; timelinefilt = []; modefilt = [];
%     return
% else
%     xfilt = nxfilt;
%     yfilt = nyfilt;
%     thetafilt = nthetafilt;
%     selected = nselected;
% end
% end
% % **********************************************************





% (4) REMOVES BORDERS OF THE ARENAS
% -------------------------------------------------------------------------

if borderFilt == 1

margin = 3; % mm to remove from the borders of the arena
limys = [0 + margin, 140 - margin];
limxs = [0 + margin, 39.5 - margin];


% Creates a matrix indicating which points to remove
outs = yfilt < limys(1);
outs(yfilt > limys(2)) = 1;
outs(xfilt < limxs(1)) = 1;
outs(xfilt > limxs(2)) = 1;

% Removes the borders from the data
xfilt(outs) = NaN;
yfilt(outs) = NaN;
thetafilt(outs) = NaN;

clear info distance margin limys limxs outs

end % -------------------------------------------------




% Wraps theta information to 0-360 interval
thetafilt = wrapTo360(thetafilt); 


end





















