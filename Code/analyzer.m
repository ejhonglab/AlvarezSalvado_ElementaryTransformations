 
function results = analyzer (data)

% THIS IS JUST A SHORTCUT FUNCTION
%
% results = analyzer (data)
%
%   data -> Structure array with flies' data.
%
% Function that calls a series of standard analysis functions. 


tic
removers = []; results = data;


% SELECT FILTERS
% Make sure to select (1) or unselect (0) the filters to apply:
velFilt = 1;       % removes trials with less total movement than 25 mm
freqFilt = 1;      % filters X and Y coordinates and Theta under 2.5 Hz
alignFilt = 0;     % warps data to fit actual stimulus experienced
borderFilt = 0;    % removes 3 mm around the arena, to exclude trajectories along the borders
% ------------------------------------


for fly = 1:length(data) % Iterates through flies
    
    % FILTERS THE DATA applying the filters previously selected
    [results(fly).xfilt, results(fly).yfilt, results(fly).thetafilt,...
        results(fly).selected] = ...
        filtering (results(fly),velFilt,freqFilt,alignFilt,borderFilt);
    
    % IF WE ARE DEALING WITH RANDOM STIMULI, USE THE FOLLOWING LINE INSTEAD
    % OF THE PREVIOUS
%     [results(fly).xfilt, results(fly).yfilt, results(fly).thetafilt, results(fly).stimulusfilt, results(fly).timefilt, results(fly).timelinefilt] =...
%         filtering (results(fly),1,1,1,0);
    % --------------------------------
    
    % Checks that there is still some data remaining after filtering and reports
    if isempty(results(fly).xfilt) || size(results(fly).xfilt,2) < 5 
        removers = [removers, fly];
        display(['Fly #', num2str(fly), ' discarded after filtering for lack of activity']);
        continue
    end
end

results(removers) = []; % Removes flies that were left without data due to lack of activity
toc

% ANALYZES THE RESULTING DATA with analysis functions written specifically
% for this type of structure. You can add or remove at will.
results = a_vels (results); toc % Calculates all kinds of velocities
results = a_warp (results); toc % Warps the data to match exact time of odor encounter













    
    
    
    
    
    
    
    
