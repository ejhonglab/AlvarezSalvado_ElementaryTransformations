
function out = a_ydist (resdata, opt_bins_number)

% out = a_ydist (resdata)
%
%   resdata -> Structure array with flies' data.
%
%   opt_bins_number -> (Optional). Number of bins to do the histogram.
%
% Calculates the Y distribution of flies' positions. It returns the same
% input structure with a new field containing the Y distribution (expressed
% in proportion of total trial time). It also returns a "ybins" field with
% the centers of each bin used.


% Decides number of bins (default or provided)
if nargin == 1
    binsnum = 70;
else
    binsnum = opt_bins_number;
end

% Copies data to modify and return as output
out = resdata; clear resdata

% Calculates bins limits
measure = 140/binsnum;
bins = ((measure/2):measure:140)';



for fly = 1:length(out) % Iterates through each fly

%     [~,yyy] = filtering(out(fly),1,1,1,1); % Filters borders out
    
    out(fly).ydist = []; % Resets variable to overwrite previous distributions
    out(fly).ybins = bins; % Saves the bins' centers
    
    % Calculates the Y distributions of each trial of each fly
    for trial = 1:size(out(fly).yfilt,2)
        
        out(fly).ydist(:,trial) = hist(out(fly).yfilt(:,trial),bins); % Calculates the distribution
        out(fly).ydist(:,trial) =...
            out(fly).ydist(:,trial)/nansum(out(fly).ydist(:,trial)); % Normalizes to total time
           
    end
        
end

    
    
    
    
    
    
    
    
    
    