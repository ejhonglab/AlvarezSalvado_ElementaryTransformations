 
function out = a_thetahist (data, opt_periods)

% This program calculates the orientation distributions of movements of
% flies, and returns the same input structure with:
%   1) the vectors of the distribution for each trial of each fly contained
%      in the "data" structure, and
%   2) the orientation bins (in radians) to represent the distributions. 
%
% The program will automatically divide trials in three periods if there is
% any stimulus (i.e. "prethetahist", "durthetahist" and "posthetahist"), or
% will analyze whole trials for blank trials. Alternatively, one can input
% a vector containing four time points (as the optional, second input argument)
% and the program will analyze the three periods defined between them.
%
% The histograms returned are expressed in proportion of total trial time
% spent traveling in each orientation. To normalize, the program doesn't 
% consider the NaN values.
%
% This function also optionally calculates the average vectors in two
% different forms (see below).


% IF YOU WANT TO OBTAIN THE AVERAGE VECTORS IN RETURN.
aveve = 0; % Set this to 1 if you want dir+length. Set to 2 if you want cartesian coordinates


tic
out = data; clear data; % Copies input data to change and return it as output
edges = 1.25:2.5:360; % Defines the bins to distribute the orientation values
oribins = edges' * (pi/180); % Edges in radians


% Takes periods of time if given
if nargin == 2
    pe = opt_periods;
end


for i = 1:length(out) % Iterates through each fly in the data
    
    % FILTERS OUT ALL ARENA BORDERS (3 mm)
    [~,~,out(i).thetafiltnb] = filtering(out(i),1,1,1,1);
    out(i) = a_warp(out(i), {'thetafiltnb'});
    thetafilt = out(i).thetafiltnbw;
    
    % Defines time periods if they are not given
    stind = find(out(i).stimulus(:,1) ~= 0);
    if nargin == 1
        stind = find(out(i).stimulus(:,1) ~= 0);
        if ~isempty(stind)
            pe = [1, stind(1), stind(end), size(thetafilt,1)]; % Periods that take all data
        elseif isempty(stind)
            pe = [1, size(thetafilt,1)]; % Period that takes all data
        end
    end
    
    % Iterates through each trial and collects orientation histograms
    for trial = 1:size(thetafilt,2)
        
        vector = thetafilt(:,trial); % Takes orientations vector of the current trial
        
        % In case there is a stimulus --------------------------------------
        if ~isempty(stind) 
            nanes = length(find(isnan(vector(pe(1):pe(2)))==1)); % Calculates the amount of NaN values to not count them in the normalization
            out(i).prethetahist(:,trial) = hist(vector(pe(1):pe(2)),edges)'./(pe(2)-pe(1)+1-nanes); % Counts and normalizes
            nanes = length(find(isnan(vector(pe(2):pe(3)))==1));
            out(i).durthetahist(:,trial) = hist(vector(pe(2):pe(3)),edges)'./(pe(3)-pe(2)+1-nanes);
            nanes = length(find(isnan(vector(pe(3):pe(4)))==1));
            out(i).posthetahist(:,trial) = hist(vector(pe(3):pe(4)),edges)'./(pe(4)-pe(3)+1-nanes);
            
            % Also collects average vectors if wanted (see line 24)
            if aveve == 1
                out(i).preavec(trial,:) = vectoraverage (oribins,out(i).prethetahist(:,trial));
                out(i).duravec(trial,:) = vectoraverage (oribins,out(i).durthetahist(:,trial));
                out(i).posavec(trial,:) = vectoraverage (oribins,out(i).posthetahist(:,trial));
            elseif aveve == 2
                [~,out(i).preavec(trial,:)] = vectoraverage (oribins,out(i).prethetahist(:,trial));
                [~,out(i).duravec(trial,:)] = vectoraverage (oribins,out(i).durthetahist(:,trial));
                [~,out(i).posavec(trial,:)] = vectoraverage (oribins,out(i).posthetahist(:,trial));
            end
            
        % In case there is NO stimulus ------------------------------------
        elseif isempty(stind) 
            
            if numel(pe) == 2 % With only one period
                nanes = length(find(isnan(vector(pe(1):pe(2)))==1));
                out(i).thetahist(:,trial) = hist(vector(pe(1):pe(2)),edges)'./(pe(2)-pe(1)+1-nanes);
                % Also collects average vectors
                if aveve == 1
                    out(i).avec(trial,:) = vectoraverage (oribins,out(i).thetahist(:,trial));
                elseif aveve == 2
                    [~,out(i).avec(trial,:)] = vectoraverage (oribins,out(i).thetahist(:,trial));
                end
            
            elseif numel(pe) == 4 % With three periods 
                nanes = length(find(isnan(vector(pe(1):pe(2)))==1));
                out(i).prethetahist(:,trial) = hist(vector(pe(1):pe(2)),edges)'./(pe(2)-pe(1)+1-nanes);
                nanes = length(find(isnan(vector(pe(2):pe(3)))==1));
                out(i).durthetahist(:,trial) = hist(vector(pe(2):pe(3)),edges)'./(pe(3)-pe(2)+1-nanes);
                nanes = length(find(isnan(vector(pe(3):pe(4)))==1));
                out(i).posthetahist(:,trial) = hist(vector(pe(3):pe(4)),edges)'./(pe(4)-pe(3)+1-nanes);
                % Also collects average vectors
                if aveve == 1
                    out(i).preavec(trial,:) = vectoraverage (oribins,out(i).prethetahist(:,trial));
                    out(i).duravec(trial,:) = vectoraverage (oribins,out(i).durthetahist(:,trial));
                    out(i).posavec(trial,:) = vectoraverage (oribins,out(i).posthetahist(:,trial));
                elseif aveve == 2
                    [~,out(i).preavec(trial,:)] = vectoraverage (oribins,out(i).prethetahist(:,trial));
                    [~,out(i).duravec(trial,:)] = vectoraverage (oribins,out(i).durthetahist(:,trial));
                    [~,out(i).posavec(trial,:)] = vectoraverage (oribins,out(i).posthetahist(:,trial));
                end
            end
        end
    end
    
    % Saves the bins info IN RADIANS
    out(i).orientbins = edges' * (pi/180);
    
end
toc
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
