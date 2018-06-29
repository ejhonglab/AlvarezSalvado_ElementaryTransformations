
function out = a_selectdata (resdata, option, opt_parameters, opt_suffix)

% This program performs a "selection" or "filtering" process in one or more
% gait parameters (if none are defined, it uses a default set of
% parameters; see below).
%
% It returns a new field/s named as the orginial ones plus a suffix (that can
% be optionally defined). This allows to do one or more filterings keeping
% all the intermediate steps. When used in this way, make sure to input the
% right field names with the right suffix from one selection to the next.
% 
% Options to select data:
%
%   'noaftertop' - Excludes the period after a fly reaches the top of the arena
%                  (during the stimulus period).
%   'onsetdown' - Selects trials where flies were at the downwind end of the arena
%                 at the stimulus onset.
%   'noborders' - Excludes moments where the fly is near the borders.
%                 Border size can be manually defined.
%   'nosideborders' - Excludes moments where the fly is near the SIDE borders
%                     of the arenas (maintains top and bottom borders).
%                     Border size can be manually defined.
%   'onlyupwind' - Excludes the whole trial if the fly spends more than 10%
%                  of the odor period moving downwind.
%   'minimum25mm' - Removes the trials in which flies didn't move more than
%                   2.5 cm during the odor period.


            
out = resdata; clear resdata

% Sets defaults and the list of parameters to filter
if nargin == 2
    if isfield(out,'pmovew')
        parameters = {'pmovew','vmovew','vymovew','angvturnw','curvaturew','pturnw'};
    elseif ~isfield(out,'pmovew')
        parameters = {'pmove','vmove','vymove','angvturn','curvature','pturn'};
    end
%     suffix = 's';
elseif nargin == 3
    parameters = opt_parameters;
%     suffix = 's';
elseif nargin == 4
    if ~isempty(opt_parameters)
        parameters = opt_parameters;
    elseif isempty(opt_parameters)
        if isfield(out,'pmovew')
            parameters = {'pmovew','vmovew','vymovew','angvturnw','curvaturew','pturnw'};
        elseif ~isfield(out,'pmovew')
            parameters = {'pmove','vmove','vymove','angvturn','curvature','pturn'};
        end
    end
    suffix = opt_suffix;
end



for fly = 1:length(out) % Iterates through each fly
    
    for trial = 1:size(out(fly).xfilt,2) % Iterates through each trial
        
        stind = find(out(fly).stimulus(:,out(fly).selected(trial)) ~= 0); % finds stimulus time
        
        % ++++++++++++++++++++  CONDITIONS  +++++++++++++++++++++++++++++++
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
            
        switch option
            

            % Selects the period after a fly reaches the top of the arena (during the stimulus period)
            case 'noaftertop' 
                if nargin<4, suffix = 's'; end
                bad = find(out(fly).yfilt(stind(1):stind(end),trial) > 135,1) + stind(1); % Finds the first encounter with the top during the stimulus period
                bads = bad:1:stind(end); % Takes the rest of the odor period after the first encounter with the top
                
%                 plot(out(fly).yfilt(:,trial)); hold on; plot(stind(1),out(fly).yfilt(stind(1),trial),'o','linewidth',3); hold off;




            % Selects flies that were at the downwind end of the arena at the stimulus onset
            case 'onsetdown' 
                if nargin<4, suffix = 'd'; end
                if out(fly).yfilt(stind(1),trial) > 28
                    bads = 1:length(out(fly).yfilt(1:end-1,trial));
                else
                    bads = [];
                end
            
                
                
            % Excludes the points where the fly is not moving upwind during
            % the odor at least 90% of the time
            case 'onlyupwind'
                if nargin<4, suffix = 'u'; end
%                 bads = find(out(fly).vymove(stind(1):stind(end),trial) < 0) + stind(1); %
                downs = find(out(fly).vymovew(stind(1):stind(end),trial) < 0);
                if numel(downs) > round((stind(end)-stind(1))/10)
                    bads = 1:length(out(fly).yfilt(1:end-1,trial));
                else
                    bads = [];
                end
            
                
                
            
            % Removes the trials in which flies didn't move more than 5 cm
            % during the odor period
            case 'minimum25mm'
                if nargin<4, suffix = 'm'; end
                distance = nansum(out(fly).vmove(stind(1):stind(end),trial).*0.02);
                if distance < 25
                    bads = 1:length(out(fly).yfilt(1:end-1,trial));
                elseif distance >= 25
                    bads = [];
                end    
                
                    
                
                
            % Excludes moments where the fly is near the borders.
            case 'noborders'
               if nargin<4, suffix = 'b'; end
%                xmargin = 3; % mm to remove from the SIDE borders of the arena
%                ymargin = 10; % mm to remove from the TOP-BOTTOM borders of the arena
               % SET LIMITS TO EXCLUDE BORDERS AREAS ------------------
               limys = [3, 137]; % Natural limits 0 - 140
               limxs = [3, 36.5]; % Natural limits 0 - 39.5
               % ------------------------------------------------------
               bads = out(fly).yfilt(1:end-1,trial) < limys(1);
               bads(out(fly).yfilt(1:end-1,trial) > limys(2)) = 1;
               bads(out(fly).xfilt(1:end-1,trial) < limxs(1)) = 1;
               bads(out(fly).xfilt(1:end-1,trial) > limxs(2)) = 1;
               if nansum(bads)<1, bads=[]; end % makes empty matrix if the fly didn't touch the borders
               
               
               
               
            % Excludes moments where the fly is near the SIDE borders of the arenas (maintains top and bottom borders).
            case 'nosideborders'
               if nargin<4, suffix = 'sb'; end
               margin = 3; % mm to remove from the borders of the arena
               limxs = [0 + margin, 39.5 - margin];
               bads = out(fly).xfilt(1:end-1,trial) < limxs(1);
               bads(out(fly).xfilt(1:end-1,trial) > limxs(2)) = 1;
               if nansum(bads)<1, bads=[]; end % makes empty matrix if the fly didn't touch the borders
        
        end
        
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
            
        
        % --------------------------------------------------------
        % NANs out the selected values, on the selected parameters
        for param = 1:length(parameters) % Iterates through each parameter
            
            data = eval(['out(fly).', parameters{param}, '(:,trial);']); % Takes the corresponding vector
            
            if ~isempty(bads) % Only filters if the conditions were met
                data(bads) = nan; % NANs out the data that meets the condition
            end
            
%             eval(['out(fly).', name, '(:,trial) = data;']); % Saves to output
            eval(['out(fly).', parameters{param}, suffix, '(:,trial) = data;']); % Saves to output
        end
    end
end















