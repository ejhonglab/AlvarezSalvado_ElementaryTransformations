
function stimseparator (data, opt_workspace)

% stimseparator (data, opt_workspace)
%
%   data -> Structure array containing data with trials in each row.
%
%   opt_workspace (optional) -> Specifies the workspace where the separated
%           data will be saved. By default it is 'base', the normal workspace. If
%           this function is invoked by another function though, workspace will
%           have to be specified as 'caller'.
%
% This function takes a structure array with data with different
% stimulation modes, it detects which are those, and returns a separated
% structure containing the trials corresponding to each stimulation.

if nargin == 1 % Specifies the default workspace or takes the one provided
    workspace = 'base';
elseif nargin ==2
    workspace = opt_workspace;   
end

% detects the names of the stimuli
stims = cell (length(data),1);

for i = 1:length(data) % Creates a list with all the stimuli names
    stims{i} = data(i).mode;
end
names = unique(stims); % Takes the different stimuli names
clear stims

assignin('caller','Stimulinames',names);

for i = 1:length(names) % Iterates through each different stimulus
    
    n = 1;
    
    for trial = 1:length(data) % Checks each trial, and stores it if corresponds to the current stimulus
        switch data(trial).mode;
            case names(i);
                seleccionados(n) = data(trial);
                n = n+1;
        end
    end
    
    name = ['s',names{i}]; % Creates a name for the structure because it cannot start with numbers 
    assignin(workspace, name, seleccionados); % Saves the structure array to the workspace
    clear seleccionados
end

clear workspace name n i trial names

end