
function out = a_warp (resdata, opt_parameters)

% This program warps the traces of certain parameters (that can be provided
% optionally as a second argument) to make them fit the exact time of odor
% encounter/offset of flies.

tic
out = resdata; clear resdata

% List of parameters to warp
if nargin == 1
    parameters = {'pmove','vmove','vymove','angvturn','curvature','pturn'};
elseif nargin == 2
    parameters = opt_parameters;
end

for fly = 1:length(out) % Iterates through each fly
    
    for param = 1:length(parameters) % Iterates through each parameter and warps it
    
        for trial = 1:size(out(fly).xfilt,2) % Iterates through each trial
        
            % Calls the function "stimwarp" for each parameter warped
            data = eval(['stimwarp(out(fly).', parameters{param}, '(:,trial), out(fly).yfilt(:,trial), out(fly).stimulus(:,out(fly).selected(trial)), out(fly).mode{out(fly).selected(trial)});']);
            eval(['out(fly).', parameters{param}, 'w(:,trial) = data;']);    
            
        end
    end
end
toc     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        