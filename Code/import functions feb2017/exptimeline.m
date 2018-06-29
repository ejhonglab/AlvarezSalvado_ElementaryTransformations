
function out = exptimeline (data,opt_timelimit)

% out = exptimeline (data)
% 
%   data -> Structure array containing experimental data with trials in the
%           rows. IMPORTANT: Must contain a field with timestamps. By
%           default the field will be considered to be called 'timestamp'.
%           Change the code to use a different name.
%
%   opt_timelimit -> (Optional). Time in seconds. If specified, it eliminates
%           the data collected after the given duration. If not specified, it does nothing.
%
% This function creates a time point for each trial in 'data', which will
% be referenced to the first trial of the whole experiment. It will add the
% field to each of the rows with the corresponding time point in seconds.


if nargin == 2; limax=opt_timelimit; elseif nargin==1, limax=0; end % Sets the time limit if there is any

data = chronosort(data); % Rearranges chronologically the data

sss = strsplit(data(1).timestamp,'.'); % Takes the first timestamp
% Creates the time zero of the experiment
time0 = str2double(sss{1})*3600 + str2double(sss{2})*60 + str2double(sss{3}); 

for i = 1:length(data) % Iterates through each trial
    % Creates a time relative to the time zero of the experiment and adds
    % it to the trial
    split = strsplit(data(i).timestamp,'.'); % Takes the timestamp
    data(i).relativetime = (str2double(split{1})*3600 + str2double(split{2})*60 + str2double(split{3})) - time0;
    if limax ~= 0 % If there is any time limit provided
        if data(i).relativetime <= limax
            continue
        elseif data(i).relativetime > limax
            out = data(1:i-1);
            return
        end
    end         
end

out = data;
clear data sss time0 i split limax
end











