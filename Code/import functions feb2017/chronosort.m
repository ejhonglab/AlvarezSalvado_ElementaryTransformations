
function out = chronosort (data)

% out = chronosort (data)
%
%   data -> Structure array containing trials in each row. IMPORTANT: It
%           has to have a field with a timestamp named 'timestamp'. If it's
%           named different, it has to be changed in the text of the
%           function.
%
% This function rearranges in chronological order the rows in a structure 
% array that contains data with trials in the different rows.

field = 'timestamp'; % Change here the name of the field if necessary

% Extracts the timestamps the number of experiments putting together hours,
% minutes and seconds in a big number
for i = 1:length(data) % Creates a list with all the experiment names
    eval(['stamp = data(i).',field,';']);
    split = strsplit(stamp,'.');
    times(i,1) = str2num([split{1},split{2},split{3}]); % Stores the big number together 
end

% Sorts the time stamps
[~, indices] = sort(times);

out = data(indices);

clear field i stamp split times indices

end