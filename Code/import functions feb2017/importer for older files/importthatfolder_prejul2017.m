
function out = importthatfolder_prejul2017 (path, expnumber)

% out = importthatfolder (path)
%
%   path -> Absolute path of the folder to import (ALWAYS END WITH /)
%
% This function imports all the data files from one day (located in the same
% folder) and arranges them into experiments, stimuli and flies.


cd(path) % Goes to the desired folder
wildcard = ['*E', num2str(expnumber), '*']; % Searches particularly for the experiment desired
filelist = dir(wildcard); % Takes the list of the files containing an 'E' in their name
clear wildcard

n = 1; % Counter to place the output trials
% prog = length(filelist); % Counter for progressbar
% progressbar; % Initializes progressbar

reverseStr = ''; % Preparation for progress report

for file = 1:length(filelist); % Iterates through each file
        
    name = filelist(file).name; % Takes the filename
    
    namecontrol = strsplit(name, '_'); % Controls that the name has more than 3 elements                                  % (which means it is not an info .txt file)                 
    if length(namecontrol) < 4 % In case it has 3 elements, it skips that file because it is an info file
        continue
    end
    clear namecontrol
       
    % Opens the file and reads x, y, theta, experimental information,
    % analog and time data
    trial = openbindata_prejul2017([path, name]); % For normal stimuli
%     trial = openbindata_rand([path, name]); % For random stimuli
%     display([trial.date, ' ', trial.experiment, ' ', trial.timestamp]);
    
    
    % A BIT OF BASIC PROCESSING OF THE SIGNAL, to get rid of tracking errors
    
    % ELIMINATES THE FIRST POINT WHICH IS USUALLY AN ERROR
    % (see notes on November 25th, 2015)
    trial.x(1,:) = [];
    trial.y(1,:) = [];
    trial.theta(1,:) = [];
    
    
    % ELIMINATES TRACKING ERRORS
    bads = []; % Empty matrix to store the marks
    for i = 1:size(trial.x,2) % In case any of the conditions is met, it remembers the column
        % In case it contains zeros (tracking error)
        if ~isempty(find(trial.x(:,i)==0)) 
            bads = [bads, i];
        % In case there is a sudden change in coordinates (tracking error)
        elseif find(abs(diff(trial.x(:,i)))>=10 | abs(diff(trial.y(:,i)))>=10);
            bads = [bads, i];
        end
    end
    trial.x(:,bads) = []; trial.y(:,bads) = []; trial.theta(:,bads) = []; % Eliminates the bad tracks
    trial.fliesid = 1:4; trial.fliesid(bads) = []; % Creates a vector specifying which flies are left

    % In case there are no data left it continues importing the next file
    if isempty(trial.x)
        display ('--- Trial discarded because all flies contained a tracking error ---');
        continue
    end
    
    % FIXES ERRORS IN ORIENTATION DETECTION
for j = 1:size(trial.theta,2) % Iterates through each column
    for pair = 1:size(trial.theta,1)-1 % Iterates through each pair of points
        trial.theta(pair+1,j) = correctorientation(trial.theta(pair,j),trial.theta(pair+1,j));        
    end
    trial.theta(:,j) = correcttheta(trial.x(:,j), trial.y(:,j), trial.theta(:,j));
end

    out(n) = trial; % Stores the trial for output
    id = [trial.date, ',', trial.experiment];
    clear trial
    
    n = n+1; % Counts up
    
    percentDone = 100 * file / length(filelist);
    msg = sprintf(['  ',id,' progress: %3.1f /100'], percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));


end

eee = exptimeline(out,7200); % Rearranges in chronological order and creates a 
                        % relative timeline for the experiment. ALSO if a
                        % time limit is provided as a second argument, it
                        % will ignore all data collected later than the
                        % given time limit.
clear out

% The following stupid turn of code is actually necessary for some reason I couldn't figure out
out = eee;
clear eee

end



