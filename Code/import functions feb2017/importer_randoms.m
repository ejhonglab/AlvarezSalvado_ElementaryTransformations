
function out = importer_randoms (varargin)

% data = importer (date1, expnumber1, date2, expnumber2, ...)
%
%   date -> String with the date of the experiment (in format YYMMDD)
%
%   expnumber -> Number of the experiment you want from among the
%                experiments of that date.
%
% This function extracts an experiment from the dataset of a desired date,
% and returns the data already separated in the different stimulus and
% flies and analyzed.


% Extracts the dates and the expnumbers
n = 1;
for i = 2:2:length(varargin)
    refs(n).date = varargin{i-1};
    refs(n).expnumber = varargin{i};
    n = n+1;
end

% Imports and separates the data of the folder
for i = 1:length(refs)
    
    % Creates the path
    path = ['/Users/efren/Data/', refs(i).date, '/']; 
    
    % Imports all the trials from the specified folder and experiment number
    ex = importthatfolder (path,refs(i).expnumber);

    % Creates a new field in the structures with the stimulus mode of the
    % preceeding trial.
    ex(1).modepre = [];
    for l = 2:length(ex)
        ex(l).modepre = ex(l-1).mode;
    end

    % Separates in the different stimuli
    stimseparator (ex, 'caller')
    clear ex

    % Separates the fly for each stimulus and stores in output
    for j = 1:length(Stimulinames); % Iterates through each different stimulus
         name = ['s', Stimulinames{j}];    
         X.(name) = eval(['flyseparator (s', Stimulinames{j}, ');']);
         eval(['clear ', 's', Stimulinames{j}]);
    end
    
    % Renames the structures so they don't overwrite
    eval(['X', num2str(i), '=X;']);
    clear X
    
end


% Concatenates together the data from different experiments and stimuli
% and ANALYZES THE DATA ===============================================
for i = 1:length(Stimulinames) % Iterates through the stimuli
    Stimulinames{i} = lower(Stimulinames{i});
    stimname = ['s', Stimulinames{i}]; % Generates a name for the current stimulus
    eval([stimname, '= [];']); % Empties the existing matrix with that name (used in a previous step, not necessary anymore)
    for j = 1:length(refs) % Iterates through each separate experiments...
        eval([stimname, '= [', stimname, ', X', num2str(j), '.', stimname, '];' ]); % ...and concatenates together the parts with the current stimulus
    end

    
    current = eval(stimname); % Creates a copy of the structure array with a generic name
    eval(['clear ', stimname]);
    
    % Checks for empty structures (flies without valid data)
    marks = [];
    for ch = 1:length(current)
        if isempty(current(ch).x)
            marks = [marks; ch];
        end
    end
    current(marks) = []; % Eliminates empty structures
    if ~isempty(marks)
        display(['Flies with numbers ', marks, 'were discarded due to tracking errors']);
    end
    
    current = normaxy(current); % NORMALIZES THE COORDINATES TO 0-140 mm
    
%     out.(stimname) = current;
%     clear current

    
    % Analyzes
%     [av,r, aav] = analyzer(current,'all');
    r = analyzer(current);
    clear current
%     eval(['clear out.', stimname]);
%     
%     % Renames and organizes the results in the output structure
    out.(['res', Stimulinames{i}]) = r;
    clear r
%     out.(['avres', Stimulinames{i}]) = av;
%     out.(['aavres', Stimulinames{i}]) = aav;
%     clear r av aav
     
end

% Clears X structures
for j = 1:length(refs) % Iterates through each separate experiments...
    eval(['clear X', num2str(j)]);
end

end % End of function





