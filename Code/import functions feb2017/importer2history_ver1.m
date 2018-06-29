
function out = importer2history (varargin)

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
    if str2double(refs(i).date) >= 170701
        ex = importthatfolder (path,refs(i).expnumber);
    elseif str2double(refs(i).date) < 170701
        ex = importthatfolder_prejul2017 (path,refs(i).expnumber);
    end
    
    fprintf('\n'); % Produces a new empty line of text

    % Creates a new field in the structures with the stimulus mode of the
    % preceeding trial.
    ex(1).modepre = [];
    for l = 1:length(ex)
        if l > 1 
            ex(l).modepre = ex(l-1).mode;
        elseif l == 1
            ex(l).modepre = [];
        end
    end

    % Separates flies
    X = flyseparator(ex);
    clear ex
    
    % Renames the structures so they don't overwrite
    eval(['X', num2str(i), '=X;']);
    clear X
    
end


% Concatenates together the different experiments
out = [];
for grupo = 1:i
    out = eval(['[out, X', num2str(grupo),'];']);
end

% Checks for empty structures (flies without valid data)
marks = [];
for ch = 1:length(out)
    if isempty(out(ch).x)
        marks = [marks, ch];
    end
end
out(marks) = []; % Eliminates empty structures
if ~isempty(marks)
    display(['Flies with numbers ', num2str(marks), 'were discarded due to tracking errors']);
end


out = normaxy(out); % NORMALIZES THE COORDINATES TO 0-140 mm


% Analyzes
out = analyzer_unfilt(out);


% Clears 'X'-named structures
for j = 1:length(refs) % Iterates through each separate experiments...
    eval(['clear X', num2str(j)]);
end


end % End of function





