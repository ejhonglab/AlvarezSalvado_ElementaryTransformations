
function out = findtrialinfo (path)

% out = findtrialinfo (path)
%
%   path -> Absolute path of a file.
%
% This function gathers the experimental details of the given trial,
% looking for information both in the name of the file and in the
% accompanying info .txt file that usually goes with binary data files.
% IMPORTANT: The route to this info file is set by default. To change it
% modify the code.

% Info default path
% infopath = '/Users/efren/Data/'; % Modify if necessary

% Splits the name and takes the useful info
pathsplit = strsplit(path,'/'); % Splits the path
filename = pathsplit{end}; % Takes the last part
filenamesplit = strsplit(filename,'_'); % Splits the subparts of the name

%Stores the info
out.date = filenamesplit{end-2};
out.experiment = filenamesplit{end-1};
out.timestamp = filenamesplit{end};
out.genotype = 'information missing';
out.conditions = 'information missing';
% out.mode = lower([filenamesplit{1}]);
% switch filenamesplit{2}
%     case 'alwayson'
%         out.mode = lower([filenamesplit{1},'_windon']);
%     case 'alwaysoff'
%         out.mode = lower([filenamesplit{1},'_windoff']);
%     case 'blank'
%         out.mode = lower([filenamesplit{1},'_windoff']);
%     case '10swind'
%         out.mode = lower([filenamesplit{1},'_wind10s']);
% end

out.mode = lower(filenamesplit{1});
% out.mode = lower([filenamesplit{1},'_',filenamesplit{2}]);

% nflies = strsplit(filenamesplit{4},'n'); nflies = nflies{end};
% out.nflies = str2num(nflies);

% Looks for the rest of the info
infofile = ['info_',out.date,'_',out.experiment,'.txt']; % Constructs the name of the corresponding info file

% Constructs the path to the folder where the trial is
pathsplit(end) = []; pathsplit(1) = []; % Removes unnecessary parts of the initial path
infopath = '/';
for i = 1:length(pathsplit)
    infopath = [infopath,pathsplit{i},'/'];
end

% Opens the file and takes the rest of the info
fid = fopen([infopath,infofile]);
if fid < 3 % In case it doesn't find the info file, it reports
    display('Info file could not be found. Some information will be missing')
else
    garbage = fgetl(fid);
    garbage = fgetl(fid);
    out.genotype = fgetl(fid);
    out.conditions = fgetl(fid);
%     out.nflies = fgetl(fid);
%     out.nflies = str2num(out.nflies);
end

fclose(fid); % Closes the file

clear pathsplit filename filenamesplit infofile infopath fid garbage

end