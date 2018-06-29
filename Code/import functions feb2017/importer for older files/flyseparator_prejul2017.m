
function out = flyseparator_prejul2017 (data)

% out = flyseparator (data)
%
%   data -> Structure array with trials in the rows.
%
% This function takes a structure array containing trials of behavioral
% data in the rows, rearranges the data in chronological order and separated
% by fly, and returns a structure array with a fly in each row.
% IMPORTANT: The program assumes there are 4 flies. To change this, modify
% the first line of the program.

nflies = 4; % CHANGE FOR DIFFERENT NUMBER OF FLIES

data = chronosort(data); % The data source is rearranged in chronological order

for fly = 1:nflies % Iterates through each fly 
    
    % Prepares the structure for the fly
    mosca.fly = ['Fly',num2str(fly)];
    mosca.date = data(1).date;
    mosca.experiment = data(1).experiment;
    mosca.timestamp = cell(1); mosca.timestamp(1) = []; % Creates an empty cell to store timestamps
    mosca.genotype = data(1).genotype;
    mosca.conditions = data(1).conditions;
    mosca.mode = cell(1); mosca.mode(1) = []; % Creates an empty cell to store stimulation modes
    mosca.modepre = cell(1); mosca.modepre(1) = []; % Does the same again
    mosca.timeline = [];
    mosca.time = [];
    mosca.stimulus = [];
%     mosca.windstimulus = [];
%     mosca.analog = [];
    mosca.x = [];
    mosca.y = [];
    mosca.theta = [];
%     mosca.xfilt = [];
%     mosca.yfilt = [];
%     mosca.dx = [];
%     mosca.dy = [];
% %     mosca.heading = [];
%     mosca.headingd = [];
%     mosca.delta = [];
%     mosca.deltamm = [];
%     mosca.v = [];
%     mosca.vx = [];
%     mosca.vy = [];
%     mosca.cdeltamm = [];
%     mosca.polartheta = [];
%         
    
    for i = 1:length(data) % Iterates through each of the trials in 'data'
    
        ind = find(data(i).fliesid == fly); % Finds whether the fly exists
                                            % in the current trial and
                                            % which is its column in the
                                            % data matrices
        
        if  isempty(ind) == 0 % If the current fly is actually present in the trial
                                                % the function takes the data
            mosca.timestamp = [mosca.timestamp; data(i).timestamp];
            mosca.mode = [mosca.mode; data(i).mode];
            mosca.modepre = [mosca.modepre; data(i).modepre];
            mosca.timeline = [mosca.timeline; data(i).relativetime];
            mosca.time = [mosca.time, data(i).time];
            mosca.stimulus = [mosca.stimulus, data(i).stimulus];
%             mosca.windstimulus = [mosca.windstimulus, data(i).windstimulus];
%             mosca.analog = [mosca.analog, data(i).analog];
            mosca.x = [mosca.x, data(i).x(:,ind)];
            mosca.y = [mosca.y, data(i).y(:,ind)];
            mosca.theta = [mosca.theta, data(i).theta(:,ind)];
%             mosca.xfilt = [mosca.xfilt, data(i).xfilt(:,ind)];
%             mosca.yfilt = [mosca.yfilt, data(i).yfilt(:,ind)];
%             mosca.dx = [mosca.dx, data(i).dx(:,ind)];
%             mosca.dy = [mosca.dy, data(i).dy(:,ind)];
% %             mosca.heading = [mosca.heading, data(i).heading(:,ind)];
%             mosca.headingd = [mosca.headingd, data(i).headingd(:,ind)];
%             mosca.delta = [mosca.delta, data(i).delta(:,ind)];
%             mosca.deltamm = [mosca.deltamm, data(i).deltamm(:,ind)];
%             mosca.v = [mosca.v, data(i).v(:,ind)];
%             mosca.vx = [mosca.vx, data(i).vx(:,ind)];
%             mosca.vy = [mosca.vy, data(i).vy(:,ind)];
%             mosca.cdeltamm = [mosca.cdeltamm, data(i).cdeltamm(:,ind)];
%             mosca.polartheta = [mosca.polartheta, data(i).polartheta(:,ind)];
            
        elseif isempty(ind) % If the current fly is not present in the trial
                                                   % it skips the trial
            continue
        end 
    end
    
    out(fly) = mosca;
    clear mosca
end

end