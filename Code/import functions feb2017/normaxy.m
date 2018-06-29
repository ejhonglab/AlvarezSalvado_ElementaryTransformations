
function out = normaxy (data)

% [x,y] = normaxy (data)
%
%   x, y -> Normalized X and Y coordinates
%
%   data -> Structure array containing a field 'x' and another field 'y',
%           both containing raw coordinates in pixels.
%
% This function takes a structure array containing behavioral data from
% flies, separates each behavioral rig and converts all X and Y
% coordinates to the real size of the arenas (39.5x140 mm as of March 31st
% of 2017). This saves further conversion and allows a better overlap of
% the different trajectories, and thus more precise analyses of flies'
% positions.

% Coordinates' limits for the borders of the arenas, as of october 2017
% RIG# = [x-minimal, x-maximal, y-minimal, y-maximal];
RIG1 = [11, 306, 8.5, 1026];
RIG2 = [11, 299, 11, 1004];
RIG3 = [11, 301, 8, 1019];



% Identifies and separates each rig's data --------------------------
rig = zeros(length(data),1);
for i = 1:length(data) % Iterates through each fly
    info = strsplit(data(i).genotype,','); % Separates the pieces of information in the 'genotype' field
    riginfo = strrep(info{end}, ' ', ''); % Eliminates any spaces in the last piece of info
    riginfo = lower(riginfo); % Makes the string only lowercase
    if strcmp(riginfo,'rig1') % If the last piece of info is 'rig1'
        rig(i) = 1;
    elseif strcmp(riginfo,'rig2') % If the last piece of info is 'rig2'
        rig(i) = 2;
    else
        rig(i) = 3;
    end
end

if ~isempty(find(rig==1))
    x1 = horzcat(data(find(rig==1)).x);
    y1 = horzcat(data(find(rig==1)).y);        
end
if ~isempty(find(rig==2))
    x2 = horzcat(data(find(rig==2)).x);
    y2 = horzcat(data(find(rig==2)).y);
end
if ~isempty(find(rig==3))
    x3 = horzcat(data(find(rig==3)).x);
    y3 = horzcat(data(find(rig==3)).y);
end




% Checks if there are enough data ----------------------------------------

% Automatic check (october 12th 2017) %%%%%%%%%%%%%%%%%
errorinterval = 3; % Amount of pixels to allow as a confidence interval for
                   % automatic detection
% RIG1 ---------------------------
if ~isempty(find(rig==1))
    % Counts error interval provided
    RIG1 = [RIG1(1) + errorinterval, RIG1(2) - errorinterval,...
        RIG1(3) + errorinterval, RIG1(4) - errorinterval];
    % Checks if data are within the limits established
    if ~isempty(find(x1<=RIG1(1))) & ~isempty(find(x1>=RIG1(2)))
        if ~isempty(find(y1<=RIG1(3))) & ~isempty(find(y1>=RIG1(4)))
            disp(['There are enough data in Rig 1 (ci=',num2str(errorinterval),')'])
            decision1 = 1;
        else
            disp(['There are NOT enough data in Rig 1 (ci=',num2str(errorinterval),')'])
            decision1 = 0;
        end
    else
        disp(['There are NOT enough data in Rig 1 (ci=',num2str(errorinterval),')'])
        decision1 = 0;
    end
end

% RIG2 -----------------------------
if ~isempty(find(rig==2))
    % Counts error interval provided
    RIG2 = [RIG2(1) + errorinterval, RIG2(2) - errorinterval,...
        RIG2(3) + errorinterval, RIG2(4) - errorinterval];
    % Checks if data are within the limits established
    if ~isempty(find(x2<=RIG2(1))) & ~isempty(find(x2>=RIG2(2)))
        if ~isempty(find(y2<=RIG2(3))) & ~isempty(find(y2>=RIG2(4)))
            disp(['There are enough data in Rig 2 (ci=',num2str(errorinterval),')'])
            decision2 = 1;
        else
            disp(['There are NOT enough data in Rig 2 (ci=',num2str(errorinterval),')'])
            decision2 = 0;
        end
    else
        disp(['There are NOT enough data in Rig 2 (ci=',num2str(errorinterval),')'])
        decision2 = 0;
    end
end

% RIG3 ----------------------------
if ~isempty(find(rig==3))
    % Counts error interval provided
    RIG3 = [RIG3(1) + errorinterval, RIG3(2) - errorinterval,...
        RIG3(3) + errorinterval, RIG3(4) - errorinterval];
    % Checks if data are within the limits established
    if ~isempty(find(x3<=RIG3(1))) & ~isempty(find(x3>=RIG3(2)))
        if ~isempty(find(y3<=RIG3(3))) & ~isempty(find(y3>=RIG3(4)))
            disp(['There are enough data in Rig 3 (ci=',num2str(errorinterval),')'])
            decision3 = 1;
        else
            disp(['There are NOT enough data in Rig 3 (ci=',num2str(errorinterval),')'])
            decision3 = 0;
        end
    else
        disp(['There are NOT enough data in Rig 3 (ci=',num2str(errorinterval),')'])
        decision3 = 0;
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%



% Supervised check --------------------
if ~isempty(find(rig==1))
subplot(1,3,1);
plot(downsample(x1,10),downsample(y1,10));
title('RIG 1');
end
if ~isempty(find(rig==2))
subplot(1,3,2);
plot(downsample(x2,10),downsample(y2,10));
title('RIG 2')
end
if ~isempty(find(rig==3))
subplot(1,3,3);
plot(downsample(x3,10),downsample(y3,10));
title('RIG 3')
end

decision = input('ARE THERE ENOUGH TRAJECTORIES? (0 if NO; anything otherwise): ');
close;





% Normalizes --------------------------------------------------------------

if decision == 0 % If there are not enough data, it adds points based on real
                 % data to scale appropriately. These points are later
                 % removed.
    if ~isempty(find(rig==1)) & decision1 == 0
        x1 = [x1, [10; 307; zeros(size(x1,1)-2,1)]];
        y1 = [y1, [7.5; 1027; zeros(size(y1,1)-2,1)]];
    end
    if ~isempty(find(rig==2)) & decision2 == 0
        x2 = [x2, [10; 300; zeros(size(x2,1)-2,1)]];
        y2 = [y2, [10; 1008; zeros(size(y2,1)-2,1)]];
    end 
    if ~isempty(find(rig==3)) & decision3 == 0
        x3 = [x3, [10; 302; zeros(size(x3,1)-2,1)]];
        y3 = [y3, [7; 1020; zeros(size(y3,1)-2,1)]];
    end   
end

% Does the actual normalization
if ~isempty(find(rig==1))
x1 = norm2(x1,0,39.5);
y1 = norm2(y1,0,140);
end
if ~isempty(find(rig==2))
x2 = norm2(x2,0,39.5);
y2 = norm2(y2,0,140);
end
if ~isempty(find(rig==3))
x3 = norm2(x3,0,39.5);
y3 = norm2(y3,0,140);
end

if decision == 0 % It removes the fake points introduced
    if ~isempty(find(rig==1)) & decision1 == 0
        x1(:,end) = [];
        y1(:,end) = [];
    end
    if ~isempty(find(rig==2)) & decision2 == 0
        x2(:,end) = [];
        y2(:,end) = [];
    end
    if ~isempty(find(rig==3)) & decision2 == 0
        x3(:,end) = [];
        y3(:,end) = [];
    end
end
    


% Re-locates the coordinates to their original places ---------------------
out = data;
if ~isempty(find(rig==1))
is = find(rig==1); % Finds which flies belong to Rig1
for i = 1:numel(is) % Iterates through each Rig1 fly
    mosca = is(i); % Takes the current fly's position in the original data structure
    ntrials = size(data(mosca).x,2); % Calculates the amount of trials that there are from that fly
    % Substitutes the corresponding normalized coordinates for the originals
    out(mosca).x = x1(:,1:ntrials);
    out(mosca).y = y1(:,1:ntrials);
    % Eliminates those coordinates already returned from the normalized pool
    x1(:,1:ntrials) = [];
    y1(:,1:ntrials) = [];
end
end
% REPEATS THE SAME OPERATION JUST DONE, BUT FOR RIG2
if ~isempty(find(rig==2))
is = find(rig==2); % Finds which flies belong to Rig2
for i = 1:numel(is) % Iterates through each Rig2 fly
    mosca = is(i); % Takes the current fly's position in the original data structure
    ntrials = size(data(mosca).x,2); % Calculates the amount of trials that there are from that fly
    % Substitutes the corresponding normalized coordinates for the originals
    out(mosca).x = x2(:,1:ntrials);
    out(mosca).y = y2(:,1:ntrials);
    % Eliminates those coordinates already returned from the normalized pool
    x2(:,1:ntrials) = [];
    y2(:,1:ntrials) = [];
end
end
if ~isempty(find(rig==3))
% REPEATS THE SAME OPERATION JUST DONE, BUT FOR RIG3
is = find(rig==3); % Finds which flies belong to Rig3
for i = 1:numel(is) % Iterates through each Rig3 fly
    mosca = is(i); % Takes the current fly's position in the original data structure
    ntrials = size(data(mosca).x,2); % Calculates the amount of trials that there are from that fly
    % Substitutes the corresponding normalized coordinates for the originals
    out(mosca).x = x3(:,1:ntrials);
    out(mosca).y = y3(:,1:ntrials);
    % Eliminates those coordinates already returned from the normalized pool
    x3(:,1:ntrials) = [];
    y3(:,1:ntrials) = [];
end
end

    














