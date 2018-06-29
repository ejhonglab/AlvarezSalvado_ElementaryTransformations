
function out = a_vels (resdata)

% out = a_vels (resdata)
%
%   resdata -> Structure array with flies' data. It must contain "xfilt",
%   "yfilt" and "thetafilt" fields. If not, run first through the
%   "filtering.m" function.
%
% Calculates the velocities (linear and angular) of the flies.


out = resdata;

% Iterates through each fly
for fly = 1:length(out)
    
    % Preallocates/resets variables to make sure it overwrites
    ntrials = size(out(fly).xfilt,2); 
    nrows = size(out(fly).xfilt,1);
    out(fly).uthetafilt = zeros(nrows,ntrials); 
    out(fly).deltamm = zeros(nrows-1,ntrials); 
    out(fly).v = zeros(nrows-1,ntrials);
    out(fly).vy = zeros(nrows-1,ntrials);
    out(fly).vx = zeros(nrows-1,ntrials);
    out(fly).angv = zeros(nrows-1,ntrials);
    out(fly).pmove = zeros(nrows-1,ntrials); 
    out(fly).vmove = zeros(nrows-1,ntrials);
    out(fly).pturn = zeros(nrows-1,ntrials);
    out(fly).angvturn = zeros(nrows-1,ntrials);
    out(fly).curvature = zeros(nrows-1,ntrials);
    
    % Unwraps the orientations
    out(fly).uthetafilt = unwrap(out(fly).thetafilt*(pi/180)) *(180/pi); % unwrap(out(fly).thetafilt * (pi/180)) * (180/pi);
    
    % Calculates the distance moved in each frame (in mm)
    out(fly).deltamm = hypot(diff(out(fly).yfilt,1,1),diff(out(fly).xfilt,1,1)); % Distance in mm traveled in each frame 
%     out(fly).deltamm = hypot(diff(out(fly).yfilt),diff(out(fly).xfilt)); % USE FOR WINDTUNNEL DATA INSTEAD;
    
    
    % Calculates radial and angular velocities
    tes = []; tes = diff(out(fly).time(1:nrows,out(fly).selected)); % Takes the exact time intervals of each frame
%     tes = []; tes = diff(out(fly).timeStamp(1:nrows,:)); % USE FOR WINDTUNNEL DATA Takes the exact time intervals of each frame
    out(fly).v = out(fly).deltamm ./tes; % Radial velocity in mm/s
    out(fly).vy = diff(out(fly).yfilt,1,1) ./tes; % * cfactor; % Upwind velocity in mm/s
    out(fly).vx = diff(out(fly).xfilt,1,1) ./tes; % * cfactor; % Crosswind velocity in mm/s
%     out(fly).vy = -diff(out(fly).yfilt) ./tes; % USE FOR WINDTUNNEL DATA % Upwind velocity in mm/s
%     out(fly).vx = diff(out(fly).xfilt) ./tes; % USE FOR WINDTUNNEL DATA % Crosswind velocity in mm/s
    out(fly).angv = abs(diff(out(fly).uthetafilt)) ./tes; % UNSIGNED Angular velocity in degrees/s
%     out(fly).angv = diff(out(fly).uthetafilt) ./tes; % SIGNED Angular velocity in degrees/s
    
    % Calculates thresholded movement and turning
    out(fly).pmove = double(out(fly).v > 1); % Finds where flies move faster than 1 mm/s
%     pstop = out(fly).v < 1; % Finds where flies move slower than 1 mm/s
    out(fly).vmove = out(fly).v; out(fly).vmove(out(fly).v < 1) = NaN; % Ground speed thresholded over 1 mm/s
    out(fly).vymove = out(fly).vy; out(fly).vymove(out(fly).v < 1) = NaN; % Upwind velocity thresholded over 1 mm/s
    out(fly).angvturn = out(fly).angv; out(fly).angvturn(out(fly).v < 1) = NaN; % Angular velocity thresholded over 10 deg/s
    out(fly).curvature = out(fly).angvturn ./ out(fly).vmove; % Thresholded angular velocity OVER thresholded ground speed
    out(fly).pturn = double(out(fly).curvature > 20); % Finds where flies turn faster than 10 deg/s
end
    
    
    
    