
function val = correcttheta (x, y, theta)

% out = correcttheta (theta,heading)
%
%   theta -> Vector with the orientations of the fly obtained by LabView.
%
%   heading -> Vector with the orientations of the fly calculated based on
%              its movements.
%
% This functions compares the theta vector (and all its possible symmetric
% variations) with the heading vector, and returns the one that best
% correlates with it. 


% Filters the trajectories
[b, a] = butter(2,2.5/25);
xfilt = filtfilt(b,a,x);
yfilt = filtfilt(b,a,y); 
heading = atan2d(diff(yfilt,1,1),diff(xfilt,1,1));
negs = heading < 0;
heading(negs) = heading(negs) + 360; % Converts the values to a 0:360 scale
clear b a xfilt yfilt negs

% Creates the symmetric variations of theta
theta1 = theta;

theta2 = 180 - theta; % Symetric with respect to Y axis
theta2(theta2<0) = theta2(theta2<0) + 360;

theta3 = 360 - theta; % Symmetric with respect to X axis

theta4 = 360 - theta2; % Inverted in both axes

% Makes the lengths of the vectors equal
if length(heading) < length(theta)
    thet1 = theta1; thet1(end) = [];
    thet2 = theta2; thet2(end) = [];
    thet3 = theta3; thet3(end) = [];
    thet4 = theta4; thet4(end) = [];
else
    thet1 = theta1;
    thet2 = theta2;
    thet3 = theta3;
    thet4 = theta4;
end


% Calculates the correlations between all theta vectors and the heading
corrs = [crosscorr(thet1,heading,50), crosscorr(thet2,heading,50), ...
    crosscorr(thet3,heading,50), crosscorr(thet4,heading,50)];
[~,indmax] = max(max(corrs));
eval(['eee = theta', num2str(indmax), ';']);
val = eee;


% In case any value exceeds 360 or is below 0, it is corrected again
val = wrapTo360(val);



clear eee corrs indmax theta1 theta2 theta3 theta4 theta
end












