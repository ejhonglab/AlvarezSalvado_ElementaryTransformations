
function [dist, bins] = accumorientations (orients, values, edges)

% [dist, bins] = acumorientations (orients, values, edges)
%
%   orients -> Vector with the orientations from the sample.
%
%   values -> Vector with the values of the sample corresponding exactly 
%             to each of the orientations in 'orients'.
%
%   edges -> Limits of the bins in which to distribute the sample. This is
%            an inferior limit of the group and a superior limit too.
%            [Example: the edges 0,1,2 would produce TWO bins, one between
%            0 and 1, and another bin between 1 and 2].
%
% This function behaves as a custom/home made version of accumarray. It
% will take a sample of orientations associated with a sample of different
% data. It will group the orientations in bins according to the 'edges'
% provided, and finally will sum together the values of data associated to
% orientations in the same bin. It returns the summed values (in 'dist') and
% the middle points of the edges (to use them in representation with polar
% plot)(in 'bins').

% In case orients is longer, it will cut out the last point because it will
% usually be only one point bigger
if length(orients) > length(values)
    display('WARNING! ORIENTS VECTOR IS BIGGER THAN VALUES!');
    display('THE LAST VALUE WILL BE CUT OUT TO MATCH');
    orients(end) = [];
end

for i = 1:length(edges)-1 % Iterates through each bin
    
    % Takes the initial and end points of the bin
    a = edges(i);
    b = edges(i+1);
    
    % Looks for coincidences in the orientations vector and takes the
    % corresponding values
    matches = values(orients >= a & orients < b);
    
    % Controls in case all the matches have NaN value
    isitnan = isnan(matches);
    if isempty(find(isitnan==0)) % In case there are only NaN values, it makes it 0.
        matches = 0;
    end
    
    % Places the numbers in the output
    bins(i,1) = mean([a,b])*(2*pi/360); % Orientation in RADIANS
    dist(i,1) = nansum(matches);
    
end
clear orients values edges a b i matches
end