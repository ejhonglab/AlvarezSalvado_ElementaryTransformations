                                 
function [out, indices] = stimwarp (datavector, ypositionsvector, stimvector, stimode)

% This program takes a reference to a given trial in a behavioral oneflydata structure
% and generates indices to warp the coordinates back to a standard stimulus
% based on the PID measurements.

% If the stimulus is blank, returns unchanged data and original indices
if strcmp(stimode,'blank') || isempty(find(stimvector ~= 0,1))
    out = datavector;
    indices = 1:numel(datavector);
    return
end

stind = find(stimvector ~= 0); % Finds stimulus timing
onset = stind(1);
% Iterates through the different stimuli types and sets the corresponding delay
switch stimode
    case '10s', delay = 48; % Delay to traverse the arena starting upwind
    case '2s', delay = 48;
    case 'sweep', [~,delay] = getwindspeed(1:stind(end)-stind(1)+1,'sweep'); % Calls a function to calculate (see at the bottom of the code)
    case 'isweep', [~,delay] = getwindspeed(1:stind(end)-stind(1)+1,'isweep'); % Calls a function to calculate (see at the bottom of the code)
    case 'walkfar', delay = 74;
    otherwise, delay = 48;
end

wspeed = delay / 140; % Calculates the wind speed in samples/mm
offset = stind(end); % Sets the maximal amount of data to warp

yfrag = ypositionsvector(onset:offset); % Takes the stimulus part of the Y positions vector
% vels = wspeed(1:length(yfrag));

% Calculates the delay for each sample, multiplying the position of the fly
% by the wind speed
iii = round((140-yfrag) .* wspeed');
iii = sum([(1:length(iii)); iii'])'; % Adds to the delays a vector with a unit series.

% Extends the indices to make them applicable to the whole trial
indices = [ ((1:onset)+(iii(1)-1))';  iii+onset;  (iii(end)+onset+1:length(datavector))' ];

% WARPS and saves to output
out = zeros(length(datavector),1); out(out==0) = NaN; % Preallocates creating a NAN vector.
out(1:length(indices)) = datavector(indices); % Takes the data using the warped indices.
                                              % Data will be missing from the end of the trial.

end % end of function





% Subfunction to calculate the specific wind speed for some stimuli in
% which it changes over the length of the stimulus (namely frequency
% sweeps). These observations were made on PID measurements, and the fits
% for isweep stimuli are quite tailor made. This process can definitely be 
% simplified but it will be less precise.

function [windspeed,delay] = getwindspeed (samplesafteronset, stimode)

switch stimode
    
    case 'sweep'
        fit = [0.0144293162539767,60.3482464855815];
        delay = samplesafteronset*fit(1) + fit(2); % gets total delay to traverse arena
        windspeed = delay / 140; % Calculates the wind speed in samples/
        
    case 'isweep'
        for i = 1:length(samplesafteronset)
            if samplesafteronset(i) < 473
                delay(i) = 78;
                windspeed(i) = delay(i) / 140;
            elseif samplesafteronset(i) >= 473 & samplesafteronset(i) <= 743 
                fit = [-0.0185730649648175,86.8476517754868]; % Fit coefficients for the Reverse Sweep stimulus
                delay(i) = samplesafteronset(i)*fit(1) + fit(2);
                windspeed(i) = delay(i) / 140;
            elseif samplesafteronset(i) > 743
                delay(i) = 73;
                windspeed(i) = delay(i) / 140;
            end
        end
end
end % End of subfunction
    















