
function out = correctorientation (a,b)
        
        dif = b-a; % Calculates the difference
        
        if abs(dif) > 220 % In case it is larger than a 200 degrees threshold
            if dif > 0 % In case the next point is larger...
                out = b - 360; % ...it subtracts 360 degrees
            elseif dif < 0 % In case the next point is smaller...
                out = b + 360; % ...it adds 360 degrees
            end 
            
        elseif abs(dif) > 45 % In case it is larger than a 45 degrees threshold
            if dif > 0 % In case the next point is larger...
                out = b - 180; % ...it subtracts 180 degrees
            elseif dif < 0 % In case the next point is smaller...
                out = b + 180; % ...it adds 180 degrees
            end
            
        else
            out = b;
            
        end
        
        % In case the final result exceeds 360 or is negative, corrects it 
        out = wrapTo360(out);
        
        clear dif a b
end