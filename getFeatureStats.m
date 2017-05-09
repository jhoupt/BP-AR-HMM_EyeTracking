function allMatrix = getFeatureStats(positionSeq, speedSeq, accSeq, stateSeq, samppersec)
allMatrix = [];
startState = 1;
currentState = stateSeq(1);

for i = 2:length(stateSeq)
    if stateSeq(i) ~= currentState || i== length(stateSeq)
     
        if i== length(stateSeq)
            endState = i;
        else
            endState = i-1;
        end
        
        duration = (1+endState-startState) / samppersec;
        stateDomain = startState:endState;
                
        ED = sqrt(sum((positionSeq(1:2,startState) - positionSeq(1:2,endState) ).^2)  );
        
        % Trajectory Length
        if( duration > 1/samppersec ) 
            dx = diff(positionSeq(1,stateDomain));
            dy = diff(positionSeq(2,stateDomain));
        elseif (endState-1) > 0 && (endState+1) < length(stateSeq)
            dx = mean(positionSeq(1,(endState-1):(endState+1)));
            dy = mean(positionSeq(2,(endState-1):(endState+1)));
        
        else
            dx =0; dy=0;
        end
        totalDist = sum(sqrt(dx.^2 + dy.^2));
        if duration > 3/samppersec
            [EM, EV] = simplepca(positionSeq(1:2,stateDomain)');
            
            % Larsson Dispersion
            lD = EV(1)/EV(2);
            
            % Larsson Consistency
            lPCD = ED/EV(2);
            
            % Larsson Positional Displacement
            lPD = ED/totalDist;
        else
            lD = NaN;
            lPCD = NaN;
            lPD = NaN;
        end
        
        % Larsson Spatial Range
        lPR = sqrt(sum( (max(positionSeq(1:2,stateDomain)') ...
                         - min(positionSeq(1:2,stateDomain)')).^2)) ;
        
        % Maximum speed reached
        maxSpeed = max( speedSeq(stateDomain) );

        % Maximum speed reached
        maxAcc = max( accSeq(stateDomain) );
        
        
        % Maximum dispersion
        maxDist = max(max(dist(positionSeq(1:2,stateDomain))));
        
        % MarkEye Unidirectionality
        %winpx = [ceil((endState-startState)/2), endState-startState];
        winp = ceil([50 100] ./ (1000/samppersec)); % window size in samples.
        if min(winp) <= 2
            winp = winp + (3-min(winp));
        end
        
        if endState-startState < winp(2)
            r = NaN;
%         elseif winpx(1)==2
%             [result,r] = pca_window(positionSeq(1:2,stateDomain),winpx(2));
        else
            [result,r] = pca_window(positionSeq(1:2,stateDomain),winp);
        end
        
        allMatrix = [allMatrix; [currentState, duration, totalDist, lD, lPCD, lPD, lPR, maxSpeed, maxDist, mean(r), maxAcc] ] ;
        currentState = stateSeq(i);
        startState = i;
    end
end
