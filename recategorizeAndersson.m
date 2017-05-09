
tempMu = zeros([nStates,3,nRetainedSamples]);

figure(state)
hold off 
for state = 1:(nStates-1)%[5 6 7];

 
    for j=1:nRetainedSamples
        idx = retainedSamples(j);
        tempMu(state,:,j) = CH.Psi(idx).theta(state).mu;
    end
    subplot(1,3,1)
    hold on
    histogram(tempMu(state, 1,:), 'Normalization', 'Probability');
    subplot(1,3,2)
    hold on
    histogram(tempMu(state, 2,:), 'Normalization', 'Probability');
    subplot(1,3,3)
    hold on
    histogram(tempMu(state, 3,:), 'Normalization', 'Probability');
end

% Recode Andersson polar input
if strcmp(parsename, 'POLAR')
    % Let Category 5 have mu(1) ~ 100
    state1 = 5;
    for idx = retainedSamples
        for state2=6:11
            fromMu = CHPOLAR.Psi(idx).theta(state2).mu(1);
            if fromMu < 105 && fromMu > 95
                CHPOLAR.Psi(idx) = swapStates(state1, state2, CHPOLAR.Psi(idx));
            end
        end
    end
    
    % Let Category 6 have all the mu =0
    state1 = 6;
    lBound = [ -4E-3; -4E-3; -3E-3 ] ;
    uBound = [  4E-3;  4E-3;  3E-3 ] ;
    for idx = retainedSamples
        for state2=7:11
            fromMu = CHPOLAR.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CHPOLAR.Psi(idx) = swapStates(state1, state2, CHPOLAR.Psi(idx));
            end
        end
    end
    
    % Let Category 7 have all the mu1 = 11.5
    state1 = 7;
    lBound = [ 11.2; 10.9; -1E-3 ] ;
    uBound = [ 12.0; 11.4;  1E-3 ] ;
    for idx = retainedSamples
        for state2=8:11
            fromMu = CHPOLAR.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CHPOLAR.Psi(idx) = swapStates(state1, state2, CHPOLAR.Psi(idx));
            end
        end
    end
    
    
    % Let Category 8 have all the mu1 = 345
    state1 = 8;
    lBound = [ 330; 42; -1E-3 ] ;
    uBound = [ 360; 45;  1E-3 ] ;
    for idx = retainedSamples
        for state2=7:11
            fromMu = CHPOLAR.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CHPOLAR.Psi(idx) = swapStates(state1, state2, CHPOLAR.Psi(idx));
            end
        end
    end
    
    
    
    
    % Let Category 9 have all the mu1 = 390
    state1 = 9;
    lBound = [ 350; 245; -1E-3 ] ;
    uBound = [ 420; 280;  1E-3 ] ;
    for idx = retainedSamples
        for state2=7:11
            fromMu = CHPOLAR.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CHPOLAR.Psi(idx) = swapStates(state1, state2, CHPOLAR.Psi(idx));
            end
        end
    end
    
    
    
    
    % Let Category 10 have all the mu1 = 9.5
    state1 = 10;
    lBound = [  9; 210; -1E-3 ] ;
    uBound = [ 10; 245;  1E-3 ] ;
    for idx = retainedSamples
        for state2=7:11
            fromMu = CHPOLAR.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CHPOLAR.Psi(idx) = swapStates(state1, state2, CHPOLAR.Psi(idx));
            end
        end
    end
    
    % Let Category 11 have all the mu3 = 1
    state1 = 11;
    lBound = [ -1E-3; -1E-3; 0.9 ] ;
    uBound = [  1E-3;  1E-3; 1.1 ] ;
    for idx = retainedSamples
        for state2=7:11
            fromMu = CHPOLAR.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CHPOLAR.Psi(idx) = swapStates(state1, state2, CHPOLAR.Psi(idx));
            end
        end
    end
    
    CH= CHPOLAR;
    
    
    
% Recode Andersson v/a input    
elseif strcmp(parsename, 'VA')
    % Let Category 1 have all the mu1 = 10.75
    state1 = 1;
    lBound = [  10.7; .013; -1E-3 ] ;
    uBound = [ 11;    .014;  1E-3 ] ;
    for idx = retainedSamples
        for state2=2:9
            fromMu = CH.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CH.Psi(idx) = swapStates(state1, state2, CH.Psi(idx));
            end
        end
    end
    % Let Category 2 have all the mu1 = 32
    state1 = 2;
    lBound = [ 32; .0192; -1E-3 ] ;
    uBound = [ 33; .0198;  1E-3 ] ;
    for idx = retainedSamples
        for state2=3:9
            fromMu = CH.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CH.Psi(idx) = swapStates(state1, state2, CH.Psi(idx));
            end
        end
    end    
    % Let Category 3 have all the mu1 = 82
    state1 = 3;
    lBound = [ 80; .041; -1E-3 ] ;
    uBound = [ 84; .043;  1E-3 ] ;
    for idx = retainedSamples
        for state2=4:9
            fromMu = CH.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CH.Psi(idx) = swapStates(state1, state2, CH.Psi(idx));
            end
        end
    end    
    % Let Category 4 have all the mu1 = 375
    state1 = 4;
    lBound = [ 365; .08; -1E-3 ] ;
    uBound = [ 380; .095;  1E-3 ] ;
    for idx = retainedSamples
        for state2=5:9
            fromMu = CH.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CH.Psi(idx) = swapStates(state1, state2, CH.Psi(idx));
            end
        end
    end    

   
    % Let Category 5 have all the mu1 = 1000
    state1 = 5;
    lBound = [ 800; .35;  -1 ] ;
    uBound = [1300; .41; .5 ] ;
    for idx = retainedSamples       
        for state2=6:9
            fromMu = CH.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CH.Psi(idx) = swapStates(state1, state2, CH.Psi(idx));
            end
        end
        fromMu = CH.Psi(idx).theta(state1).mu;
        if any(fromMu > uBound) || any(fromMu < lBound)
            CH.Psi(idx) = swapStates(state1, state1+1, CH.Psi(idx));
        end         
    end
    
    
    % Let Category 6 have all the mu3 = 1
    state1 = 6;
    lBound = [-.5; 1; .5 ] ;
    uBound = [ .5;  3; 2 ] ;
    for idx = retainedSamples
        for state2=5:9
            fromMu = CH.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CH.Psi(idx) = swapStates(state1, state2, CH.Psi(idx));
            end
        end
        fromMu = CH.Psi(idx).theta(state1).mu;
        if any(fromMu > uBound) || any(fromMu < lBound)
            CH.Psi(idx) = swapStates(state1, state1+1, CH.Psi(idx));
        end        
    end


    
    
    % Let Category 7 have all the mu1 = 0
    state1 = 7;
    lBound = [ -1E-3;-1E-3;  .2 ] ;
    uBound = [  1E-3; .3; 1 ] ;
    for idx = retainedSamples
        for state2=8:9
            fromMu = CH.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CH.Psi(idx) = swapStates(state1, state2, CH.Psi(idx));
            end
        end
        fromMu = CH.Psi(idx).theta(state1).mu;
        if any(fromMu > uBound) || any(fromMu < lBound)
            CH.Psi(idx) = swapStates(state1, state1+1, CH.Psi(idx));
        end 
    end 
 
    % Let Category 8 have all the mu3 <1
    state1 = 8;
    lBound = [ -1.5;-2;  -1 ] ;
    uBound = [  1.5; 2; .2 ] ;
    for idx = retainedSamples
        for state2=9
            fromMu = CH.Psi(idx).theta(state2).mu;
            if all(fromMu < uBound) && all(fromMu > lBound)
                CH.Psi(idx) = swapStates(state1, state2, CH.Psi(idx));
            end
        end

        fromMu = CH.Psi(idx).theta(state1).mu;
        if any(fromMu > uBound) || any(fromMu < lBound)
            CH.Psi(idx) = swapStates(state1, state1+1, CH.Psi(idx));
        end        
    end
end 