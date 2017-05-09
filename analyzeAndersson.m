%% Script to anlayze Andersson data with BP-HMM
runNew = false;
addpath(genpath('/Users/jhoupt/MATLAB/NPBayesHMM'));
addpath(genpath('/Users/jhoupt/Documents/git/saliency'));

%parsename = 'POLAR';
%parsename = 'L';
%parsename = 'S';
%parsename = 'V';
parsename = 'VA';


dataname = 'Andersson';
%dataname = 'DIEM';

%% Load Data

switch dataname
    case 'Andersson'
        [dataV, dataS, dataP, dataVA, dataPOLAR] = loadAndersson();
        samppersec = 500;
        standardStates = 1:6;
        lastDatePolar = '27-Mar-2017';
    case 'DIEM'
        [dataV, dataS, dataP, dataVA, dataPOLAR] = loadDIEM();
        samppersec = 30;
        standardStates = 1:3;
        lastDatePolar = '27-Mar-2017';
end

%% Rund BP-AR-HMM (or load from file) %%

% Run Location analysis
if strcmp(parsename, 'L')
    if runNew
        modelL = {'bpM.gamma', 1, 'bpM.c', .25, 'hmmM.kappa', 25};
        algL   = {'Niter', 5000, 'HMM.doSampleHypers',1,'BP.doSampleMass',1,'BP.doSampleConc',1};
        % Start out with just one feature for all objects
        initL  = {'F.nTotal', 10};
        CHL = runBPHMM( dataL, modelL, {1, 1}, algL, initL );
        save(strcat(dataname, '-CHL', date, '.mat'), 'CHL');
    else
        load(strcat(dataname, '-CHL', '25-Mar-2017.mat'));
    end
    
    CH = CHL;
    data = dataL;
end


% Run Speed analysis 
if strcmp(parsename, 'S')
    if runNew
        modelS = {'bpM.gamma', 1, 'bpM.c', .25, 'hmmM.kappa', 25};
        algS   = {'Niter', 5000, 'HMM.doSampleHypers',1,'BP.doSampleMass',1,'BP.doSampleConc',1};
        % Start out with just one feature for all objects
        initS  = {'F.nTotal', 10};
        CHS = runBPHMM( dataS, modelS, {1, 1}, algS, initS );
        save(strcat(dataname, '-CHS', date, '.mat'), 'CHS');
    else
        load(strcat(dataname, '-CHS', '22-Dec-2016.mat'))
    end    
    CH = CHS;
    data = dataS;
end


% Run Velocity analysis
if strcmp(parsename, 'V')
    if runNew
        modelV = {'bpM.gamma', 1, 'bpM.c', .25, 'hmmM.kappa', 25};
        algV   = {'Niter', 5000, 'HMM.doSampleHypers',1,'BP.doSampleMass',1,'BP.doSampleConc',1};
        initV  = {'F.nTotal', 10};
        CHV = runBPHMM( dataV, modelV, {1, 1}, algV, initV );
        save(strcat(dataname, '-CHV', date, '.mat'), 'CHS');
    else
        load(strcat(dataname, '-CHV', '22-Dec-2016.mat'))
    end
    CH = CHV;
    data = dataV;
end


% Run Speed/Acceleration analysis
if strcmp(parsename,'VA')
    if runNew
        modelVA = {'bpM.gamma', 1, 'bpM.c', .25, 'hmmM.kappa', 25};
        algVA   = {'Niter', 5000, 'HMM.doSampleHypers',1,'BP.doSampleMass',1,'BP.doSampleConc',1};
        initVA  = {'F.nTotal', 10};
        CHVA = runBPHMM( dataVA, modelVA, {1, 1}, algVA, initVA );
        save(strcat(dataname, '-CHVA', date, '.mat'), 'CHVA');
    else
        load(strcat(dataname, '-CHVA', '29-Dec-2016.mat'))
    end
    
    CH = CHVA;
    data = dataVA;
end


% Run Polar Speed analysis
if strcmp(parsename,'POLAR')
    if runNew
        modelPOLAR = {'bpM.gamma', 1, 'bpM.c', .25, 'hmmM.kappa', 25};
        algPOLAR   = {'Niter', 5000, 'HMM.doSampleHypers',1,'BP.doSampleMass',1,'BP.doSampleConc',1};
        initPOLAR  = {'F.nTotal', 10};
        CHPOLAR = runBPHMM( dataPOLAR, modelPOLAR, {1, 1}, algPOLAR, initPOLAR );
        save(strcat(dataname, '-CHPOLAR', date, '.mat'), 'CHPOLAR');
    else
        load(strcat(dataname, '-CHPOLAR', lastDatePolar, '.mat'))
    end
    
    CH = CHPOLAR;
    data = dataPOLAR;
end


%% Remove Burn in (first half of samples)
nSamples = length(CH.Psi)-2;
retainedSamples = (nSamples/2+3):(nSamples+2);
nRetainedSamples = length(retainedSamples);

%% Calculate posterior distribution over number of states
X = estimateStateProbability(CH, retainedSamples);
nStates = size(X,2);
% for i=1:nRetainedSamples
%     X(i,:) = sort(X(i,:));
% end
% Mean number of features used
mean(sum(X>0,2))

% Distribution over category probabilities
meanandquantiles = zeros(3,nStates);
meanandquantiles(1,:) = mean(X,1);
for i=1:nStates
    meanandquantiles(2,i) = quantile(X(:,i),.975);
    meanandquantiles(3,i) = quantile(X(:,i),.025);
end
csvwrite(strcat(dataname, '-DistributionCategoryProbability_', parsename, '.csv'), ...
    meanandquantiles);

%% Comparison between Standard Features and BP-AR-HMM 
nStandardStates = length(standardStates);
confusionProbability= zeros([nStandardStates nStates nRetainedSamples]);
for j=1:nRetainedSamples
    idx = retainedSamples(j);
    for sj = 1:data.N
        seqsj = (dataP.aggTs(sj)+1):dataP.aggTs(sj+1);
        tagseq = dataP.Xdata(3,seqsj);
        for st = 1:nStates
            for i=1:nStandardStates
                confusionProbability(i,st,j) = confusionProbability(i,st,j) +...
                    sum( tagseq==standardStates(i) & (CH.Psi(idx).stateSeq(sj).z == st));
            end
        end
    end
end
for j=1:(nRetainedSamples)
    confusionProbability(:,:,j) = confusionProbability(:,:,j) ./...
        repmat(sum(confusionProbability(:,:,j),2), [1 nStates]);
end
confusionProbability = mean(confusionProbability, 3);
csvwrite(strcat(dataname, '-ConfusionProbability_', parsename, '.csv'), ...
   confusionProbability);


%% Comparison between MarkEye Features and our feature
confusionProbabilityMarkEye = zeros([6 nStates]);
for sj = 1:data.N
    seqsj = (dataP.aggTs(sj)+1):dataP.aggTs(sj+1);
    tagseq = dataP.Xdata(4,seqsj);
    for st = 1:nStates
         for i=1:6
             confusionProbabilityMarkEye(i,st) = confusionProbabilityMarkEye(i,st) +...
               sum( tagseq==(i-1) & (CH.Psi(idx).stateSeq(sj).z == st));
         end
    end
end
confusionProbabilityMarkEye = confusionProbabilityMarkEye ./ ...
    repmat(sum(confusionProbabilityMarkEye,2), [1 nStates]);
csvwrite(strcat('ConfusionProbabilityMarkEye_', parsename, '.csv'), ...
    confusionProbabilityMarkEye);


%% Comparison between Andersson Feature and MarkEye Feature
confusionProbabilityAM = zeros([6 nStates]);
for sj = 1:data.N
    seqsj = (dataP.aggTs(sj)+1):dataP.aggTs(sj+1);
    tagseq = dataP.Xdata(3,seqsj);
    for j = 1:6
         for i=1:6
             confusionProbabilityAM(i,j) = confusionProbabilityAM(i,j) +...
               sum( tagseq==i & dataP.Xdata(4,seqsj)==(j-1)) ;
         end
    end
end
confusionProbabilityAM = confusionProbabilityAM ./ repmat(sum(confusionProbabilityAM,2), [1 nStates]);

%% Standard Set Stats
standardMatrix = [];
standardDurations = cell([1,nStandardStates]);
standardD = cell([1,nStandardStates]);
standardPCD = cell([1,nStandardStates]);
standardDists = cell([1,nStandardStates]);
standardPD = cell([1,nStandardStates]);
standardPR = cell([1,nStandardStates]);

standardMaxSpeed = cell([1,nStandardStates]);
standardMaxDists = cell([1,nStandardStates]);
standardMaxAcc = cell([1,nStandardStates]);

standardPCA = cell([1,nStandardStates]);


for j=1:size(dataP.Xdata,2)
    if j==1 
        beginState = j;
        currentState = dataP.Xdata(3,j);
    elseif ( dataP.Xdata(3,j) ~= currentState ) || any(dataP.aggTs==(j+1))
        endState = j-1;
        
        %Duration
        duration = (1+endState-beginState) / samppersec;
        standardDurations{currentState} = [ standardDurations{currentState} duration]; 

        
        ED = sqrt(sum((dataP.Xdata(1:2,beginState) - dataP.Xdata(1:2,endState)).^2));
        
        % Trajectory Length
        dx = diff(dataP.Xdata(1,beginState:endState) );
        dy = diff(dataP.Xdata(2,beginState:endState) );
        totalDist = sum(sqrt(dx.^2 + dy.^2));
        
        standardDists{currentState} = [ standardDists{currentState} totalDist];
        
        if duration > 3/samppersec
            [EM, EV] = simplepca(dataP.Xdata(1:2,beginState:endState)');
            
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

        standardD{currentState} = [standardD{currentState} lD];
        standardPCD{currentState} = [standardPCD{currentState} lPCD];
        standardPD{currentState} = [standardPD{currentState} lPD];
        
        % Larsson Spatial Range
        lPR = sqrt(sum( (max(dataP.Xdata(1:2,beginState:endState)')- min(dataP.Xdata(1:2,beginState:endState)')).^2)) ;
        standardPR{currentState} = [standardPR{currentState} lPR];
        
        % Maximum speed reached
        maxSpeed = max( dataS.Xdata(1,beginState:endState) );
        standardMaxSpeed{currentState} = [ standardMaxSpeed{currentState} maxSpeed];
        
        % Maximum acceleration reached
        maxAcc = max( dataVA.Xdata(2,beginState:endState) );
        standardMaxAcc{currentState} = [ standardMaxAcc{currentState} maxAcc];
        
        % Maximum dispersion
%         maxDist = max(sqrt(sum((dataP.Xdata(1:2,beginState:endState) -...
%             repmat(dataP.Xdata(1:2,beginState), [1,1+endState-beginState])).^2)));
        maxDist = max(max(dist(dataP.Xdata(1:2,beginState:endState))));
        standardMaxDists{currentState} = [ standardMaxDists{currentState}  maxDist ];
        
        % MarkEye Unidirectionality
        %winpx = [ceil((endState-beginState)/2), endState-beginState];
        winp = ceil([50 100] ./ (1000/samppersec)); % window size in samples.
        if min(winp) <= 2
            winp = winp + (3-min(winp));
        end
        
        if endState-beginState < winp(2)
            r = NaN;
%         elseif winpx(1)==2
%             [result,r] = pca_window(dataP.Xdata(1:2,beginState:endState),winpx(2));
        else
            [result,r] = pca_window(dataP.Xdata(1:2,beginState:endState),winp);
        end
        standardPCA{currentState} = [ standardPCA{currentState}  mean(r) ];
        
        standardMatrix = [standardMatrix; [currentState, duration, totalDist, lD, lPCD, lPD, lPR, maxSpeed, maxDist, mean(r), maxAcc] ] ;
        currentState = dataP.Xdata(3,j);
        beginState = j;
        
    end
end
csvwrite(strcat(dataname, '-StandardStats.csv'), standardMatrix);



%% BP-AR-HMM Stats
allMatrix = [];
allMatrixC = [];

%collapseTo1 = [2];
%collapseTo2 = [1 3];
%collapseTo3 = [4 5 7 8 9];
if strcmp(dataname, 'DIEM')
%     collapseTo1 = [1 2 6];
%     collapseTo2 = [3 5 7 8 9];
%     collapseTo3 = [4];
%      collapseTo1 = [2];
%      collapseTo2 = [1 3 6 8 9];
%      collapseTo3 = [5 7 9];
     collapseTo1 = [1 2];
     collapseTo2 = [3 4 5];
     collapseTo3 = [6 7];

elseif strcmp(dataname, 'Andersson')
%     collapseTo1 = [1 2 3];  
%     collapseTo2 = [4 7 10]; 
%     collapseTo3 = [5 8 9];  
%     collapseTo1 = [1 3 7 9 10];
%     collapseTo2 = [2];
%     collapseTo3 = [4 5 8];
      collapseTo1 = [1];
      collapseTo2 = [2 3 5];
      collapseTo3 = [4];

end


for sj = 1:dataP.N
    seqsj = (dataP.aggTs(sj)+1):(dataP.aggTs(sj+1));
    pos = dataP.Xdata(1:2,seqsj);
    spd = dataS.Xdata(1,seqsj);
    acc = dataVA.Xdata(2,seqsj);
    
    
    stt1Mat = zeros(nRetainedSamples, length(CH.Psi(1).stateSeq(sj).z));
    stt2Mat = zeros(nRetainedSamples, length(CH.Psi(1).stateSeq(sj).z));
    for j =1:nRetainedSamples
        idx=retainedSamples(j);
        CH.Psi(idx).stateSeq(sj).z(2,:) = CH.Psi(idx).stateSeq(sj).z(1,:);
        
        for i=collapseTo1
            CH.Psi(idx).stateSeq(sj).z(2,CH.Psi(idx).stateSeq(sj).z(1,:)==i) = 1;
        end
        for i=collapseTo2
            CH.Psi(idx).stateSeq(sj).z(2,CH.Psi(idx).stateSeq(sj).z(1,:)==i) = 2;
        end
        for i=collapseTo3
            CH.Psi(idx).stateSeq(sj).z(2,CH.Psi(idx).stateSeq(sj).z(1,:)==i) = 3;
        end
        
        
        stt1Mat(j,:) = CH.Psi(idx).stateSeq(sj).z(1,:);
        stt2Mat(j,:) = CH.Psi(idx).stateSeq(sj).z(2,:);
    end
    stt1 = mode(stt1Mat, 1);
    stt2 = mode(stt2Mat, 1);
    
    sjMat1 = getFeatureStats(pos, spd, acc, stt1, samppersec);
    sjMat2 = getFeatureStats(pos, spd, acc, stt2, samppersec);
    
    allMatrix = [allMatrix; sjMat1];
    allMatrixC = [allMatrixC; sjMat2];
    
end

csvwrite(strcat(dataname, '-bparhmmStats_', parsename, '.csv'), allMatrix);
csvwrite(strcat(dataname, '-bparhmmStatsCollapsed_', parsename, '.csv'), allMatrixC);

%%
%dat = allMatrix;
dat = allMatrixC;


durations = cell([1,nStates]);
totalDists = cell([1,nStates]);
lDispersions = cell([1,nStates]);
lConsistencies = cell([1,nStates]);
lDisplacements = cell([1,nStates]);
lSpatialRanges = cell([1,nStates]);
maxSpeeds = cell([1,nStates]);
maxDists = cell([1,nStates]);
meDispersions = cell([1,nStates]);

for i = 1:nStates
    durations{i} = dat(dat(:,1)== i,2);
    totalDists{i} = dat(dat(:,1)== i,3);
    lDispersions{i} = dat(dat(:,1)== i,4);
    lConsistencies{i} = dat(dat(:,1)== i,5);
    lDisplacements{i} = dat(dat(:,1)== i,6);
    lSpatialRanges{i} = dat(dat(:,1)== i,7);
    maxSpeeds{i} = dat(dat(:,1)== i,8);
    maxDists{i} = dat(dat(:,1)== i,9);
    meDispersions{i} = dat(dat(:,1)== i,10);
end

%% Outlier Detection

outlier = cell([1,nStates]);
for k =1:nStates
    if ~isempty(totalDists{k})
        outlier{k} = (maxSpeeds{k} > 3*iqr(maxSpeeds{k})) | ...
            (maxDists{k}  > 3*iqr(maxDists{k})) | ...
            (totalDists{k}  > 3*iqr(totalDists{k}));% | ...
            %(durations{k} == 1/samppersec);
%         outlier{k} = (maxSpeeds{k} > quantile(maxSpeeds{k}, .99)) | ...
%             (maxDists{k}  > quantile(maxDists{k}, .99)) | ...
%             (totalDists{k}  > quantile(totalDists{k}, .99)) | ...
%             (durations{k} == 1/samppersec);
    end
end


standardOutlier = cell([1,nStandardStates]);
for k =1:nStandardStates
    standardOutlier{k} = ...
        (standardMaxSpeed{k} > quantile(standardMaxSpeed{k}, .99)) | ...
        (standardMaxDists{k} > quantile(standardMaxDists{k}, .99)) | ...
        (   standardDists{k} > quantile(   standardDists{k}, .99)) ;
    standardOutlier{k} = ...
        (standardMaxSpeed{k} > iqr(standardMaxSpeed{k})) | ...
        (standardMaxDists{k} > iqr(standardMaxDists{k})) | ...
        (   standardDists{k} > iqr(   standardDists{k})) ;
end


%% Plotting 

plotStates = [1 2 3];
nPlot = length(plotStates);
figure

if ( nPlot > 3) 
    figrows = ceil(nPlot/3);
    figcols = 3;
else
    figrows = 1;
    figcols = nPlot;
end

% Maximum Distance
for i = 1:nPlot
    k = plotStates(i);
    subplot( figrows,figcols, i);
    if size(maxDists{k},1) > 0
        histogram(maxDists{k}(~outlier{k}), 100, 'Normalization', 'probability')
        title(strcat('State ', int2str(k)));
    end
    xlabel('Maximum Distance')
    %axis([0 15 0 0.22])
end

figure
hold on
for i = 1:nPlot
    k = plotStates(i);
    if size(maxDists{k},1) > 0
        histogram(maxDists{k}(~outlier{k}), 'BinWidth', .1, 'Normalization', 'probability')
    end
    %axis([0 15 0 0.22])
    xlabel('Maximum Distance')
end

figure
% Maximum Speed
for i = 1:nPlot
    k = plotStates(i);
    subplot( figrows,figcols, i);
    if size(maxSpeeds{k},1) > 0
        histogram(maxSpeeds{k}(~outlier{k}), 100, 'Normalization', 'probability')
        %histogram(maxSpeeds{k}, 100, 'Normalization', 'probability')
        title(strcat('State ', int2str(k)));
    end
    xlabel('Maximum Speed')
    %axis([0 15 0 0.22])
end

figure
hold on
for i = 1:nPlot
    k = plotStates(i);
    if size(maxSpeeds{k},1) > 0
        histogram(maxSpeeds{k}(~outlier{k}), 'BinWidth', 5, 'Normalization', 'probability')
    end
    xlabel('Maximum Speed')
    %axis([0 15 0 0.22])
end

% 2-D histogram
figure
for i = 1:nPlot
    k = plotStates(i);
    subplot( figrows,figcols, i);
    n = hist3([maxSpeeds{k}(~outlier{k}), maxDists{k}(~outlier{k})], ...
        [50 50]);
    n1 = n';
    n1(size(n,1) + 1, size(n,2) + 1) = 0;
    xn = linspace(min(maxSpeeds{k}(~outlier{k})),max(maxSpeeds{k}(~outlier{k})),size(n,1)+1);
    yn = linspace(min(maxDists{k}(~outlier{k})),max(maxDists{k}(~outlier{k})),size(n,2)+1);
    pcolor(xn, yn, n1)
    xlabel('Max Speed'); ylabel('Maximum Distance');
    %axis([1 20 1 20])
end



figure
for i = 1:nPlot
    k = plotStates(i);
    subplot( figrows,figcols, i);
    if sum(~isnan(meDispersions{k}(~outlier{k}))) > 5
        n = hist3([meDispersions{k}(~outlier{k}), maxDists{k}(~outlier{k})], ...
            [50 50]);
        n1 = n';
        n1(size(n,1) + 1, size(n,2) + 1) = 0;
        xn = linspace(min(meDispersions{k}(~outlier{k})),max(meDispersions{k}(~outlier{k})),size(n,1)+1);
        yn = linspace(min(maxDists{k}(~outlier{k})),max(maxDists{k}(~outlier{k})),size(n,2)+1);
        pcolor(xn, yn, n1)
        xlabel('MarkEye Dispersion'); ylabel('Maximum Distance');
    end
    %axis([1 20 1 20])
end

%%
oldCode = allMatrix(:,1);
%% Recode based on MarkEye cutoff
mecutoff = 0.02;
useCollapse = true;
oneFigure = false;

close all;
if useCollapse
    allMatrixC( allMatrixC(:,1)==4, 1) = 1;
    allMatrixC( allMatrixC(:,1)==1 & allMatrixC(:,10) > mecutoff) = 4;
    
    %colsB = ['b' 'r' 'k' 'g' 'r' 'k' 'k' 'k' 'k' 'k'];
    colsB = ['b' 'r' 'k' 'r' 'k' 'k'];
    bounds = int64([0; 500*cumsum(allMatrixC(:,2))]);
    
else
    allMatrix(:, 1) = oldCode;
    allMatrix( (allMatrix(:,1)==1 | allMatrix(:,1)==2 ) ...
        & allMatrix(:,10) > mecutoff) = 10;
    close all;
    
    bounds = int64([0; 500*cumsum(allMatrix(:,2))]);
    colsB = ['b' 'c' 'm' 'r' 'r' 'k' 'k' 'k' 'k' 'g'];
end

switch dataname
    case 'Andersson'
        colsA = ['b' 'r' 'm' 'g' 'k' 'k' ];
    case 'DIEM'
        colsA = ['k' 'b' 'r'];
end
i = 1;

for sj = 1:10% 43:53%23:33
    
    
    if oneFigure
        figure
        sjSeq = (1+dataP.aggTs(sj)):dataP.aggTs(sj+1);
        boundsA = [sjSeq(1)-1 sjSeq( diff(dataP.Xdata(3,sjSeq)) ~=0) sjSeq(end)];
        
        plot(dataP.Xdata(1,sjSeq), dataP.Xdata(2,sjSeq), 'k', 'LineWidth', 3);
        hold on
        j = 1;
        while boundsA(j) < dataP.aggTs(sj+1)
            ftSeq = (boundsA(j)+1):boundsA(j+1);
            plot(dataP.Xdata(1,ftSeq), dataP.Xdata(2,ftSeq), ...
                'Color', colsA(dataP.Xdata(3,ftSeq(1)+1) ), 'LineWidth', 3 );
            j = j +1;
        end
        
        plot(dataP.Xdata(1,sjSeq), dataP.Xdata(2,sjSeq), 'k');
        while bounds(i) < dataP.aggTs(sj+1)
            ftSeq = (bounds(i)+1):bounds(i+1);
            if useCollapse
                plot(dataP.Xdata(1,ftSeq), dataP.Xdata(2,ftSeq), 'Color', colsB(allMatrixC(i,1)));
            else
                plot(dataP.Xdata(1,ftSeq), dataP.Xdata(2,ftSeq), 'Color', colsB(allMatrix(i,1)));
            end
            i = i +1;
        end
        
    else
        figure
        
        sjSeq = (1+dataP.aggTs(sj)):dataP.aggTs(sj+1);
        x = dataP.Xdata(1,:);
        y = dataP.Xdata(2,:);
        x(x < -15) = NaN;
        y(y < -6) = NaN;
        
        subplot(1,2,1)
        plot(x(sjSeq), y(sjSeq), 'k');
        hold on
        while bounds(i) < dataP.aggTs(sj+1)
            ftSeq = (bounds(i)+1):bounds(i+1);
            if useCollapse
                plot(x(ftSeq), y(ftSeq), 'Color', colsB(allMatrixC(i,1)));
            else
                plot(x(ftSeq), y(ftSeq), 'Color', colsB(allMatrix(i,1)));
            end
            i = i +1;
        end
        axis([-15 15 -10 10]);
        
        
        boundsA = [sjSeq(1)-1 sjSeq( diff(dataP.Xdata(3,sjSeq)) ~=0) sjSeq(end)];
        subplot(1,2,2)
        plot(x(sjSeq), y(sjSeq), 'k');
        hold on
        j = 1;
        while boundsA(j) < dataP.aggTs(sj+1)
            ftSeq = (boundsA(j)+1):boundsA(j+1);
            plot(dataP.Xdata(1,ftSeq), dataP.Xdata(2,ftSeq), 'Color', colsA(dataP.Xdata(3,ftSeq(1)+1) ) );
            j = j +1;
        end
        axis([-15 15 -10 10]);
    end
    
end

%% Andersson Data Vis
figure
for k = [1:4,6]
    subplot(2,3,k)
    histogram(standardMaxDists{k}, 100, 'Normalization', 'probability')
    %axis([0 50 0 0.22])
end
figure
for k = [1:4,6]
    subplot(2,3,k)
    histogram(standardDists{k}, 100, 'Normalization', 'probability')
    %axis([0 50 0 0.22])
end


figure
for k = [1:4,6]
    subplot(2,3,k)
    n = hist3([standardMaxSpeed{k}; standardDists{k}]', ...
        [50 50]);
    n1 = n';
    n1(size(n,1) + 1, size(n,2) + 1) = 0;
    xn = linspace(min(standardMaxSpeed{k}(~standardOutlier{k})), ...
                  max(standardMaxSpeed{k}(~standardOutlier{k})),size(n,1)+1);
    yn = linspace(min(standardMaxDists{k}(~standardOutlier{k})), ...
                  max(standardMaxDists{k}(~standardOutlier{k})),size(n,2)+1);
    pcolor(xn, yn, n1)
    xlabel('Max Speed'); ylabel('Dispersion');
    %axis([1 20 1 20])
end


figure
for k = [1:4,6]
    subplot(2,3,k)
    n = hist3([standardPCA{k}; standardDists{k}]', ...
        [50 50]);
    n1 = n';
    n1(size(n,1) + 1, size(n,2) + 1) = 0;
    xn = linspace(min(standardPCA{k}(~standardOutlier{k})), ...
                  max(standardPCA{k}(~standardOutlier{k})),size(n,1)+1);
    yn = linspace(min(standardDists{k}(~standardOutlier{k})), ...
                  max(standardDists{k}(~standardOutlier{k})),size(n,2)+1);
    pcolor(xn, yn, n1)
    xlabel('PCA Ratio'); ylabel('Total Distance');
    %axis([1 20 1 20])
end

xn = linspace(0, 150);
yn = linspace(0, 14);
na = hist3([standardMaxSpeed{1}; standardDists{1}]', {xn, yn});
nb = hist3([standardMaxSpeed{4}; standardDists{4}]', {xn, yn});
n1 = (na + nb)';
n1(size(na,1) + 1, size(na,2) + 1) = 0;

n1a = na';
n1a(size(na,1) + 1, size(na,2) + 1) = 0;
n1b = nb';
n1b(size(na,1) + 1, size(na,2) + 1) = 0;

xb = linspace(0,150,101);
yb = linspace(0,14, 101);

figure
h = pcolor(xb, yb, n1);

figure
subplot(1,2,1)
pcolor(xb, yb, n1a);
subplot(1,2,2)
pcolor(xb, yb, n1b);




