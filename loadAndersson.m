function [dataV, dataS, dataP, dataVA, dataPOLAR] = loadAndersson()
basedir =  '/Users/jhoupt/Documents/git/EyeMovementDetectorEvaluation/annotated_data/';
dirlist = {'dots' 'images' 'videos'};

dataV = SeqData();
dataS = SeqData();
dataP = SeqData();
dataVA = SeqData();
dataPOLAR = SeqData();

for d = 1:size(dirlist,2)
    matFilelist = dir(strcat(basedir, char(dirlist(d)), '/*RA.mat'));
    expdir = strcat(basedir, char(dirlist(d)), '/');
    for f = 1:size(matFilelist,1)
        load(strcat(expdir, matFilelist(f).name));
        
        if ETdata.sampFreq ~= 500
            matFilelist(f).name
        end
        
        
        pxperdeg = ETdata.screenRes(1) / ( 2*360 * atan((ETdata.screenDim(1)/2) / ETdata.viewDist) / (2*pi) );
        center = ETdata.screenRes ./ 2;
        samppersec = ETdata.sampFreq;
        
        tag = matFilelist(f).name(1:(size(matFilelist(f).name,2)-4));
        
        
        markEyeParse = runMarkEye(ETdata);
        
        blink = (ETdata.pos(:,6)==5)';
        if sum(ETdata.pos(:,1) < 0)
            matFilelist(f).name
        end
        
        xy = ETdata.pos(:,4:5)';
        % Convert to degrees with 0,0 at center
        xy = (xy - repmat(center', [1 size(xy,2)]))./pxperdeg;
            
        blink(xy(2,:)< -11.8) = 1;
        blink(xy(2,:)> 11.8) = 1;
        blink(xy(1,:)< -15.8) = 1;
        blink(xy(1,:)> 15.8) = 1;
        
        ETdata.pos(blink,6)=5;
        sig = 1.5;
        [dx, dy, ddx, ddy, dr, dt] = calcSmoothPathDerivative(xy, blink, samppersec, sig);
        
        
        dataS = dataS.addSeq([sqrt(dx.^2 + dy.^2); blink], tag);
        dataV = dataV.addSeq([dx; dy; blink], tag);
        dataVA = dataVA.addSeq([sqrt(dx.^2 + dy.^2); sqrt(ddx.^2 + ddy.^2); blink], tag);
        dataPOLAR = dataPOLAR.addSeq([abs(dr); abs(dt); blink], tag);    
        dataP = dataP.addSeq([xy; ETdata.pos(:,6)'; markEyeParse], tag);
    end
end
end
