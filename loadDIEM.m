function [dataV, dataS, dataP, dataVA, dataPOLAR] = loadDIEM()

basedir =  '/Users/jhoupt/Documents/Research/ParseGaze/Data/';
dirlist = {'50_people_brooklyn_1280x720' '50_people_london_1280x720' 'game_trailer_lego_indiana_jones_1280x720' 'game_trailer_wrath_lich_king_shortened_subtitles_1280x548' 'nigella_chocolate_pears_1280x712' 'pingpong_long_shot_960x720' 'sport_wimbledon_baltacha_1280x704' 'sport_wimbledon_federer_final_1280x704' 'sport_wimbledon_murray_1280x704' 'tv_the_simpsons_860x528'};


screenRes = [1280 960];
screenDim = .508 * [1280 960] / sqrt(960^2 + 1280^2);
viewDist = 0.9;
pxperdeg = screenRes(1) / ( 2*360 * atan((screenDim(1)/2) / viewDist) / (2*pi) );
center = screenRes ./ 2;
samppersec = 30;

dataV = SeqData();
dataS = SeqData();
dataP = SeqData();
dataVA = SeqData();
%dataL = SeqData();
dataPOLAR = SeqData();

fncheck = '^diem\d+s\d+\.asc';
sj = '(\d)s(\d\d)';
for k=1:length(dirlist)
    dirname = char(strcat(basedir, dirlist(k), '/event_data/'));
    allfiles = dir(dirname);
    
    %for i = 1:(length(allfiles)-109)
    l = 1;
    m = 1;
    while l <= 10
        fname = allfiles(m).name;
        sname = regexp(fname, fncheck, 'match');
        if (size(sname) > 0)
            tag = fname;%regexp(allfiles(l).name, sj, 'match');
            %X = csvread(strcat(dirname, fname), 1, 1)';
            % File Columns:  Frame; Left X; Left Y; Left Pupil Dilation; Left Event; Right X; Right Y; Right Dilation; Right Event
            
            X = load(strcat(dirname, fname));
            xy = X(:,2:3)';
            % Convert to degrees with 0,0 at center
            xy = (xy - repmat(center', [1 size(xy,2)]))./pxperdeg;
            
            eyelinkParse = X(:,5);
            blink = (eyelinkParse <= 0)'; %(X(:,2) <= 0 | X(:,3) <= 0)';
            
            % VERIFY !!! %
            blink(xy(2,:)< -9.54) = 1;
            blink(xy(2,:)> 9.54) = 1;
            blink(xy(1,:)< -12.72) = 1;
            blink(xy(1,:)> 12.72) = 1;
            
            eyelinkParse(blink) = 0;
            % Start counting at 1 rather than 0
            eyelinkParse = eyelinkParse + 1;
            
            sig =1.;
            [dx, dy, ddx, ddy, dr, dt] = calcSmoothPathDerivative(xy, blink, samppersec, sig);
            
            dataS = dataS.addSeq([sqrt(dx.^2 + dy.^2); blink], tag);
            dataV = dataV.addSeq([dx; dy; blink], tag);
            dataVA = dataVA.addSeq([sqrt(dx.^2 + dy.^2); sqrt(ddx.^2 + ddy.^2); blink], tag);
            
            dataPOLAR = dataPOLAR.addSeq([abs(dr); abs(dt); blink], tag);
            
            dataP = dataP.addSeq([xy; eyelinkParse'], tag);
            %dataL = dataL.addSeq([xy; blink], tag);
            
            l = l+1;
        end
        m = m+1;
    end
end