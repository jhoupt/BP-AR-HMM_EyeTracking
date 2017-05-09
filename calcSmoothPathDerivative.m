function [dx, dy, ddx, ddy, dr, dt] = calcSmoothPathDerivative(xy, blink, samppersec, sig)

% Convert to polar
rt = zeros(size(xy));
rt(1,:) = sqrt(xy(1,:).^2+ xy(2,:).^2);
rt(2,:) = atan2( xy(2,:), xy(1,:));

rt(:,blink) = NaN;
if any(any(isnan(rt)))
    rt(1,:) = fillMissing(rt(1,:));
    rt(2,:) = fillMissing(rt(2,:));
end

xy(:,blink) = repmat([-20;-20], [1, sum(blink)]);


% begin from Itti
hw = floor(sig * sqrt(-2 * log(0.001))); % half filter width
fil = -hw:1:hw; fil = exp(-fil.*fil/(2*sig*sig));
fil = fil / sum(fil);  % normalized Gaussian filter
filV = convCleanMod(fil, [-1 0 1]); % do a sobel on this
filA = convCleanMod(fil, [1 -2 1]);

dx = convClean(xy(1,:), filV);
dx(1:hw) = 0; dx(length(dx)-hw+1:length(dx))=0;
dy = convClean(xy(2, :), filV);
dy(1:hw) = 0; dy(length(dy)-hw+1:length(dy))=0;
% end from Itti

dr = convClean(rt(1,:), filV);
drm = dr;
drm(dr > pi) = [2 * pi - dr(dr>pi)];
drm(dr < -pi) = [2 * pi + dr(dr<-pi)];
dr = drm;
dr(1:hw) = 0; dr(length(dr)-hw+1:length(dr))=0;


dt = convClean(rt(2, :), filV);
dt(1:hw) = 0; dt(length(dt)-hw+1:length(dt))=0;


ddx = convClean(xy(1,:), filA);
ddx(1:hw) = 0; ddx(length(ddx)-hw+1:length(ddx))=0;
ddy = convClean(xy(2, :), filA);
ddy(1:hw) = 0; ddy(length(ddy)-hw+1:length(ddy))=0;

dx(blink) = 0;
dy(blink) = 0;

dr(blink) = 0;
dt(blink) = 0;

dx = dx * samppersec;
dy = dy * samppersec;
ddx = ddx * samppersec^2;
ddy = ddy * samppersec^2;

dr = dr * samppersec;
dt = dt * samppersec;
end