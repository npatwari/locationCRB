sigmaOverN   = 1.7;   % unitless
sigma_d      = 6.1e-9 * 3e8; % meters (from measurements in journal article)
sideLength   = 20;    % meters
gridSizeList = 4:15;  % # sensors along each side
radioRange   = 10;    % meters

%| Set up the coordinates of the reference devices
xRef = sideLength.*[0 0 1 1];
yRef = sideLength.*[0 1 1 0];

% Set up the coordinates of the blind devices
for i = 1:length(gridSizeList)
    disp(i);
    gridSize = gridSizeList(i);
    delta    = 1/(gridSize-1);
    coords   = sideLength.*(0:delta:1);
    xMatrix  = ones(gridSize,1)*coords;
    yMatrix  = xMatrix';
    total    = gridSize^2;
    xBlind   = [xMatrix(2:gridSize-1), ...
            xMatrix(gridSize+1:total-gridSize), ...
            xMatrix(total-gridSize+2:total-1)];
    yBlind   = [yMatrix(2:gridSize-1), ...
            yMatrix(gridSize+1:total-gridSize), ...
            yMatrix(total-gridSize+2:total-1)];
    blinds   = length(xBlind);
    
    % Put the coordinates together
    x = [xBlind, xRef];
    y = [yBlind, yRef];
    
    %| Call the CRB calculation routine.
    measMade = ones(total);
    [std_RSS_Inf]   = calcCRBGivenLocsAny('R', x, y, blinds, total, sigmaOverN);
    RMS_RSS_Inf(i)  = sqrt(sum(std_RSS_Inf.^2)/blinds);
    [std_TOA_Inf]   = calcCRBGivenLocsAny('T', x, y, blinds, total, sigma_d);    
    RMS_TOA_Inf(i)  = sqrt(sum(std_TOA_Inf.^2)/blinds);
    [std_RSS_range] = calcCRBGivenLocsAny('R', x, y, blinds, total, sigmaOverN, radioRange);
    RMS_RSS_range(i)  = sqrt(sum(std_RSS_range.^2)/blinds);
    [std_TOA_range] = calcCRBGivenLocsAny('T', x, y, blinds, total, sigma_d, radioRange);    
    RMS_TOA_range(i)  = sqrt(sum(std_TOA_range.^2)/blinds);
end

clf;
h1 = plot(gridSizeList, RMS_RSS_Inf, 'b-o', gridSizeList, RMS_RSS_range, 'b--o', ...
    gridSizeList, RMS_TOA_Inf, 'r-v', gridSizeList, RMS_TOA_range, 'r--v');
set(h1,'MarkerSize',10)
set(h1,'LineWidth',2)
set(gca,'FontSize',16)
legend('RSS', sprintf('RSS r=%dm',radioRange),...
    'TOA', sprintf('TOA r=%dm',radioRange))
set(gca,'FontSize',20)
set(h1(1),'MarkerFaceColor','b')
set(h1(2),'MarkerFaceColor','b')
set(h1(3),'MarkerFaceColor','r')
set(h1(4),'MarkerFaceColor','r')
xlabel('Number of Devices per Side')
ylabel('Lower Bound on RMS Localization Error (m)')
set(gca,'ylim',[0 3])
set(gca,'xlim',[4 15])
grid;

%print -depsc plotCRB_Eg_Indoor.eps