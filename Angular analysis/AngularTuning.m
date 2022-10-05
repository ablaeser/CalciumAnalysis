function [ dirBinResponse, orBinResponse, resultant ] = AngularTuning( thetaBin, thetaDirTest, thetaOrTest, testActivity, Nboot, show ) % , CI

% Bin the test activity
[ dirBinResponse, orBinResponse ] = AngularResponse( thetaDirTest, thetaOrTest, testActivity, thetaBin, false); % , true 
% Suppress negative response means 
negDirVal = dirBinResponse.mean < 0;
negOrVal = orBinResponse.mean < 0;
if any(negDirVal)
    %fprintf('\nSupressing %i negative directional response means', sum(negDirVal));
    dirBinResponse.mean(negDirVal) = 0;
    dirBinResponse.std(negDirVal) = 0;
    dirBinResponse.sem(negDirVal) = 0;
end
if any(negOrVal)
    %fprintf('\nSupressing %i negative orientational response means', sum(negOrVal));
    orBinResponse.mean(negOrVal) = 0;
end
%fprintf('\n');

resultant = struct('dir',[], 'or',[]);
% Directional tuning 
dirResult = AngularResultant( thetaBin.dir.cent(1:thetaBin.dir.Nbin)', dirBinResponse.mean', false ); % cent has a redundant entry for plotting purposes
resultant.dir = struct('norm',dirResult(1), 'normCI',[NaN,NaN], 'angle',dirResult(2), 'angleCI',[NaN,NaN], 'mag',dirResult(3), 'magCI',[NaN,NaN], 'width',NaN);
% Orientational tuning 
orResult = AngularResultant( thetaBin.or.cent, orBinResponse.mean', false );
resultant.or = struct('norm',orResult(1), 'normCI',[NaN,NaN], 'angle',orResult(2), 'angleCI',[NaN,NaN], 'mag',orResult(3), 'magCI',[NaN,NaN], 'width',NaN);

% Bootstrap confidence 95% intervals (orientation commented out to save time since we don't use it)
if Nboot > 0
    Ntest = size( testActivity, 1 );
    tic
    dirBootResult = nan(Nboot,3); %orBootResult = nan(Nboot,3);
    parfor n = 1:Nboot
        bootInd = randsample( 1:Ntest, Ntest, true );
        [ bootDirResponse,  ~] = AngularResponse( thetaDirTest(bootInd), thetaOrTest(bootInd), testActivity(bootInd), thetaBin, false ); % , false  bootOrResponse
        negDirVal = bootDirResponse.mean < 0;
        bootDirResponse.mean(negDirVal) = 0;
        dirBootResult(n,:) = AngularResultant( thetaBin.dir.cent(1:thetaBin.dir.Nbin)', bootDirResponse.mean', false );
        %negOrVal = bootOrResponse.mean < 0;
        %bootOrResponse.mean(negOrVal) = 0;   
        %orBootResult(n,:) = AngularResultant( thetaBin.or.cent, bootOrResponse.mean', false );
    end
    toc
    lowerDirLim = prctile( dirBootResult, 5, 1 ); 
    upperDirLim = prctile( dirBootResult, 95, 1 );
    resultant.dir.normCI = [lowerDirLim(1), upperDirLim(1)]; % CI.dir.norm = [lowerDirLim(1), upperDirLim(1)];
    resultant.dir.magCI = [lowerDirLim(3), upperDirLim(3)];
    % Need to be account for wrap around effect for percentiles of radians
    binWidth = 30;
    angleBinEdge = 0:binWidth:360;
    angleBinCent = angleBinEdge(2:end)-binWidth/2;
    Nbin = numel(angleBinCent);
    bootAngle = Rad2PosDeg(dirBootResult(:,2));
    % find the bin opposite to the most populated bin 
    angleBinCounts = histcounts(bootAngle, angleBinEdge);
    [~, maxCountsBin] = max(angleBinCounts);  % maxCounts
    maxBinCent = angleBinCent(maxCountsBin);
    if maxBinCent <= 180
        minBootAngle = maxBinCent + 180;
    else
        minBootAngle = maxBinCent - 180;
    end
    subAngle = bootAngle - minBootAngle + 360;
    subAngle(subAngle > 360) = subAngle(subAngle > 360)-360;
    %subAngle(subAngle >= 0 & subAngle < 90) = subAngle(subAngle >= 0 & subAngle < 90) + 360; % add 360 to first quadrant to avoid the seam where 0 and 360 are the same
    %bootAngle(bootAngle >= 0 & bootAngle < 90) = bootAngle(bootAngle >= 0 & bootAngle < 90) + 360; % add 360 to first quadrant to avoid the seam where 0 and 360 are the same
    %{
    alphaVal = 0.4;
    figure;
    subplot(1,2,1);
    for n = 1:Nboot
        polarplot(deg2rad(bootAngle(n))*[1,1], [0,1], 'color',[0,0,0,alphaVal]); hold on;
    end
    polarplot(deg2rad(minBootAngle)*[1,1], [0,1], 'r', 'LineWidth',1.5)
    title('Original data');
    subplot(1,2,2);
    for n = 1:Nboot
        polarplot(deg2rad(subAngle(n))*[1,1], [0,1], 'color',[0,0,0,alphaVal]); hold on;
    end
    title(sprintf('Rotated by %2.0f deg', minBootAngle));
    %}
    resultant.dir.angleCI = deg2rad(prctile(subAngle, [5,95])+minBootAngle);  % Rad2PosDeg(dirBootResult(:,2))
    resultant.dir.angleCIwidth = rad2deg(abs(diff(resultant.dir.angleCI)));
end

% Show the confidence interval for directional case
if show
    figure;
    for n = 1:Nboot
        polarplot(dirBootResult(n,2)*[1,1], [0,1], 'color',[0,0,0,0.2]); hold on;
    end
    polarplot(dirResult(1,2)*[1,1], [0,1], 'color','b', 'LineWidth',1.5); %hold on;
    polarplot(resultant.dir.angleCI(1)*[1,1], [0,1], 'color','r', 'LineWidth',1.5); 
    polarplot(resultant.dir.angleCI(2)*[1,1], [0,1], 'color','r', 'LineWidth',1.5); 
    title(sprintf('Conf Interval Width = %2.1f deg', resultant.dir.angleCIwidth))
    pause;
end
end