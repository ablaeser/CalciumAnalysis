%% Fit test range responses to cosine tuning curves
thetaRange = linspace(-pi, pi, 1000)'; 
tuningModel = @(b,x)( b(1) + b(2)*cos(b(4)*x - b(3)) ); % beta = [offset, tuning sharpness, peak angle, frequency ]
figure('WindowState','maximized');
tic;
cosFitResults = cell(1,Nexpt);
for x = 6
    for r = flip(1:Nclick(x))
        % Fit to a cosine tuning curve
        cosFit = fitnlm( testStruct(x).dir, testStruct(x).activity(:,r), tuningModel, [0,1,0,1] );
        cosFitResults{x}(r) = struct( 'angle',cosFit.Coefficients.Estimate(3), 'offset',cosFit.Coefficients.Estimate(1), 'factor',cosFit.Coefficients.Estimate(2), 'frequency',cosFit.Coefficients.Estimate(4),...
            'Pangle',cosFit.Coefficients.pValue(3), 'Poffset',cosFit.Coefficients.pValue(1), 'Pfactor',cosFit.Coefficients.pValue(2), 'Pfrequency',cosFit.Coefficients.pValue(4), ...
            'centPred',predict(cosFit, thetaBin.dir.cent) );
        toc;   
        % {
        cosAct = predict( cosFit, thetaDirStretch{x} );
        sp(1) = subplot(1,3,1);
        plot( testStruct(x).dir, testStruct(x).activity(:,r), '.' ); hold on;
        plot( thetaRange, predict(cosFit, thetaRange), 'k' );
        plot( thetaBin.dir.cent, dirBinResponse{x}(r).mean ); hold off;
        title(sprintf('angle = %2.1f, offset = %2.1f, factor = %2.1f, frequency = %2.1f', cosFitResults{x}(r).angle, cosFitResults{x}(r).offset, cosFitResults{x}(r).factor, cosFitResults{x}(r).frequency ));
        axis square;

        sp(2) = subplot(1,3,2);
        plot( thetaDirStretch{x}, merged(x).deconvolved(:,r), 'b.' ); hold on;
        plot( thetaDirStretch{x}, cosAct, 'k.' ); hold off;

        title('All Data');
        axis square;
        linkaxes(sp,'xy');
        xlim([-pi,pi]); ylim([0,0.06]);

        subplot(1,3,3);
        polarplot( thetaBin.dir.cent, dirBinResponse{x}(r).mean, 'b' ); hold on;
        polarplot( thetaBin.dir.cent, cosFitResults{x}(r).centPred, 'k' ); hold off;
        title('Test Range Data');
        legend('Data','Fit');
        pause;
        %}
    end
end





%% Compare model fit parameters to DSI results
figure;
for x = 6
    tempResults = [testRes{x}.dir]; %sort([tempResults.norm], 'asc
    tempFitParams = [cosFitResults{x}];
    absAngleDiff = abs([tempResults.angle] - [tempFitParams.angle]);
    for r = 1:Nclick(x) %
        plot( testRes{x}(r).dir.norm, absAngleDiff, 'k.' ); hold on;
    end
end
axis square;
xlabel('DSI');
ylabel('Absolute difference in angle');

LW = 1.5; LC = [0,0,0,1];
figure('WindowState','maximized');
for r = 1:Nclick(x)
    subtightplot(1,2,1,opt{:});   
    polarplot( testRes{x}(r).dir.angle*[1,1], testRes{x}(r).dir.norm*[0,1], 'color',LC, 'LineWidth',LW );  hold on;
    
    subtightplot(1,2,2,opt{:});  
    polarplot( cosFitResults{x}(r).angle*[1,1], cosFitResults{x}(r).factor*[0,1], 'color',LC, 'LineWidth',LW );  hold on;
    %pause;
end
subtightplot(1,2,1,opt{:});  title('Directional Selectivity Index');
subtightplot(1,2,2,opt{:});  title('Cosine Fit');

%%
tuningModel = @(b,x)( b(1) + b(2)*cos(b(4)*x - b(3)) ); % beta = [offset, tuning sharpness, peak angle, frequency ]
figure('WindowState','maximized');
Ncycle = 20;
thetaRange = linspace(-pi, pi, 1000)'; 

cosFitResults = cell(1,Nexpt); thetaMSEcheck = cell(1,Nexpt); 
for x = 6
    NmergedFrames = size(merged(x).exp,1);
    % Define training data
    trainWave = square( (2*pi*Ncycle)*(0:NmergedFrames-1)/(NmergedFrames-1), 0 );
    plot( trainWave ); xlabel('Frame'); set(gca,'Ytick',[0,1], 'box','off');
    ylim([-1.01, 1.01] );
    trainFrames = find( trainWave > 0 );
    checkFrames = find( trainWave < 0 );
    trainActivity = merged(x).deconvolved(trainFrames, :);
    trainDir = thetaDirStretch{x}(trainFrames);
    checkActivity = merged(x).deconvolved(checkFrames, :);
    checkDir = thetaDirStretch{x}(checkFrames);
    Ntrain = numel(trainFrames); Ncheck = numel(checkFrames);
    thetaMSEcheck{x} = nan(Ncheck, Nclick(x));
    thetaVectorMethod = nan(Ncheck, Nclick(x));
    tuneCompMat = nan(1, Nclick(x), 2);
    tic;
    for r = flip(1:Nclick(x))
        % Fit to a cosine tuning curve
        cosFit = fitnlm( testStruct(x).dir, testStruct(x).activity(:,r), tuningModel, [0,1,0,1] );
        cosFitResults{x}(r) = struct( 'angle',cosFit.Coefficients.Estimate(3), 'offset',cosFit.Coefficients.Estimate(1), 'factor',cosFit.Coefficients.Estimate(2), 'frequency',cosFit.Coefficients.Estimate(4),...
            'Pangle',cosFit.Coefficients.pValue(3), 'Poffset',cosFit.Coefficients.pValue(1), 'Pfactor',cosFit.Coefficients.pValue(2), 'Pfrequency',cosFit.Coefficients.pValue(4));
        cosAct = predict( cosFit, thetaDirStretch{x} );
        %tuneCompMat(1,r,:) = [cos(cosFitResults{x}(r).angle), sin(cosFitResults{x}(r).angle)];
        toc;   
        %{
        rangePredict = predict( cosFit, thetaRange );
        for z = 1:Ncheck
            %checkPredict = predict( cosFit, thetaRange );
            %[mse, mseMinInd] = min( (checkActivity(z,r) - rangePredict).^2, [], 1 );
            %thetaMSEcheck{x}(z,r) = thetaRange(mseMinInd);
            thetaVectorMat{x}(z,r) = 
        end
        %}
        %{
        sp(1) = subplot(2,2,1);
        plot( testStruct(x).dir, testStruct(x).activity(:,r), '.' ); hold on;
        plot( thetaRange, predict(cosFit, thetaRange), 'k' );
        plot( thetaBin.dir.cent, dirBinResponse{x}(r).mean ); hold off;
        title(sprintf('angle = %2.1f, offset = %2.1f, factor = %2.1f, frequency = %2.1f', cosFitResults{x}(r).angle, cosFitResults{x}(r).offset, cosFitResults{x}(r).factor, cosFitResults{x}(r).frequency ));
        axis square;

        sp(2) = subplot(2,2,3);
        plot( thetaDirStretch{x}, merged(x).deconvolved(:,r), 'b.' ); hold on;
        plot( thetaDirStretch{x}, cosAct, 'k.' ); hold off;

        title('All Data');
        axis square;
        linkaxes(sp,'xy');
        xlim([-pi,pi]); ylim([0,0.04]);

        subplot(2,2,2);
        polarplot( thetaBin.dir.cent, dirBinResponse{x}(r).mean, 'b' ); hold on;
        polarplot( thetaBin.dir.cent, predict(cosFit, thetaBin.dir.cent), 'k' ); hold off;
        title('Test Range Data');
        legend('Data','Fit');
        
        subplot(2,2,4);
        %plot( checkDir, checkActivity(:,r), 'b.'); hold on;
        plot( checkDir, vectorMethodTheta, '.', 'color', [0,0,0,0.1] );
        %plot( checkDir, mean(thetaMSEcheck{x},2, 'omitnan'), 'k.' ); hold off;
        axis square;
        xlim([-pi,pi]); ylim([-pi,pi])

        %pause;
        %}
    end
    % Find a subset of well tuned units covering the range without too much redundancy
    
    
    tuneCompMat = repmat( tuneCompMat, Ncheck, 1, 1 );
    vectorMethodMat = squeeze( sum( tuneCompMat.*repmat(checkActivity, 1, 1, 2), 2 ) );
    vectorMethodTheta = atan( vectorMethodMat(:,2) ./ vectorMethodMat(:,1) );
end

close all; figure;

for z = 1:100
    %line( [0,vectorMethodMat(z,1)], [0,vectorMethodMat(z,2)], 'color','r' )
    line( [0,vectorMethodMat(z,1)], [0,vectorMethodMat(z,2)], 'color','r' )
    axis square;
    pause(0.2)
end
