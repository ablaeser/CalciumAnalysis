% For each ROI, get the tuning angle from the GLM. Gather stretch and
% activity around this angle (+/- 30?) and fit the stimulus-response curve
% to a sigmoid
dTheta = pi/6;
modelSig = @(b,x)( b(1)./(1+exp(-(x-b(3))/b(2))) + b(4) ); % beta = [factor, growth rate, threshold, constant response offset]
sigmFit = cell(Nclick(x),1); 
fitResult = repmat( struct('R2',NaN, 'thresh',NaN, 'pThresh',NaN, 'rate',NaN, 'pRate',NaN, 'factor',NaN, 'pFactor',NaN, 'offset',NaN, 'pOffset',NaN ), 1, Nclick(x) );

close all
figure('Units','normalized','OuterPosition',[0,0,1,1],'color','w','PaperOrientation','landscape');
for x = 1
    for r = find( ~isnan([DAresult{x}.angle]) )
                                           
        
        subplot(1,3,1); cla;
        polarplot( thetaData, magData, 'k.' ); hold on;
        polarplot( thetaData(refInd), magData(refInd), 'b.')
        polarplot( refAngle*[1,1], [0,10], 'k', 'LineWidth',2 )     
        polarplot( thetaData(oppInd), magData(oppInd), 'r.')
        polarplot( oppAngle*[1,1], [0,10], 'k-.', 'LineWidth',2 )
        
        % Sigmoidal fitting
        fitInd = refInd; % [refInd; oppInd];
        subplot(1,3,2); cla;
        % {
        % Fit this data directly to a sigmoid (unbinned)
        % Determine initial parameters for sigmoidal fit
        betaInit(4) = 0; % Constant offset response
        betaInit(1) = max(Rda{x}(fitInd,r)); % Saturation response 
        betaInit(2) = DAresult{x}(r).mag; % Response growth rate
        betaInit(3) = median(magData(fitInd)); % Threshold stimulus
        [ sigmFit{r}, fitResult(r) ] = FitSigmoid( magData(fitInd), Rda{x}(fitInd,r), betaInit );
        %sigmFit{r} = fitnlm( magData(fitInd), Rda{x}(fitInd,r), modelSig, betaInit);  
        plot( magData(fitInd), Rda{x}(fitInd,r), 'b.' ); hold on; axis square;
        xRange = linspace(0,max(magData(fitInd)))';  %(0:0.01:groupStimMean(end))';
        plot( xRange, predict(sigmFit{r}, xRange), 'k' ); 
        title( sprintf('R^2adj = %2.3f, Threshold = %2.5f (p=%1.3f), Rate = %2.5f (p=%1.3f)', fitResult(r).R2, fitResult(r).thresh, fitResult(r).pThresh, fitResult(r).rate, fitResult(r).pRate ) );
        %}        
        ylim([0,0.05]);
        % {
        subplot(1,3,3);
        [ groupStimMean, groupResponseMean, groupResponseStd, groupSize, fitResult, fitSumm, sigmFit{r} ] = StimResponse( magData(fitInd), Rda{x}(fitInd,r), 0, 10, 0, 'show',false, 'fit',false );
        errorbar( groupStimMean, groupResponseMean, groupResponseStd./sqrt(groupSize), 'LineStyle','none' ); hold on; axis square;
        xRange = linspace(0,groupStimMean(end))';  %(0:0.01:groupStimMean(end))';
        plot( xRange, predict(sigmFit{r}{1}, xRange), 'k' ); 
        ylim([0,0.05]);
        %}
        title( sprintf('R^2adj = %2.3f, Threshold = %2.5f (p=%1.3f), Rate = %2.5f (p=%1.3f)', fitResult(r).R2, fitResult(r).thresh, fitResult(r).pThresh, fitResult(r).rate, fitResult(r).pRate ) );
        pause;
    end
end
