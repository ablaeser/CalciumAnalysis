Ntop = 5;
minDur = 10;
minIso = 10;
fitWindow = [10,10];

waveType = 'spont';
modelStep = @(b,x)( b(1)./(1+exp(-(x-b(3))/b(2))) + b(4)./(1+exp((x-b(6))/b(5))) + b(7)); %  
betaInit(1) = 1; % Saturation response
betaInit(2) = 0.1; % onset time constant
betaInit(3) = 0; % Onset time
betaInit(4) = 1; % offset saturations
betaInit(5) = 0.2; % offset time constant
betaInit(6) = 3; % Offset time
betaInit(7) = 0; % Constant term
fitOpts = statset('Display','off', 'MaxIter',1000);
for x = 37 %xWave
    zCat = [fluor{x}(expt(x).preRuns).z];
    zCat = vertcat(zCat.ROI);
    Tcat = vertcat(Tscan{x}{expt(x).preRuns});
    for f = 1:Nfiber(x)
        %plot( fluor{x}(1).z.ROI(:,fiber{x}(f).ROI) );
        
        % Find the Ntop brightest spontaneous (non-loco) events meeting duration/isolation criteria
        %fiberSAFE = {SAFE{x}{fiber{x}(f).ROI}}; %{fluorEvent{x}{run}{fiber{x}(f).ROI}};
        %fiberEventDur = [fiberEvent]
        zTop = nan(Ntop, fiber{x}(f).Nroi);
        topSAFE = repmat( struct('runScan',NaN, 'T',NaN, 'Tstart',NaN, 'Tstop',NaN, 'Tpeak',NaN, 'dur',NaN, 'zPeak',NaN, 'dFFpeak',NaN, 'iso',[]), Ntop, fiber{x}(f).Nroi);
        for r = 1:fiber{x}(f).Nroi %fiber{x}(f).ROI
            tempSAFE = SAFE{x}{fiber{x}(f).ROI(r)};
            tempSAFEiso = vertcat(tempSAFE.iso);
            tempSAFEiso = tempSAFEiso(:,1)';
            tempSAFE = tempSAFE(tempSAFEiso >= minIso & [tempSAFE.dur] > minDur);
            [zSafe,zSort] = sort([tempSAFE.zPeak], 'descend');
            Ntemp = min(numel(tempSAFE), Ntop);
            topSAFE(1:Ntemp,r) = tempSAFE(zSort(1:Ntemp))';
        end
        [~,rTrigger] = max( median(reshape([topSAFE.zPeak], Ntop, fiber{x}(f).Nroi), 1, 'omitnan') ); % use the overall brightest ROI for trigger
        topDur = reshape([topSAFE.dur], Ntop, fiber{x}(f).Nroi);
        Ttrigger = [topSAFE(:,rTrigger).Tstart];
        Ntrigger = sum(~isnan(Ttrigger));
        for n = 1:Ntrigger
            tempScans = find(Tcat >= Ttrigger(n)-fitWindow(1) & Tcat < Ttrigger(n)+fitWindow(2)+topDur(n,rTrigger) );
            Ttemp = Tcat(tempScans);
            Ttemp = Ttemp - Ttrigger(n);
            for roi = fiber{x}(f).ROI
                tempTrace = zCat(tempScans,roi);
                %Nnan = sum(isnan(tempTrace));
                tempTrace = (tempTrace - median( tempTrace(Ttemp < 0), 'omitnan')); % ,
                tempTrace = tempTrace/max(tempTrace, [], 1, 'omitnan');
                %fiber(f).(waveType).trace{k,roi} = tempTrace;
                plot(Ttemp, tempTrace ); hold on; % zCat(tempScans,roi)
                %{
                % Fit peri-onset signal to step function
                betaInit(6) = topDur(n,rTrigger); % Offset time    
                fitResult = struct('R2',NaN, 'Tonset',NaN, 'pTonset',NaN, 'tauOn',NaN, 'pTauOn',NaN, 'Toffset',NaN, 'pToffset',NaN, 'tauOff',NaN, 'pTauOff',NaN);
                try
                    tempFit = fitnlm( Ttemp, tempTrace, modelStep, betaInit, 'Options',fitOpts ); % modelSigmoid
                    %fiber(f).(waveType).traceFit{k,roi} = fitnlm( Twindow, tempTrace, modelStep, betaInit, 'Options',fitOpts ); % modelSigmoid
                    %[x, R2] = lsqcurvefit(modelStep, betaInit, Twindow, tempTrace, betaLower, betaUpper)
                    %{
                    fitResult.R2 = fiber(f).(waveType).traceFit{k,roi}.Rsquared.Adjusted; % goodness of fit
                    fitResult.Tonset = fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(3);
                    fitResult.pTonset = fiber(f).(waveType).traceFit{k,roi}.Coefficients.pValue(3);
                    fitResult.tauOn = fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(2);
                    fitResult.pTauOn = fiber(f).(waveType).traceFit{k,roi}.Coefficients.pValue(2);
                    fitResult.Toffset = fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(6);
                    fitResult.pToffset = fiber(f).(waveType).traceFit{k,roi}.Coefficients.pValue(6);
                    fitResult.tauOff = fiber(f).(waveType).traceFit{k,roi}.Coefficients.Estimate(5);
                    fitResult.pTauOff = fiber(f).(waveType).traceFit{k,roi}.Coefficients.pValue(5);
                    %}
                    plot(Ttemp, predict(tempFit, Ttemp), 'k--' ); %hold on;
                    pause
                catch
                    fprintf('\nSkipped fitting - no data available');
                end
                %}
            end
            
            title(sprintf('[x, f, n] = [%i, %i, %i]', x, f, n));
            pause; cla;
        end
    end
end

%%

fitOpts = statset('Display','off', 'MaxIter',1000, 'Robust',[], 'RobustWgtFun','bisquare');
tempFit = fitnlm( Twindow, tempTrace, modelStep, betaInit, 'Options',fitOpts ); % fiber(f).(waveType).traceFit{k,roi}

betaTest = betaInit; betaTest(5) = 1;
testData = modelStep(betaInit, Twindow);
betaLower = [0, 0, -2, 0, 0, 0, -1]; %lower bounds on [Saturation response, onset time constant, Onset time, offset saturations, offset time constant, Offset time, Constant term]
betaUpper = [2, 10, 2, 2, 10, 100, 1];
[lsfit, R2] = lsqcurvefit(modelStep, betaInit, Twindow, tempTrace, betaLower, betaUpper ); % , resid, exitFlag

figure;
plot(Twindow, tempTrace); hold on;
plot(Twindow, modelStep(lsfit, Twindow) );
