%% Use GLM to assess contribution of different variables
locoDeform_vol_pred = cell(1,Nexpt); locoDeform_vol_resp = cell(1,Nexpt); locoDeform_vol_opts = cell(1,Nexpt); locoDeform_vol_result = cell(1,Nexpt); locoDeform_vol_summary = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
%rLocoPreFit = cell(1,Nexpt);
GLMname = 'LocoDeform_3D_60s_planes'; % 'LocoDeform_3D_60s_scale_planes'; %  'LocoDeform_3D_60s_accelMag'; %
GLMrate = 15.49/30;
for x = xPresent %x3Dcsd
    % GLMparallel options
    locoDeform_vol_opts{x}.name = sprintf('%s_%s', expt(x).name, GLMname); %strcat(expt(x).name, , '_locoDeform_volglm');
    %locoDeform_vol_opts{x}.show = true;
    locoDeform_vol_opts{x}.rShow = NaN;
    locoDeform_vol_opts{x}.figDir = ''; % figDir;
    locoDeform_vol_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    locoDeform_vol_opts{x}.standardize = true;
    locoDeform_vol_opts{x}.trainFrac = 0.75; % 1; %
    locoDeform_vol_opts{x}.Ncycle = 20;
    locoDeform_vol_opts{x}.distribution = 'gaussian'; % 'poisson'; %
    locoDeform_vol_opts{x}.CVfold = 10;
    locoDeform_vol_opts{x}.nlamda = 1000;
    locoDeform_vol_opts{x}.maxit = 5*10^5;
    locoDeform_vol_opts{x}.minDev = 0.05;
    locoDeform_vol_opts{x}.minDevFrac = 0.1;
    locoDeform_vol_opts{x}.maxP = 0.05;
    locoDeform_vol_opts{x}.Nshuff = 0;  Nshuff = locoDeform_vol_opts{x}.Nshuff;
    locoDeform_vol_opts{x}.minShuff = 15;
    locoDeform_vol_opts{x}.window = [-60,60]; % [0,0]; % [-0.5, 0.5]; %
    locoDeform_vol_opts{x}.lopo = true; %false; %
    locoDeform_vol_opts{x}.frameRate = GLMrate; %expt(x).scanRate;
    locoDeform_vol_opts{x}.binSize = expt(x).scanRate/GLMrate;
    locoDeform_vol_opts{x}.minShuffFrame = round( locoDeform_vol_opts{x}.frameRate*locoDeform_vol_opts{x}.minShuff );
    windowFrame = [ceil(locoDeform_vol_opts{x}.window(1)*locoDeform_vol_opts{x}.frameRate), floor(locoDeform_vol_opts{x}.window(2)*locoDeform_vol_opts{x}.frameRate)]; % floor(locoDeform_vol_opts{x}.window*locoDeform_vol_opts{x}.frameRate);
    locoDeform_vol_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    locoDeform_vol_opts{x}.maxShift = max( abs(windowFrame) );
    locoDeform_vol_opts{x}.Nshift = numel( locoDeform_vol_opts{x}.shiftFrame );  %Nshift = locoDeform_volOpts(x).Nshift;
    locoDeform_vol_opts{x}.lags = locoDeform_vol_opts{x}.shiftFrame/locoDeform_vol_opts{x}.frameRate;
    
    %{
    % Define predictors
    tempVelocityCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).Vdown), locoDeform_vol_opts{x}.binSize );
    tempAccelCat = BinDownMean( abs(vertcat(loco{x}(expt(x).preRuns).Adown)), locoDeform_vol_opts{x}.binSize ); %[NaN; diff(tempVelocityCat)];
    %tempSpeedCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).speedDown), locoDeform_vol_opts{x}.binSize );
    tempStateCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).stateDown), locoDeform_vol_opts{x}.binSize );
    locoDeform_vol_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); % LocoDeformPred(x).data = []; %LocoDeformPred(x).data = LocoDeformCat.T; LocoDeformPred(x).name = {'Time'};
    locoDeform_vol_pred{x}.data = [tempVelocityCat, tempAccelCat, tempStateCat]; % tempSpeedCat,,
    locoDeform_vol_pred{x}.name = {'Velocity','|Accel|','State'}; % ,'Speed',
    locoDeform_vol_pred{x}.N = size(locoDeform_vol_pred{x}.data,2);
    for p = flip(1:locoDeform_vol_pred{x}.N), locoDeform_vol_pred{x}.lopo.name{p} = ['No ',locoDeform_vol_pred{x}.name{p}]; end
    % Set up leave-one-family-out
    locoDeform_vol_pred{x}.fam.col = {1:2, 3}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9};
    locoDeform_vol_pred{x}.fam.N = numel(locoDeform_vol_pred{x}.fam.col);
    locoDeform_vol_pred{x}.fam.name = {'Kinematics', 'State'}; %{'All'};%  'Onset Time',
    
    % Define planar responses
    tempTransMag = vertcat(deform{x}(expt(x).preRuns).transMag);
    tempTransMag = BinDownMean(tempTransMag(:,segParams{x}.zProj), locoDeform_vol_opts{x}.binSize ); %
    transNames = sprintfc('|Trans| %i', segParams{x}.zProj);
    tempScaleMag = vertcat(deform{x}(expt(x).preRuns).scaleMag);
    tempScaleMag = BinDownMean(tempScaleMag(:,segParams{x}.zProj), locoDeform_vol_opts{x}.binSize ); %
    scaleNames = sprintfc('|Scale| %i', segParams{x}.zProj);
    tempShearMag = vertcat(deform{x}(expt(x).preRuns).shearMag);
    tempShearMag = BinDownMean(tempShearMag(:,segParams{x}.zProj), locoDeform_vol_opts{x}.binSize ); %
    shearNames = sprintfc('|Shear| %i', segParams{x}.zProj);
    tempShift = vertcat(deform{x}(expt(x).preRuns).shiftZ);
    tempShift = BinDownMean(tempShift(:,intersect(segParams{x}.zProj, 4:expt(x).Nplane-3)), locoDeform_vol_opts{x}.binSize ); %
    shiftNames = sprintfc('|Shift| %i', intersect(segParams{x}.zProj, 4:expt(x).Nplane-3));
    
    %{
    % Define responses (averaging across planes)
    tempTrans = vertcat(deform{x}(expt(x).preRuns).transMag);
    tempTrans = BinDownMean( mean( tempTrans(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_vol_opts{x}.binSize ); % (:,segParams{x}.zProj)
    tempTransSpd = vertcat(deform{x}(expt(x).preRuns).DtransMag);
    tempTransSpd = BinDownMean( mean( tempTransSpd(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_vol_opts{x}.binSize ); % vertcat(deform{x}(expt(x).preRuns).DtransMag)
    tempScaleMag = vertcat(deform{x}(expt(x).preRuns).scaleMag);
    tempScaleMag = BinDownMean( mean( tempScaleMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_vol_opts{x}.binSize );
    tempStretchMag = vertcat(deform{x}(expt(x).preRuns).stretchMag);
    tempStretchMag = BinDownMean( mean( tempStretchMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_vol_opts{x}.binSize );
    tempShearMag = vertcat(deform{x}(expt(x).preRuns).shearMag);
    tempShearMag = BinDownMean( mean( tempShearMag(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_vol_opts{x}.binSize );
    tempShearRate = vertcat(deform{x}(expt(x).preRuns).DshearMag);
    tempShearRate = BinDownMean( mean( tempShearRate(:,segParams{x}.zProj), 2, 'omitnan'), locoDeform_vol_opts{x}.binSize );
    tempShift = vertcat(deform{x}(expt(x).preRuns).shiftZ);
    tempShiftMean =  BinDownMean( mean( tempShift(:,4:end-3), 2, 'omitnan' ), locoDeform_vol_opts{x}.binSize );
    locoDeform_vol_resp{x}.data = [tempTrans, tempScaleMag, tempShearMag, tempTransSpd, tempStretchMag, tempShearRate, tempShiftMean]; %  , tempDshearMag tempExp, tempComp,  , tempCompressMean
    locoDeform_vol_resp{x}.data = normalize(locoDeform_vol_resp{x}.data, 1);
    locoDeform_vol_resp{x}.data(abs(locoDeform_vol_resp{x}.data) > 5) = NaN; % suppress extreme outliers
    locoDeform_vol_resp{x}.name = {'|Translation|', '|Scale|', '|Shear|', 'Trans. Speed', '|Stretch|', 'Shear Rate', 'Z Shift'}; % , 'Z Shift', 'Compression' , 'Z Compression'
    locoDeform_vol_resp{x}.N = size(locoDeform_vol_resp{x}.data,2);
    %}
    locoDeform_vol_resp{x}.data = [tempTransMag, tempScaleMag, tempShearMag, tempShift]; % [tempTrans, tempScaleMag, tempShearMag, tempTransSpd, tempStretchMag, tempShearRate, tempShiftMean]; 
    %subplot(2,1,1); imagesc(normalize(locoDeform_vol_resp{x}.data, 1)'); %colorbar; impixelinfo
    %outlierInd = find(abs(normalize(locoDeform_vol_resp{x}.data, 1)) > 5);
    %locoDeform_vol_resp{x}.data(outlierInd) = NaN; % suppress extreme outliers
    %subplot(2,1,2); imagesc(normalize(locoDeform_vol_resp{x}.data, 1)'); colorbar; impixelinfo
    locoDeform_vol_resp{x}.name = [transNames, scaleNames, shearNames, shiftNames]; %[, shiftNames]; %{'|Scale|','Z Shift'}; %{'|Translation|', '|Scale|', '|Shear|', 'Trans. Speed', '|Stretch|', 'Shear Rate', 'Z Shift'}; % , 'Z Shift', 'Compression' , 'Z Compression'
    locoDeform_vol_resp{x}.N = size(locoDeform_vol_resp{x}.data,2);
    
    % Remove scans with missing data
    nanFrame = find(any(isnan([locoDeform_vol_pred{x}.data, locoDeform_vol_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) );
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    locoDeform_vol_pred{x}.data(nanFrame,:) = []; locoDeform_vol_resp{x}.data(nanFrame,:) = [];
    %}
    % Run the GLM
    locoDeform_vol_opts{x}.load = true; % false; %
    locoDeform_vol_opts{x}.saveRoot = expt(x).dir; % sprintf('%s%s', expt(x).dir, locoDeform_vol_opts{x}.name  ); % ''; %
    [locoDeform_vol_result{x}, locoDeform_vol_summary{x}, locoDeform_vol_opts{x}, locoDeform_vol_pred{x}, locoDeform_vol_resp{x}] = GLMparallel(locoDeform_vol_pred{x}, locoDeform_vol_resp{x}, locoDeform_vol_opts{x}); % LocoDeform_result{x}, LocoDeform_summary(x), ~, LocoDeformPred(x), LocoDeformResp(x)
    %ViewGLM(LocoDeform_pred{x}, LocoDeform_resp{x}, locoDeform_vol_opts{x}, LocoDeform_result{x}, LocoDeform_summary{x}); %GLMresultFig =
    %pause;
end

%% View GLM results
for x = xPresent
    locoDeform_vol_opts{x}.rShow = 3; %1:locoDeform_vol_resp{x}.N; % 1:LocoDeform_resp{x}.N; %NaN;
    locoDeform_vol_opts{x}.xVar = 'Time';
    ViewGLM(locoDeform_vol_pred{x}, locoDeform_vol_resp{x}, locoDeform_vol_opts{x}, locoDeform_vol_result{x}, locoDeform_vol_summary{x}); %GLMresultFig =
end

%%
devSummMat = cell(1,Nexpt);
GLMresultFig = figure('WindowState','maximized', 'color','w');
for x = xPresent
    cla;
    devSummMat{x} = [locoDeform_vol_summary{x}.dev; locoDeform_vol_summary{x}.lopo.dev; locoDeform_vol_summary{x}.lofo.dev];
    % {
    imagesc( devSummMat{x}' );
    set(gca,'XtickLabel',  [{'All'}; locoDeform_vol_summary{xPresent(1)}.lopo.name(:); locoDeform_vol_summary{xPresent(1)}.lofo.name(:)], ...
        'Ytick',1:locoDeform_vol_resp{xPresent(1)}.N , 'YtickLabel',locoDeform_vol_resp{xPresent(1)}.name, 'TickDir','out', 'TickLength',[0.003,0]);
    title(sprintf('%s: GLM Fit Summary', locoDeform_vol_opts{x}.name), 'Interpreter','none');
    xtickangle(30);
    axis image;
    CB = colorbar; CB.Label.String = 'Deviance Explained';
    caxis([0,0.2]);
    impixelinfo;
    pause;
    %}
end

devSummCat = cat(3, devSummMat{:});
devSummMed = median(devSummCat, 3, 'omitnan' );

locoDeform_vol_GLMdevFig = figure('WindowState','maximized', 'color','w');
imagesc( devSummMed' );
set(gca,'XtickLabel',  [{'All'}; locoDeform_vol_summary{x}.lopo.name(:); locoDeform_vol_summary{x}.lofo.name(:)], ...
    'YtickLabel',locoDeform_vol_resp{x}.name, 'TickDir','out', 'TickLength',[0.003,0]);
title('GLM Fit Summary (All 3D Data)', 'Interpreter','none');
xtickangle(30);
axis image;
CB = colorbar; CB.Label.String = 'Median Deviance Explained';
figPath = sprintf('%s%s_DevSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
%print( locoDeform_vol_GLMdevFig, figPath, '-dtiff');
impixelinfo;

%% Show relative explanatory value from submodels
exptColor = distinguishable_colors(Npresent);
locoDeform_vol_relExpFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_vol_resp{xPresent(1)}.N;
%tiledlayout(1, Ncol);
opt = {[0.05,0.04], [0.06,0.02], [0.04,0.02]};
Ndrops = locoDeform_vol_pred{xPresent(1)}.N+locoDeform_vol_pred{xPresent(1)}.fam.N;
for col = 1:Ncol
    tempDevFracMat = nan( Npresent, Ndrops );
    k = 0;
    for x = xPresent
        k = k + 1;
        if locoDeform_vol_result{x}(col).dev > locoDeform_vol_opts{x}.minDev
            tempDevFracMat(k,:) = [locoDeform_vol_result{x}(col).lopo.devFrac, locoDeform_vol_result{x}(col).lofo.devFrac];
        end
    end
    subtightplot(1,Ncol,col,opt{:}); %nexttile
    JitterPlot( 1-tempDevFracMat, 0.45 ); hold on;
    line(0.5+[0,Ndrops], [0,0], 'color','k');
    set(gca, 'Xtick',1:Ndrops); %, 'Xticklabel'
    ylim([-0.4, 1.1]); xlim(0.5+[0,Ndrops])
    %colororder(exptColor)
    axis square;
    if col == 1, ylabel('Relative Explanatory Value'); end
    title( locoDeform_vol_resp{xPresent(1)}.name{col} );
    set(gca, 'XtickLabel', [locoDeform_vol_pred{xPresent(1)}.name, locoDeform_vol_pred{xPresent(1)}.fam.name])  % [locoDeform_vol_summary{x}.lopo.name(:); locoDeform_vol_summary{x}.lofo.name(:)]
    xtickangle(45)
end
figPath = sprintf('%s%s_ExpValSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
print( locoDeform_vol_relExpFig, figPath, '-dtiff');

%% Show time filters per response
exptColor = distinguishable_colors(Npresent);
locoDeform_vol_GLMcoeffFig = figure('WindowState','maximized', 'color','w');
Ncol = locoDeform_vol_resp{xPresent(1)}.N;   Nrow = locoDeform_vol_pred{xPresent(1)}.N;
tiledlayout(Nrow, Ncol);
k = 0;
for row = 1:Nrow
    for col = 1:Ncol
        nexttile
        tempCoeffMat = zeros( numel(locoDeform_vol_opts{xPresent(1)}.lags), Npresent );
        k = 0;
        for x = xPresent
            k = k + 1;
            if locoDeform_vol_result{x}(col).dev > locoDeform_vol_opts{x}.minDev
                tempCoeffMat(:,k) = locoDeform_vol_result{x}(col).coeff(:,row);
            end
        end
        plot( locoDeform_vol_opts{x}.lags, tempCoeffMat ); hold on;
        colororder(exptColor)
        axis square;
        if col == 1, ylabel(sprintf('%s Coefficient', locoDeform_vol_pred{xPresent(1)}.name{row} )); end
        if row == 1, title( locoDeform_vol_resp{xPresent(1)}.name{col} ); end
        if row == Nrow, xlabel('Response Lag (s)'); end
        %pause;
    end
end

figPath = sprintf('%s%s_CoeffSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
print( locoDeform_vol_GLMcoeffFig, figPath, '-dtiff');

%% Summarize results for LocoDeform_3D_60s_planes version
tempSumm = [locoDeform_vol_summary{xPresent}];
planeScaleDev = cell2padmat( {tempSumm.dev} );
%{
JitterPlot( planeScaleDev ); hold on;
xlabel('Experiment'); ylabel('Dev explained');
set(gca, 'Xtick',1:Npresent);
line([0,numel(tempSumm)+1], 0.05*[1,1], 'linestyle','--', 'color','r');
%}

transFrac = []; scaleFrac = []; shearFrac = []; shiftFrac = [];
for x = xPresent
    transFrac = [transFrac; 1-locoDeform_vol_summary{x}.lopo.devFrac(:,intersect(locoDeform_vol_summary{x}.rGood, find(contains( locoDeform_vol_resp{1}.name, 'Trans'))))'];
    scaleFrac = [scaleFrac; 1-locoDeform_vol_summary{x}.lopo.devFrac(:,intersect(locoDeform_vol_summary{x}.rGood, find(contains( locoDeform_vol_resp{1}.name, 'Scale'))))'];
    shearFrac = [shearFrac; 1-locoDeform_vol_summary{x}.lopo.devFrac(:,intersect(locoDeform_vol_summary{x}.rGood, find(contains( locoDeform_vol_resp{1}.name, 'Shear'))))'];
    shiftFrac = [shiftFrac; 1-locoDeform_vol_summary{x}.lopo.devFrac(:,intersect(locoDeform_vol_summary{x}.rGood, find(contains( locoDeform_vol_resp{1}.name, 'Shift'))))'];
end
JitterPlot( shiftFrac ); hold on;


devFrac_mean = [mean(transFrac, 1); mean(scaleFrac, 1); mean(shearFrac, 1); mean(shiftFrac, 1)];
devFrac_sem = [SEM(transFrac, 1); SEM(scaleFrac, 1); SEM(shearFrac, 1); SEM(shiftFrac, 1)];

locoDeform_vol_GLMcoeffFig = figure('WindowState','maximized', 'color','w');
ngroups = size(devFrac_mean, 1);
nbars = size(devFrac_mean, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
bar(devFrac_mean);
legend({'','',''}, 'AutoUpdate','off', 'box','off') % 'Velocity','Acceleration','State'
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars); hold on; % Calculate the width for each bar group
    errorbar(x, devFrac_mean(:,i), devFrac_sem(:,i), '.', 'color','k');
end
ylim([0,0.2]); xlim([0.5,5.5]);
set(gca,'Xtick',1:4, 'XtickLabel',[], 'Ytick',0:0.1:0.2, 'YtickLabel',[], 'TickDir','out', 'TickLength',[TL,0]);

%% Look for depth dependence in planar loco->deform GLM
quadModel = @(b,x)b(1) + b(2)*x + b(3)*x.^2; % + b(3)*(x^.2)
quadOpts = statset('nlinfit');
quadOpts.RobustWgtFun = 'bisquare';
getPlane = @(y)( str2double(y(9:end) ) );
locoDeform_vol_GLMcoeffFig = figure('WindowState','maximized', 'color','w');
beta_COM = cell(1,Nexpt); beta_AUC = cell(1,Nexpt);
respVar = ["Trans","Scale","Shear","Shift"];
minFitPlanes = 8;

for x = xPresent %(sum(planeScaleDev >= locoDeform_vol_opts{xPresent(1)}.minDev, 1) >= 10)
    negDelay = find(locoDeform_vol_opts{x}.shiftFrame < -1);
    colPlane = cellfun(getPlane, locoDeform_vol_resp{x}.name);
    goodCol = find(locoDeform_vol_summary{x}.dev >= locoDeform_vol_opts{x}.minDev); %
    beta_COM{x} = nan(3, locoDeform_vol_pred{x}.N, 4, 2); beta_AUC{x} = nan(3, locoDeform_vol_pred{x}.N, 4, 2); % coeff x predictor x resp x p-val
    for resp = 1:numel(respVar) % 4 %
        respCol = find(contains( locoDeform_vol_resp{x}.name, respVar(resp))); % 'Scale'
        goodRespCol = intersect(respCol, goodCol);
        if numel(goodRespCol) > minFitPlanes
            goodRespPlane = colPlane(goodRespCol);
            tempCoeffMat = permute(cat(3, locoDeform_vol_result{x}(goodRespCol).coeff), [3,1,2]); % goodShiftCol
            coeffMat = nan(expt(x).Nplane, locoDeform_vol_opts{x}.Nshift, locoDeform_vol_pred{x}.N);
            coeffMat(goodRespPlane,:,:) = tempCoeffMat;
            coeffMat = flip(coeffMat, 1); % flip so superficial planes are on top
            coeffMat = abs(coeffMat);
            posCoeff = coeffMat;  
            posCoeff(:,negDelay,:) = 0; % supress negative-delay coefficients
            
            Xtick = linspace(1, locoDeform_vol_opts{x}.Nshift, 3);
            XtickLabel = linspace(locoDeform_vol_opts{x}.window(1), locoDeform_vol_opts{x}.window(2), 3);
            set(gcf, 'Name',expt(x).name);
            for p = flip(1:locoDeform_vol_pred{x}.N)
                coeffCOM = sum((1:locoDeform_vol_opts{x}.Nshift).*posCoeff(:,:,p), 2)./sum(posCoeff(:,:,p), 2); % center of mass          
                subplot(2,locoDeform_vol_pred{x}.N, p)
                imagesc( coeffMat(:,:,p) ); hold on;
                plot(coeffCOM, 1:expt(x).Nplane, 'k.');
                set(gca, 'xtick',Xtick, 'xtickLabel',XtickLabel, 'Ytick',[])
                axis square;
                xlabel('Delay (s)'); ylabel(sprintf('%s Plane', respVar(resp)));
                if p == 1, set(gca, 'Ytick',[1,expt(x).Nplane], 'YtickLabel',["Superficial","Deep"]); end
                % Fit center of mass to quadratic
                beta0_COM = [median(coeffCOM, 'omitnan'),0,0]; %
                tempFit_COM = fitnlm((1:expt(x).Nplane)', coeffCOM, quadModel, beta0_COM, 'Options',quadOpts);
                pZ_COM = tempFit_COM.Coefficients.pValue;
                % If a non-constant term is significant, use those coefficients, otherwise use the median
                if any(pZ_COM(2:3) < 0.05)
                    tempBeta_COM = tempFit_COM.Coefficients.Estimate; % beta_COM{x}(:,p)
                    beta_COM{x}(:,p,resp,1) = tempBeta_COM;
                else
                    beta_COM{x}(:,p,resp,1) = beta0_COM';
                end
                beta_COM{x}(:,p,resp,2) = pZ_COM;

                title(sprintf('%s: p = [%1.3f, %1.3f, %1.3f]', locoDeform_vol_pred{x}.name{p}, pZ_COM) ) % beta = [%1.2f, %1.2f, %1.2f],  beta_COM{x},
                plot( polyval(flip(beta_COM{x}(:,p,resp,1)),(1:expt(x).Nplane)' ), 1:expt(x).Nplane, 'linestyle','--' ) % beta_COM{x}(:,p)
                
                % Fit magnitude to quadratic
                coeffAUC = sum(posCoeff(:,:,p), 2); % , 'omitnan'  coeffMat
                beta0_AUC = [median(coeffAUC, 'omitnan'),0,0]; %
                tempFit_AUC = fitnlm((1:expt(x).Nplane)', coeffAUC, quadModel, beta0_AUC, 'Options',quadOpts);
                pZ_AUC = tempFit_AUC.Coefficients.pValue;
                % If a non-constant term is significant, use those coefficients, otherwise use the median
                if any(pZ_AUC(2:3) < 0.05)
                    tempBeta_AUC = tempFit_AUC.Coefficients.Estimate; % beta_AUC{x}(:,p)
                    beta_AUC{x}(:,p,resp,1) = tempBeta_AUC;
                else
                    beta_AUC{x}(:,p,resp,1) = beta0_AUC';
                end
                beta_AUC{x}(:,p,resp,2) = pZ_AUC;

                subplot(2,locoDeform_vol_pred{x}.N, p+locoDeform_vol_pred{x}.N)
                plot(1:expt(x).Nplane, coeffAUC, '.'); hold on;
                plot(1:expt(x).Nplane, polyval(flip(beta_AUC{x}(:,p,resp,1)),(1:expt(x).Nplane)'), 'linestyle','--' )
                title(sprintf('p = [%1.3f, %1.3f, %1.3f]', pZ_AUC) ) % beta = [%1.2f, %1.2f, %1.2f],  beta_AUC{x}, 
                xlabel('Plane'); ylabel(sprintf('%s AUC', respVar(resp))); %ylabel('AUC');
                %pause;
            end
            impixelinfo;
            %pause;  clf;
        else
            fprintf('\nx = %i: %s has %i well-fit planes', x, respVar(resp), numel(goodRespCol));
        end
    end
end

beta_AUC_pool = cat( 5, beta_AUC{xPresent} );
beta_COM_pool = cat( 5, beta_COM{xPresent} );
zeroShift = find(locoDeform_vol_opts{x}.shiftFrame == 0);
transDelay = GLMrate*(squeeze(beta_COM_pool(1,3,strcmpi(respVar, 'Trans'),1,:)) - zeroShift);
scaleDelay = GLMrate*(squeeze(beta_COM_pool(1,3,strcmpi(respVar, 'Scale'),1,:)) - zeroShift);
shearDelay = GLMrate*(squeeze(beta_COM_pool(1,3,strcmpi(respVar, 'Shear'),1,:)) - zeroShift);
shiftDelay = GLMrate*(squeeze(beta_COM_pool(1,1,strcmpi(respVar, 'Shift'),1,:)) - zeroShift); % velocity is a better predictor than state for z-shift
%{
delayMat = [transDelay, scaleDelay, shearDelay, shiftDelay];
delayMean = mean(delayMat, 'omitnan');
bar( delayMean ); hold on;
errorbar(1:4, delayMean, SEM(delayMat) );
%}

% How many experiments showed siginifcant z-effect for COM / AUC
zEffect_trans_COM = nan(1, Npresent);
zEffect_trans_COM(min(squeeze(beta_COM_pool(2:3, 3, strcmpi(respVar, 'Trans'), 2,:))) < 0.05) = 1;
zEffect_trans_COM(min(squeeze(beta_COM_pool(2:3, 3, strcmpi(respVar, 'Trans'), 2,:))) >= 0.05) = 0;
fprintf('\nFound evidence of COM z-effect for state -> translation in %i of %i experiments', sum(zEffect_trans_COM == 1), sum(~isnan(zEffect_trans_COM)) );

zEffect_trans_AUC = nan(1, Npresent);
zEffect_trans_AUC(min(squeeze(beta_AUC_pool(2:3, 3, strcmpi(respVar, 'Trans'), 2,:))) < 0.05) = 1;
zEffect_trans_AUC(min(squeeze(beta_AUC_pool(2:3, 3, strcmpi(respVar, 'Trans'), 2,:))) >= 0.05) = 0;
fprintf('\nFound evidence of AUC z-effect for state -> translation in %i of %i experiments', sum(zEffect_trans_AUC == 1), sum(~isnan(zEffect_trans_AUC)) )

zEffect_scale_COM = nan(1, Npresent);
zEffect_scale_COM(min(squeeze(beta_COM_pool(2:3, 3, strcmpi(respVar, 'Scale'), 2,:))) < 0.05) = 1;
zEffect_scale_COM(min(squeeze(beta_COM_pool(2:3, 3, strcmpi(respVar, 'Scale'), 2,:))) >= 0.05) = 0;
fprintf('\nFound evidence of COM z-effect for state -> scaling in %i of %i experiments', sum(zEffect_scale_COM == 1), sum(~isnan(zEffect_scale_COM)) );

zEffect_scale_AUC = nan(1, Npresent);
zEffect_scale_AUC(min(squeeze(beta_AUC_pool(2:3, 3, strcmpi(respVar, 'Scale'), 2,:))) < 0.05) = 1;
zEffect_scale_AUC(min(squeeze(beta_AUC_pool(2:3, 3, strcmpi(respVar, 'Scale'), 2,:))) >= 0.05) = 0;
fprintf('\nFound evidence of AUC z-effect for state -> scaling in %i of %i experiments', sum(zEffect_scale_AUC == 1), sum(~isnan(zEffect_scale_AUC)) )

zEffect_shear_COM = nan(1, Npresent);
zEffect_shear_COM(min(squeeze(beta_COM_pool(2:3, 3, strcmpi(respVar, 'Shear'), 2,:))) < 0.05) = 1;
zEffect_shear_COM(min(squeeze(beta_COM_pool(2:3, 3, strcmpi(respVar, 'Shear'), 2,:))) >= 0.05) = 0;
fprintf('\nFound evidence of COM z-effect for state -> shearing in %i of %i experiments', sum(zEffect_shear_COM == 1), sum(~isnan(zEffect_shear_COM)) );

zEffect_shear_AUC = nan(1, Npresent);
zEffect_shear_AUC(min(squeeze(beta_AUC_pool(2:3, 3, strcmpi(respVar, 'Shear'), 2,:))) < 0.05) = 1;
zEffect_shear_AUC(min(squeeze(beta_AUC_pool(2:3, 3, strcmpi(respVar, 'Shear'), 2,:))) >= 0.05) = 0;
fprintf('\nFound evidence of AUC z-effect for state -> scaling in %i of %i experiments', sum(zEffect_shear_AUC == 1), sum(~isnan(zEffect_shear_AUC)) )

zEffect_shift_COM = nan(1, Npresent);
zEffect_shift_COM(min(squeeze(beta_COM_pool(2:3, 1, strcmpi(respVar, 'Shift'), 2,:))) < 0.05) = 1;
zEffect_shift_COM(min(squeeze(beta_COM_pool(2:3, 1, strcmpi(respVar, 'Shift'), 2,:))) >= 0.05) = 0;
fprintf('\nFound evidence of COM z-effect for velocity -> z-shift in %i of %i experiments', sum(zEffect_shift_COM == 1), sum(~isnan(zEffect_shift_COM)) );

zEffect_shift_AUC = nan(1, Npresent);
zEffect_shift_AUC(min(squeeze(beta_AUC_pool(2:3, 1, strcmpi(respVar, 'Shift'), 2,:))) < 0.05) = 1;
zEffect_shift_AUC(min(squeeze(beta_AUC_pool(2:3, 1, strcmpi(respVar, 'Shift'), 2,:))) >= 0.05) = 0;
fprintf('\nFound evidence of AUC z-effect for velocity -> z-shift in %i of %i experiments', sum(zEffect_shift_AUC == 1), sum(~isnan(zEffect_shift_AUC)) );

