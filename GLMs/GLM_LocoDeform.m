%% Use GLM to assess contribution of different variables
locoDeform_pred = cell(1,Nexpt); locoDeform_resp = cell(1,Nexpt); locoDeform_opts = cell(1,Nexpt); locoDeform_result = cell(1,Nexpt); locoDeform_summary = cell(1,Nexpt);
figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
%rLocoPreFit = cell(1,Nexpt);
GLMname = 'LocoDeformGLM_2D';
GLMrate = 15.49/2;
bindownmean = @(x,n)( nanmean(reshape([x(:); nan(mod(-numel(x),n),1)],n,[]))' ); % bin vector, padding with NaNs if necessary, then take the mean of each bin
for x = 10 %xPresent %x3Dcsd
    % GLMparallel options
    locoDeform_opts{x}.name = sprintf('%s_%s', expt(x).name, GLMname); %strcat(expt(x).name, , '_locoDeformglm');
    %locoDeform_opts{x}.show = true;
    locoDeform_opts{x}.rShow = 1:5;
    locoDeform_opts{x}.figDir = ''; % figDir;
    locoDeform_opts{x}.alpha = 0.01;  % The regularization parameter, default is 0.01
    locoDeform_opts{x}.standardize = true; 
    locoDeform_opts{x}.trainFrac = 0.75; % 1; %
    locoDeform_opts{x}.Ncycle = 20;
    locoDeform_opts{x}.distribution = 'gaussian'; % 'poisson'; %  
    locoDeform_opts{x}.CVfold = 10;
    LocoDeform.opts(x).nlamda = 1000;
    LocoDeform.opts(x).maxit = 5*10^5;
    locoDeform_opts{x}.minDev = 0.1;  minDev = locoDeform_opts{x}.minDev;
    locoDeform_opts{x}.minDevFrac = 0.1;
    locoDeform_opts{x}.maxP = 0.05;
    locoDeform_opts{x}.Nshuff = 0;  Nshuff = locoDeform_opts{x}.Nshuff;
    locoDeform_opts{x}.minShuff = 15; 
    locoDeform_opts{x}.window = [-5,5]; % [0,0]; % [-0.5, 0.5]; % 
    locoDeform_opts{x}.lopo = true; %false; %
    locoDeform_opts{x}.frameRate = GLMrate; %expt(x).scanRate; 
    locoDeform_opts{x}.binSize = expt(x).scanRate/GLMrate;
    locoDeform_opts{x}.minShuffFrame = round( locoDeform_opts{x}.frameRate*locoDeform_opts{x}.minShuff );
    windowFrame = round(locoDeform_opts{x}.window*locoDeform_opts{x}.frameRate);
    locoDeform_opts{x}.shiftFrame = windowFrame(1):windowFrame(2);
    locoDeform_opts{x}.maxShift = max( abs(windowFrame) );
    locoDeform_opts{x}.Nshift = numel( locoDeform_opts{x}.shiftFrame );  %Nshift = locoDeformOpts(x).Nshift;
    locoDeform_opts{x}.lags = locoDeform_opts{x}.shiftFrame/locoDeform_opts{x}.frameRate;

    % Concatenate and bin input variables 
    if ~isnan(expt(x).csd),  preCSDruns = 1:expt(x).csd-1;  else,  preCSDruns = 1:expt(x).Nruns;  end
    tempVelocityCat = bindownmean( vertcat(loco{x}(preCSDruns).Vdown), locoDeform_opts{x}.binSize );
    tempAccelCat = bindownmean( vertcat(loco{x}(preCSDruns).Adown), locoDeform_opts{x}.binSize ); %[NaN; diff(tempVelocityCat)];
    tempSpeedCat = bindownmean( vertcat(loco{x}(preCSDruns).speedDown), locoDeform_opts{x}.binSize );
    tempStateCat = bindownmean( vertcat(loco{x}(preCSDruns).stateDown), locoDeform_opts{x}.binSize );
    tempTransVel = bindownmean( abs( mean( vertcat(deform{x}(preCSDruns).DtransMag), 2, 'omitnan') ), locoDeform_opts{x}.binSize );
    tempScaleMag = bindownmean( mean( vertcat(deform{x}(preCSDruns).scaleMag), 2, 'omitnan'), locoDeform_opts{x}.binSize );
    tempStretchMag = bindownmean( mean( vertcat(deform{x}(preCSDruns).stretchMag), 2, 'omitnan'), locoDeform_opts{x}.binSize );
    tempShearMag = bindownmean( mean( vertcat(deform{x}(preCSDruns).shearMag), 2, 'omitnan'), locoDeform_opts{x}.binSize );
    tempDshearMag = bindownmean( mean( vertcat(deform{x}(preCSDruns).DshearMag), 2, 'omitnan'), locoDeform_opts{x}.binSize ); %circshift(tempAccelCat, 5); %
    if expt(x).Nplane > 1    
        tempShift = vertcat(deform{x}(preCSDruns).shiftZ);
        tempShiftMean =  bindownmean( mean( tempShift(:,4:end-3), 2, 'omitnan' ), locoDeform_opts{x}.binSize );
        tempCompress = vertcat(deform{x}(preCSDruns).compressZ);
        tempCompressMean = bindownmean( mean( tempCompress(:,4:end-2), 2, 'omitnan' ), locoDeform_opts{x}.binSize );
    end

    % Define predictors
    LocoDeform_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); % LocoDeformPred(x).data = []; %LocoDeformPred(x).data = LocoDeformCat.T; LocoDeformPred(x).name = {'Time'};
    LocoDeform_pred{x}.data = [tempVelocityCat, tempAccelCat, tempStateCat]; % tempSpeedCat,, 
    LocoDeform_pred{x}.name = {'Velocity','Acceleration','State'}; % ,'Speed',
    LocoDeform_pred{x}.N = size(LocoDeform_pred{x}(r).data,2);
    for p = flip(1:LocoDeform_pred{x}(r).N), LocoDeform_pred{x}(r).lopo.name{p} = ['No ',LocoDeform_pred{x}(r).name{p}]; end
    % Set up leave-one-family-out
    LocoDeform_pred{x}.fam.col = {1:2, 3}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    LocoDeform_pred{x}.fam.N = numel(LocoDeform_pred{x}.fam.col); 
    LocoDeform_pred{x}.fam.name = {'No Kinematics', 'No State'}; %{'All'};%  'Onset Time',

    % Define response
    LocoDeform_resp{x}.data = [tempTransVel, tempScaleMag, tempStretchMag, tempShearMag, tempDshearMag]; % tempScaleMag; % , tempShiftMean, tempCompressMean
    LocoDeform_resp{x}.name = {'Trans Speed','|Scale|','|Stretch|','|Shear|', 'd|Shear|/dt'}; % , 'Z Shift', 'Compression' 
    LocoDeform_resp{x}.N = size(LocoDeform_resp{x}.data,2); 
    
    % Remove scans with missing data 
    nanFrame = find(any(isnan([LocoDeform_pred{x}.data, LocoDeform_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    LocoDeform_pred{x}.data(nanFrame,:) = []; LocoDeform_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    locoDeform_opts{x}.load = false; %  true; %
    locoDeform_opts{x}.saveRoot = sprintf('%s%s', expt(x).dir  ); %''; % , locoDeform_opts{x}.name

    [LocoDeform_result{x}, LocoDeform_summary{x}, ~, LocoDeform_pred{x}, LocoDeform_resp{x}] = GLMparallel(LocoDeform_pred{x}, LocoDeform_resp{x}, locoDeform_opts{x}); % LocoDeform_result{x}, LocoDeform_summary(x), ~, LocoDeformPred(x), LocoDeformResp(x)
    %ViewGLM(LocoDeform_pred{x}, LocoDeform_resp{x}, locoDeform_opts{x}, LocoDeform_result{x}, LocoDeform_summary{x}); %GLMresultFig = 
    %pause;
end

%% View GLM results
for x = xPresent
    locoDeform_opts{x}.rShow = 1:5;
    ViewGLM(LocoDeform_pred{x}, LocoDeform_resp{x}, locoDeform_opts{x}, LocoDeform_result{x}, LocoDeform_summary{x}); %GLMresultFig = 
end
%%
devSummMat = cell(1,Nexpt);
GLMresultFig = figure('WindowState','maximized', 'color','w');
for x = xPresent
    cla;
    devSummMat{x} = [LocoDeform_summary{x}.dev; LocoDeform_summary{x}.lopo.dev; LocoDeform_summary{x}.lofo.dev];
    % {
    imagesc( devSummMat{x}' );
    set(gca,'XtickLabel',  [{'All'}; LocoDeform_pred{xPresent(1)}.lopo.name(:); LocoDeform_pred{xPresent(1)}.fam.name(:)], ...
        'Ytick',1:LocoDeform_resp{xPresent(1)}.N , 'YtickLabel',LocoDeform_resp{xPresent(1)}.name, 'TickDir','out', 'TickLength',[0.003,0]);
    title(sprintf('%s: GLM Fit Summary', locoDeform_opts{x}.name), 'Interpreter','none');
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

GLMdevFig = figure('WindowState','maximized', 'color','w');
imagesc( devSummMed' );
set(gca,'XtickLabel',  [{'All'}; LocoDeform_summary{x}.lopo.name(:); LocoDeform_summary{x}.lofo.name(:)], ...
    'YtickLabel',LocoDeform_resp{x}.name, 'TickDir','out', 'TickLength',[0.003,0]);
title('GLM Fit Summary (All 3D Data)', 'Interpreter','none');
xtickangle(30);
axis image;
CB = colorbar; CB.Label.String = 'Median Deviance Explained';
figPath = sprintf('%s%s_DevSummary.tif', figDir, GLMname );
fprintf('\nSaving %s\n', figPath);
print( GLMdevFig, figPath, '-dtiff');
impixelinfo;
