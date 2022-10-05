%% Use GLM to assess contribution of different variables
locoDeform_vol_pred = cell(1,Nexpt); locoDeform_vol_resp = cell(1,Nexpt); locoDeform_vol_opts = cell(1,Nexpt); locoDeform_vol_result = cell(1,Nexpt); locoDeform_vol_summary = cell(1,Nexpt);
%figDir = 'D:\MATLAB\LevyLab\Figures\3D\GLM\'; mkdir( figDir )
%rLocoPreFit = cell(1,Nexpt);
GLMname = 'LocoDeform_3D_planar';
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

    % Define predictors
    tempVelocityCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).Vdown), locoDeform_vol_opts{x}.binSize );
    tempAccelCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).Adown), locoDeform_vol_opts{x}.binSize ); %[NaN; diff(tempVelocityCat)];
    tempSpeedCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).speedDown), locoDeform_vol_opts{x}.binSize );
    tempStateCat = BinDownMean( vertcat(loco{x}(expt(x).preRuns).stateDown), locoDeform_vol_opts{x}.binSize );
    locoDeform_vol_pred{x} = struct('data',[], 'name',[], 'N',NaN, 'TB',[], 'lopo',[], 'fam',[]); % LocoDeformPred(x).data = []; %LocoDeformPred(x).data = LocoDeformCat.T; LocoDeformPred(x).name = {'Time'};
    locoDeform_vol_pred{x}.data = [tempVelocityCat, tempAccelCat, tempStateCat]; % tempSpeedCat,, 
    locoDeform_vol_pred{x}.name = {'Velocity','Acceleration','State'}; % ,'Speed',
    locoDeform_vol_pred{x}.N = size(locoDeform_vol_pred{x}(r).data,2);
    for p = flip(1:locoDeform_vol_pred{x}(r).N), locoDeform_vol_pred{x}(r).lopo.name{p} = ['No ',locoDeform_vol_pred{x}(r).name{p}]; end
    % Set up leave-one-family-out
    locoDeform_vol_pred{x}.fam.col = {1:2, 3}; %{1:2, 3:4, 5:6, 7:8, 9:10, 11:12};  % {1:12};%{1, 2:3, 4:5, 6:7, 8, 9}; 
    locoDeform_vol_pred{x}.fam.N = numel(locoDeform_vol_pred{x}.fam.col); 
    locoDeform_vol_pred{x}.fam.name = {'Kinematics', 'State'}; %{'All'};%  'Onset Time',
    
    % Define responses
    tempTrans = vertcat(deform{x}(expt(x).preRuns).transMag);
    tempTrans = BinDownMean( tempTrans(:,segParams{x}.zProj) , locoDeform_vol_opts{x}.binSize ); % (:,segParams{x}.zProj) 
    transName = sprintfc('Trans_Z%i', segParams{x}.zProj);
    tempTransSpd = vertcat(deform{x}(expt(x).preRuns).DtransMag);
    tempTransSpd = BinDownMean( tempTransSpd(:,segParams{x}.zProj) , locoDeform_vol_opts{x}.binSize );
    transSpdName = sprintfc('TransSpd_Z%i', segParams{x}.zProj);
    tempScaleMag = vertcat(deform{x}(expt(x).preRuns).scaleMag);
    tempScaleMag = BinDownMean( tempScaleMag(:,segParams{x}.zProj) , locoDeform_vol_opts{x}.binSize );
    scaleName = sprintfc('Scale_Z%i', segParams{x}.zProj);
    tempStretchMag = vertcat(deform{x}(expt(x).preRuns).stretchMag);
    tempStretchMag = BinDownMean( tempStretchMag(:,segParams{x}.zProj) , locoDeform_vol_opts{x}.binSize );
    stretchName = sprintfc('Stretch_Z%i', segParams{x}.zProj);
    tempShearMag = vertcat(deform{x}(expt(x).preRuns).shearMag);
    tempShearMag = BinDownMean( tempShearMag(:,segParams{x}.zProj) , locoDeform_vol_opts{x}.binSize );
    shearName = sprintfc('Shear_Z%i', segParams{x}.zProj);
    tempShearRate = vertcat(deform{x}(expt(x).preRuns).DshearMag);
    tempShearRate = BinDownMean( tempShearRate(:,segParams{x}.zProj) , locoDeform_vol_opts{x}.binSize );
    shearRateName = sprintfc('ShearRate_Z%i', segParams{x}.zProj);
    tempShift = vertcat(deform{x}(expt(x).preRuns).shiftZ);
    shiftPlanes = intersect(4:expt(x).Nplane-3, segParams{x}.zProj);
    tempShiftMean =  BinDownMean( tempShift(:,shiftPlanes) , locoDeform_vol_opts{x}.binSize );
    shiftName = sprintfc('Shift_Z%i', shiftPlanes);

    locoDeform_vol_resp{x}.data = [tempTrans, tempTransSpd, tempScaleMag, tempStretchMag, tempShearMag, tempShearRate, tempShiftMean];
    %locoDeform_vol_resp{x}.data = normalize(locoDeform_vol_resp{x}.data, 1);
    %locoDeform_vol_resp{x}.data(abs(locoDeform_vol_resp{x}.data) > 5) = NaN; % suppress extreme outliers
    locoDeform_vol_resp{x}.name = [transName, transSpdName, scaleName, stretchName, shearName, shearRateName, shiftName]; %{'|Translation|','Trans. Speed','|Scale|','|Stretch|','|Shear|', 'Shear Rate','Z Shift'}; % , 'Z Shift', 'Compression' , 'Z Compression' 
    locoDeform_vol_resp{x}.N = size(locoDeform_vol_resp{x}.data,2); 
    
    % Remove scans with missing data 
    nanFrame = find(any(isnan([locoDeform_vol_pred{x}.data, locoDeform_vol_resp{x}.data]),2)); % find( isnan(sum(pred(x).data,2)) ); 
    fprintf('\nRemoving %i NaN-containing frames', numel(nanFrame));
    locoDeform_vol_pred{x}.data(nanFrame,:) = []; locoDeform_vol_resp{x}.data(nanFrame,:) = [];

    % Run the GLM
    locoDeform_vol_opts{x}.load = true; %false; %  
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

