function [ theta, stretchAngStats, stretchDirPref, stretchOrPref, testRange, dirResponse, dirPref, orResponse, orPref ] = AngularSensitivity( MLstretch, APstretch, affActivity, varargin )
%OrientationSensitivity  

IP = inputParser;
addRequired( IP, 'MLstretch', @isnumeric )
addRequired( IP, 'MLactivity', @isnumeric )
addRequired( IP, 'activity', @isnumeric )
%addRequired( IP, 'deform', @isstruct )
%addRequired( IP, 'activity', @isstruct )
addParameter( IP, 'show', false, @checkShow ) % 
addParameter( IP, 'dTheta', pi/6, @isnumeric )
addParameter( IP, 'dPct', 10, @isnumeric )
addParameter( IP, 'TestPct', [60,90], @isnumeric )
addParameter( IP, 'rot2stretch', false, @islogical );
addParameter( IP, 'Nshuff', 1000, @isnumeric )
addParameter( IP, 'pMax', 0.05, @isnumeric )
parse( IP, MLstretch, APstretch, affActivity, varargin{:} ); 
show = IP.Results.show;
dTheta = IP.Results.dTheta;
dPct = IP.Results.dPct;
testPct = IP.Results.TestPct;
Nshuff = IP.Results.Nshuff;
pNormMax = IP.Results.pMax;
rot2stretch = IP.Results.rot2stretch;

% Get the angular distribution of stretch data
[ theta, thetaDirData, thetaOrData, magData, ~, ~, stretchAngStats ] = AngularDist( MLstretch, APstretch, 'dPct', dPct, 'show',false, 'dThetaDir', dTheta, 'dThetaOr', dTheta/2 ); % 

% Define test range: what magnitude of stretch is consistently observed across all angles?
[~,testPctInd(2)] = min( abs(stretchAngStats.dir.pctRange - testPct(2)) );  [~,testPctInd(1)] = min( abs(stretchAngStats.dir.pctRange - testPct(1)) );
[testMin, ~] = min( stretchAngStats.dir.pct(:,testPctInd(1)) );  % testIndT
[testMax, ~] = min( stretchAngStats.dir.pct(:,testPctInd(2)) ); % test range extends from minimum 80th percentile value to minimum of 100th percentile
testRange = [testMin, testMax];  %[testMin, stretchAngStats.dir.pct(testIndT,end)]; % test range extends from minimum 80th percentile value to the maximum value at that angle

% Quantify the extent to which the observed stretch is concentrated along a specific direction/orientation 
stretchDirPref = AngularPreference( theta.dir(1:theta.Ndir-1), stretchAngStats.dir.pct(1:theta.Ndir-1,testPctInd(1)), pi, 0, false );
stretchOrPref  = AngularPreference( theta.or(1:theta.Nor), stretchAngStats.or.pct(1:theta.Nor,testPctInd(1)), pi/2, 0, false );

% Rotate coord system so that angle of maximum stretch orientation = 0
if rot2stretch
    referenceAngle = stretchOrPref.max.angle;
else
    referenceAngle = 0; 
end
thetaDirRot = thetaDirData - referenceAngle;
thetaDirRot(thetaDirRot < -pi) = thetaDirRot(thetaDirRot < -pi) + 2*pi; thetaDirRot(thetaDirRot >= pi) = thetaDirRot(thetaDirRot >= pi) - 2*pi; % direction: [-pi, pi)
thetaOrRot = thetaOrData - referenceAngle;
thetaOrRot(thetaOrRot < -pi/2) = thetaOrRot(thetaOrRot < -pi/2) + pi; thetaOrRot(thetaOrRot > pi/2) = thetaOrRot(thetaOrRot > pi/2) - pi; % orientation: [-pi, pi]
%{
figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
subplot(2,3,1);
polarhistogram(thetaOrData); title( sprintf('Unrotated Orientation (ref angle = %2.2f)', referenceAngle ));
%polarplot( thetaOrData, magData, '.' ); 
subplot(2,3,2);
polarhistogram(thetaOrData-referenceAngle); title('Rotated, but not wrapped');
%polarplot( thetaOrData-referenceAngle, magData, '.' ); 
subplot(2,3,3);
polarhistogram(thetaOrRot); title('Rotated, and wrapped');
%polarplot( thetaOrRot, magData, '.');
subplot(2,3,4);
polarhistogram(thetaDirData); title('Unrotated Direction');
%polarplot( thetaOrData, magData, '.' ); 
subplot(2,3,5);
polarhistogram(thetaDirData-referenceAngle);
%polarplot( thetaOrData-referenceAngle, magData, '.' ); 
subplot(2,3,6);
polarhistogram(thetaDirRot);
%}



% How do afferents respond, directionally, to stretches in the specific test-range defined above?
for t = flip(2:theta.Ndir-1)   
    dirTestInd{t} = find( magData >= testRange(1) & magData <= testRange(2) & thetaDirRot >= theta.dirBin(t) & thetaDirRot < theta.dirBin(t+1) ); %#ok<AGROW>
    dirResponse.N(t) = numel( dirTestInd{t} );
    dirResponse.mean(t,:) = mean( affActivity(dirTestInd{t},:), 1 );
    dirResponse.std(t,:) = std( affActivity(dirTestInd{t},:), 0, 1 );
    dirResponse.sem(t,:) = dirResponse.std(t,:)/sqrt(dirResponse.N(t));
end
dirTestInd{1} = find( magData >= testRange(1) & magData <= testRange(2) & (thetaDirRot >= theta.dirBin(1) | thetaDirRot < theta.dirBin(2)) );
dirResponse.N(1) = numel( dirTestInd{1} );
dirResponse.mean(1,:) = mean( affActivity(dirTestInd{1},:), 1 );
dirResponse.std(1,:) = std( affActivity(dirTestInd{1},:), [], 1 );
dirResponse.sem(1,:) = dirResponse.std(1,:)/sqrt(dirResponse.N(1));
% Close the circles for plotting purposes
dirResponse.mean(theta.Ndir,:) = dirResponse.mean(1,:);
dirResponse.std(theta.Ndir,:) = dirResponse.std(1,:);
dirResponse.sem(theta.Ndir,:) = dirResponse.sem(1,:);

% How do afferents respond, orientationally, to stretches in the specific test-range defined above?
for t = flip(1:theta.Nor)   
    orTestInd{t} = find( magData >= testRange(1) & magData <= testRange(2) & thetaOrRot >= theta.orBin(t) & thetaOrRot < theta.orBin(t+1) );
    orResponse.N(t) = numel( orTestInd{t} );
    orResponse.mean(t,:) = mean( affActivity(orTestInd{t},:), 1 );
    orResponse.std(t,:) = std( affActivity(orTestInd{t},:), [], 1 );
    orResponse.sem(t,:) = orResponse.std(t,:)/sqrt(orResponse.N(t));
end

% Calculate test-range orientation selectivity for each afferent
Nroi = size(affActivity,2);
fprintf('\nCalculating angular preferences for afferents...  '); tic
for r = flip(1:Nroi) %1:Nroi
    dirPref(r) = AngularPreference( theta.dir(1:theta.Ndir-1), dirResponse.mean(1:theta.Ndir-1,r), pi, Nshuff, false ); %#ok<AGROW>
    %pause; close all
    orPref(r) = AngularPreference( theta.or(1:theta.Nor), orResponse.mean(1:theta.Nor,r), pi/2, Nshuff, false ); %#ok<AGROW>
    %pause; close all
end
toc
%Show the results (optional)
if ~strcmpi(show,'none')
    if any(strcmpi(show, {'all','indv'}))
        opt = {[0.10,0.05], [0.10,0.1], [0.03,0.03]};  % {[vert, horz], [bottom, top], [left, right] } 
        figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
        subtightplot(2,2,1,opt{:})
        h(1) = polarplot( theta.dir, stretchAngStats.dir.mean, 'k', 'LineWidth',1 ); hold on;%pause;
        h(2) = polarplot( theta.dir, stretchAngStats.dir.pct(:,testPctInd(1)), 'b' );
        polarplot( stretchDirPref.max.angle, stretchDirPref.max.mag, 'bx' )
        polarplot( stretchDirPref.orth.angle, stretchDirPref.orth.mag, 'bo' )
        h(3) = polarplot( theta.dir, stretchAngStats.dir.pct(:,testPctInd(2)), 'g' );
        h(4) = polarplot( theta.dir, repmat(testRange(1),theta.Ndir,1), 'k--' );
        h(5) = polarplot( [stretchDirPref.angle; stretchDirPref.angle], [0; stretchDirPref.mag], 'k', 'LineWidth',2 );
        polarplot( theta.dir, repmat(testRange(2), theta.Ndir, 1), 'k--' );
        legend(h, {'Mean','80th','100th','Test','Resultant' }, 'Location','NorthWest', 'AutoUpdate','off');
        title( sprintf('Directional Distribution of Stretch: Norm Mag = %1.3f,  Max-Orth Index = %1.3f', stretchDirPref.norm, stretchDirPref.index  ) ); %title('Directional Distribution of Stretch');

        subtightplot(2,2,2,opt{:})
        h(1) = polarplot( theta.or, stretchAngStats.or.mean, 'k', 'LineWidth',1 ); hold on;%pause;
        h(2) = polarplot( theta.or, stretchAngStats.or.pct(:,testPctInd(1)), 'b' );
        polarplot( stretchOrPref.max.angle, stretchOrPref.max.mag, 'bx' )
        polarplot( stretchOrPref.orth.angle, stretchOrPref.orth.mag, 'bo' )
        h(3) = polarplot( theta.or, stretchAngStats.or.pct(:,testPctInd(2)), 'g' );
        h(4) = polarplot( theta.or, repmat(testRange(1),theta.Nor,1), 'k--' );
        h(5) = polarplot( [stretchOrPref.angle; stretchOrPref.angle], [0; stretchOrPref.mag], 'k', 'LineWidth',2 );
        thetalim([-90,90]);
        legend(h, {'Mean','80th','100th','Test','Resultant' }, 'Location','NorthWest', 'AutoUpdate','off');
        polarplot( theta.or, repmat(testRange(2),theta.Nor,1), 'k--' );
        title( sprintf('Orientational: Norm Mag = %1.3f,  Max-Orth Index = %1.3f', stretchOrPref.norm, stretchOrPref.index  ) );
        
        for r = 1:Nroi
            subtightplot(2,2,3,opt{:}); cla;
            polarplot( theta.dir, dirResponse.mean(:,r), 'k', 'LineWidth',2 ); hold on;
            polarplot( theta.dir, dirResponse.mean(:,r)-dirResponse.sem(:,r), 'k--', 'LineWidth',1 ); 
            polarplot( theta.dir, dirResponse.mean(:,r)+dirResponse.sem(:,r), 'k:', 'LineWidth',1 ); 
            polarplot( [dirPref(r).angle; dirPref(r).angle], [0; dirPref(r).mag], 'k', 'LineWidth',2 );
            title( sprintf('r = %i: Afferent Response to Stretch Direction (rotated). Norm mag = %1.2f (p=%1.4f)', r, dirPref(r).norm, dirPref(r).pNorm ) ); 
            
            subtightplot(2,2,4,opt{:}); cla;
            polarplot( theta.or, orResponse.mean(:,r), 'k', 'LineWidth',2 ); hold on;
            polarplot( theta.or, orResponse.mean(:,r)-orResponse.sem(:,r), 'k--', 'LineWidth',1 ); 
            polarplot( theta.or, orResponse.mean(:,r)+orResponse.sem(:,r), 'k:', 'LineWidth',1 ); 
            polarplot( [orPref(r).angle; orPref(r).angle], [0; orPref(r).mag], 'k', 'LineWidth',2 );
            title( sprintf('Afferent Response to Stretch Orientation (rotated). Norm mag = %1.2f (p=%1.4f)', orPref(r).norm, orPref(r).pNorm ) ); 
            thetalim([-90,90]);
            pause;  
        end
    end

    % Summarize results
    if any(strcmpi(show, {'all','summ'}))
        opt = {[0.10,0.05], [0.10,0.1], [0.06,0.05]};  % {[vert, horz], [bottom, top], [left, right] } 
        figure('Units','normalized', 'OuterPosition',[0,0,1,1], 'Color','w', 'PaperOrientation','landscape');
        subtightplot(2,2,1,opt{:})
        histogram( [dirPref.norm], 12 ); 
        xlabel('Normalized Magnitude'); ylabel('Frequency'); title('Directional');
        box off; axis square;
        
        subtightplot(2,2,3,opt{:})
        oriTheta = [dirPref.angle];
        polarhistogram( oriTheta( [dirPref.pNorm] < pNormMax ), 12 ) 
        title( sprintf('Distribution of Preferred Direction (p < %1.2f)', pNormMax) );
        
        subtightplot(2,2,2,opt{:})
        histogram( [orPref.norm], 12 ); 
        xlabel('Normalized Magnitude'); ylabel('Frequency'); title('Orientational');
        box off;  axis square;
        subtightplot(2,2,4,opt{:})
        oriTheta = [orPref.angle];
        polarhistogram( oriTheta( [orPref.pNorm] < pNormMax ), 12 ) 
        title( sprintf('Distribution of Preferred Orientation (p < %1.2f)', pNormMax) );
    end
end
end