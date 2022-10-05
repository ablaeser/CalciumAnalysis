function [event, Nevent, eventRate, eventRaster] = DetectEvents(T, zData, magData, varargin) % expt, 
% Parse inputs
IP = inputParser;
%addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'T', @isnumeric )
addRequired( IP, 'zData', @isnumeric )
addOptional( IP, 'magData', [], @isnumeric )
addParameter( IP, 'minZ', 1, @isnumeric )
addParameter( IP, 'minDur', 0, @isnumeric )
addParameter( IP, 'minMag', NaN, @isnumeric )
addParameter( IP, 'show', false, @islogical )
parse( IP, T, zData, magData, varargin{:} ); % expt,
magData = IP.Results.magData; 
minZ =IP.Results.minZ;
minDur = IP.Results.minDur;
minMag = IP.Results.minMag;
show = IP.Results.show;
dT = (T(2)-T(1));
scanRate = 1/dT;
minScan = round(minDur*scanRate);
magToggle = ~isempty(magData);
if ~magToggle, putMagMax = NaN; end
%if isnan(minDFF) && dFFtoggle, error('
[Nscan, Nroi] = size(zData);
event = cell(1,Nroi);
Nevent = nan(1,Nroi);
eventRaster = false(Nscan, Nroi);

% Binarize the z data
binaryData = zData > minZ;
%subplot(2,1,1); imagesc(zData); caxis([-1,minZ]);
%subplot(2,1,2); imagesc(binaryData);

% Find events for each ROI
%tic;
for roi = 1:Nroi
    putEvent = regionprops( binaryData(:,roi), zData(:,roi), 'Area', 'PixelIdxList', 'MaxIntensity');  % 
    %putDur = [putEvent.Area]/scanRate;
    %putEvent(putDur < minDur) = []; 
    putEvent([putEvent.Area] < minScan) = [];
    Nput = numel(putEvent);
    E = 0;
    if Nput > 0
        for e = 1:Nput
            putRunScan = putEvent(e).PixelIdxList;
            Tput = T(putRunScan); 
            if magToggle % Check that the event satisfies minimum magnitude requirement
                putMag = magData(putRunScan,roi);
                [putMagMax, peakScan] = max(putMag);
            else
                putZ = zData(putRunScan,roi);
                [~, peakScan] = max(putZ);
            end
            if ~magToggle || putMagMax > minMag
                E = E+1; 
                event{roi}(E).runScan = putRunScan; % scans within the run
                event{roi}(E).T = Tput;
                event{roi}(E).Tstart = Tput(1);
                event{roi}(E).Tstop = Tput(end);
                event{roi}(E).Tpeak = Tput(peakScan);
                event{roi}(E).dur = putEvent(e).Area/scanRate;
                event{roi}(E).zPeak = putEvent(e).MaxIntensity;
                event{roi}(E).magPeak = putMagMax;  
                event{roi}(E).magAUC = dT*sum(putMag);
                eventRaster(event{roi}(E).runScan,roi) = true; 
            end
        end
    end
    Nevent(roi) = E;   
    % Determine isolation of each event
    if Nevent(roi) > 0
        tempStart = [event{roi}.Tstart];
        tempStop = [event{roi}.Tstop];
        tempIso =  [event{roi}(1).Tstart-T(1), tempStart(2:end)-tempStop(1:end-1), T(end)-event{roi}(end).Tstop];
        for e = 1:Nevent(roi), event{roi}(e).iso = [tempIso(e), tempIso(e+1)]; end
    end
end
eventRate = Nevent/(60*(T(end)-T(1)));
%toc
if show
    figure;
    subplot(4,1,1);
    imagesc( zData' );
    caxis([-1,3]);
    
    subplot(4,1,2);
    imagesc( binaryData' );
    title('Binarized');
    
    subplot(4,1,3);
    imagesc( eventRaster' );
    title('Event Raster');
    impixelinfo;
    
    subplot(4,1,4);
    for roi = 1:Nroi
        cla;
        plot(T, zData(:,roi)); hold on;
        line([T(1),T(end)], minZ*[1,1], 'color','r', 'lineStyle','--');
        axis tight;
        for e = 1:Nevent(roi)
            plot(event{roi}(e).T, zData(event{roi}(e).runScan,roi), 'k' );
        end
        ylabel('Z-score'); xlabel('Time (s)');
        title( sprintf('ROI %i', roi) );
        pause;
    end
end
end