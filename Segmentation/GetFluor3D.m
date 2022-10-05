function fluor = GetFluor3D(sbxInfo, sbxPath, varargin) % mouse, exptDate expt, run
% Get data from registered 3D calcium imaging experiments

IP = inputParser;
addRequired( IP, 'sbxInfo', @isstruct )
addParameter( IP, 'path', '', @ischar)
addParameter( IP, 'chan', 'green', @ischar)
addParameter( IP, 'save', true, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, sbxInfo, sbxPath, varargin{:} );  %  expt, run
sbxPath = IP.Results.path;
chan = IP.Results.chan;
saveToggle = IP.Results.save;
overwrite = IP.Results.overwrite;

fSep = '\';
usePMTind = DeterminePMT(chan, sbxInfo);
if isempty(sbxPath)
    sbxPath = sbxInfo.path;
end
[fDir, fName, ~] = fileparts(sbxPath); 
savePath = strcat( fDir, fSep, fName, '_fluor.mat' );
tic;
if ~exist(savePath, 'file') || overwrite
    fprintf('Extracting fluor signals from %s', sbxPath );
    sbxData = readSBX(sbxPath, sbxInfo, 1, sbxInfo.totScan, usePMTind, []);  % sbxInfo.totScan vs sbxInfo.Nscan
    %
    fluor = struct('F',[], 'Fo',[], 'dFF',[]);
    fluor.F = struct('vol',nan(sbxInfo.totScan,1), 'plane',nan(sbxInfo.totScan, sbxInfo.Nplane) ); 
    fluor.F = struct('vol',nan(sbxInfo.totScan,1), 'plane',nan(sbxInfo.totScan, sbxInfo.Nplane) ); 
    fluor.dFF = struct('vol',nan(sbxInfo.totScan,1), 'plane',nan(sbxInfo.totScan, sbxInfo.Nplane) ); 
    % Average over whole volume
    fluor.F.vol = squeeze( mean( mean( mean( sbxData, 3), 2 ), 1) ); 
    fluor.Fo.vol = prctile(fluor.F.vol, 10);
    fluor.dFF.vol = (fluor.F.vol - fluor.Fo.vol)/fluor.Fo.vol;
    fluor.z.vol = zscore( fluor.dFF.vol, 0, 1 );
    % Average by plane
    fluor.F.plane = squeeze( mean( mean( sbxData, 2 ), 1) )'; 
    fluor.Fo.plane = prctile(fluor.F.plane, 10, 1);
    fluor.dFF.plane = (fluor.F.plane - fluor.Fo.plane)./fluor.Fo.plane;
    fluor.z.plane = zscore( fluor.dFF.plane, 0, 1 );
    
    % Calculate statistics about each plane
    fluor.stat.mean = squeeze( mean( mean( mean(sbxData, 1), 2), 4) )';
    fluor.stat.std = squeeze( std( mean( mean(sbxData, 1), 2),0, 4) )';

    if saveToggle
        fprintf('\nSaving %s... ', savePath);
        save( savePath, 'fluor' ); % 
    end
else
    fprintf('\nLoading %s... ', savePath);
    load( savePath );
end
toc;
end

