function [ROI, varargout] = WriteROIproj(expt, catInfo, ROI, varargin)
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'catInfo', @isstruct )
addRequired( IP, 'ROI', @isstruct )
addOptional( IP, 'axon', [], @isstruct )
addParameter( IP, 'edges', [0,0,0,0], @isnumeric )
addParameter( IP, 'buffer', 10*ones(1,4), @isnumeric )
addParameter( IP, 'thickness', 1, @isnumeric)
addParameter( IP, 'rSet', 1:expt.Nroi, @isnumeric)
addParameter( IP, 'overwrite', false, @islogical)
parse( IP, expt, catInfo, ROI, varargin{:} );
axon = IP.Results.axon;

segEdges = IP.Results.edges;
edgeBuffer = IP.Results.buffer;
zThick = IP.Results.thickness;
rSet = IP.Results.rSet;
Nset = numel(rSet);
overwrite = IP.Results.overwrite;
projDir = [expt.dir, 'ROI\']; mkdir( projDir );
if numel(edgeBuffer) == 1, edgeBuffer = edgeBuffer*[1,1,1,1]; end % edgeBuffer = [left, right, top, bottom]
if nargout > 1, crop_proj_out = cell(1,expt.Nroi); end
Nroi = numel(ROI);
Naxon = numel(axon);
% Generate projections over relevant planes
%{
allCent = vertcat(ROI.cent ); %unique([ROI.zUse]);
zUse = unique( round( allCent(:,3) ) )';
allZproj = cell(1,expt.Nplane);
tic
for z = zUse
    fprintf('\nz = %i (thickness = %i)', z, zThick);
    zProj = z-zThick:z+zThick;
    zProj(zProj<1) = []; zProj(zProj>expt.Nplane) = [];
    allZproj{z} = pipe.zproj(expt.sbx, 1, expt.totScan-1, 1, zProj, 'mtype','.sbx_interp', 'registration',false);
    toc
end
%}
if expt.Nplane > 1
    zCurrent = NaN;
    if isempty(axon)
        w = waitbar(0, 'Generating ROI z-projection tifs');
        for r = rSet% 1:expt.Nroi
            ROI(r).projPath = sprintf('%s%s_ROI_%02i.tif', projDir, expt.name, r );
            tic;
            if ~exist(ROI(r).projPath, 'file') || overwrite
                % Perform z-projection if it doesn't already exist (taking advantage of the fact that ROI are sorted by zCent)
                fprintf('\nROI %i:  ', r); % ,  expt.sbx
                zCent = round(ROI(r).cent(1,3));
                if zCent ~= zCurrent
                    zCurrent = zCent;
                    zProj = zCent-zThick:zCent+zThick;
                    zProj(zProj < 1) = []; 
                    zProj(zProj > expt.Nplane) = [];
                    fprintf(' projecting  z = %i (thickness = %i)', zCent, zThick);
                    ROI_proj = pipe.zproj(expt.sbx, 1, expt.totScan-1, 1, zProj, 'mtype','.sbx_interp', 'registration',false); % last interpolated scan is bad
                end

                % Determine range for cropping and z-projection
                ROI(r).crop(1) = max(floor(ROI(r).box_x(1))-edgeBuffer(1), segEdges(1)+1); % left edge
                ROI(r).crop(2) = min(ceil(ROI(r).box_x(2))+edgeBuffer(2), expt.Ncol-segEdges(2)); % right edge
                ROI(r).crop(3) = max(floor(ROI(r).box_y(1))-edgeBuffer(3), segEdges(3)+1); % top edge
                ROI(r).crop(4) = min(ceil(ROI(r).box_y(2))+edgeBuffer(4), expt.Nrow-segEdges(4)); % bottom edge
                crop_proj = ROI_proj(ROI(r).crop(3):ROI(r).crop(4), ROI(r).crop(1):ROI(r).crop(2),:); % Crop projection
                % Write tiff
                if ~isempty(crop_proj)      
                    %ROI(r).projPath = sprintf('%s%s_ROI_%02i.tif', projDir, expt.name, r );
                    fprintf('\nWriting %s... ', ROI(r).projPath );
                    pipe.io.writeTiff(uint16(crop_proj), ROI(r).projPath); 
                else
                    fprintf('\nWriting %s failed, crop_proj is empty... ', ROI(r).projPath );
                end   
            else
                fprintf('\nLoading %s...  ', ROI(r).projPath );
                crop_proj = double(loadtiff(ROI(r).projPath)); 
            end

            % Save mean, std dev and max projections
            ROI(r).cropMean = uint16(mean(crop_proj, 3));
            ROI(r).cropStd = uint16(std(crop_proj, 0, 3));
            ROI(r).cropMax = uint16(max(crop_proj, [], 3));

            if nargout > 1,  crop_proj_out{r} = crop_proj;  end
            waitbar(r/expt.Nroi, w)
            toc
        end
    else
        w = waitbar(0, 'Writing axon projection tifs');
        for a = 1:Naxon
            axon(a).projPath = sprintf('%s%s_axon_%02i.tif', projDir, expt.name, a );
            if ~exist(axon(a).projPath, 'file') || overwrite
                axonROIcent = vertcat( ROI(axon(a).ROI).cent);
                zCent = round( median( axonROIcent(:,3) ) );
                if zCent ~= zCurrent
                    fprintf(' projecting  z = %i (thickness = %i)', zCent, zThick);
                    zProj = zCent-zThick:zCent+zThick;
                    zProj(zProj<1) = []; zProj(zProj>expt.Nplane) = [];
                    axon_proj = pipe.zproj(expt.sbx, 1, expt.totScan-1, 1, zProj, 'mtype','.sbx_interp', 'registration',false); % last interpolated scan is bad
                end
                tempCrop(1) = max(floor(axon(a).crop(1))-edgeBuffer(1), segEdges(1)+1); % left edge
                tempCrop(2) = min(ceil(axon(a).crop(2))+edgeBuffer(2), expt.Ncol-segEdges(2)); % right edge
                tempCrop(3) = max(floor(axon(a).crop(3))-edgeBuffer(3), segEdges(3)+1); % top edge
                tempCrop(4) = min(ceil(axon(a).crop(4))+edgeBuffer(4), expt.Nrow-segEdges(4)); % bottom edge
                crop_proj = axon_proj(tempCrop(3):tempCrop(4), tempCrop(1):tempCrop(2),:); % Crop projection
                % Write tiff
                if ~isempty(crop_proj)      
                    fprintf('\nWriting %s... ', axon(a).projPath );
                    pipe.io.writeTiff(uint16(crop_proj), axon(a).projPath); 
                else
                    fprintf('\nWriting %s failed, crop_proj is empty... ', axonProjPath );
                end   
                waitbar(a/Naxon, w)
                toc 
            else
                fprintf('\nLoading %s...  ', axon(a).projPath );
                crop_proj = double(loadtiff(axon(a).projPath)); 
            end
            % Save mean, std dev and max projections
            if nargout > 1,  crop_proj_out{a} = crop_proj;  end
            if nargout > 2
                axon(a).cropMean = uint16(mean(crop_proj, 3));
                axon(a).cropStd = uint16(std(crop_proj, 0, 3));
                axon(a).cropMax = uint16(max(crop_proj, [], 3));
            end
        end
    end
else
    %D:\2photon\DL72\170614_DL72\Ztifs\DL72_170614_affZ01_green.tif
    movPath = sprintf('%sZtifs\\%s_affZ01_green.tif', expt.dir, expt.name ); %sprintf('%s%s_zproj.tif', expt.dir, expt.name); % exist(movPath, 'file')
    if ~exist(movPath, 'file')
        ZtifDir = sprintf('%sZtifs\\',expt.dir);
        affMov = WriteSbxPlaneTif(expt.sbx, catInfo, 1, 'dir',ZtifDir, 'name',expt.name, 'type','aff', 'binT',10, 'verbose',true, 'chan',1 );
        %error('%s does not exist!', movPath);  
    else
        fprintf('\nLoading  %s...', movPath); tic;
        affMov = pipe.io.read_tiff(movPath); toc
    end

    w = waitbar(0, 'Writing ROI projection tifs');
    if isempty(axon)
        affMean = mean(affMov, 3);
        %affStd = std(double(affMov), 0, 3);
        affMax = max(affMov, [], 3);
        w = waitbar(0, 'Writing ROI projection tifs');
        for r = rSet%1:expt.Nroi
            ROI(r).projPath = sprintf('%s%s_ROI_%02i.tif', projDir, expt.name, r );
            tic;
            if ~exist(ROI(r).projPath, 'file') || overwrite
                % Perform z-projection if it doesn't already exist (taking advantage of the fact that ROI are sorted by zCent)
                fprintf('\nROI %i:  ', r); % ,  expt.sbx
                % Determine range for cropping and z-projection
                ROI(r).crop(1) = max(floor(ROI(r).box_x(1))-edgeBuffer(1), segEdges(1)+1); % left edge
                ROI(r).crop(2) = min(ceil(ROI(r).box_x(2))+edgeBuffer(2), expt.Ncol-segEdges(2)); % right edge
                ROI(r).crop(3) = max(floor(ROI(r).box_y(1))-edgeBuffer(3), segEdges(3)+1); % top edge
                ROI(r).crop(4) = min(ceil(ROI(r).box_y(2))+edgeBuffer(4), expt.Nrow-segEdges(4)); % bottom edge
                crop_proj = affMov(ROI(r).crop(3):ROI(r).crop(4), ROI(r).crop(1):ROI(r).crop(2),:); % Crop projection
                % Write tiff
                if ~isempty(crop_proj)      
                    ROI(r).projPath = sprintf('%s%s_ROI_%03i.tif', projDir, expt.name, r );
                    fprintf('\nWriting %s... ', ROI(r).projPath );
                    pipe.io.writeTiff(uint16(crop_proj), ROI(r).projPath); 
                else
                    fprintf('\nWriting %s failed, crop_proj is empty... ', ROI(r).projPath );
                end   
            else
                fprintf('\nLoading %s...  ', ROI(r).projPath );
                crop_proj = double(loadtiff(ROI(r).projPath)); 
            end

            % Save mean, std dev and max projections
            % {
            ROI(r).cropMean = affMean(ROI(r).crop(3):ROI(r).crop(4), ROI(r).crop(1):ROI(r).crop(2)); % uint16(mean(crop_proj, 3));
            %ROI(r).cropStd = affStd(ROI(r).crop(3):ROI(r).crop(4), ROI(r).crop(1):ROI(r).crop(2)); %uint16(std(double(crop_proj), 0, 3));
            ROI(r).cropMax = affMax(ROI(r).crop(3):ROI(r).crop(4), ROI(r).crop(1):ROI(r).crop(2)); %uint16(max(crop_proj, [], 3));
            %}
            if nargout > 1,  crop_proj_out{r} = crop_proj;  end
            waitbar(r/expt.Nroi, w)
            toc
        end
    else
        for a = 1:Naxon
            axon(a).projPath = sprintf('%s%s_axon%02i.tif', projDir, expt.name, a );
            if ~exist(axon(a).projPath, 'file') || overwrite
                %axonROIcent = vertcat( ROI(axon(a).ROI).cent);              
                tempCrop(1) = max(floor(axon(a).crop(1))-edgeBuffer(1), segEdges(1)+1); % left edge
                tempCrop(2) = min(ceil(axon(a).crop(2))+edgeBuffer(2), expt.Ncol-segEdges(2)); % right edge
                tempCrop(3) = max(floor(axon(a).crop(3))-edgeBuffer(3), segEdges(3)+1); % top edge
                tempCrop(4) = min(ceil(axon(a).crop(4))+edgeBuffer(4), expt.Nrow-segEdges(4)); % bottom edge
                crop_proj = affMov(tempCrop(3):tempCrop(4), tempCrop(1):tempCrop(2),:); % Crop projection
                % Write tiff
                if ~isempty(crop_proj)      
                    fprintf('\nWriting %s... ', axon(a).projPath );
                    pipe.io.writeTiff(uint16(crop_proj), axon(a).projPath); 
                else
                    fprintf('\nWriting %s failed, crop_proj is empty... ', axonProjPath );
                end   
                waitbar(a/Naxon, w)
                toc 
            elseif nargout > 0
                fprintf('\nLoading %s...  ', axon(a).projPath );
                crop_proj = double(loadtiff(axon(a).projPath)); 
            end
            % Save mean, std dev and max projections
            if nargout > 1,  crop_proj_out{a} = crop_proj;  end
            if nargout > 2
                axon(a).cropMean = uint16(mean(crop_proj, 3));
                %axon(a).cropStd = uint16(std(crop_proj, 0, 3));
                axon(a).cropMax = uint16(max(crop_proj, [], 3));
            end
        end
    end
end
delete(w);
if nargout > 1,  varargout{1} = crop_proj_out;  end
end