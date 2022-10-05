function InterpTemporalCat(sbxInfo, varargin) 

IP = inputParser;
addRequired( IP, 'sbxInfo', @isstruct )
addParameter( IP, 'pmt', -1, @isnumeric )
addParameter( IP, 'edges', [0,0,0,0], @isnumeric )
addParameter( IP, 'chunkSize',15);
parse( IP, sbxInfo, varargin{:} );
pmt = IP.Results.pmt;
edges = IP.Results.edges;
chunkSize = IP.Results.chunkSize;

Nx = sbxInfo.sz(1); 
Ny = sbxInfo.sz(2);

% set up paths
affPath = sprintf('%s%s.sbxreg', sbxInfo.dir, sbxInfo.exptName );  % sprintf('%s%s.sbx_affine', sbxInfo.dir, sbxInfo.exptName ); 
interpPath = sprintf('%s%s.sbx_interp', sbxInfo.dir, sbxInfo.exptName );

%crop edges of video
Xrange = edges(3)+1:Nx-edges(4); 
Yrange = edges(1)+1:Ny-edges(2); 
%make the vectors for interpolation (same for each plane)
vx = 1:numel(Xrange);
vy = 1:numel(Yrange);

[chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, sbxInfo.totScan-1, sbxInfo.totScan-1, 'size',chunkSize); % prevent going beyond totScan
w = waitbar(0,'Starting interpolation...'); %start the waitbar clock
tic
rw = SbxWriter(interpPath, sbxInfo, '.sbx_interp', true); % rw = pipe.io.RegWriter(interpPath, sbxInfo, '.sbx_interp', true);
for chunk = 1:Nchunk % 
    % Load and crop data chunk to be interpolated
    chunk_data = double(readSBX(affPath, sbxInfo, chunkLims(chunk,1), chunkLength(chunk)+1, pmt, []));
    chunk_crop = chunk_data(:,Xrange,Yrange,:,:);
    vt = 1:chunkLength(chunk)+1;
    % Interpolate data from each channel and plane
    if chunk ~= Nchunk
        chunk_interp = zeros(sbxInfo.nchan, Nx, Ny, sbxInfo.Nplane, chunkLength(chunk));
    else
        chunk_interp = zeros(sbxInfo.nchan, Nx, Ny, sbxInfo.Nplane, chunkLength(chunk)+1);
    end
    for chan = 1:sbxInfo.nchan
        for z = 1:sbxInfo.Nplane % parfor is slower
            %{
            %load the whole chunk. Have to load one extra volume to account for the interpolation shift
            %M = double(readSBX(affPath, sbxInfo, chunkLims(chunk,1), chunkLength(chunk)+1, pmt, z));
            %M = double(pipe.io.sbxRead(affPath, chunk_idx(c), ChunkLength(c)+1, pmt, z));

            %crop it based on the same range above
            M = M(Xrange,Yrange, :);
            if size(M,3) == chunkLength(chunk)
                M(:,:,chunkLength(chunk)+1) = zeros(size(M,1),size(M,2));
            end
            %}
            M = squeeze(chunk_crop(chan, :,:,z,:));
            % create the gridded interpolant based on the registered movie and the interpolation vectors
            F = griddedInterpolant({vx,vy,vt}, M);

            % calculate the partial shift for each z-plane. Earlier planes are shifted forward more than later ones
            s = (sbxInfo.Nplane-z)/sbxInfo.Nplane;
             
            qt = linspace(1+s, chunkLength(chunk)+1+s, chunkLength(chunk)+1); % %make vector of shifted timepoints to evaluate the interpolant at

            %evaluate the interpolant at the same XY grid points, but shifted timepoints (qt).
            M_s = F({vx,vy,qt});
            if chunk ~= Nchunk
                chunk_interp(chan,Xrange,Yrange,z,:) = M_s(:,:,1:chunkLength(chunk));
            else
                chunk_interp(chan,Xrange,Yrange,z,:) = M_s(:,:,1:chunkLength(chunk)+1);
            end
        end 
    end
    %write to the binary file
    chunk_interp = reshape(chunk_interp, sbxInfo.nchan, Nx, Ny,[]);
    rw.write(uint16(chunk_interp));
    
    %estimate how much time is left and update the waitbar.
    t_elapse = minutes(seconds(toc));
    t_remain = minutes(seconds(toc*((Nchunk/chunk)-1)));
    msg = sprintf('Interpolating... %.f minutes elapsed \n about %.f minutes remaining',t_elapse,t_remain);
    waitbar(chunk/Nchunk,w,msg);
end
rw.delete;
close(w);
%clear M F M_i M_reg chunk_interp;
end