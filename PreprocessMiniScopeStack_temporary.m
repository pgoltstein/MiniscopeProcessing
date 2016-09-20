function [ImData, RawImData] = PreprocessMiniScopeStack_temporary
% Pre-processes raw miniscope data
    fprintf('\nPre-processing data in: %s\n', pwd);
    
    %% Pre-processing settings
    ProcessingSettings.AuxSampFreq = 1000; % Hz
    ProcessingSettings.DoRegistration = false; % True / False
    ProcessingSettings.ImRegTemplateLength = 10; % seconds
    
    %% Load Experimental data

    % Load aux recorder data
    AuxFile = dir(['.' filesep '*.lvd']);
    fprintf('\nLoading Aux-data: %s ... ', AuxFile(1).name);
    AuxData = load_lvd( AuxFile(1).name );
    fprintf('done\n');
    
    % Find frame times
    FrameCh = AuxData(10,:);
    FrameOnsets = find( diff(FrameCh>2)>0.1 );
    InterFrameTime = mean(FrameOnsets(2:end-2)-FrameOnsets(1:end-3));
    fprintf('Inter frame time = %5.1f ms\n', InterFrameTime);
    ProcessingSettings.SamplingFreq = ProcessingSettings.AuxSampFreq/InterFrameTime;
    ProcessingSettings.FrameOnsets = FrameOnsets;
    fprintf('Estimated frame rate = %4.1f Hz\n', ProcessingSettings.SamplingFreq);
    FrameTimes = FrameOnsets + ( 0.5 * InterFrameTime );
    fprintf('Nr of frames detected in aux-data: %d\n', length(FrameTimes));
    save('ProcessingSettings.mat','ProcessingSettings');
    
    % Load stimulus data if stimulusfile is present
    StimFile = dir([pwd filesep 'Stimuli*.mat']);
    if ~isempty(StimFile)
        fprintf('\nLoading stimulus data: %s ... ', StimFile(1).name);
        StimData = load( StimFile(1).name );
        X = LOCAL_LoadStimulusFrames( AuxData, StimData, ProcessingSettings, FrameOnsets );
        fprintf('done\n');
    else
        StimData = [];
        fprintf('\nNo stimulus data found.\n');
    end
    
    % Load raw imaging stacks
    ImFile = dir(['.' filesep '*.tif']);
    fprintf('\nLoading imaging stack: %s\n', ImFile(1).name);
    warning('off'); % Tiff class gives stupid warnings because of missing fields
    
    TiffInfo = Tiff(ImFile(1).name);
    ImInfo = imfinfo(ImFile(1).name);
    nFrames = length(ImInfo);
    xRes = ImInfo(1).Width;
    yRes = ImInfo(1).Height;
    fprintf('Image dimensions: x=%d, y=%d, t=%d\n', xRes, yRes, nFrames);
    ImData = zeros( yRes, xRes, nFrames, 'uint16' );
    fprintf('Loading frame: %6d',0);
    for f = 1:nFrames
        fprintf('\b\b\b\b\b\b%6d',f);
        TiffInfo.setDirectory(f);
        ImData(:,:,f) = TiffInfo.read();
    end
    fprintf(' ... done\n');
    warning('on'); % Aaaand warnings back on..
    RawImData = ImData;
    
    I = LOCAL_AutoScale( mean(ImData,3), 'uint16' );
    imwrite(I,'RawAverage.tiff','tiff');

    %% Perform image registration
    if ProcessingSettings.DoRegistration
        fprintf('\nPerforming image registration\n');

        % Calulate template to register to
        fprintf('Making template\n');
        TemplateData = ImData( :, :, 1:round(ProcessingSettings.ImRegTemplateLength*ProcessingSettings.SamplingFreq) );
        RegOutput = LOCAL_ImageRegistration( TemplateData, mean(TemplateData,3) );
        Template = mean( LOCAL_shift_data( TemplateData, -1*RegOutput.y, -1*RegOutput.x ), 3 );

        % Get registration parameters
        fprintf('Calculating registration parameters\n');
        RegOutput = LOCAL_ImageRegistration( ImData, Template );

        % Register data from channel 2
        fprintf('Shifting frames ... ');
        ImData = LOCAL_shift_data( ImData, -1*RegOutput.y, -1*RegOutput.x );
        save('ImageRegistrationParameters.mat','RegOutput');    
        fprintf('done\n');

        I = LOCAL_AutoScale( mean(ImData,3), 'uint16' );
        imwrite(I,'RegisteredAverage.tiff','tiff');
    else
        fprintf('\nSkipping image registration\n');
    end
    
    %% Calculate response maps
    if ~isempty(StimData)
        fprintf('\nCalculating response maps\n');
        nStimuli = length(StimData.StimSettings.Stimulus);
        LOCAL_ResponseMap( ImData, 1:nStimuli, 'AllStimuli', 0, ...
            ProcessingSettings.SamplingFreq, [1 2 4 8], X, StimData );
        
        % Retinotopy
        if isfield(StimData.StimSettings,'NumXpatches')
            xPatch = []; yPatch = [];
            for s = 1:nStimuli
                xPatch(s) = StimData.StimSettings.Stimulus(s).X;
                yPatch(s) = StimData.StimSettings.Stimulus(s).Y;
            end
            LOCAL_ResponseMap( ImData, xPatch, 'Azimuth', 0, ...
                ProcessingSettings.SamplingFreq, [1 2 4 8], X, StimData );
            LOCAL_ResponseMap( ImData, yPatch, 'Elevation', 0, ...
                ProcessingSettings.SamplingFreq, [1 2 4 8], X, StimData );
        end
    end
    
    fprintf('\nFinished.\n');
end

function I = LOCAL_AutoScale( I, Type )
    
    % assume double
    I = double(I);
    I = (I-min(I(:))) ./ (max(I(:))-min(I(:)));

    if strcmpi(Type,'uint8')
        I = uint8( I .* (2^8) );
    elseif strcmpi(Type,'int8')
        I = int8( (I.*(2.^8)) - (2.^7) );
    elseif strcmpi(Type,'uint16')
        I = uint16( I .* (2^16) );
    elseif strcmpi(Type,'int16')
        I = int16( (I.*(2.^16)) - (2.^15) );
    end
    
end

function data = LOCAL_shift_data( data, shift_x, shift_y )
    median_data=median(double(data(1,1,:)));
    dimensions = size(data);
    for ind=1:size(data,3)
        tmp_im=uint16(median_data)*ones(size(data,1),size(data,2),'uint16');
        tmp_im(max(1,-round(shift_x(ind))+1):min(dimensions(1),dimensions(1)-round(shift_x(ind))),max(1,-round(shift_y(ind))+1):min(dimensions(2),dimensions(2)-round(shift_y(ind)))) = ...
            data(max(1,round(shift_x(ind))+1):min(dimensions(1)+round(shift_x(ind)),dimensions(1)),max(1,round(shift_y(ind))+1):min(dimensions(2)+round(shift_y(ind)),dimensions(2)),ind);
        data(:,:,ind)=tmp_im;
    end
end

function C = LOCAL_ColorMap( N, Name )
    figure;
    CM = colormap( Name );
    close;
    IX = round(linspace(1,64,N));
    for i = 1:length(IX)
        C{i} = CM(IX(i),:);
    end
end

function RegOutput = LOCAL_ImageRegistration( data, TMPL )
    RegOutput = struct;
    
    % Get number of frames in stack
    nFrames = size(data,3);
    
    % Get fft of template
%     TMPLfft = fft2( TMPL );
    TMPL = TMPL ./ max(TMPL(:));
    AverageFilter = fspecial('average',31);
    TMPLfft = fft2( imfilter( TMPL, AverageFilter, 'replicate' ) );

    % Loop stack and realign data
    DT = clock;
    tic;
    fprintf('Registering %d frames\n',nFrames);
    fprintf('Time started: %02.0f:%02.0f:%02.0f \n', DT(4), DT(5), DT(6) );
    fprintf('Registering frame: %d',0);
    for f = 1:nFrames
        fprintf('\b\b\b\b\b\b%6d',f);

        % Register image
%         Output = dftregistration( TMPLfft, ...
%             fft2( double(data(:,:,f)) ), 1 );
        ImFr = imfilter( double(data(:,:,f)), AverageFilter, 'replicate' );
        ImFr = ImFr ./ max(ImFr(:));
        Output = dftregistration( TMPLfft, ...
            fft2( ImFr ), 1 );
        
        % Get realignment parameters
        RegOutput.Performance(f) = Output(1);
        RegOutput.PhaseCorr(f) = Output(2);
        RegOutput.y(f) = Output(3);
        RegOutput.x(f) = Output(4);
        
    end
    fprintf('\n');
    DT = clock;
    fprintf( 'Time finished: %02.0f:%02.0f:%02.0f (%1.0f seconds)\n', DT(4), DT(5), DT(6), toc );
    
end

function X = LOCAL_LoadStimulusFrames( AuxData, StimData, ProcessingSettings, FrameOnsets )
    
    % Find framenumber of stimuli
    StimCh = AuxData(4,:);
    StimOn = StimCh>0.85;
    StimOnsets = find(diff(StimOn==1)>0);
    StimOffsets = find(diff(StimOn==1)<0);
    
    % Get stimchannel values
    for Nr = 1:length(StimOnsets)
        StimChannel(Nr) = mean( StimCh(StimOnsets(Nr):StimOffsets(Nr)) );
    end
    
    % Find stimulus id's
    AllIds = cell2mat(StimData.StimSettings.RandStims);
    for Nr = 1:length(StimOnsets)
        StimId(Nr) = AllIds(Nr);
    end
    
    % Find frames corresponding to stimulus presentations
    for t = 1:length(StimOnsets)
        [~,StimOnsetFrames(t)] = min( abs(FrameOnsets-StimOnsets(t)) );
        [~,StimOffsetFrames(t)] = min( abs(FrameOnsets-StimOffsets(t)) );
    end
        
    % Finds period of no-stimulus
    NoStimOn = StimCh>0.4 & StimCh<0.6;
    NoStimOnsets = find(diff(NoStimOn==1)>0);
    NoStimOffsets = find(diff(NoStimOn==1)<0);
    if length(NoStimOffsets) < length(NoStimOnsets)
        NoStimOffsets = [NoStimOffsets length(NoStimOn)];
    end
    
    % Remove periods that are too short to be useful (<1 s)
    NoStimLength = NoStimOffsets-NoStimOnsets;
    NoStimOnsets( NoStimLength < (ProcessingSettings.AuxSampFreq*1) ) = [];
    NoStimOffsets( NoStimLength < (ProcessingSettings.AuxSampFreq*1) ) = [];
    
    % Find frames corresponding to no-stimulus periods
    NoStimOnsetFrames = [];
    NoStimOffsetFrames = [];
    for t = 1:length(NoStimOnsets)
        [~,NoStimOnsetFrames(t)] = min( abs(FrameOnsets-NoStimOnsets(t)) );
        [~,NoStimOffsetFrames(t)] = min( abs(FrameOnsets-NoStimOffsets(t)) );
    end
    
    % Output struct
    X.onFrames = StimOnsetFrames;
    X.offFrames = StimOffsetFrames;
    X.onFramesNoStim = NoStimOnsetFrames;
    X.offFramesNoStim = NoStimOffsetFrames;
    X.stimChannel = StimChannel;
    X.stimId = StimId;
end

function LOCAL_ResponseMap( Data, StimINDX, WriteName, GroupStimuli, SamplingFreq, Scaling, X, StimData )

    % Settings and filters
    AvgFilter = fspecial('disk', 1 );
    
    % Load average images
    if exist('RegisteredAverage.tiff','file')
        Iavg = imread('RegisteredAverage.tiff');
    else
        Iavg = imread('RawAverage.tiff');
    end
    Iavg = im2double(Iavg);
    Iavg = imadjust(Iavg);

    % Get stimulus information
    StimDur = StimData.StimSettings.StimDur;	
    ITIlength = StimData.StimSettings.ITIlength;
    
    % Re-organize stimuli according to the list provided in StimINDX
    SON = cell(1,max(StimINDX));
    SOF = cell(1,max(StimINDX));
    StimLst = cell(1,length(StimINDX));
    for s = 1:length(StimINDX)
        if StimINDX(s) > 0
            SON{StimINDX(s)} = [ SON{StimINDX(s)} X.onFrames( X.stimId == s ) ];
            SOF{StimINDX(s)} = [ SOF{StimINDX(s)} X.offFrames( X.stimId == s ) ];
            StimLst{StimINDX(s)} = [StimLst{StimINDX(s)} s];
        end
    end

    % Number of stimuli
    NumStims = length(SON);
%     for s = 1:NumStims
%         disp(['Stim ' num2str(s) ': ' num2str(StimLst{s})]);
%     end
    
    % Get baseline map
    BSmap = zeros( size(Data,1), size(Data,2), length(SON{NumStims})*NumStims);
    cnt = 0;
    for s = 1:NumStims
        for t = 1:length(SON{s})
            Start = SON{s}(t)-round(0.7*ITIlength*SamplingFreq);
            Stop = SON{s}(t);
            Ibs = sum(Data(:,:,Start:Stop),3)./length(Start:Stop);
            cnt = cnt + 1;
            BSmap(:,:,cnt) = imfilter(Ibs, AvgFilter);
        end
    end
    BSmap = mean(BSmap,3);
        
    % Get response maps
    dFoFmaps = cell(1,NumStims);
    for s = 1:NumStims
        dFoFmaps{s} = zeros( size(Data,1), size(Data,2), length(SON{s}) );
        for t = 1:length(SON{s})
            Start = SON{s}(t)+1;
            Stop = SOF{s}(t)+round((0.2*StimDur*SamplingFreq));
            if Stop <= size(Data,3)
                I = sum(Data(:,:,Start:Stop),3)./length(Start:Stop);
                I = imfilter(I, AvgFilter);
                dFoFmaps{s}(:,:,t) = (I-BSmap)./BSmap;
            else
                disp(['Missing data: stim=' num2str(s) ', trial=' num2str(t)]);
                dFoFmaps{s}(:,:,t) = NaN;
            end
        end
    end

    % Get average map per stimulus
    if GroupStimuli > 1
        NumStims = round(NumStims/GroupStimuli);
        dFoFavgMap = zeros( size(Data,1), size(Data,2), NumStims );
        for s = 1:NumStims
            GrId = [];
            I = zeros( size(Data,1), size(Data,2), GroupStimuli );
            for r = 1:GroupStimuli
                GrId(r) = s+((r-1)*NumStims);
                I(:,:,r) = nanmean(dFoFmaps{s+((r-1)*NumStims)},3);
            end
            dFoFavgMap(:,:,s) = mean(I,3);
%             disp(['Stim-gr ' num2str(s) ': ' num2str(GrId)]);
        end
    elseif GroupStimuli < -1
        NumStims = abs(round(NumStims/GroupStimuli));
        dFoFavgMap = zeros( size(Data,1), size(Data,2), NumStims );
        for s = 1:NumStims
            GrId = [];
            I = zeros( size(Data,1), size(Data,2), abs(GroupStimuli) );
            for r = 1:abs(GroupStimuli)
                GrId(r) = r+((s-1)*abs(GroupStimuli));
                I(:,:,r) = nanmean(dFoFmaps{r+((s-1)*abs(GroupStimuli))},3);
            end
            dFoFavgMap(:,:,s) = mean(I,3);
%             disp(['Stim-gr ' num2str(s) ': ' num2str(GrId)]);
        end
    else
        dFoFavgMap = zeros( size(Data,1), size(Data,2), NumStims );
        for s = 1:NumStims
            dFoFavgMap(:,:,s) = nanmean(dFoFmaps{s},3);
        end
    end
    StimulusAngles = ((1:NumStims)/NumStims)*360;
    
    % Calculate max stimulus, amplitude and tuning width (circular variance
    % of sorted map)
    [ResponseAmplitude,PrefStim] = max( dFoFavgMap, [], 3 );
    PrefStim(PrefStim<1) = 1; PrefStim(PrefStim>360) = 360;

    SortedResponseMap = sort(dFoFavgMap,3,'descend');
    SortedResponseMap(SortedResponseMap<0) = 0;
    Phasor = zeros(size(dFoFavgMap,1),size(dFoFavgMap,2),NumStims);
    for s = 1:NumStims
        Phasor(:,:,s) = SortedResponseMap(:,:,s) * exp( 1i * (StimulusAngles(s)/360) *2 *pi );
    end
    ResLen = abs( sum(Phasor,3) ) ./ sum(abs(Phasor),3);
    CircVar = 1-ResLen;
    
    % Scaling
    fprintf('HLSmap: %s (Ampl=', WriteName);
    for Sc = 1:length(Scaling)
        ResponseAmplitudeX = ResponseAmplitude * Scaling(Sc);
        ScalingStr = ['x' num2str(Scaling(Sc))];

        % Create HLS map
        if NumStims > 2
            C = LOCAL_ColorMap( NumStims, 'jet' );
        else
            C{1} = [0 0.8 1];
            C{2} = [1 0 0.4];
        end
        HLSmap = zeros( size(dFoFavgMap,1), size(dFoFavgMap,2), 3 );
        HLS_I_map = zeros( size(dFoFavgMap,1), size(dFoFavgMap,2), 3 );
        for y = 1:size(dFoFavgMap,1)
            for x = 1:size(dFoFavgMap,2)
                for c = 1:3
                    if ~isnan(PrefStim(y,x)) && ResponseAmplitudeX(y,x) > 0
                        HLSmap(y,x,c) = (C{PrefStim(y,x)}(c) .* ResponseAmplitudeX(y,x) .* ResLen(y,x)) ...
                                                           + (ResponseAmplitudeX(y,x).*CircVar(y,x));
                        HLS_I_map(y,x,c) = HLSmap(y,x,c) .* Iavg(y,x);
                    end
                end
            end
        end

        % Clip maps
        HLSmap(HLSmap<0) = 0;
        HLSmap(HLSmap>1) = 1;
        HLS_I_map(HLS_I_map<0) = 0;
        HLS_I_map(HLS_I_map>1) = 1;
        ResponseAmplitudeX(ResponseAmplitudeX<0) = 0;
        ResponseAmplitudeX(ResponseAmplitudeX>1) = 1;

        % Add color index
        [yRes,xRes,~] = size(HLS_I_map);
        xRes = xRes - 100;
        for s = 1:NumStims
            xRange = round(50+(((s*(xRes/NumStims)) - (xRes/(NumStims*2))):(s*(xRes/NumStims))));
            for c = 1:3
                HLSmap(yRes:yRes+10,xRange,c) = C{s}(c);
                HLS_I_map(yRes:yRes+10,xRange,c) = C{s}(c);
            end
        end

        % Display maps and write to tiff
        figure;
        imshow(HLSmap);
        title(['HLSmap ' WriteName ' (Ampl ' ScalingStr ')']);

        figure;
        imshow(ResponseAmplitudeX);
        title(['dFoFmap ' WriteName ' (Ampl ' ScalingStr ')']);

        figure;
        imshow(HLS_I_map);
        title(['HLSmap_Scaled_ ' WriteName ' (Ampl ' ScalingStr ')']);

        fprintf('%s',ScalingStr);
        if Sc<length(Scaling); fprintf(','); else fprintf(')\n'); end

        % Save data
        if ~exist('.\MAPS','dir')
            mkdir('MAPS');
        end
        imwrite(HLSmap,['.' filesep 'MAPS' filesep 'HLSmap_' WriteName '_' ScalingStr '.tiff'],'tiff');
        imwrite(HLS_I_map,['.' filesep 'MAPS' filesep 'HLSmapSc_' WriteName '_' ScalingStr '.tiff'],'tiff');
        imwrite(ResponseAmplitudeX,['.' filesep 'MAPS' filesep 'dFoFmap_' WriteName '_' ScalingStr '.tiff'],'tiff');    
    end
    
    save(['.' filesep 'MAPS' filesep 'AvgMaps_' WriteName '.mat'], 'dFoFavgMap');
end



