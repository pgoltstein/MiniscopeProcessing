function ImData = PreprocessMiniScopeStack_temporary
% Pre-processes raw miniscope data
    fprintf('\nPre-processing data in: %s\n', pwd);
    
    %% Pre-processing settings
    ProcessingSettings.DoRegistration = true; % True / False
    ProcessingSettings.ImRegTemplateLength = 10; % 10 seconds
    ProcessingSettings.ImReg_BG_FoF = 31; % 31 px; Background subtraction for image registration
    ProcessingSettings.ImReg_FG_FoF = 5; % 5 px; Foreground blur for image registration
    
    ProcessingSettings.DoBackgroundSubtraction = true; % True / False
    ProcessingSettings.BackgroundSubtractionFoF = 17; % 17 px; Real background subtraction
    
    %% Load Experimental data

    % Load aux recorder data
    AuxFile = dir(['.' filesep '*.lvd']);
    fprintf('\nLoading Aux-data: %s ... ', AuxFile(1).name);
    [AuxData, ProcessingSettings.AuxSampFreq] = LoadLvdFile( AuxFile(1).name );
    fprintf('done\n');
    fprintf('Aux sampling frequency = %5.1f ms\n', ProcessingSettings.AuxSampFreq);
        
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
        X = LoadStimulusFrames( AuxData, StimData, ProcessingSettings, FrameOnsets );
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
    I = AutoScaleImage( mean(ImData,3), 'uint16' );
    imwrite(I,'RawAverage.tiff','tiff');

    %% Perform image registration
    if ProcessingSettings.DoRegistration
        fprintf('\nPerforming image registration\n');

        % Calulate template to register to
        fprintf('Making template\n');
        TemplateData = ImData( :, :, 1:round(ProcessingSettings.ImRegTemplateLength*ProcessingSettings.SamplingFreq) );
        RegOutput = ImageRegistration( TemplateData, mean(TemplateData,3), ProcessingSettings.ImReg_BG_FoF, ProcessingSettings.ImReg_FG_FoF );
        Template = mean( ShiftImagingData( TemplateData, -1*RegOutput.y, -1*RegOutput.x ), 3 );

        % Get registration parameters
        fprintf('Calculating registration parameters\n');
        RegOutput = ImageRegistration( ImData, Template, ProcessingSettings.ImReg_BG_FoF, ProcessingSettings.ImReg_FG_FoF );

        % Register data from channel 2
        fprintf('Shifting frames ... ');
        ImData = ShiftImagingData( ImData, -1*RegOutput.y, -1*RegOutput.x );
        save('ImageRegistrationParameters.mat','RegOutput','Template');    
        fprintf('done\n');

        I = AutoScaleImage( mean(ImData,3), 'uint16' );
        imwrite(I,'RegisteredAverage.tiff','tiff');
    else
        fprintf('\nSkipping image registration\n');
    end
    
    %% Perform background subtraction
    if ProcessingSettings.DoBackgroundSubtraction
        fprintf('Background subtraction (FoF=%dpx); frame %6d',ProcessingSettings.BackgroundSubtractionFoF,0);
        BGfilter = fspecial('average',ProcessingSettings.BackgroundSubtractionFoF);
        for f = 1:nFrames
            fprintf('\b\b\b\b\b\b%6d',f);
            BG = imfilter( ImData(:,:,f), BGfilter, 'replicate' );
            ImData(:,:,f) = ImData(:,:,f)-BG;
        end
        fprintf(' ... done\n');

        I = AutoScaleImage( mean(ImData,3), 'uint16' );
        imwrite(I,'BackgroundSubtractedAverage.tiff','tiff');
    end    
    
    %% Calculate response maps
    if ~isempty(StimData)
        fprintf('\nCalculating response maps\n');
        nStimuli = length(StimData.StimSettings.Stimulus);

        % Retinotopy
        if isfield(StimData.StimSettings,'NumXpatches')
            xPatch = []; yPatch = [];
            for s = 1:nStimuli
                xPatch(s) = StimData.StimSettings.Stimulus(s).X;
                yPatch(s) = StimData.StimSettings.Stimulus(s).Y;
            end
            MakeResponseMap( ImData, xPatch, 'Azimuth', 0, ...
                ProcessingSettings.SamplingFreq, [1 2 4 8], X, StimData );
            MakeResponseMap( ImData, yPatch, 'Elevation', 0, ...
                ProcessingSettings.SamplingFreq, [1 2 4 8], X, StimData );
        
        % Gratings
        elseif isfield(StimData.StimSettings,'Angles')
            MakeResponseMap( ImData, 1:nStimuli, 'Directions', 0, ...
                ProcessingSettings.SamplingFreq, [1 2], X, StimData );
            MakeResponseMap( ImData, 1:nStimuli, 'Orientations', 2, ...
                ProcessingSettings.SamplingFreq, [1 2], X, StimData );
        
        % Just all stimuli
        else
            MakeResponseMap( ImData, 1:nStimuli, 'AllStimuli', 0, ...
                ProcessingSettings.SamplingFreq, [1 2], X, StimData );
        end
        
    end
    
    fprintf('\nFinished.\n');
end


