function MakeResponseMap( Data, StimINDX, WriteName, GroupStimuli, SamplingFreq, Scaling, X, StimData )

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
            C = GetColorMap( NumStims, 'jet' );
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
