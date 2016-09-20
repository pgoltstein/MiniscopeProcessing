function X = LoadStimulusFrames( AuxData, StimData, ProcessingSettings, FrameOnsets )
    
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