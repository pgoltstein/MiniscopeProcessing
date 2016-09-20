function [AuxData, AuxScanRate, NumChannels, TimeStamp, InputRange, Status] = LoadLvdFile(FileName)

    % Output variables
    Status = 0;
    AuxData = [];
    AuxScanRate = [];

    % try to open file, otherwise quit
    [FileIndex, Message] = fopen(FileName, 'r', 'ieee-be');

    if FileIndex == -1
        disp('There was a problem reading the following file:');
        disp(FileName);
        disp(Message);
        Status = -1;
        return;
    end

    % Get header information

    AuxScanRate = double(fread(FileIndex, 1, 'double'));
    NumChannels = fread(FileIndex, 1, 'double');
    TimeStamp = fread(FileIndex, 1, 'double');
    InputRange = fread(FileIndex, 1, 'double');

    % Read auxdata
    ReadSamples = [NumChannels inf];
    AuxData = fread(FileIndex, ReadSamples, 'double');

    % Close file
    fclose(FileIndex);

end
