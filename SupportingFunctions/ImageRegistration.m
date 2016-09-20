function RegOutput = ImageRegistration( data, TMPL )
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