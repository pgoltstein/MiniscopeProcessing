function RegOutput = ImageRegistration( data, TMPL, BG_FoF, FG_FoF )

    % Output variable
    RegOutput = struct;
    
    % Prepare filters
    BGfilter = fspecial('average',BG_FoF);
    FGfilter = fspecial('average',FG_FoF);
    
    % Get number of frames in stack
    nFrames = size(data,3);
    
    % Convert to double
    TMPL = double(TMPL);
    
    % Subtract background
    BG = imfilter( TMPL, BGfilter, 'replicate' );
    TMPL = TMPL-BG;
    
    % Blur foreground
    TMPL = imfilter( TMPL, FGfilter, 'replicate' );
    
    % FFT of template
    TMPLfft = fft2( TMPL );

    % Loop stack and realign data
    DT = clock;
    tic;
    fprintf('Registering %d frames (BG_FoF=%d; FG_FoF=%d)\n',nFrames,BG_FoF,FG_FoF);
    fprintf('Time started: %02.0f:%02.0f:%02.0f \n', DT(4), DT(5), DT(6) );
    fprintf('Registering frame: %6d',0);
    for f = 1:nFrames
        fprintf('\b\b\b\b\b\b%6d',f);

        % Convert to double, subtract background, blur foreground
        IMAGE = double(data(:,:,f));
        BG = imfilter( IMAGE, BGfilter, 'replicate' );
        IMAGE = IMAGE-BG;
        IMAGE = imfilter( IMAGE, FGfilter, 'replicate' );

        % Register image
        Output = dftregistration( TMPLfft, fft2( IMAGE ), 1 );
                
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