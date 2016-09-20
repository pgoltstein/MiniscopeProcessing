function I = AutoScaleImage( I, Type )
    
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
