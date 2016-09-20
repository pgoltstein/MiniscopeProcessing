function C = GetColorMap( N, Name )
    figure;
    CM = colormap( Name );
    close;
    IX = round(linspace(1,64,N));
    for i = 1:length(IX)
        C{i} = CM(IX(i),:);
    end
end
