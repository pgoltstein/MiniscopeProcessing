function FileSize = GetFileSize(FileIndex)
   
    % go to the end of the file
    fseek(FileIndex, 0, 'eof');
    FileSize = ftell(FileIndex);
    
end

