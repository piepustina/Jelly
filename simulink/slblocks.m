function blkStruct = slblocks
% This function specifies that the library 'SoftRobotics'
% should appear in the Library Browser with the 
% name 'Jelly'

    Browser.Library = 'SoftRobotics';
    % Browser.Library is the name of the library

    Browser.Name = 'Jelly';
    % Browser.Name is the library name that appears
    % in the Library Browser

    blkStruct.Browser = Browser;