function blkStruct = slblocks
% This function specifies that the library 'SoftRobotics'
% should appear in the Library Browser with the 
% name 'Soft Robotics Toolbox'

    Browser.Library = 'SoftRobotics';
    % Browser.Library is the name of the library

    Browser.Name = 'Soft Robotics Toolbox';
    % Browser.Name is the library name that appears
    % in the Library Browser

    blkStruct.Browser = Browser;