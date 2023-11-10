##########################################################################
The toolbox must be initialized calling "initToolbox.m".
This should also add automatically the Simulink Blocks in the Library Browser.
If this is not the case, in Simulink, open the Library Browser, right click on the window and press "Refresh library browser".
You should see the "Soft Robotics Toolbox" 

##########################################################################
To define a new body type of body in the chain perform the following steps:
1) Create a new folder in the classes directory witht the name of the body
2) Define a body class that inherits from "classes/Body" and override any method you need
3) Define a joint class that inherits from "classes/Joint" and ovverride any method you need
4) Create a folder named "+model_functions" and add there an implementation of all the methods called in "bodyUpdate" and "jointUpdate" methods of "classes/Body" and "classes/Joint" respectively
   Note that you can ovverride also "bodyUpdate" and "jointUpdate" if you want.
5) Modify TreeFactory to add the new body and joint