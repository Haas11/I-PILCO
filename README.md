# I-PILCO

Adaptation of the PILCO algorithm for learning variable impedance contact tasks.

Two scenarios:

1) Stiffness learning along predefind trajectory

2) Stiffness + trajectory learning

Both scenarios require 2 files: a main run file and a settings file.
(A simulation scheme is provided when rollouts are not performed in the real world).


Additionally required:
- MATLAB 

- SIMULINK

- Robotics Toolbox by Peter Corke   (v9.10)
http://petercorke.com/Robotics_Toolbox.html

- PILCO Toolbox by Rasmussen & Deisenroth   (v0.9)
http://mlg.eng.cam.ac.uk/pilco/

- SDE Tools
https://github.com/horchler/SDETools

- 

Optional:
- Double Precision Toolbox

- The Lightspeed MATLAB Toolbox   (v2.7)
http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/

- GP Training w/ Input Noise Library (non-sparse)
http://mlg.eng.cam.ac.uk/?portfolio=andrew-mchutchon
