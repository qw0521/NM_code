This code was written by Zhili Hu and adjusted by Wen Qin, both works in NUAA as in 2024.
Email: zhili.hu@nuaa.edu.cn

It contains two runnable MATLAB/Octave files, "droplet_climbs_the_film.m" and "Droplets_form_films.m".

File droplet_climbs_the_film.m is based on the phase-field model by Borcia, R. et al. and the original paper is, Static and dynamic contact angles – A phase field modelling. Eur. Phys. J. Spec. Top. 166, 127-131 (2009).
We use it to simulate the equilibrium shape profile of a nanodroplet that sits at an atomic step. The droplet can either sit still at the step or climb up the step, depending on its relative size and wettabilities. The default parameter set is for a Pt droplet that sits at the edge of a PtSex monolayer, on a MoS2 substrate. Users could modify the code for other scenarios.

File Droplets_form_films.m is based on our previous code written for paper He, Y. et al. Nat. Commun. 11, 57 (2020), and the phase field model in it is based on Karma, Phys. Rev. Lett. 81, 4444-4447 (1998).
We used it to simulate the morphological evolution in our experiments, where Pt nanodroplets drive the formation of earthworm like PtSex nanoribbons, which further stitch into complete films. Users could modify this code for their own purposes.

Details of the models can be found in the Methods section in our paper.

Other information:
1. Systems requirements
This code runs in either MATLAB or Octave. We have only tested it in Windows OS, but it should also work under Mac and Linux, as long as MATLAB/Octave is installed.
The code currently has only one (alpha) version. We have tested it in Octave version 9.1.0 in Windows 11, and MATLAB version 9.11.0.1809720 (R2021b) Update 1.
It doesn't require any non-standard hardware.

2. Installation
Users need to make sure a MATLAB/Octave is installed in advance.
After that, users simply need to unpackage the code with an omittable time cost and run it inside MATLAB/Octave.

3. Demo
The intact files themselves contains a well adjusted parameter set. By running it, users are expected to obtain figures similar to Fig.3 in our work.
It would take around 0.5~1 hours to finish the test on a "normal" desktop computer.

4. Instruction
The intact code already contains well adjusted data for running. Users just need to open "droplet_climbs_the_film.m" or "Droplets_form_films.m", click the "run" button and select the local folder in a MATLAB/Octave to run it, to obtain figures similar to Fig.3 in our work.
For users who would modify the parameter set for other purposes, they could pay attention to lines marked with "USER-ATTENTION". As for users who would modify the code structure for other purposes, they may need to read the comments in the code as well as the model of the code in the Methods section in our work.
