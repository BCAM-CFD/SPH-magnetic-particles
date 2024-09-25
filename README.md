# SPH-magnetic-particles
<p>
<img src="https://github.com/BCAM-CFD/SPH-magnetic-particles/blob/main/freq1_6_crit_omega.gif" height="350rm" align="right">

 Code for simulating suspensions with magnetic particles using the SPH  method.




----------------------- DESCRIPTION ------------------------

 Code for simulating suspensions with magnetic particles using the SPH
 method. The  following magnetic  interactions have  been implemented
 
    1. Dipole-dipole  forces  between  particles.   
    2. Dipole-dipole  torques between particles.
    3. Magnetic torque  on individual particles.
    
 In the following  references, the  simplified  model,  which only  considers
 dipole-dipole  forces  between  particles and  assumes  the  magnetic
 moment  is  instantaneously  aligned  with the  magnetic  field,  is
 used:
  - Emanuele  Rossi;   Jose  A   Ruiz-López;  A   Vázquez-Quesada;  M
     Ellero.   Dynamics    and   rheology    of   a    suspension   of
     super-paramagnetic chains  under the  combined effect of  a shear
     flow  and  a  rotating   magnetic  field.  Soft  Matter,  17(24),
     6006-6019. 2021.
  - Adolfo Vázquez-Quesada;  Thomas Franke; Marco Ellero.  Theory and
     simulation of the  dynamics, deformation, and breakup  of a chain
     of   superparamagnetic   beads    under   a   rotating   magnetic
     field. Physics of Fluids. 29, pp. 032006, 2017.

 This code is  based on the original MCF code  developed by Xin Bian.
 The  current version  has  been developed  in collaboration  between
 - Marco Ellero,  leader of the  CFD Modelling and Simulation  group at
   BCAM (Basque Center  for Applied Mathematics) in  Bilbao, Spain.
 - Adolfo Vazquez-Quesada from  the Department of Fundamental Physics
   at UNED, in Madrid, Spain.

 Developers:
     Xin Bian.
     Adolfo Vazquez-Quesada.

 Contact: a.vazquez-quesada@fisfun.uned.es, mellero@bcamath.org
 
--------------------------------------------------------------------

-------------------------- INSTALLATION -------------------------------

After download MCF/mcf/ folder,
I suppose you have already installed PPM library
in a proper path,
which needs to be given when you configure mcf client.
I also suppose you have installed makedepf90,
which is to check routines dependency of mcf/src/*.F90.

1) run script
 
  ./clean.sh

to clean all old files, which may be machine dependant.

2) run script

  ./bootstrap.sh

to make use of GNU autoconf + automake,
which generates configure script.

3) run configuration.

./configure --prefix=$HOME/MCF/mcf/mcf_install/ FC=ifort MPIFC=mpif90 LDFLAGS=-L$HOME/MCF/ppm/local/lib/ FCFLAGS=-I$HOME/MCF/ppm/local/include/ MAKEDEPF90=$HOME/MCF/ppm/local/bin/makedepf90

to generate Makefile, which is used to compile mcf code.

LDFLAGS: indicate ppm library object files.
FCFLAGS: indicate ppm library header files.


4) run compiler

  make -j 8

to compile mcf code,
-j 8 specify to use 8 processors to accelerate compiling.

5) run installation

  make install

to install mcf executable binary at ...../mcf/mcf_install/ folder.

-------------- USE -----------------------------------------------

Three input files are required to launch a simulation. Examples of
these input files can be found in the 'mcf_config' directory. The
details of the inputs are explained within those files.

-------------------------------------------------------------------
