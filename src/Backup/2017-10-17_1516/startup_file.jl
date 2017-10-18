##########################################################################
##### 10/17/2017: Modified by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>                                
##### 10/17/2017: Previsouly modified by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>                                       
##### 10/17/2017: Created by Miura, Hirotaka <Hirotaka.Miura@ny.frb.org>
##########################################################################   
##### Description: 
##### 	- Program to set values at startup.
##### Modifications:
#####	10/17/2017: 
#####		- Duplicated from /home/rcehxm10/Hiro/02_Project/05_Marco/01_Program
#####		- Localize ENV["JULIA_PKGDIR"].
##########################################################################
##### For testing: Add current working dirrectory to LOAD_PATH.
push!(LOAD_PATH,".")
##### Set environmental variable to point to shared package repository.
ENV["JULIA_PKGDIR"]="/home/rcehxm10/Hiro/02_Project/06_Gizem/01_LSMSM/01_Program/lib";

