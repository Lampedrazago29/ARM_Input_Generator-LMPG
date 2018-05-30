from pymol import cmd,stored 
bg_color white 
load 1u19.ARM.pdb 
set auto_zoom, off 
hide cartoon, all 
hide spheres, all 
hide nb_spheres, all 
hide sticks, all 
# Chromophore 
select resn RET 
create CHROMOPHORE, sele 
show sticks, CHROMOPHORE 
color tv_green, CHROMOPHORE 
#Linker amino acid 
select resi 296 
create linkerAA, sele 
show sticks, linkerAA 
color blue, linkerAA 
#Main counterion 
select resi 113 
create mainCounter, sele 
show sticks, mainCounter 
color lightblue, mainCounter 
#CL 
select resn CL 
create CL, sele 
show spheres, CL 
color red, CL 
#NA 
select resn NA 
create NA, sele 
show spheres, NA 
color blue, NA 
#Chain 
select resn HOH 
create WaterHOH, sele 
show nb_spheres, WaterHOH 
select all and not CHROMOPHORE and not CL and not NA and not WaterHOH and not mainCounter and not linkerAA 
create mainChain, sele 
show cartoon, mainChain 
color gray, mainChain 
set_bond stick_radius, 0.35, all 
set sphere_scale, 0.5 
set sphere_quality, 2 
set cartoon_transparency, 0.6 
set ray_trace_mode, 3 
set ray_shadows, 0 
set antialias, 2 
rotate x, -90 
rotate y, -90 
ray 1200,1200 
png 1u19.ARM.png 
