#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {orthographic
  right -12.04*x up 12.55*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}

#declare simple = finish {phong 0.7}
#declare pale = finish {ambient 0.5 diffuse 0.85 roughness 0.001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 specular 0.5 }
#declare jmol = finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.05 diffuse 0.3 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.01 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.070;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

cylinder {< -8.10,  -7.77, -16.00>, < -8.10,  -7.77,   0.00>, Rcell pigment {Black}}
cylinder {<  7.90,  -7.77, -16.00>, <  7.90,  -7.77,  -0.00>, Rcell pigment {Black}}
cylinder {<  7.90,   8.23, -16.00>, <  7.90,   8.23,  -0.00>, Rcell pigment {Black}}
cylinder {< -8.10,   8.23, -16.00>, < -8.10,   8.23,   0.00>, Rcell pigment {Black}}
cylinder {< -8.10,  -7.77, -16.00>, <  7.90,  -7.77, -16.00>, Rcell pigment {Black}}
cylinder {< -8.10,  -7.77,   0.00>, <  7.90,  -7.77,  -0.00>, Rcell pigment {Black}}
cylinder {< -8.10,   8.23,   0.00>, <  7.90,   8.23,  -0.00>, Rcell pigment {Black}}
cylinder {< -8.10,   8.23, -16.00>, <  7.90,   8.23, -16.00>, Rcell pigment {Black}}
cylinder {< -8.10,  -7.77, -16.00>, < -8.10,   8.23, -16.00>, Rcell pigment {Black}}
cylinder {< -8.10,  -7.77,   0.00>, < -8.10,   8.23,   0.00>, Rcell pigment {Black}}
cylinder {<  7.90,  -7.77,  -0.00>, <  7.90,   8.23,  -0.00>, Rcell pigment {Black}}
cylinder {<  7.90,  -7.77, -16.00>, <  7.90,   8.23, -16.00>, Rcell pigment {Black}}
atom(<  1.36,   1.14,  -8.93>, 1.48, rgb <0.67, 0.36, 0.95>, 0.0, ase2) // #0 
atom(< -1.01,   3.41,  -9.25>, 1.48, rgb <0.67, 0.36, 0.95>, 0.0, ase2) // #1 
atom(< -2.06,   0.44,  -8.28>, 1.48, rgb <0.67, 0.36, 0.95>, 0.0, ase2) // #2 
atom(< -2.83,  -2.63,  -7.34>, 1.48, rgb <0.67, 0.36, 0.95>, 0.0, ase2) // #3 
atom(<  3.52,  -1.29,  -8.42>, 1.48, rgb <0.67, 0.36, 0.95>, 0.0, ase2) // #4 
atom(<  0.37,  -2.09,  -7.86>, 1.48, rgb <0.67, 0.36, 0.95>, 0.0, ase2) // #5 
