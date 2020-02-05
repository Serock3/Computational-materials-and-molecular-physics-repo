#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {orthographic
  right -24.96*x up 26.01*y
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

cylinder {<  9.28,   7.11, -19.26>, < 11.89,  -3.08,  -7.21>, Rcell pigment {Black}}
cylinder {< -4.02,  12.27, -12.02>, < -1.41,   2.07,   0.03>, Rcell pigment {Black}}
cylinder {<-12.52,   1.07, -19.65>, < -9.91,  -9.13,  -7.60>, Rcell pigment {Black}}
cylinder {<  0.78,  -4.09, -26.89>, <  3.39, -14.29, -14.84>, Rcell pigment {Black}}
cylinder {<  9.28,   7.11, -19.26>, < -4.02,  12.27, -12.02>, Rcell pigment {Black}}
cylinder {< 11.89,  -3.08,  -7.21>, < -1.41,   2.07,   0.03>, Rcell pigment {Black}}
cylinder {<  3.39, -14.29, -14.84>, < -9.91,  -9.13,  -7.60>, Rcell pigment {Black}}
cylinder {<  0.78,  -4.09, -26.89>, <-12.52,   1.07, -19.65>, Rcell pigment {Black}}
cylinder {<  9.28,   7.11, -19.26>, <  0.78,  -4.09, -26.89>, Rcell pigment {Black}}
cylinder {< 11.89,  -3.08,  -7.21>, <  3.39, -14.29, -14.84>, Rcell pigment {Black}}
cylinder {< -1.41,   2.07,   0.03>, < -9.91,  -9.13,  -7.60>, Rcell pigment {Black}}
cylinder {< -4.02,  12.27, -12.02>, <-12.52,   1.07, -19.65>, Rcell pigment {Black}}
atom(< -0.11,   3.14, -12.96>, 1.48, rgb <0.76, 1.00, 0.00>, 0.0, ase2) // #0 
atom(<  1.14,   0.50, -14.81>, 1.48, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #1 
atom(< -2.14,   1.59, -15.17>, 1.48, rgb <0.76, 1.00, 0.00>, 0.0, ase2) // #2 
atom(< -1.48,  -1.72, -15.06>, 1.48, rgb <0.76, 1.00, 0.00>, 0.0, ase2) // #3 
atom(<  0.88,  -2.23, -12.67>, 1.48, rgb <0.76, 1.00, 0.00>, 0.0, ase2) // #4 
atom(<  1.75,   0.77, -11.42>, 1.48, rgb <0.76, 1.00, 0.00>, 0.0, ase2) // #5 
atom(< -1.58,   0.11, -12.10>, 1.48, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #6 
