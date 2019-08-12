% Generated on: 190810
% Last modification : 190812
% Author: Suwon Lee from Seoul National University
%
% This code runs the examples considered in Farouki(2019).
% https://linkinghub.elsevier.com/retrieve/pii/S0747717119300197
% 
% User can use spacePH handle class to allocate a new Pythagorean-hodograph
% (PH) curve.
% - The parameters p_i, p_f, t_i, t_f, S should be handed over to allocate
% a new instance of spacePH.
% -- p_i: initial position of the curve
% -- p_f: final position of the curve
% -- t_i: initial unit tangent of the curve
% -- t_f: final unit tangent of the curve
% -- S  : arc length of the curve
% - A point on the PH curve for parameter xi (in [0,1]) can be evaluated
% using spacePH.curve() method.
%
% User can use paintPH handle class to visualize the PH curve in 3D space.
% - The paintPH instance is constructed with input arguments of spacePH
% instance, and figureNumber.
% - The PH curve can be visualized using paintPH.visualize() method.

%%
clear all; clc; close all;

p_i = [0 0 0];
p_f = [1 0 0];
t_i = [cos(0.65*pi),sin(0.65*pi)*cos(-0.25*pi),sin(0.65*pi)*sin(-0.25*pi)];
t_f = [cos(0.25*pi),sin(0.25*pi)*cos(0.25*pi),sin(0.25*pi)*sin(0.25*pi)];
%% Example 1, Fig. 2.
for i = 1:5
  aa(i) = spacePH(p_i,p_f,t_i,t_f,1+0.1*i);
  bb(i) = paintPH(aa(i),1);
  bb(i).visualize;
  hold on;
  set(gca,'color','none');
  view([-156.2736 43.8201]);
end
%% Example 1, Fig. 3.
close all; clc;
A1 = spacePH(p_i,p_f,t_i,t_f,1.25);
A2 = spacePH(p_i,p_f,t_i,t_f,1.25);
A1.psi0 = 0;
A1.psi2 = 0;
A2.psi0 = 0;
A2.psi2 = 0;
B1 = paintPH(A1,1);
B2 = paintPH(A2,2);
B1.visualize;
B2.visualize;
for i = 1:7
  A1.psi0 = pi/4*i;
  A2.psi2 = pi/4*i;
  B1.visualize;
  B2.visualize;
end
figure(1);
set(gca,'color','none');
view([-156.2736 43.8201]);
figure(2);
set(gca,'color','none');
view([-156.2736 43.8201]);