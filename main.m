clear all; clc; close all;

p_i = [0 0 0];
p_f = [1 0 0];
% ti = [-1 0 -0.2];
% tf = [1 -1 -0.5];
t_i = [cos(0.65*pi),sin(0.65*pi)*cos(-0.25*pi),sin(0.65*pi)*sin(-0.25*pi)];
t_f = [cos(0.25*pi),sin(0.25*pi)*cos(0.25*pi),sin(0.25*pi)*sin(0.25*pi)];
%%
for i = 1:5
  aa(i) = spacePH(p_i,p_f,t_i,t_f,1+0.1*i);
  aa(i) = aa(i).updateTransients;
  bb(i) = paintPH(aa(i));
  hold on;
  set(gca,'color','none');
  view([-156.2736 43.8201]);
end