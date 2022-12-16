format shortg
format compact
close all
clear all
clc

% uncertain state space system
% must be linear in the parameters
zeta = ureal('zeta',0.5,'PlusMinus',1);
eta = ureal('eta',0,'PlusMinus',1);
A = [-0.6, 4+zeta+0.3*eta; -4, eta];
B = [0, 0; 1.5, 0];
C = [0.3*eta, -1.2];
D = [0, 1];
usys = uss(A,B,C,D);

[Poly,legacy] = uss2pol(usys)

legacy.Sys
Poly.Sys
