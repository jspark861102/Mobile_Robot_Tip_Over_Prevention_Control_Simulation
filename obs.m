clear all
close all
clc

M = 10;
m1 = 2;
m2 = 0.1;
k1 = 10000;
k2 = 100;

Mmatrix = diag([M m1 m2]);
Kmatrix = [k1 -k1 0; -k1 k1+k2 -k2; 0 -k2 k2];
Fmatrix = [1;0;0];

A = -inv(Mmatrix)*Kmatrix;
B = inv(Mmatrix)*Fmatrix;
C = [1 0 0; 0 0 1];
% C = [0 0 1];


rank(ctrb(A,B))
rank(obsv(A,C))

obsv(A,C)
