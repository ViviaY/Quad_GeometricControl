clc
clear all
close all

kr=0.6;

kw=1;

r=1;

A=[zeros(3,3) r*ones(3,3); -kr*eye(3,3), -kw*eye(3,3)];

At=0.5*(A+transpose(A));

max(eig(At))

