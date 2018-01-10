clc
clear
close all
clearvars -global

addpath('ginputc');
addpath('timetic-1.0');
addpath('mesh_based_warping');

aa = imread('image/a.jpg');
bb = imread('image/b.jpg');

[VCA,VCB,C,xyAB] = mesh_based_warping(aa,bb);