% This demo implements the algorithm proposed in:
% Cui et al., "On semiparametric clutter estimation for ship detection in synthetic aperture radar images." 
% in IEEE transactions on geoscience and remote sensing 51, no. 5 (2013): 3170-3180.

%%
clear;
close all;

%% load sample image
load('./data/radarsat2-tj.mat')

%% do multi-looking to increase SCR
I_multilook = f_multilooking(I, 2, 2);

%% detection using semiparametric algorithm (KDE and Gaussian copula)
% NOTE: This will first ask you to drag a rectangle of sea clutter. After
% that, double click on the boundary of the rectangle to proceed.
[I_prob, I_bw] = f_semiparametric(I_multilook, 15, 10, 0.0005);

%% show some results
figure, imshow(I_multilook,[]);
figure, imshow(I_prob);
figure, imshow(I_bw);