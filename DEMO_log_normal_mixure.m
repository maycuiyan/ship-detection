% This demo implements the algorithm proposed in:
% Cui et al., “CFAR ship detection in SAR images based on lognormal mixture models,” 
% in Proc. 3rd IEEE Int. Asia-Pacific Conf. on Synthetic Aperture Radar, pp. 1–3, IEEE, Seoul, South Korea (2011).

%%
clear;
close all;

%% load sample image
load('radarsat2-tj.mat')

%% do multi-looking to increase SCR
I_multilook = f_multilooking(I, 2, 2);

%% detection using lognormal mixture model
[I_prob, I_bw] = f_lognormal_mixture(I_multilook, 15, 10, 3, 0.0005);

%% show some results
figure, imshow(I_multilook,[]);
figure, imshow(I_prob);
figure, imshow(I_bw);