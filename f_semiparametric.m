function [I_prob, I_bw] = f_semiparametric(I, r_c, r_g, Pf)
    % Ship detection based on KDE and Gaussian copula
    % INPUT:
        % I:        SAR intensity image
        % r_c:      radius of the reference window, e.g., r_c = 15
        % r_g:      radius of the guard area, e.g., r_g = 10 (NB 1 < r_g < r_c)
        % Pf:       false alarm rate, e.g., Pf = 0.0005
    % OUTPUT:
        % I_prob:   image same size of I, pixel value is probability
        % I_bw:     binary image of detected ships, I_bw = I > 1-Pf
    
    % transform the image into log-domain
    I_log = log(I);
    
    % user needs to select a homogeneous sea area
    figure, imshow(I, []);
    title('Please drag a rectangle of sea clutter and double click to proceed')
    h = imrect;
    pos = wait(h); 
    close;
    xmin = pos(1);
    ymin = pos(2);
    width = pos(3);
    height = pos(4);
    
    % estimate correlation structure of sea clutter
    I_train = I_log(ymin:(ymin+height),xmin:(xmin+width)); % train samples 
    N = numel(I_train);
    [bandwidth,~,~,~]=f_kde(I_train(:));
    [mrows_t,ncols_t] = size(I_train);
    I_unif = I_train;
    for i = 1:mrows_t
        for j = 1:ncols_t
            I_unif(i,j) = sum(normcdf(I_train(i,j),I_train(:),bandwidth))/N;
        end
    end

    U1 = zeros(floor(mrows_t/3),floor(ncols_t/3));
    U2 = zeros(floor(mrows_t/3),floor(ncols_t/3));
    U3 = zeros(floor(mrows_t/3),floor(ncols_t/3));
    U4 = zeros(floor(mrows_t/3),floor(ncols_t/3));
    U5 = zeros(floor(mrows_t/3),floor(ncols_t/3));
    U6 = zeros(floor(mrows_t/3),floor(ncols_t/3));
    U7 = zeros(floor(mrows_t/3),floor(ncols_t/3));
    U8 = zeros(floor(mrows_t/3),floor(ncols_t/3));
    U9 = zeros(floor(mrows_t/3),floor(ncols_t/3));
    for i = 2:3:mrows_t-1
        for j = 2:3:ncols_t-1
            tmp = I_unif(i-1:i+1,j-1:j+1);
            U1((i-2)/3+1,(j-2)/3+1) = tmp(2,2); % central pixel
            U2((i-2)/3+1,(j-2)/3+1) = tmp(1,1);
            U3((i-2)/3+1,(j-2)/3+1) = tmp(1,2);
            U4((i-2)/3+1,(j-2)/3+1) = tmp(1,3);
            U5((i-2)/3+1,(j-2)/3+1) = tmp(2,3);
            U6((i-2)/3+1,(j-2)/3+1) = tmp(3,3);
            U7((i-2)/3+1,(j-2)/3+1) = tmp(3,2);
            U8((i-2)/3+1,(j-2)/3+1) = tmp(3,1);
            U9((i-2)/3+1,(j-2)/3+1) = tmp(2,1);
        end
    end
    U1 = U1(:);
    U2 = U2(:);
    U3 = U3(:);
    U4 = U4(:);
    U5 = U5(:);
    U6 = U6(:);
    U7 = U7(:);
    U8 = U8(:);
    U9 = U9(:);
    RHOHAT = copulafit('Gaussian',[U1,U2,U3,U4,U5,U6,U7,U8,U9]);
    
    % pad array to avoid boundary issues
    I_log = padarray(I_log, [r_c, r_c], 'symmetric');
    
    % dimension of the image
    [mrows,ncols] = size(I_log);
    
    % setup the detection window
    ref_win = ones(2*r_c+1,2*r_c+1);
    ref_win(r_c-r_g+1:r_c+r_g+1,r_c-r_g+1:r_c+r_g+1)=0;
    index = find(ref_win==1);
    
    % bandwidth for KDE
    N0 = numel(index);
    hCV = bandwidth*(N/N0)^(1/5);
    
    % generate test statistic image
    Z = zeros(size(I_log));
    T_MQD = zeros(size(I_log));
    T_SPD = zeros(size(I_log));
    h = waitbar(0, 'Processing');
    for i=(r_c+1):(mrows-r_c)
        waitbar(i/mrows);
        for j=(r_c+1):(ncols-r_c)
            x = I_log(i-r_c:i+r_c,j-r_c:j+r_c);
            x = x(index);
            x = x(:);
            T_SPD(i,j) = sum(normcdf(I_log(i,j),x,hCV))/N0; 
            Z(i,j) = norminv(abs(T_SPD(i,j)-eps),0,1);
        end
    end
    close(h)
    z = zeros(9,1);
    for i=(r_c+1):(mrows-r_c)
        for j=(r_c+1):(ncols-r_c)
            z(1) = Z(i,j);
            z(2) = Z(i-1,j-1);
            z(3) = Z(i-1,j);
            z(4) = Z(i-1,j+1);
            z(5) = Z(i,j+1);
            z(6) = Z(i+1,j+1);
            z(7) = Z(i+1,j);
            z(8) = Z(i+1,j-1);
            z(9) = Z(i,j-1);
            if sum(z)<0
                T_MQD(i,j) = 0;
            else
                T_MQD(i,j) = z.'*(RHOHAT\z);
            end
        end
    end
    
    % obtain probability image and binary image
    I_prob = chi2cdf(T_MQD, 9);
    I_bw = I_prob > 1-Pf;
    
    % get the original size of image
    I_prob = I_prob((r_c+1):(mrows-r_c), (r_c+1):(ncols-r_c));
    I_bw = I_bw((r_c+1):(mrows-r_c), (r_c+1):(ncols-r_c));
end