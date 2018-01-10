function [I_prob, I_bw] = f_lognormal_mixture(I, r_c, r_g, K, Pf)
    % Ship detection based on lognormal mixture models
    % INPUT:
        % I:        SAR intensity image
        % r_c:      radius of the reference window, e.g., r_c = 15
        % r_g:      radius of the guard area, e.g., r_g = 10 (NB r_g < r_c)
        % K:        number of components, e.g., K = 3
        % Pf:       false alarm rate, e.g., Pf = 0.0005
    % OUTPUT:
        % I_prob:   image same size of I, pixel value is probability
        % I_bw:     binary image of detected ships, I_bw = I > 1-Pf
    
    % transform image into log-domain
    I_log = log(I);

    % pad array to avoid boundary issues
    I_log = padarray(I_log, [r_c, r_c], 'symmetric');
        
    % dimension of the image
    [mrows,ncols] = size(I_log);
    
    % setup the detection window
    ref_win = ones(2*r_c+1,2*r_c+1); 
    ref_win(r_c-r_g+1:r_c+r_g+1,r_c-r_g+1:r_c+r_g+1)=0;
    index = (ref_win==1);
    
    % generate probability image
    I_prob = zeros(mrows, ncols);
    h = waitbar(0, 'Processing');
    for i=(r_c+1):(mrows-r_c)
        waitbar(i/mrows);
        for j=(r_c+1):(ncols-r_c)
            temp = I_log(i-r_c:i+r_c,j-r_c:j+r_c);
            x = temp(index);
            if j==(r_c+1)
                [mu,variance,weight,~] = fMyEM(x,K,20,0.01);
            else
                [mu,variance,weight,~] = fMyEM(x,K,20,0.01,mu,variance,weight);
            end
            I_prob(i,j) = sum(weight.*normcdf(I_log(i,j),mu,sqrt(variance))); 
        end
    end
    close(h);
    
    % generate binary image
    I_bw = I_prob > 1-Pf;
    
    % get the original size of image
    I_prob = I_prob((r_c+1):(mrows-r_c), (r_c+1):(ncols-r_c));
    I_bw = I_bw((r_c+1):(mrows-r_c), (r_c+1):(ncols-r_c));
end
    
function [mu,variance,p,iterNum] = fMyEM(X,numOfComponent,maxIter,tol,muInit,varianceInit,pInit)
    % Customized version of Gaussian mixture model
    % This one is supposed to be faster than matlab's own version
    
    X = X(:);
    N = numel(X); % number of samples
    K = numOfComponent; % number of components

    if nargin==4
        % set initial values
        temp = randperm(N);
        mu = X(temp(1:K));mu = mu(:);
        clear temp;
        variance = var(X) * ones(K,1);
        p = (1/K) * ones(K,1);
    end
    if nargin==5
        % set initial values
        mu = muInit(:);
        variance = var(X) * ones(K,1);
        p = (1/K) * ones(K,1);
    end
    if nargin==6
        % set initial values
        mu = muInit(:);
        variance = varianceInit(:);
        p = (1/K) * ones(K,1);
    end
    if nargin==7
        mu = muInit(:);
        variance = varianceInit(:);
        p = pInit(:);
    end

    for i=1:maxIter
        G = repmat(p',N,1) ./ sqrt(repmat(variance',N,1)) .* exp(-((repmat(X,1,K)-repmat(mu',N,1)).^2) ./ (2*repmat(variance',N,1)));
        g = sum(G,2);g=g(:);
        G = G./ repmat(g,1,K);
        M = sum(G,1); M = M(:); M = M + 1e-6; M = M/sum(M) * N;
        muNEW = (G'* X) ./ M;
        varianceNEW = sum(G.*((repmat(X,1,K)-repmat(muNEW',N,1)).^2),1)'./M + 1e-6;
        pNEW = M/N;
        if (norm(mu-muNEW,1)/norm(mu,1)<tol && norm(varianceNEW-variance,1)/norm(variance,1)<tol && norm(pNEW-p,1)/norm(p,1)<tol)
            break;
        end
        mu = muNEW;
        variance = varianceNEW;
        p = pNEW;
    end
    iterNum = i;
    mu = muNEW;
    variance = varianceNEW;
    p = pNEW;
end
    