function I_multilook = f_multilooking(I,az_look,rg_look)

N_az = floor(size(I,1)/az_look);
N_rg = floor(size(I,2)/rg_look);
I_multilook = zeros(N_az,N_rg);
for i = 1:N_az
%     i
    for j = 1:N_rg
        I_multilook(i,j) = mean2(I((i-1)*az_look+1:i*az_look,(j-1)*rg_look+1:j*rg_look));
    end
end