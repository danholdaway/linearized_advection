close all
clear
clc

evals
jac
xmin = 1;
xmax = 640;

x = 1:length(maxeval_ppm_lin);

figure
set(gcf,'position',[125 34 1154 885])

subplot(3,3,1)
plot(x,maxeval_ppm_lin)
title('PPM Lin limiter')
ylabel('\lambda')

subplot(3,3,2)
plot(x,maxeval_ppm_cwl)
title('PPM CW limiter')
ylabel('\lambda')

subplot(3,3,3)
plot(x,maxeval_ppm_csl)
title('PPM CS limiter')
ylabel('\lambda')

subplot(3,3,4)
plot(x,maxeval_ppm_unl)
title('PPM no limiter')
ylabel('\lambda')

subplot(3,3,5)
plot(x,maxeval_3rd)
title('Third order')
ylabel('\lambda')

subplot(3,3,6)
plot(x,maxeval_3rdul)
title('3rd UL')
ylabel('\lambda')

subplot(3,3,7)
plot(x,maxeval_3rdrl)
title('3rd RUL')
ylabel('\lambda')

subplot(3,3,8)
plot(x,maxeval_slice)
title('SLICE')
ylabel('\lambda')

subplot(3,3,9)
plot(x,maxeval_slicebs)
title('SLICEBS')
ylabel('\lambda')

count_ppm_lin = 0;
count_ppm_cw = 0;
count_ppm_cs = 0;
count_3rd_ul = 0;
count_3rd_rl = 0;

for i = 1:length(maxeval_slicebs)
    
    if maxeval_ppm_lin(i)
        count_ppm_lin = count_ppm_lin + 1;
    end
    
    if maxeval_ppm_cwl(i)
        count_ppm_cw = count_ppm_cw + 1;
    end
    
    if maxeval_ppm_csl(i)
        count_ppm_cs = count_ppm_cs + 1;
    end
    
    if maxeval_3rdul(i)
        count_3rd_ul = count_3rd_ul + 1;
    end
    
    if maxeval_3rdrl(i)
        count_3rd_rl = count_3rd_rl + 1;
    end
    
    
end

count_ppm_lin
count_ppm_cw
count_ppm_cs
count_3rd_ul
count_3rd_rl

mean(maxeval_ppm_lin)
mean(maxeval_ppm_cwl)
mean(maxeval_ppm_csl)
mean(maxeval_3rdul)
mean(maxeval_3rdrl)