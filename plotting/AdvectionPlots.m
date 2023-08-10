close all
clear
clc
mydir = pwd;

output

N = 64;
x = 0:1/(N-1):1;

ymin = -0.2;
ymax = 1.1;

yminp = -0.2e-5;
ymaxp = 1.2e-5;

figure
set(gcf,'position',[125 34 1154 885])
subplot(3,4,1)
plot(x,q_free)
hold on
plot(x,q1st_free,'r')
ylim([ymin ymax])
title('1st order FD')
ylabel('q')

subplot(3,4,2)
plot(x,q_free)
hold on
plot(x,q2nd_free,'r')
ylim([ymin ymax])
title('2nd order FD')

subplot(3,4,3)
plot(x,q_free)
hold on
plot(x,q3rd_free,'r')
ylim([ymin ymax])
title('3rd order FD')

subplot(3,4,4)
plot(x,q_free)
hold on
plot(x,qppm_unl_free,'r')
ylim([ymin ymax])
title('PPM unlimited')

subplot(3,4,5)
plot(x,q_free)
hold on
plot(x,qppm_lin_free,'r')
ylim([ymin ymax])
title('PPM Lin limiter')
ylabel('q')

subplot(3,4,6)
plot(x,q_free)
hold on
plot(x,qppm_cwl_free,'r')
ylim([ymin ymax])
title('PPM Collela Woodward limiter')

subplot(3,4,7)
plot(x,q_free)
hold on
plot(x,qppm_csl_free,'r')
ylim([ymin ymax])
title('PPM Collela Sekora limiter')

subplot(3,4,8)
plot(x,q_free)
hold on
plot(x,q3rdrl_free,'r')
ylim([ymin ymax])
title('3rd order with relaxed limiter')

subplot(3,4,9)
plot(x,q_free)
hold on
plot(x,q3rdul_free,'r')
ylim([ymin ymax])
title('3rd order with universal limiter')
xlabel('x')
ylabel('q')

subplot(3,4,10)
plot(x,q_free)
hold on
plot(x,qslice_free,'r')
ylim([ymin ymax])
title('SLICE unlimited')
xlabel('x')

subplot(3,4,11)
plot(x,q_free)
hold on
plot(x,qslicebs_free,'r')
ylim([ymin ymax])
title('SLICE with Bermejo Staniforth limiter')
xlabel('x')

subplot(3,4,12)
plot(x,q_free)
hold on
plot(x,qspec_free,'r')
ylim([ymin ymax])
title('Spectral method')
xlabel('x')


figure
set(gcf,'position',[125 34 1154 885])
subplot(3,4,1)
plot(x,q_repl)
hold on
plot(x,q1st_repl,'r')
ylim([ymin ymax])
title('1st order FD')
ylabel('q')

subplot(3,4,2)
plot(x,q_repl)
hold on
plot(x,q2nd_repl,'r')
ylim([ymin ymax])
title('2nd order FD')

subplot(3,4,3)
plot(x,q_repl)
hold on
plot(x,q3rd_repl,'r')
ylim([ymin ymax])
title('3rd order FD')

subplot(3,4,4)
plot(x,q_repl)
hold on
plot(x,qppm_unl_repl,'r')
ylim([ymin ymax])
title('PPM unlimited')

subplot(3,4,5)
plot(x,q_repl)
hold on
plot(x,qppm_lin_repl,'r')
ylim([ymin ymax])
title('PPM Lin limiter')
ylabel('q')

subplot(3,4,6)
plot(x,q_repl)
hold on
plot(x,qppm_cwl_repl,'r')
ylim([ymin ymax])
title('PPM Collela Woodward limiter')

subplot(3,4,7)
plot(x,q_repl)
hold on
plot(x,qppm_csl_repl,'r')
ylim([ymin ymax])
title('PPM Collela Sekora limiter')

subplot(3,4,8)
plot(x,q_repl)
hold on
plot(x,q3rdrl_repl,'r')
ylim([ymin ymax])
title('3rd order with relaxed limiter')

subplot(3,4,9)
plot(x,q_repl)
hold on
plot(x,q3rdul_repl,'r')
ylim([ymin ymax])
title('3rd order with universal limiter')
xlabel('x')
ylabel('q')

subplot(3,4,10)
plot(x,q_repl)
hold on
plot(x,qslice_repl,'r')
ylim([ymin ymax])
title('SLICE unlimited')
xlabel('x')

subplot(3,4,11)
plot(x,q_repl)
hold on
plot(x,qslicebs_repl,'r')
ylim([ymin ymax])
title('SLICE with Bermejo Staniforth limiter')
xlabel('x')

subplot(3,4,12)
plot(x,q_repl)
hold on
plot(x,qspec_repl,'r')
ylim([ymin ymax])
title('Spectral method')
xlabel('x')


figure
set(gcf,'position',[1341 53 1154 885])

subplot(3,4,1)
plot(x,q_repl-q_free,'k')
hold on
plot(x,q1st_repl-q1st_free,'b')
plot(x,q1st_tlm,'r--')
%ylim([yminp ymaxp])
title('1st order FD')
ylabel('q^prime')

subplot(3,4,2)
plot(x,q_repl-q_free,'k')
hold on
plot(x,q2nd_repl-q2nd_free,'b')
plot(x,q2nd_tlm,'r--')
%ylim([yminp ymaxp])
title('2nd order FD')

subplot(3,4,3)
plot(x,q_repl-q_free,'k')
hold on
plot(x,q3rd_repl-q3rd_free,'b')
plot(x,q3rd_tlm,'r--')
%ylim([yminp ymaxp])
title('3rd order FD')

subplot(3,4,4)
plot(x,q_repl-q_free,'k')
hold on
plot(x,qppm_unl_repl-qppm_unl_free,'b')
plot(x,qppm_unl_tlm,'r--')
%ylim([yminp ymaxp])
title('PPM unlimited')

subplot(3,4,5)
plot(x,q_repl-q_free,'k')
hold on
plot(x,qppm_lin_repl-qppm_lin_free,'b')
plot(x,qppm_lin_tlm,'r--')
%ylim([yminp ymaxp])
title('PPM Lin limiter')
ylabel('q^prime')

subplot(3,4,6)
plot(x,q_repl-q_free,'k')
hold on
plot(x,qppm_cwl_repl-qppm_cwl_free,'b')
plot(x,qppm_cwl_tlm,'r--')
%ylim([yminp ymaxp])
title('PPM Collela Woodward limiter')

subplot(3,4,7)
plot(x,q_repl-q_free,'k')
hold on
plot(x,qppm_csl_repl-qppm_csl_free,'b')
plot(x,qppm_csl_tlm,'r--')
%ylim([yminp ymaxp])
title('PPM Collela Sekora limiter')

subplot(3,4,8)
plot(x,q_repl-q_free,'k')
hold on
plot(x,q3rdrl_repl-q3rdrl_free,'b')
plot(x,q3rdrl_tlm,'r--')
%ylim([yminp ymaxp])
title('3rd order with relaxed limiter')

subplot(3,4,9)
plot(x,q_repl-q_free,'k')
hold on
plot(x,q3rdul_repl-q3rdul_free,'b')
plot(x,q3rdul_tlm,'r--')
%ylim([yminp ymaxp])
title('3rd order with universal limiter')
xlabel('x')
ylabel('q^prime')

subplot(3,4,10)
plot(x,q_repl-q_free,'k')
hold on
plot(x,qslice_repl-qslice_free,'b')
plot(x,qslice_tlm,'r--')
%ylim([yminp ymaxp])
title('SLICE unlimited')
xlabel('x')

subplot(3,4,11)
plot(x,q_repl-q_free,'k')
hold on
plot(x,qslicebs_repl-qslicebs_free,'b')
plot(x,qslicebs_tlm,'r--')
%ylim([yminp ymaxp])
title('SLICE with Bermejo Staniforth limiter')
xlabel('x')

subplot(3,4,12)
plot(x,q_repl-q_free,'k')
hold on
plot(x,qspec_repl-qspec_free,'b')
plot(x,qspec_tlm,'r--')
%ylim([yminp ymaxp])
title('Spectral method')
xlabel('x')


















cd(mydir)