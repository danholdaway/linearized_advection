close all
clear
clc
mydir = pwd;

output

ymin = -0.2;
ymax = 1.1;

figure
set(gcf,'position',[46 590 1154 289])

subplot(1,4,1)
plot(q_free)
hold on
plot(qsl_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('Semi-Lagrangian')

subplot(1,4,2)
plot(q_free)
hold on
plot(qslice_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('SLICE')

subplot(1,4,3)
plot(q_free)
hold on
plot(qslbs_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('Semi-Lagrangian with Bermejo Staniforth')

subplot(1,4,4)
plot(q_free)
hold on
plot(qslicebs_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('SLICE with Bermejo Staniforth')






figure
set(gcf,'position',[46 149 1154 289])

subplot(1,4,1)
plot(q_repl)
hold on
plot(qsl_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('Semi-Lagrangian')

subplot(1,4,2)
plot(q_repl)
hold on
plot(qslice_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('SLICE')

subplot(1,4,3)
plot(q_repl)
hold on
plot(qslbs_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('Semi-Lagrangian with Bermejo Staniforth')

subplot(1,4,4)
plot(q_repl)
hold on
plot(qslicebs_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('SLICE with Bermejo Staniforth')




figure
set(gcf,'position',[1300 604 1154 289])

subplot(1,4,1)
plot(q_tlm)
hold on
plot(qsl_tlm,'r')
xlim([0 64])
title('Semi-Lagrangian')

subplot(1,4,2)
plot(q_tlm)
hold on
plot(qslice_tlm,'r')
xlim([0 64])
title('SLICE')

subplot(1,4,3)
plot(q_tlm)
hold on
plot(qslbs_tlm,'r')
xlim([0 64])
title('Semi-Lagrangian with Bermejo Staniforth')

subplot(1,4,4)
plot(q_tlm)
hold on
plot(qslicebs_tlm,'r')
xlim([0 64])
title('SLICE with Bermejo Staniforth')













cd(mydir)