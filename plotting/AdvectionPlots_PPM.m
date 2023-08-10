close all
clear
clc
mydir = pwd;


output

ymin = -0.2;
ymax = 1.1;

figure
set(gcf,'position',[125 34 1154 885])

subplot(2,4,1)
plot(q_free)
hold on
plot(qppm_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM 4 With Limiter')

subplot(2,4,2)
plot(q_free)
hold on
plot(qppm1_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM No Limiter')

subplot(2,4,3)
plot(q_free)
hold on
plot(qppm2_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CW Part 1')

subplot(2,4,4)
plot(q_free)
hold on
plot(qppm3_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CW Part 2')

subplot(2,4,5)
plot(q_free)
hold on
plot(qppm4_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CW Full')

subplot(2,4,6)
plot(q_free)
hold on
plot(qppm4a_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CW Full - FV Version')

subplot(2,4,7)
plot(q_free)
hold on
plot(qppm5_free,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CS Full')




figure
set(gcf,'position',[125 34 1154 885])

subplot(2,4,1)
plot(q_repl)
hold on
plot(qppm_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM 4 With Limiter')

subplot(2,4,2)
plot(q_repl)
hold on
plot(qppm1_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM No Limiter')

subplot(2,4,3)
plot(q_repl)
hold on
plot(qppm2_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CW Part 1')

subplot(2,4,4)
plot(q_repl)
hold on
plot(qppm3_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CW Part 2')

subplot(2,4,5)
plot(q_repl)
hold on
plot(qppm4_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CW Full')

subplot(2,4,6)
plot(q_repl)
hold on
plot(qppm4a_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CW Full - FV Version')

subplot(2,4,7)
plot(q_repl)
hold on
plot(qppm5_repl,'r')
ylim([ymin ymax])
xlim([0 64])
title('PPM CS Full')




figure
set(gcf,'position',[1376 48 1154 885])
subplot(2,4,1)
plot(q_tlm)
hold on
% plot(qppm_tlm,'r')
plot(qppm_repl-qppm_free,'k--')
title('PPM 4 With Limiter')
xlim([0 64])

subplot(2,4,2)
plot(q_tlm)
hold on
plot(qppm1_tlm,'r')
plot(qppm1_repl-qppm1_free,'k--')
title('PPM No Limiter')
xlim([0 64])

subplot(2,4,3)
plot(q_tlm)
hold on
plot(qppm2_tlm,'r')
plot(qppm2_repl-qppm2_free,'k--')
title('PPM CW Part 1')
xlim([0 64])

subplot(2,4,4)
plot(q_tlm)
hold on
plot(qppm3_tlm,'r')
plot(qppm3_repl-qppm3_free,'k--')
title('PPM CW Part 2')
xlim([0 64])

subplot(2,4,5)
plot(q_tlm)
hold on
plot(qppm4_tlm,'r')
plot(qppm4_repl-qppm4_free,'k--')
title('PPM CW Full')
xlim([0 64])

subplot(2,4,6)
plot(q_tlm)
hold on
plot(qppm4a_tlm,'r')
plot(qppm4a_repl-qppm4a_free,'k--')
title('PPM CW Full - FV Version')
xlim([0 64])

subplot(2,4,7)
plot(q_tlm)
hold on
plot(qppm5_tlm,'r')
plot(qppm5_repl-qppm5_free,'k--')
title('PPM CS Full')
xlim([0 64])

























cd(mydir)