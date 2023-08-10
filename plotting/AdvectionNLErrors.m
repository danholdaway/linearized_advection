close all
clear
clc

output_step

N = 64;

e1st = q1st_free - q_free;
e2nd = q2nd_free - q_free;
e3rd = q3rd_free - q_free;
e3rdul = q3rdul_free - q_free;
e3rdrl = q3rdrl_free - q_free;
eppm_unl = qppm_unl_free - q_free;
eppm_lin = qppm_lin_free - q_free;
eppm_cw = qppm_cwl_free - q_free;
eppm_cs = qppm_csl_free - q_free;
eslice = qslice_free - q_free;
eslicebs = qslicebs_free - q_free;
espec = qspec_free - q_free;

l2_1st = sqrt(sum(e1st.^2)/N)
l2_2nd = sqrt(sum(e2nd.^2)/N)
l2_3rd = sqrt(sum(e3rd.^2)/N)
l2_3rdul = sqrt(sum(e3rdrl.^2)/N)
l2_3rdrl = sqrt(sum(e3rdrl.^2)/N)
l2_ppm_unl = sqrt(sum(eppm_unl.^2)/N)
l2_ppm_lin = sqrt(sum(eppm_lin.^2)/N)
l2_ppm_cw = sqrt(sum(eppm_cw.^2)/N)
l2_ppm_cs = sqrt(sum(eppm_cs.^2)/N)
l2_slice = sqrt(sum(eslice.^2)/N)
l2_slicebs = sqrt(sum(eslicebs.^2)/N)
l2_spec = sqrt(sum(espec.^2)/N)

plot(qspec_free)
hold on
plot(qppm_unl_free,'r')