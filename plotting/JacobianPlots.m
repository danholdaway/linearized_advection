close all
clear
clc
mydir = pwd;

jac

x = 1:64;

%Normal matrix test
J3rd_norm = J3rd*J3rd' - J3rd'*J3rd;
J3rdul_norm = J3rdul*J3rdul' - J3rdul'*J3rdul;
J3rdrl_norm = J3rdrl*J3rdrl' - J3rdrl'*J3rdrl;
Jlin_norm = Jlin*Jlin' - Jlin'*Jlin;
Jppmunl_norm = Jppmunl*Jppmunl' - Jppmunl'*Jppmunl;
Jcsl_norm = Jcsl*Jcsl' - Jcsl'*Jcsl;
Jcwl_norm = Jcwl*Jcwl' - Jcwl'*Jcwl;
Jslice_norm = Jslice*Jslice' - Jslice'*Jslice;
Jslicebs_norm = Jslicebs*Jslicebs' - Jslicebs'*Jslicebs;

disp(max(J3rd_norm(:)))
disp(max(J3rdul_norm(:)))
disp(max(J3rdrl_norm(:)))
disp(max(Jppmunl_norm(:)))
disp(max(Jlin_norm(:)))
disp(max(Jcsl_norm(:)))
disp(max(Jcwl_norm(:)))
disp(max(Jslice_norm(:)))
disp(max(Jslicebs_norm(:)))

sad

[U1,eig1] = eig(J3rdul);
eig1 = diag(eig1);
[U2,eig2] = eig(J3rd);
eig2 = diag(eig2);


%Sort eigenvalues
[ieig1, ind] = sort(imag(eig1),'descend');
reig1 = real(eig1(ind));
U1 = U1(:,ind);

[reig2, ind] = sort(real(eig2),'descend');
ieig2 = imag(eig2(ind));
U2 = U2(:,ind);

% figure
% scatter(reig1,ieig1)
% hold on
% scatter(reig2,ieig2,'rx')
% 
% box on
% legend('PPM Lin Limiter','3rd Order Universal Lim','Location','NorthWest')
% xlabel('Re(\lambda)')
% ylabel('Im(\lambda)')
% xlim([0.75 1.05])
% 
% close


for i = 1:2:63

    %PPM SCHEME
    figure(3)
    set(gcf,'position',[228 553 1051 366])
    subplot(2,2,[1,3])
    scatter(reig1,ieig1)
    hold on
    scatter(reig1(i),ieig1(i),'rx')
    scatter(reig1(end-i+1),ieig1(end-i+1),'gx')
    hold off
    box on
%     legend('PPM Lin Limiter','3rd Order Universal Lim','Location','NorthWest')
    xlabel('Re(\lambda)')
    ylabel('Im(\lambda)')
    
    subplot(2,2,2)
    plot(x,real(U1(:,i)),'b')
    hold on
    plot(x,imag(U1(:,i)),'r')
    hold off
    
    subplot(2,2,4)
    plot(x,real(U1(:,end-i+1)),'b')
    hold on
    plot(x,imag(U1(:,end-i+1)),'r')
    hold off


    %Third Order universal limiter scheme
    figure(4)
    set(gcf,'position',[228 85 1051 366])
    subplot(2,2,[1,3])
    scatter(reig2,ieig2)
    hold on
    scatter(reig2(i),ieig2(i),'rx')
    scatter(reig2(i+1),ieig2(i+1),'gx')
    hold off
    box on
%     legend('PPM Lin Limiter','3rd Order Universal Lim','Location','NorthWest')
    xlabel('Re(\lambda)')
    ylabel('Im(\lambda)')

    subplot(2,2,2)
    plot(x,real(U2(:,i)),'b')
    hold on
    plot(x,imag(U2(:,i)),'r')
    hold off

    subplot(2,2,4)
    plot(x,real(U2(:,i+1)),'b')
    hold on
    plot(x,imag(U2(:,i+1)),'r')
    hold off

    pause
    
end





















cd(mydir)