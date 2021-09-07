%% Powermeter measure
% Measurement wavelength was fixed to 514 nm (arbitrary)

%% 0806 Measure without blackcorrection
p_blk = 598.003 % nW

p_rand_test2 = 848.263 % nW
p_rand_R_test2 = 238.250
p_rand_G_test2 = 624.224
sum_test2 = 862.4740

p_rand_test3 = 1036
p_rand_R_test3 = 471.076
p_rand_G_test3 = 573.109
sum_test3 = 1044.23

p_rand_test6 = 927.322
p_rand_R_test6 = 312.058
p_rand_G_test6 = 624.625
sum_test6 = 936.683

p_rand_set = [p_rand_test2 p_rand_test3 p_rand_test6];
p_rand_R_set = [p_rand_R_test2 p_rand_R_test3 p_rand_R_test6];
p_rand_G_set = [p_rand_G_test2 p_rand_G_test3 p_rand_G_test6];
p_rand_sum_set = p_rand_R_set + p_rand_G_set;

figure(3); hold on;
plot(p_rand_set,p_rand_sum_set,'r.','MarkerSize',13)
plot([1 1400],[1 1400],'k-')
xlabel('Power meter measure (nW) (Single)')
ylabel('Power meter measure (nW) (Sum)')
% xlim([800 1050])
% ylim([800 1050])
legend('Test color with 2 primaries (N=3)','45-deg line','location','southeast')

%% 0810 measure with black correction
p_blk = 610.528 % nW

p_rand_set = [1375 1562 1451] - p_blk;
p_rand_R_set = [803.896 1022 876.498] - p_blk;
p_rand_G_set = [1163 1127 1164] - p_blk;

p_rand_sum_set = p_rand_R_set + p_rand_G_set;

figure(3); hold on;
plot(p_rand_set,p_rand_sum_set,'r.','MarkerSize',13)
plot([1 1400],[1 1400],'k-')
xlabel('Power meter measure (nW) (Single)')
ylabel('Power meter measure (nW) (Sum)')
% xlim([800 1050])
% ylim([800 1050])
legend('Test color with 2 primaries (N=3)','45-deg line','location','southeast')


fw_rand_G = [fw_rand_G_test2 fw_rand_G_test3 fw_rand_G_test6];

fw_blks = ones(201,3).*fw_blk; 

w = [380:2:780];
fw_rand_blk = fw_rand -fw_blks;
fw_rand_R_blk = fw_rand_R - fw_blks;
fw_rand_G_blk = fw_rand_G - fw_blks;
fw_rand_sum = fw_rand_R_blk + fw_rand_G_blk;

figure(1); hold on;
plot(w,fw_rand_blk(:,3),'k-');
plot(w,fw_rand_sum(:,3),'r--');

