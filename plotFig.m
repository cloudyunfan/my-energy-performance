load('VarN_MAC(UP0-6,NH1-1-9)(P1_x0.9)(NL1-1-9)noCHnoE.mat');
N_all = 2:2:18;
plot(N_all, EE_MAP_t(1,:), N_all, EE_RAP_t(1,:), '--');
grid;
axis([2, 18, 0, 0.025]);
title('energy efficiency of CSMA/CA and TDMA');
xlabel('Number of nodes in WBAN');
ylabel('energy efficiency');
legend('TDMA', 'CSMA/CA');
