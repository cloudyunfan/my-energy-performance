load('VarN_MAC(UP0-6,NH1-1-9)(P1_x0.9)(NL1-1-9)(Pgb0.4)(Pbg0.4)no.mat');
N_all = 2:2:18;
plot(N_all, Pktloss_rate_MAP(1,:), N_all, Pktloss_rate_RAP(1,:), '--');
grid;
title('successful transmission probability of CSMA/CA and TDMA');
xlabel('Number of nodes in WBAN');
ylabel('successful transmission probability');
legend('TDMA', 'CSMA/CA');
