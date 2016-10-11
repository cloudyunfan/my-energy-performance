function [ pl,ps,outcome,CHN_sta_f,Succ_TX_flag,channelslt] = pktsendTDMA( CHNbefore_leng,CHNafter_leng,CHN_sta_ini,slotNO,Pu,Pd,channelslt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
%CHNbefore_leng : the former slots length.
%CHNafter_leng : the latter slots length.
%slotNO : NO. of allocated slots of the node.
%CHN_sta_ini : temperal variable to record every INITIAL state in a superframe
%Pu : the  transition probability from the bad state to the good state.(up:0->1)
%Pd : the  transition probability from the good state to the bad state.(down:1->0)
%Succ_TX_time: record the time when node send pkt successfully
%ind_SF: index of this superframe
%output:
%CHN_sta_f : the  INITIAL state  after the pktsend.
%ps : the NO. of successful packets
%pl : the NO. of lossed packets.
%outcome :  last slot state of current node
%CHN_sta_f : the last state after the superframe.
%Succ_TX_time： record the time when node send pkt successfully after thie superframe
%*******与pktsend的区别是多返回一个平均成功发包间隔Interval_avg***********
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% channel state lasting time
global statelast
%%
ps = 0;
pl = 0;
Succ_TX_flag = zeros(1,slotNO);
CHN_sta = CHN_sta_ini; % CHN_sta is a temperal variable updating every loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel state is updating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the channel state before current transmission，使用马尔科夫链来计算当前信道状态

for c = 1:floor( (channelslt + CHNbefore_leng)/statelast ) 
        if CHN_sta == 1
            CHN_sta = randsrc(1,1,[0 1;Pd 1-Pd]); %%%%%% channel model
        else
            CHN_sta = randsrc(1,1,[0 1;1-Pu Pu]); %%%%%% using Markov chain
        end
end
channelslt = mod(channelslt + CHNbefore_leng, statelast);
% transmitting  slotNO packets
for d = 1:slotNO  % each node i transmit slotNO slots
    if floor( (channelslt + slotNO)/statelast ) > 0
        if CHN_sta == 1
            CHN_sta = randsrc(1,1,[0 1;Pd 1-Pd]); %%%%%% channel model
        else
            CHN_sta = randsrc(1,1,[0 1;1-Pu Pu]); %%%%%% using Markov chain
        end
        channelslt = mod(channelslt + slotNO, statelast);
    end
        if CHN_sta == 0
            pl = pl+1; %%%%% calculate the NO. of lossed packets
        else
            ps = ps+1; %%%%% calculate the NO. of successful packets
            Succ_TX_flag(d) = 1;  %对成功发包的时隙进行标记
        end
end
% outcome of last slot of current node,记录发送包时的信道状态（0：不好；1：好）
outcome = CHN_sta;
% update the channel state after transmission
for e = 1:CHNafter_leng 
        if CHN_sta == 1
            CHN_sta = randsrc(1,1,[0 1;Pd 1-Pd]); %%%%%% channel model
        else
            CHN_sta = randsrc(1,1,[0 1;1-Pu Pu]); %%%%%% using Markov chain
        end
end
% the final channel state of current node
CHN_sta_f = CHN_sta;
end

