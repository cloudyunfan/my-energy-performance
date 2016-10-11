function [ pl,ps,outcome,CHN_sta_f,channelslt ] = pktsend( CHNbefore_leng,CHNafter_leng,CHN_sta_ini,slotNO,Pu,Pd,channelslt )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
%CHNbefore_leng is the former slots length.
%CHNafter_leng is the latter slots length.
%slotNO is NO. of allocated slots of the node.
%CHN_sta_ini is temperal variable to record every INITIAL state in a superframe
%Pu is the  transition probability from the bad state to the good state.(up:0->1)
%Pd is the  transition probability from the good state to the bad state.(down:1->0)
%output:
%CHN_sta_f is the  INITIAL state  after the pktsend.
%ps is the NO. of successful packets
%pl is the NO. of lossed packets.
%outcome is  last slot state of current node
%CHN_sta_f is the last state after the superframe.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% channel state lasting time
global statelast
%%
ps = 0;
pl = 0;
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

