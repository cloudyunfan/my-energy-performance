%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ����802.15.6����µĶ�̬MACʱ϶���䣬���ѡ��ڵ����
%         Author:Ljg
%         Date:2015/03/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

%yf energy threshold,energy for transmission and CCA
global E_th
global N Tslot Data_rate TB Pbg Pgb CW CWmin CWmax UP UPnode Pkt_len E_TX E_CCA Tsim %isMAP lambdaE Emax Bmax isRAP
global channelslot statelast
%% ����MDP���
% load MDP_result(UP0-6,NH4)(Em20-Bm20,B0.05-E0.05)(NL3-3-18).mat %*******************************%
%------------802.15.6��ص�ȫ�ֲ���---------
UP = 0:7;  %8�����ȼ�
CWmin = [16,16,8,8,4,4,2,1];
CWmax = [64,32,32,16,16,8,8,4];

%----------------------���ò���-----------------------------------------------------
Tsim = 200; %NO. of superframes simulated
Tslot = 1;  % slot length (ms)
Pkt_len = 512; %packet length, unit is bit
Data_rate = 51.2; % transmission rate (kilo bits per second)
% yf initialize energy
E_th = 3;
E_CCA = E_th;   %�ŵ�������ĵ�����,���͡����ܡ���������1:1:1%*******************************%
E_TX = E_th;       %�������ݰ���Ҫ������

TB = 2000; %len_TDMA + len_RAP <=255
statelast = 10;
act = 2;
%yf omit the MAP
M = 25;   %MAP��ѯ�ʵ�ʱ϶���� M = 7;
T_block = 25;  %MAP��ÿһ�����ʱ϶��
len_MAP = M*T_block;  %MAP�ĳ���%*******************************%
len_RAP = TB-len_MAP; %RAP�׶ι̶���100��ʱ϶%*******************************%
UPH = 6;
UPN = 0;
NL = 2:2:18; %3:3:18
%% ------------------------------------------------------------------------
for indE = 1:length(NL)%   �������ȼ������
%     lambdaE = E_rate(indE);   
    NLnode = NL(indE)/2;
%   -------------���ýڵ�����ȼ�----------------------------
    UPclass = [UPH,UPN];%0,,6,7
    NH = NL(indE)-NLnode;
    N_UP = [NH,NLnode];  %ÿһ�����ȼ��ڵ�ĸ���
    UPnode = [];
    for up=1:length(UPclass)
       node = UPclass(up)*ones(1,N_UP(up)); 
       UPnode = [UPnode node];
    end
    N = length(UPnode);
    %��ʼ����������
    for n=1:N
        CW(n) = CWmin(find(UP==UPnode(n)));  %��ʼ��CWΪ�ڵ��Ӧ���ȼ���CWmin
    end

    %-----�ŵ�ģ��ʹ�������Ʒ������ŵ�״̬��ģ�������ŵ�״̬ת�Ƹ���---------
    Pbg = 0.4*ones(1,N);
    Pgb = 0.4*ones(1,N);   %���Ϊ�����ŵ��������ŵ��㶨ΪGOOD����
    channelslot = zeros(1, N);
    %------------------��ʼ����غ����ݻ�����----------------
    E_buff = (E_th)*ones(1,N); % ��ʼ�����ڵ�����״̬Ϊ0 yf��Ϊ����������    
     %----------------�������-----------------------------
    RAP_CHN_Sta = ones(1,N); % initial channel assumed to be GOOD. temperal variable to record every INITIAL state in a superframe
    TDMA_CHN_Sta = ones(1,N);    % initial channel assumed to be GOOD
    last_TX_time = ones(1,N); 
    CSMA_Sta_Pre = zeros(1,N);   % initial CSMA state assumed to be all 0 (0:initialization;1:backoff counter;2:sending packets)
    Def_Time_Pre = (-1)*ones(1,N); % initial deferred time -1
    ReTX_time_pre = zeros(1,N);  % ��ǽڵ��ش�����
    Succ_TX_time = zeros(Tsim*TB,N);   %��¼�ɹ������ʱ��
    Req = zeros(1,N); %�Ƿ��ڷ����������ݰ�����ʼ��Ϊ0���������������ݰ�    

    %--------------һ��Ҫͳ�ƵĽ��-------------------------------
    PL_RAP_sp = zeros(Tsim,N);  %������
    PL_MAP_sp = zeros(Tsim,N);    
    Colli_RAP_sp = zeros(Tsim,N);
    PS_RAP_sp = zeros(Tsim,N);   %�ɹ�����İ���
    PS_MAP_sp = zeros(Tsim,N);
    
    Swait = waitbar(0,'�������');   %���ý�����
    %last_TX_time_RAP = last_TX_time;
    for j = 1: Tsim
        %--------------------RAP�׶�ʹ��ʱ϶CSMA/CAʱ϶���䷽ʽ,���нڵ��������׶�------
        last_TX_time_RAP = last_TX_time-(j-1)*TB;
%         tic
        [ReTX_time_pre,Def_Time_Pre,CSMA_Sta_Pre,PL_RAP,PS_RAP,Colli_RAP,Succ_TX_time_RAP] = slotCSMACA_unsat_newnoE(len_RAP,CSMA_Sta_Pre,Def_Time_Pre,RAP_CHN_Sta,ReTX_time_pre,CW,last_TX_time_RAP,E_buff);%ELE_RAP,���������ģ�,E_buff,E_flow�������������һ�¸� CW����ȫ�ִ���ȥ
%         toc
        
        PL_RAP_sp(j,:) = PL_RAP;
        PS_RAP_sp(j,:) = PS_RAP;
        Colli_RAP_sp(j,:) = Colli_RAP;
        %ELE_RAP_sp(j,:) = ELE_RAP;
        for n=1:N
            %�������һ�γɹ�������ʱ��
            if( ~isempty(Succ_TX_time_RAP{n}) )
                %���³ɹ�����ʱ���¼
                ind_TX_RAP = Succ_TX_time_RAP{n} + (j-1)*TB;  %recover the real index
                last_TX_time(n) =  ind_TX_RAP(end);
                Succ_TX_time(ind_TX_RAP,n) = 1;
                %last_TX_time_RAP(n) = last_TX_time(n)-(j-1)*TB;
            end  
        end
        
        %--------------MAP �׶Σ�ʹ��TDMA��ʽ����ʱ϶--------------;               
         start = (j-1)*TB + len_RAP + 1; 
         TDMA_sift = 0;   %ƫ����  
         indMAP = find(ones(1,N)==1); %���нڵ㶼����MAP
         indPoll = getPollNode(indMAP,M);  %ȷ������poll�Ľڵ�
         %hist_MAP(j,1:length(indPoll))=indPoll;
         for poll =1:length(indPoll)   %�������и����ȼ��ڵ�ľ�����Ϊ     
             ind_node_poll = indPoll(poll); %ȡ�±�
            %----------scheduled slots---------------
            CHNafter_leng = 0;
            CHNbefore_leng = start + TDMA_sift - last_TX_time(ind_node_poll);
            last_TX_time_MAP = last_TX_time(n) - CHNbefore_leng - TDMA_sift;
            [PL_td,PS_td,lastout(ind_node_poll),TDMA_CHN_Sta(ind_node_poll),Succ_TX_time_td,channelslot(ind_node_poll)] = pktsendTDMA( CHNbefore_leng,CHNafter_leng,TDMA_CHN_Sta((ind_node_poll)),T_block,Pbg((ind_node_poll)),Pgb((ind_node_poll)),channelslot(ind_node_poll));
            if(~isempty(Succ_TX_time_td))
                %recover the real index                        
                ind_TX_MAP = Succ_TX_time_td + start + TDMA_sift;
                last_TX_time(ind_node_poll) = ind_TX_MAP(end);
                Succ_TX_time(ind_TX_MAP,ind_node_poll) = 1;
            end
           %---------------------����ͳ�Ʊ���------------------------
            TDMA_sift = TDMA_sift + T_block; %����ÿ���ڵ���䵽��ʱ϶������ƫ����                    
            %ELE_MAP_sp(j,ind_node_poll) = ELE_MAP_sp(j,ind_node_poll) + ELE_MAP;  %���ĵ��������� 
            PL_MAP_sp(j,ind_node_poll) = PL_td;
            PS_MAP_sp(j,ind_node_poll) = PS_td;                                              
         end 
         
        %--------------������ʾ-------------------------------
         str = ['�������', num2str(j*100/Tsim), '%'];     
         waitbar(j/Tsim,Swait,str);
    end
    close(Swait);

        %--------------yf���ŵ�������-------------------------------

    for n=1:N
        ind_Intv = find(Succ_TX_time(:,n)==1);
        Intv(n) = mean( diff(ind_Intv) );
        Slot_ulti(n) = length(ind_Intv)*100/length(Succ_TX_time(:,n));  %�ŵ�������
    end
    
    %-------------ͳ��ͨ�Ų�����Ҫ�Ľ��-------------------------------
    
    for up=1:length(UPclass)
        indUP = find(UPnode==UPclass(up));
        
        %EH_total(up,indE) = mean( sum( EH_sp(:,indUP) ) );  %�ܲɼ���������      
        %ELE_RAP_t(up,indE) = mean( sum( ELE_RAP_sp(:,indUP) ) );
        PS_RAP_total(up,indE) = mean( sum( PS_RAP_sp(:,indUP) ) );      %�ܵ�RAP�׶η��͵ĳ�֡����ȡ���нڵ��ƽ����   
        PS_MAP_total(up,indE) = mean( sum( PS_MAP_sp(:,indUP) ) );
        PL_RAP_total(up,indE) = mean( sum( PL_RAP_sp(:,indUP) ) );      %�ܵ�RAP�׶η��͵ĳ�֡����ȡ���нڵ��ƽ����   
        PL_MAP_total(up,indE) = mean( sum( PL_MAP_sp(:,indUP) ) );
        
        Colli_t(up,indE) = mean( sum( Colli_RAP_sp(:,indUP) ) );  %�ܳ�ͻ����ȡ���нڵ��ƽ����                
        Interval_avg(up,indE) = mean( Intv(indUP) );  %ƽ���ɹ�������� ,ȥ���ڵ��ƽ����
        Ulit_rate(up,indE) = mean( Slot_ulti(indUP) ); %ƽ���ŵ�������
%         Pktloss = PL_t./(PL_t+PS_t);
        Pktloss_rate(up,indE) = mean( sum( PL_RAP_sp(:,indUP)+PL_MAP_sp(:,indUP) )./sum( PS_RAP_sp(:,indUP)+PS_MAP_sp(:,indUP) ) );   %������ͬһ���ȼ��Ľڵ�(RAP + MAP)��ƽ�������ʱ�������        
        Pktloss_rate_RAP(up,indE) = mean( sum( PL_RAP_sp(:,indUP) )./sum( PS_RAP_sp(:,indUP) + PL_RAP_sp(:,indUP) ) ) ;   %������ͬһ���ȼ��Ľڵ�RAP��ƽ�������ʱ�������
        Pktloss_rate_MAP(up,indE) = mean( sum( PL_MAP_sp(:,indUP) )./sum( PS_MAP_sp(:,indUP) + PL_MAP_sp(:,indUP) ) ) ;   %������ͬһ���ȼ��Ľڵ��ƽ�������ʱ�������
                                                           
    end
%       %-----------------ͳ������WBAN�Ľ��---------------------------   

%         Act_time(indE) = mean(actTime_sp);
        %EH_WBAN(indE) = sum( sum(EH_sp) )/N;
        %yf
        %ELE_of_WBAN(indE) = sum( sum(ELE_RAP_sp) )/N;
        
        Interval_WBAN(indE) = mean(Intv);
        %yf
        %ELE_WBAN(indE) = sum( sum(ELE_RAP_sp) )/N;
        %yf
        PS_WBAN(indE) = sum( sum(PS_RAP_sp) )/N;       
        Colli_WBAN(indE) = sum( sum(Colli_RAP_sp) )/N;
        Pktloss_WBAN(indE) = sum(sum( PL_RAP_sp+PL_MAP_sp ))/sum(sum( PS_RAP_sp+PS_MAP_sp ));
        Pktloss_WBAN_RAP(indE) = sum(sum( PL_RAP_sp ))/sum(sum( PS_RAP_sp ));
        Pktloss_WBAN_RAP(indE) = sum(sum( PL_MAP_sp ))/sum(sum( PS_MAP_sp ));
%     disp(['indE NumUP lambdaE: ',num2str([indE N lambdaE])]) ; % yf һ������ģ�ȥ��lambdaE
      disp(['indE NumUP: ',num2str([indE N])]) 
end
disp('saturation VaringN simulation done!')
save('VarN_MAC(UP0-6,NH1-1-9)(P1_x0.9)(NL1-1-9)(Pgb0.4)(Pbg0.4)no.mat');