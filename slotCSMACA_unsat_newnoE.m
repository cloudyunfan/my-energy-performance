function [ReTX_time,backoff_after,CSMA_sta,ELE_ex,pl_t,ps_t,PL_colli,TX_time] = slotCSMACA_unsat_newnoE(rap_length,CSMA_sta,def_time_pre,last_CHN_sta,ReTX_times_pre,CW,last_TX_time,E_buff) %act,,B_buff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ELE_ex(��������),E_buff,E_flow�����������
%%% CSMA/CA transmission under unsaturation condition
% Input:
%     2.rap_length: duration of RAP
%     3.CSMA_sta: last CSMA state in pre superframe
%     4.def_time_pre: backoff counter before superframe(-1: node has not obtained any slot before or just send packet successfully)
%     5.last_CHN_sta: Channel state (0: BAD,1: GOOD)
%     6.ReTX_times_pre: last retransmission times at the end of RAP
%     7.CW: contention window of every node
%     8.last_TX_time�� last time sending packet successfully
%     10:E_buff:����������
%     11:e_flow:������
% Output:
%     1.ReTX_time: remaining retransmission times after this superframe
%     2.backoff_after: remaining backoff counter after this superframe
%     3.CSMA_sta:  last CSMA state in this superframe(0:initialization;1:backoff counter;2:sending packets)
%     4.tdma_CHN_sta:last Channel state(0:
%     idle,1:busy)%*******************************%
%     5.pl_t: number of lost packets
%     6.ps_t:number of packets sent succesiffuly 
%     7.PL_colli: collision times
%     8.ELE_ex: total energy exhost
%     9.TX_time: record time when node send packet successfully

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters initialization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Data_rate Pkt_len Tslot Pbg Pgb CWmin CWmax UP UPnode E_TX E_CCA channelslot statelast %isMAP isRAP 
%yf probability of arrive one unit energy in each slot
global P1_x
%-----------����---------------------------------------------------------
%yf
P1_x = 0.6;
SIFS = 0.5;      %���֡����0 yf �е����⣬һ�����0.5
T_pkt = ceil( Pkt_len/(Data_rate*Tslot) );  %time to send a packet and receive ACK,unit slot
M = 4;        %����ش�����
N = length(E_buff);  %��ȡʹ��CSMA/CA�Ľڵ���

ReTX_time = ReTX_times_pre;  %��ȡ��һ��֡�е��ش�������¼
TX_time_rap = zeros(rap_length,N); % the time when node send pkt successfully

Backoff_time = zeros(1,N); % temporal variable to record backoff time
CHN = zeros(1,rap_length); % ��ʼ��ÿ��ʱ϶���ŵ�״��Ϊ����
ELE_ex = zeros(1,N); % ��¼����֡���ĵ�������
backoff_after = zeros(1,N); % deferred backoff should be output,defoult 0
backoff_lock = zeros(1,N);  % flag if lock the backoff counter
TX_finish_time = ones(1,N);
TX_ready = zeros(1,N); %to see how many nodes want to TX at a specific slot
numIdleCHN = zeros(1,N); %record number of successive idle CHN when backoff counter is locked
no_use = zeros(1,N);
pl_t = zeros(1,N);
ps_t = zeros(1,N);
PL_colli = zeros(1,N); % pkt loss caused by collisions
%yf������ʼ��
EH_sp = zeros(rap_length,N);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization is over, protocol is started.
% to check every time point during CAP period there are four states for
% nodes. When instant time exceed CAP length, protocol stops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
t = 1;
while ( t<=rap_length ) 
    CHNflag = ones(1,N);  %��ʼ���ŵ���ǣ�1��æ��0�����У�
    isCCA = zeros(1,N); %yf�ж��ǲ��ǵ�һ�ν���
    %% ��һ�鿴ÿ���ڵ㣬�޸����ǵ�״̬
    %����Ƿ�Ҫ�����˱ܼ�����
    for n=1:N 
         %-------------CCA----------------%
%          if ( isFirst==0 &&  Backoff_time(n) == 0 && CSMA_sta(n)==2 ) %���ǵ�һ�����˱ܼ�����Ϊ0���Ͳ���CCA��
%              break;
%          else
%          if (Backoff_time(n) == 1) %yf
           %  if(E_buff(n)>=(E_CCA+E_TX) ) %&&isRAP(n)==1
              %�����ŵ�״̬���CCA
              ELE_ex(n) = ELE_ex(n) + E_CCA;  %���ĵ������ۻ�����
              E_buff(n) = E_buff(n) - E_CCA;   %���ĵ�����
              if(CHN(t)==0) %yf���ʱ϶���ŵ�״��Ϊ����
                   CHNflag(n) = 0;                 %����ŵ�Ϊ���� 
              end
              if(backoff_lock(n)==1)           
                   if( CHNflag(n)==0 && CSMA_sta(n)~=2 )  %�ڵ㴦�ڷ�������״̬ʱ�����п����ŵ��ļ���
                      numIdleCHN(n) = numIdleCHN(n) + 1;  %count the number of successive idle CHN during backoff counter is locked
                   else
                      numIdleCHN(n) = 0;    %reset the number
                   end
                   if(numIdleCHN(n)>=SIFS && (rap_length-t-T_pkt)>=0 )  %channel has been idle for SIFS time and remaining time is enough to send pkt
                       backoff_lock(n) = 0;    %unlock the backoff counter
                       numIdleCHN(n) = 0;    %reset the number
                   end
              end   
              isCCA(n) = 1;
         %    end
             
%          else %yf
%              
%              if(E_buff(n)>=(E_CCA) )%&&isRAP(n)==1
%               %�����ŵ�״̬���CCA
%               ELE_ex(n) = ELE_ex(n) + E_CCA;  %���ĵ������ۻ�����
%               E_buff(n) = E_buff(n) - E_CCA;   %���ĵ�����
%               if(CHN(t)==0) %yf���ʱ϶���ŵ�״��Ϊ����
%                    CHNflag(n) = 0;                 %����ŵ�Ϊ���� 
%               end
%               if(backoff_lock(n)==1)           
%                    if( CHNflag(n)==0 && CSMA_sta(n)~=2 )  %�ڵ㴦�ڷ�������״̬ʱ�����п����ŵ��ļ���
%                       numIdleCHN(n) = numIdleCHN(n) + 1;  %count the number of successive idle CHN during backoff counter is locked
%                    else
%                       numIdleCHN(n) = 0;    %reset the number
%                    end
%                    if(numIdleCHN(n)>=SIFS && (rap_length-t-T_pkt)>=0 )  %channel has been idle for SIFS time and remaining time is enough to send pkt
%                        backoff_lock(n) = 0;    %unlock the backoff counter
%                        numIdleCHN(n) = 0;    %reset the number
%                    end
%               end          
%              end
%              
%          end %yf
%          end %yf
    end
    %���ڵ�״̬
    for n=1:N  
        if( (rap_length-t-T_pkt)>=0 ) % &&isRAP(n)==1  &&E_buff(n)>=(E_TX)
            %%----- case 0:initialize the backoff counter
            if( CSMA_sta(n)==0 ) % 0 state for backoff                    
                if def_time_pre(n) == -1 %�ڵ�û�л�ȡ���κ�RAP�׶ε�ʱ϶                    
                    Backoff_time(n) = randint(1,1,[1,CW(n)]);  
                    backoff_lock(n) = 1;  %lock the backoff counter after reset it
                else % nodes with defer ��һ��֡�е�backoff time��ʣ��
                    Backoff_time(n) = def_time_pre(n); 
                    def_time_pre(n) = -1;
                end % yf�Ѿ�ȷ��def_time_pre =1
                CSMA_sta(n) = 1;
            end

            %%-----------------case 1:�����˱�------------------
            if( CSMA_sta(n)==1&&backoff_lock(n)==0&&isCCA(n)== 1 ) 
                if(Backoff_time(n) == 0) %yf�˱ܵ�1����һ��ʱ϶��0�����뷢��״̬,���뷢��״̬ Backoff_time(n) == 0
                      %һ��ʱ϶������ɵĶ�����1��CCA+�Ѱ����͵�����㣻2��ͨ�����߽����ݰ����ͳ�ȥ�ͽ��ܶ�Ӧ��ACK�������ʱ϶����������
                      if( CHNflag(n)==0) % channel is idle,energy is enough,there has pkt to send                                
                        % now it can be TX_ready���Ƚ��ڵ�����Ϊ׼����������״̬
%                           Backoff_time(n) = 0;%yf�˱ܼ�����Ϊ0�����뷢��״̬
                          TX_ready(n) = 1;    %׼������һ��ʱ϶��ʼռ���ŵ���������                              
                          CHN( t+1:t+T_pkt ) = 1;  %set the CHN of following slot busy                              
                          TX_finish_time(n) =  t + T_pkt;  %���۷��ͳɹ������Ҫ������������֪�������ڼ�ڵ�ʼ�մ��ڷ��ͽ׶�
                          CSMA_sta(n) = 2; % set the state of node to 2
                      else %channel is busy yf,�����⣬һ���
                         backoff_lock(n) = 1;   %CHN is busy ,lock the backoff counter
                         numIdleCHN(n) = 0;    %reset the number of slot with idle CHN
                      end %end CHN                          
                else %�����˱�                         
                     if( CHNflag(n)==0 ) % �ŵ����в�ִ���˱�ʱ���1
                         Backoff_time(n) = Backoff_time(n) - 1;  %�˱�ʱ���1                    
                     else
                         backoff_lock(n) = 1;   %CHN is busy ,lock the backoff counter
                         numIdleCHN(n) = 0;    %reset the number of slot with idle CHN
                     end
                end 
            end

             % --case 2:check the nodes are sending packets--------------                 
            if( CSMA_sta(n)==2 )   
                 if(t == TX_finish_time(n))  %finish TX
                     CSMA_sta(n) = 0;
                 end                                           
            end
        else
            if( (rap_length - t - T_pkt)<0 )
                %�������ʣ���ʱ���޷���ɷ��Ͱ�������,���˱�ʱ�������ȵ���һ��֡������
                backoff_after(n) = Backoff_time(n); %yf,next def_time_pre
                backoff_lock(n) = 1;   %if remaining time is not enough, lock backoff counter
                numIdleCHN(n) = 0;    %reset the number of slot with idle CHN
                CSMA_sta(n) = 0;      %���س�ʼ״̬
            end
        end %end if

    end  %end for
    %% check how many nodes want to send packet at this moment
%     �˲�ͬʱ϶�д��������Ͱ���״̬��ʱ�ڵ�����������1�����ͻ�����͵����ݽ�ʧ�ܣ�������Է������ݡ�
    ind_TX = find(TX_ready==1);
    if( length(ind_TX)>1 )
        %%----------------׼������ʱ���ʧ���������һ���ش�---------------
%         disp('Collision!')
        PL_colli(ind_TX) = PL_colli(ind_TX) + 1;         %��ͻ������1
        ReTX_time(ind_TX) = ReTX_time(ind_TX) + 1;       %�ش�������1 
        ELE_ex(ind_TX) = ELE_ex(ind_TX) + E_TX;%���ĵ������ۻ�����
        E_buff(ind_TX) = E_buff(ind_TX) - E_TX; %���ĵ�����
        %%-----------------�޸ľ�������----------------------------------
        for n=1:length(ind_TX)
            n1 = ind_TX(n);
            if( ReTX_time(n1)>M )           %�ﵽ����ش�������������ǰ��
                pl_t(n1) = pl_t(n1) + 1;    %����������1
                ReTX_time(n1) = 0;           %�ش�������0
                CW(n1) = CWmin(find(UP==UPnode(n1)));  %�޸ľ�������Ϊ�ڵ��Ӧ���ȼ���CWmin 
            else
                %�޸ľ�������ż�����ش����ڼӱ�
                if ( ReTX_time(n1)>0 && mod(ReTX_time(n1),2)==0 )
                    CW(n1) = min(2*CW(n1),CWmax(find(UP==UPnode(n1))));
                end                                    
            end
        end               
    else
        if( length(ind_TX)==1 )
            %%---------�ڵ�ind_TX���Է�������------------------------
 
            CHNb_leng = t + 1 - last_TX_time(ind_TX); %�������һ�η������ݰ������ڵ�ʱ��
            % TX one packet, the channel state when the pkt is finished is
            % recorded and used  3.last slot state of current node 4.the last state after the superframe.
            [ PL_cap,PS_cap,last_CHN_sta(ind_TX), no_use(ind_TX), channelslot(ind_TX)] = pktsend( CHNb_leng,0,last_CHN_sta(ind_TX),1,Pbg(ind_TX),Pgb(ind_TX), channelslot(ind_TX)); 
%             disp(['node ',num2str(ind_TX),' send Pkt ',num2str(PS_cap),' successfully!']);
           
            %%--------------------�޸ķ������------------------------
            CW(n) = CWmin( find(UP==UPnode(n)) );  %���ͳɹ������þ�������
            %������ͳɹ����ڵ�ǰʱ϶�����
            if(PS_cap>0)                                
                last_TX_time(ind_TX) = t + T_pkt;   %��ǰʱ϶�ɹ�������һ����,�������ɵ��Ǹ�ʱ϶
%                 TX_time{ind_TX} = [TX_time{ind_TX},last_TX_time(ind_TX)];
                TX_time_rap(last_TX_time(ind_TX),ind_TX) = 1;
            end
            ELE_ex(ind_TX) = ELE_ex(ind_TX) + E_TX; %���ĵ������ۻ����� 
            E_buff(ind_TX) = E_buff(ind_TX) - E_TX; %���ĵ�����           
             %yf
             pl_t(ind_TX) = pl_t(ind_TX) + PL_cap;  %��pktsend�׶���Ϊ�ŵ�����ԭ�򶪰���Ҳ�����
             ps_t(ind_TX) = ps_t(ind_TX) + PS_cap; 
        end
    end
    TX_ready(ind_TX) = 0;   %���ñ�־λ 
    
    %yfÿ��ʱ϶�������ı�
   % [e_flow,E_buff] = buff_update_new(E_buff);
   % EH_sp(t,:) = e_flow;
    %yf����ʱ϶
    t = t + 1;
    %yf����first
 %   if isFirst == 1
 %   isFirst = 0;
 %   end;
end %end while

TX_time = cell(1,N);
for n=1:N
    ind_TX = find( TX_time_rap(:,n)==1 );
    if( ~isempty(ind_TX) )
        TX_time{1,n} = ind_TX;
    end
  %  E_flow = sum( EH_sp(:,n) );%yfÿ��ʱ϶���������
end
%% update channel condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the CSMA/CA and the TX is over
% calculate the channel state after CAP, i.e., the begining of TDMA, for
% each node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
 for n=1:N
    for g=1:floor( (channelslot(n) + floor (rap_length) - floor (last_TX_time(n)))/statelast )
        if last_CHN_sta(n) == 1
            last_CHN_sta(n) = randsrc(1,1,[0 1;Pgb(n) 1-Pgb(n)]); %%%%%% channel model
        else
            last_CHN_sta(n) = randsrc(1,1,[0 1;1-Pbg(n) Pbg(n)]); %%%%%% using Markov chain
        end
    end
    channelslot(n) = mod( (channelslot(n) + floor (rap_length) - floor (last_TX_time(n)) ), statelast);
      %yf
      %last_CHN_sta(n) = 1;
end % end if
% t5 = toc
    
end % function end

