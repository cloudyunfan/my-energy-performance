function [ReTX_time,backoff_after,CSMA_sta,pl_t,ps_t,PL_colli,TX_time] = slotCSMACA_unsat_newnoE(rap_length,CSMA_sta,def_time_pre,last_CHN_sta,ReTX_times_pre,CW,last_TX_time,E_buff) %act,,B_buff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ELE_ex(倒数第四),E_buff,E_flow（最后两个）
%%% CSMA/CA transmission under unsaturation condition
% Input:
%     2.rap_length: duration of RAP
%     3.CSMA_sta: last CSMA state in pre superframe
%     4.def_time_pre: backoff counter before superframe(-1: node has not obtained any slot before or just send packet successfully)
%     5.last_CHN_sta: Channel state (0: idle,1:busy)
%     6.ReTX_times_pre: last retransmission times at the end of RAP
%     7.CW: contention window of every node
%     8.last_TX_time： last time sending packet successfully
%     10:E_buff:能量缓存区
%     11:e_flow:能量流
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

global Data_rate Pkt_len Tslot Pbg Pgb CWmin CWmax UP UPnode E_TX E_CCA channelslot statelast%isMAP isRAP 
%yf probability of arrive one unit energy in each slot
global P1_x
%-----------参数---------------------------------------------------------
%yf
P1_x = 0.6;
SIFS = 0.5;      %最短帧间间隔0 yf 有点问题，一会儿改0.5
T_pkt = ceil( Pkt_len/(Data_rate*Tslot) );  %time to send a packet and receive ACK,unit slot
M = 4;        %最大重传次数
N = length(E_buff);  %获取使用CSMA/CA的节点数

ReTX_time = ReTX_times_pre;  %读取上一超帧中的重传次数记录
TX_time_rap = zeros(rap_length,N); % the time when node send pkt successfully

Backoff_time = zeros(1,N); % temporal variable to record backoff time
CHN = zeros(1,rap_length); % 初始化每个时隙的信道状况为空闲
ELE_ex = zeros(1,N); % 记录本超帧消耗的总能量
backoff_after = zeros(1,N); % deferred backoff should be output,defoult 0
backoff_lock = zeros(1,N);  % flag if lock the backoff counter
TX_finish_time = ones(1,N);
TX_ready = zeros(1,N); %to see how many nodes want to TX at a specific slot
numIdleCHN = zeros(1,N); %record number of successive idle CHN when backoff counter is locked
no_use = zeros(1,N);
pl_t = zeros(1,N);
ps_t = zeros(1,N);
PL_colli = zeros(1,N); % pkt loss caused by collisions
%yf能量初始化
EH_sp = zeros(rap_length,N);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization is over, protocol is started.
% to check every time point during CAP period there are four states for
% nodes. When instant time exceed CAP length, protocol stops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
t = 1;
while ( t<=rap_length ) 
    CHNflag = ones(1,N);  %初始化信道标记（1：忙；0：空闲）
    isCCA = zeros(1,N); %yf判断是不是第一次进入
    %% 逐一查看每个节点，修改他们的状态
    %检查是否要解锁退避计数器
    for n=1:N 
         %-------------CCA----------------%
%          if ( isFirst==0 &&  Backoff_time(n) == 0 && CSMA_sta(n)==2 ) %不是第一次且退避计数器为0，就不用CCA了
%              break;
%          else
%          if (Backoff_time(n) == 1) %yf
           %  if(E_buff(n)>=(E_CCA+E_TX) ) %&&isRAP(n)==1
              %进行信道状态检查CCA
              ELE_ex(n) = ELE_ex(n) + E_CCA;  %消耗的能量累积记下
              E_buff(n) = E_buff(n) - E_CCA;   %消耗掉能量
              if(CHN(t)==0) %yf如果时隙的信道状况为空闲
                   CHNflag(n) = 0;                 %标记信道为空闲 
              end
              if(backoff_lock(n)==1)           
                   if( CHNflag(n)==0 && CSMA_sta(n)~=2 )  %节点处于发送数据状态时不进行空闲信道的计数
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
%               %进行信道状态检查CCA
%               ELE_ex(n) = ELE_ex(n) + E_CCA;  %消耗的能量累积记下
%               E_buff(n) = E_buff(n) - E_CCA;   %消耗掉能量
%               if(CHN(t)==0) %yf如果时隙的信道状况为空闲
%                    CHNflag(n) = 0;                 %标记信道为空闲 
%               end
%               if(backoff_lock(n)==1)           
%                    if( CHNflag(n)==0 && CSMA_sta(n)~=2 )  %节点处于发送数据状态时不进行空闲信道的计数
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
    %检查节点状态
    for n=1:N  
        if( (rap_length-t-T_pkt)>=0 ) % &&isRAP(n)==1  &&E_buff(n)>=(E_TX)
            %%----- case 0:initialize the backoff counter
            if( CSMA_sta(n)==0 ) % 0 state for backoff                    
                if def_time_pre(n) == -1 %节点没有获取过任何RAP阶段的时隙                    
                    Backoff_time(n) = randint(1,1,[1,CW(n)]);  
                    backoff_lock(n) = 1;  %lock the backoff counter after reset it
                else % nodes with defer 上一超帧中的backoff time有剩余
                    Backoff_time(n) = def_time_pre(n); 
                    def_time_pre(n) = -1;
                end % yf已经确定def_time_pre =1
                CSMA_sta(n) = 1;
            end

            %%-----------------case 1:进入退避------------------
            if( CSMA_sta(n)==1&&backoff_lock(n)==0&&isCCA(n)== 1 ) 
                if(Backoff_time(n) == 0) %yf退避到1，下一个时隙到0，进入发送状态,进入发送状态 Backoff_time(n) == 0
                      %一个时隙可以完成的动作：1、CCA+把包发送到物理层；2、通过天线将数据包发送出去和接受对应的ACK（这个是时隙的整数倍）
                      if( CHNflag(n)==0) % channel is idle,energy is enough,there has pkt to send                                
                        % now it can be TX_ready，先将节点设置为准备发送数据状态
%                           Backoff_time(n) = 0;%yf退避计数器为0，进入发送状态
                          TX_ready(n) = 1;    %准备在下一个时隙开始占用信道发送数据                              
                          CHN( t+1:t+T_pkt ) = 1;  %set the CHN of following slot busy                              
                          TX_finish_time(n) =  t + T_pkt;  %无论发送成功与否，需要到发送完后才能知道，这期间节点始终处于发送阶段
                          CSMA_sta(n) = 2; % set the state of node to 2
                      else %channel is busy yf,有问题，一会儿
                         backoff_lock(n) = 1;   %CHN is busy ,lock the backoff counter
                         numIdleCHN(n) = 0;    %reset the number of slot with idle CHN
                      end %end CHN                          
                else %继续退避                         
                     if( CHNflag(n)==0 ) % 信道空闲才执行退避时间减1
                         Backoff_time(n) = Backoff_time(n) - 1;  %退避时间减1                    
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
                %如果后面剩余的时间无法完成发送包的任务,则退避时间锁定等到下一超帧再启动
                backoff_after(n) = Backoff_time(n); %yf,next def_time_pre
                backoff_lock(n) = 1;   %if remaining time is not enough, lock backoff counter
                numIdleCHN(n) = 0;    %reset the number of slot with idle CHN
                CSMA_sta(n) = 0;      %返回初始状态
            end
        end %end if

    end  %end for
    %% check how many nodes want to send packet at this moment
%     核查同时隙中处即将发送包的状态的时节点数，若大于1，则冲突，发送的数据将失败，否则可以发送数据。
    ind_TX = find(TX_ready==1);
    if( length(ind_TX)>1 )
        %%----------------准备传输时如果失败则进入下一次重传---------------
%         disp('Collision!')
        PL_colli(ind_TX) = PL_colli(ind_TX) + 1;         %冲突次数加1
        ReTX_time(ind_TX) = ReTX_time(ind_TX) + 1;       %重传次数加1 
   %     ELE_ex(ind_TX) = ELE_ex(ind_TX) + E_TX;%消耗的能量累积记下
   %     E_buff(ind_TX) = E_buff(ind_TX) - E_TX; %消耗掉能量
        %%-----------------修改竞争窗口----------------------------------
        for n=1:length(ind_TX)
            n1 = ind_TX(n);
            if( ReTX_time(n1)>M )           %达到最大重传次数，丢弃当前包
                pl_t(n1) = pl_t(n1) + 1;    %丢包次数加1
                ReTX_time(n1) = 0;           %重传次数归0
                CW(n1) = CWmin(find(UP==UPnode(n1)));  %修改竞争窗口为节点对应优先级的CWmin 
            else
                %修改竞争窗，偶数次重传窗口加倍
                if ( ReTX_time(n1)>0 && mod(ReTX_time(n1),2)==0 )
                    CW(n1) = min(2*CW(n1),CWmax(find(UP==UPnode(n1))));
                end                                    
            end
        end               
    else
        if( length(ind_TX)==1 )
            %%---------节点ind_TX可以发送数据------------------------
 
            CHNb_leng = t + 1 - last_TX_time(ind_TX); %计算从上一次发送数据包到现在的时间
            % TX one packet, the channel state when the pkt is finished is
            % recorded and used  3.last slot state of current node 4.the last state after the superframe.
            [ PL_cap,PS_cap,last_CHN_sta(ind_TX), no_use(ind_TX), channelslot(ind_TX)] = pktsend( CHNb_leng,0,last_CHN_sta(ind_TX),1,Pbg(ind_TX),Pgb(ind_TX), channelslot(ind_TX)); 
%             disp(['node ',num2str(ind_TX),' send Pkt ',num2str(PS_cap),' successfully!']);
           
            %%--------------------修改仿真变量------------------------
            CW(n) = CWmin( find(UP==UPnode(n)) );  %发送成功，重置竞争窗口
            %如果发送成功则在当前时隙做标记
            if(PS_cap>0)                                
                last_TX_time(ind_TX) = t + T_pkt;   %当前时隙成功发送了一个包,标记在完成的那个时隙
%                 TX_time{ind_TX} = [TX_time{ind_TX},last_TX_time(ind_TX)];
                TX_time_rap(last_TX_time(ind_TX),ind_TX) = 1;
            end
           % ELE_ex(ind_TX) = ELE_ex(ind_TX) + E_TX;%消耗的能量累积记下 
           % E_buff(ind_TX) = E_buff(ind_TX) - E_TX; %消耗掉能量           
             %yf
             pl_t(ind_TX) = pl_t(ind_TX) + PL_cap;  %在pktsend阶段因为信道条件原因丢包了也算进来
             ps_t(ind_TX) = ps_t(ind_TX) + PS_cap; 
        end
    end
    TX_ready(ind_TX) = 0;   %重置标志位 
    
    %yf每个时隙的能量改变
   % [e_flow,E_buff] = buff_update_new(E_buff);
   % EH_sp(t,:) = e_flow;
    %yf更新时隙
    t = t + 1;
    %yf更新first
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
  %  E_flow = sum( EH_sp(:,n) );%yf每个时隙到达的能量
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

