function indMAP = getPollNode(IndNode_MAP,M)
% Input:
%     IndNode_MAP:index of nodes which can us MAP
%     M: max poll times in MAP
% Output:
%     indMAP: index of nodes which will be poll in MAP
indMAP = [];
if(length(IndNode_MAP)<=M)
   indMAP =  IndNode_MAP;  %节点数小于M，保持原来的下标返回即可
else
    i=1;
    while (i<=M )
       ind = randint(1,1,[1,length(IndNode_MAP)]); %随机选取一个节点poll 
       %检查是否已被poll过
       flag = 0;
%        for j=1:length(indMAP)
%           if( indMAP(j)==ind)
%              flag = 1;               
%           end
%        end
       if(flag==0)
           indMAP = [indMAP ind];
           i = i + 1;
       end       
    end
end

end %end function