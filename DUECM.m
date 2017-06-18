DECM;
%queue是升序排列AFT的节点队列
[result,queue]=sort(AFT);
temp=[];
for i=1:1:N
    %%
    u=node_inprocess(queue(i));
    %计算LFT
    if queue(i)==n_exit
        LFT(queue(i))=D_G;
    else
        %寻找该节点的子节点
        son=find(C(queue(i),:)~=0);
        [num1,num2]=size(son);
        temp=[];
        for j=1:1:num2
            u2=node_inprocess(son(j));
            if u==u2
                temp(j)=AST(son(j));
            else
                temp(j)=AST(son(j))-C(queue(i),son(j));
            end
        end
            index=find(process(u).member==queue(i));
            [a,b]=size(process(u).member);
        if  index~=b   %判断该节点存在后继结点
            index=index+1;
            after=process(u).member(index);
            LFT(queue(i))=min([min(temp),AST(after)]);
        else
            LFT(queue(i))=min(temp);
        end
    end
    %%
    %计算新的f
    f2(queue(i))=f(queue(i))*(AFT(queue(i))-AST(queue(i)))/(LFT(queue(i))-AST(queue(i)));
    f2(queue(i))=max([f2(queue(i)),f_low(u)]);
    f2(queue(i))=round(f2(queue(i))*100)/100;%保留两位小数
    %%
    %更新 AET,AFT,AST
    AET(queue(i))=f(queue(i))*(LFT(queue(i))-AST(queue(i)))/f2(queue(i));
    AFT(queue(i))=LFT(queue(i));
    AST(queue(i))=AFT(queue(i))-AET(queue(i));
    AST(queue(i))=round(AST(queue(i))*10^4)/10^4;%保留四位小数
    E_duecm(queue(i))=(p_ind(u)+c_ef(u)*(f2(queue(i))^mm(u)))*W(queue(i),u)*f_max(u)/f2(queue(i));
end
    E_total2=sum(E_duecm);