DECM;
%queue����������AFT�Ľڵ����
[result,queue]=sort(AFT);
temp=[];
for i=1:1:N
    %%
    u=node_inprocess(queue(i));
    %����LFT
    if queue(i)==n_exit
        LFT(queue(i))=D_G;
    else
        %Ѱ�Ҹýڵ���ӽڵ�
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
        if  index~=b   %�жϸýڵ���ں�̽��
            index=index+1;
            after=process(u).member(index);
            LFT(queue(i))=min([min(temp),AST(after)]);
        else
            LFT(queue(i))=min(temp);
        end
    end
    %%
    %�����µ�f
    f2(queue(i))=f(queue(i))*(AFT(queue(i))-AST(queue(i)))/(LFT(queue(i))-AST(queue(i)));
    f2(queue(i))=max([f2(queue(i)),f_low(u)]);
    f2(queue(i))=round(f2(queue(i))*100)/100;%������λС��
    %%
    %���� AET,AFT,AST
    AET(queue(i))=f(queue(i))*(LFT(queue(i))-AST(queue(i)))/f2(queue(i));
    AFT(queue(i))=LFT(queue(i));
    AST(queue(i))=AFT(queue(i))-AET(queue(i));
    AST(queue(i))=round(AST(queue(i))*10^4)/10^4;%������λС��
    E_duecm(queue(i))=(p_ind(u)+c_ef(u)*(f2(queue(i))^mm(u)))*W(queue(i),u)*f_max(u)/f2(queue(i));
end
    E_total2=sum(E_duecm);