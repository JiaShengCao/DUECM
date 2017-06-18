clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=10;    %Ni 节点数目
%W     %W(i,j) 节点i在进程j执行的通信花费
%C     %Cij 节点i到j的通信花费
%P     %Pi 进程的集合;
N_p=3;  %进程的数目
C=zeros(N,N);

%%
%设置输入的图
C(1,2)=18;
C(1,3)=12;
C(1,4)=9;
C(1,5)=11;
C(1,6)=14;
C(2,8)=19;
C(2,9)=16;
C(3,7)=23;
C(4,8)=27;
C(4,9)=23;
C(5,9)=13;
C(6,8)=15;
C(7,10)=17;
C(8,10)=11;
C(9,10)=13;

W(1,1)=14;
W(1,2)=16;
W(1,3)=9;
W(2,1)=13;
W(2,2)=19;
W(2,3)=18;
W(3,1)=11;
W(3,2)=13;
W(3,3)=19;
W(4,1)=13;
W(4,2)=8;
W(4,3)=17;
W(5,1)=12;
W(5,2)=13;
W(5,3)=10;
W(6,1)=13;
W(6,2)=16;
W(6,3)=9;
W(7,1)=7;
W(7,2)=15;
W(7,3)=11;
W(8,1)=5;
W(8,2)=11;
W(8,3)=14;
W(9,1)=18;
W(9,2)=12;
W(9,3)=20;
W(10,1)=21;
W(10,2)=7;
W(10,3)=16;
%画出任务DAG
 draw_graph(C);
 for i=1:1:N
     W_ave(i)=sum(W(i,:))/N_p;
 end
 W_ave=W_ave';
 
 
 %%
 
 %从exit节点开始计算ranku
 ranku(N)=W_ave(N,1);
 for i=N-1:-1:1
     %先找节点i的下一个节点,next是下一个节点的集合
     next=find(C(i,:)~=0);
     %temp=C(i,:)+ranku([next]);
     [num1,num2]=size(next);
     for m=1:1:num2
         temp(m)=C(i,next(m))+ranku(next(m));
     end
     max_num=max(temp);
     ranku(i)=W_ave(i)+max_num;
 end
 
 %降序排列，prior为优先级从高到低的节点。
 [result,prior]=sort(ranku,'descend');
 
 %%
 %初始化enter节点
 for i=1:1:N_p
    EST(prior(1),i)=0;
    EFT(prior(1),i)=EST(prior(1),i)+W(prior(1),i);
    process(i).member=[];
 end
    EFT_min(prior(1))=min(EFT(prior(1),:));
    index=find(min(EFT(prior(1),:))==EFT(prior(1),:));

    now(index)=prior(1);%某个进程的最后一个任务
    node_inprocess(prior(1))=index;%某个节点所在的进程
    process(index).member=[process(index).member,now(index)];%某个进程的所有的任务
    
    
    %遍历每一个节点
 for i=2:1:N         
    father=find(C(:,prior(i))~=0);%该节点的父节点
    [num_father,temp]=size(father);
    
    %每个进程
    for m=1:1:N_p
        
        for k=1:1:num_father
            if m==node_inprocess(father(k))
                temp(m,k)=EFT_min(father(k));
            else
                temp(m,k)=EFT_min(father(k))+C(father(k),prior(i));
            end
        end
        
        if now(m)==0
            EST(prior(i),m)=max(temp(m,:));
            EFT(prior(i),m)=EST(prior(i),m)+W(prior(i),m);
        else
            EST(prior(i),m)=max([EFT_min(now(m)),max(temp(m,:))]);
            EFT(prior(i),m)=EST(prior(i),m)+W(prior(i),m);
        end
    end
        
        EFT_min(prior(i))=min(EFT(prior(i),:));%找每个节点最小的EFT
        node_inprocess(prior(i))=find(min(EFT(prior(i),:))==EFT(prior(i),:));%找到其所在的进程
        now(node_inprocess(prior(i)))=prior(i);
        process(node_inprocess(prior(i))).member=[process(node_inprocess(prior(i))).member,prior(i)];%某个进程的所有的任务
        
 end
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 