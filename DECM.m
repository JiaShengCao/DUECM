clear;
clc;
%DECM�㷨������HEFT
HEFT_new;
D_G=100;
n_exit=10;

LB_G=EFT_min(N);
DS_G=D_G-LB_G;
AFT=EFT_min;%ĳ�������ʵ�����ʱ��
LB=AFT;

%%
%����L��ni��,D(ni)
for i=1:1:N
    if i==1
        L(i)=1;
    end
    father=find(C(:,i)~=0);
    [num,temp]=size(father);
    l_temp=0;
    for j=1:1:num
        l_temp(j)=L(father(j));
    end
    L(i)=max(l_temp)+1;
end

for i=1:1:N
 D(i)=LB(i)+DS_G/L(n_exit)*L(i);
end



%%
%������ʼ��
p_ind(1)=0.03;
p_ind(2)=0.04;
p_ind(3)=0.07;
c_ef(1)=0.8;
c_ef(2)=0.8;
c_ef(3)=1.0;
mm(1)=2.9;
mm(2)=2.5;
mm(3)=2.5;
f_low(1)=0.22;
f_low(2)=0.21;
f_low(3)=0.29;
f_max(1)=1.0;
f_max(2)=1.0;
f_max(3)=1.0;


%%
%��ʼ����һ���ڵ�
%%%%%%%���㷨ʹ�õ�EFT2,EST2�ȣ���Ϊ�˱���workspacec�е�һ�����н���ĸ���
%%%%%%%�ҵ���һ���ڵ�ֵ���Ƶ���Լ��ܺ�
u=node_inprocess(prior(1));
ff=f_low(u):0.01:f_max(u);
[temp,num]=size(ff);
E(prior(1))=999;
 for i=1:1:N_p
    EST2(prior(1),i)=0;
    EFT2(prior(1),i)=EST2(prior(1),i)+W(prior(1),i);
 end
 EFT_min2(prior(1))=min(EFT2(prior(1),:));
for j=num:-1:1
        EFT_t(prior(1),j)=EST2(prior(1),u)+W(prior(1),u)*f_max(u)/ff(j);
        if EFT_t(prior(1),j)<=D(prior(1))
            E_t(prior(1),j)=(p_ind(u)+c_ef(u)*(ff(j)^mm(u)))*W(prior(1),u)*f_max(u)/ff(j);
            if E_t(prior(1),j)<E(prior(1))
                f(prior(1))=ff(j);
                AFT(prior(1))= EFT_t(prior(1),j);
                E(prior(1))=E_t(prior(1),j);
                EFT_min2(prior(1))=EFT_t(prior(1),j);
            end
        end
end
index2=find(min(EFT2(prior(1),:))==EFT2(prior(1),:));
now2(index2)=prior(1);%ĳ�����̵����һ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for i=2:1:N
    f(prior(i))=0;
    AFT(prior(i))=0;
    E(prior(i))=999;
    u=node_inprocess(prior(i));
    
    %%
    %ÿ��ѭ���µĽڵ�֮ǰ����Ҫ���¼���һ��֮ǰ�ڵ��EST
    father=find(C(:,prior(i))~=0);%�ýڵ�ĸ��ڵ�
    [num_father,temp]=size(father);
    
    %ÿ������
    for m=1:1:N_p
        
        for k=1:1:num_father
            if m==node_inprocess(father(k))
                temp(m,k)=EFT_min2(father(k));
            else
                temp(m,k)=EFT_min2(father(k))+C(father(k),prior(i));
            end
        end
        
        if now2(m)==0
            EST2(prior(i),m)=max(temp(m,:));
            EFT2(prior(i),m)=EST2(prior(i),m)+W(prior(i),m);
        else
            EST2(prior(i),m)=max([EFT_min2(now2(m)),max(temp(m,:))]);
            EFT2(prior(i),m)=EST2(prior(i),m)+W(prior(i),m);
        end
    end

    %%
    %ff��Ƶ�ʵļ���
    ff=f_low(u):0.01:f_max(u);
    [temp,num]=size(ff);
    
    %%%%%%%%%%%%%%%    Ѱ��ÿ���ڵ��ܺ���С��Ƶ��    %%%%%%%%%%%%%%%%%%%%%%%%
    for j=num:-1:1
        EFT_t(prior(i),j)=EST2(prior(i),u)+W(prior(i),u)*f_max(u)/ff(j);
        if EFT_t(prior(i),j)<=D(prior(i))
            E_t(prior(i),j)=(p_ind(u)+c_ef(u)*(ff(j)^mm(u)))*W(prior(i),u)*f_max(u)/ff(j);
            if E_t(prior(i),j)<E(prior(i))
                f(prior(i))=ff(j);
                AFT(prior(i))= EFT_t(prior(i),j);
                EFT_min2(prior(i))=EFT_t(prior(i),j);
                E(prior(i))=E_t(prior(i),j);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    AST(prior(i))=AFT(prior(i))-f_max(u)*W(prior(i),u)/f(prior(i));
    now2(node_inprocess(prior(i)))=prior(i);
    
end

%%
E_sum=sum(E);%�������������ܺ�
SL_G=EFT_min2(n_exit);%�������ִ�г���













