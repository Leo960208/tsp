function GATSP_PPT
%D�Ǿ������
%popsizeΪ��Ⱥ����
%geneΪֹͣ�������Ŵ�����gene��ʱ����ֹͣ,gene�ľ���ȡֵ������Ĺ�ģ�ͺķѵ�ʱ�����
%accΪ��Ӧֵ��һ����̭����ָ��,���ȡΪ1,2,3,4,����̫��
%alphaΪ��̭����ָ��,��ȡΪ0~1֮������С��,ȡ1ʱ�رձ�������,����ȡ0.8~1.0֮���ֵ
popsize=1000;
gene=800;
ps=1;
pc=0.65;%������
pm=0.2;%������
position=pos;%��������
D=P2D(position);
[R,Rlength]=geneticTSP(D,position,popsize,gene,ps,pc,pm);

%%����һ
%bestΪ���·��,bestlengthΪ·������
function [best,bestlength]=geneticTSP(D,a,n,C,ps,pc,pm)
[N,NN]=size(D);
pop=zeros(n,N);%���ڴ洢��Ⱥ
for i=1:n
    pop(i,:)=randperm(N);%������ɳ�ʼ��Ⱥ
end
best=pop(1,:);
figure(1)
subplot(1,3,1)
scatter(a(:,1),a(:,2),'ro')
pause(1)
subplot(1,3,2)
plotaiwa(a,best)
pause(1)
len=zeros(n,1);%�洢·������
fitness=zeros(n,1);%�洢��һ����Ӧֵ
counter=0;
while counter<C
    for i=1:n
        len(i,1)=myLength(D,pop(i,:));%����·������
    end
    maxlen=max(len);
    minlen=min(len);
    fitness=fit(len,maxlen,minlen);%�����һ����Ӧֵ
    rr=find(len==minlen);
    best=pop(rr(1,1),:);%�������·��
    
    [newpop]=selection(pop,fitness,ps);%ѡ��
    [newpop]=crossover(newpop,pc); %����
    [newpop]=mutation(newpop,pm); %����

    pop=newpop;
    counter=counter+1   
    figure(2)
    plot(counter,minlen,'b.');
    hold on;
    plot(counter,mean(len),'r.');
end
bestlength=myLength(D,best);
figure(1)
subplot(1,3,3)
plotaiwa(a,best)


%%������������ڽӾ���
%�������a���й�31�����е�����
%�������D������ͼ�ĸ�Ȩ�ڽӾ���
function D=P2D(a)
[c,d]=size(a);
D=zeros(c,c);
for i=1:c
    for j=i:c
        bb=(a(i,1)-a(j,1)).^2+(a(i,2)-a(j,2)).^2;
        D(i,j)=bb^(0.5);
        D(j,i)=D(i,j);
    end
end


%%�������������һ����Ӧֵ
%�����һ����Ӧֵ���ӳ���
function fitness=fit(len,maxlen,minlen)
fitness=len;
for i=1:length(len)
    fitness(i,1)=(1-((len(i,1)-minlen)/(maxlen-minlen+0.0001))).^2;
    %fitness(i,1)=1/len(i,1);
end

%ѡ���ӳ���
function [newpop]=selection(pop,fitvalue,ps)
totalfit=sum(fitvalue); %����Ӧֵ֮��
fitvalue=fitvalue/totalfit; %�������屻ѡ��ĸ���
fitvalue=cumsum(fitvalue); %�� fitvalue=[1 2 3 4]���� cumsum(fitvalue)=[1 3 6 10] 
[px,py]=size(pop);
ms=sort(rand(px,1)); %��С��������
fitin=1;
newin=1;
while newin<=px %���ؿ��巽������
    if(ms(newin))<fitvalue(fitin)
        newpop(newin,:)=pop(fitin,:);
        newin=newin+1;
    else
        fitin=fitin+1;
    end
end


%�����㷨���õ�����Goldberg��Lingle��1985�������PMX(����ƥ�佻��)
function [newpop]=crossover(pop,pc) %����
N=size(pop,1);
newpop=pop;
for i=1:2:N-1
    if(rand<pc)
        [newpop(i,:),newpop(i+1,:)]=intercross(pop(i,:),pop(i+1,:));
    end
end
function [a,b]=intercross(a,b)
L=length(a);
W=randperm(L);
W=W(1);
p=unidrnd(L-W+1);%���ѡ�񽻲淶Χ����p��p+W
for i=1:W  %����
    x=find(a==b(1,p+i-1));
    y=find(b==a(1,p+i-1));
    [a(1,p+i-1),b(1,p+i-1)]=exchange(a(1,p+i-1),b(1,p+i-1));
    [a(1,x),b(1,y)]=exchange(a(1,x),b(1,y));   
end
function [x,y]=exchange(x,y)
temp=x;
x=y;
y=temp;


%����
function [newpop]=mutation(pop,pm)
[px,py]=size(pop);
newpop=pop;
for i=1:px
    if(rand<pm)
        mpoint=randperm(py);
        newpop(i,mpoint(1):mpoint(2))=fliplr(newpop(i,mpoint(1):mpoint(2)));
    end
end

%����·��
%��·��������һ���պϵ�·���ĳ���
function len=myLength(D,p)
[N,NN]=size(D);
len=D(p(1,N),p(1,1));
for i=1:(N-1)
    len=len+D(p(1,i),p(1,i+1));
end

%%����·��ʾ��ͼ
function plotaiwa(a,R)
scatter(a(:,1),a(:,2),'x')
hold on
plot([a(R(1),1),a(R(31),1)],[a(R(1),2),a(R(31),2)])
hold on
for i=2:length(R)
    x0=a(R(i-1),1);
    y0=a(R(i-1),2);
    x1=a(R(i),1);
    y1=a(R(i),2);
    xx=[x0,x1];
    yy=[y0,y1];
    plot(xx,yy)
    hold on
end

%�й�31������������
function position=pos
position =...
   [1304,2312;
    3639,1315;
    4177,2244;
    3712,1399;
    3488,1535;
    3326,1556;
    3238,1229;
    4196,1004;
    4312,790;
    4386,570;
    3007,1970;
    2562,1756;
    2788,1491;
    2381,1676;
    1332,695;
    3715,1678;
    3918,2179;
    4061,2370;
    3780,2212;
    3676,2578;
    4029,2838;
    4263,2931;
    3429,1908;
    3507,2367;
    3394,2643;
    3439,3201;
    2935,3240;
    3140,3550;
    2545,2357;
    2778,2826;
    2370,2975;]