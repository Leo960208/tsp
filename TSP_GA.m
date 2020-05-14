function GATSP_PPT
%D是距离矩阵
%popsize为种群个数
%gene为停止代数，遗传到第gene代时程序停止,gene的具体取值视问题的规模和耗费的时间而定
%acc为适应值归一化淘汰加速指数,最好取为1,2,3,4,不宜太大
%alpha为淘汰保护指数,可取为0~1之间任意小数,取1时关闭保护功能,建议取0.8~1.0之间的值
popsize=1000;
gene=800;
ps=1;
pc=0.65;%交叉率
pm=0.2;%变异率
position=pos;%调入坐标
D=P2D(position);
[R,Rlength]=geneticTSP(D,position,popsize,gene,ps,pc,pm);

%%程序一
%best为最短路径,bestlength为路径长度
function [best,bestlength]=geneticTSP(D,a,n,C,ps,pc,pm)
[N,NN]=size(D);
pop=zeros(n,N);%用于存储种群
for i=1:n
    pop(i,:)=randperm(N);%随机生成初始种群
end
best=pop(1,:);
figure(1)
subplot(1,3,1)
scatter(a(:,1),a(:,2),'ro')
pause(1)
subplot(1,3,2)
plotaiwa(a,best)
pause(1)
len=zeros(n,1);%存储路径长度
fitness=zeros(n,1);%存储归一化适应值
counter=0;
while counter<C
    for i=1:n
        len(i,1)=myLength(D,pop(i,:));%计算路径长度
    end
    maxlen=max(len);
    minlen=min(len);
    fitness=fit(len,maxlen,minlen);%计算归一化适应值
    rr=find(len==minlen);
    best=pop(rr(1,1),:);%更新最短路径
    
    [newpop]=selection(pop,fitness,ps);%选择
    [newpop]=crossover(newpop,pc); %交叉
    [newpop]=mutation(newpop,pm); %变异

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


%%程序二：计算邻接矩阵
%输入参数a是中国31个城市的坐标
%输出参数D是无向图的赋权邻接矩阵
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


%%程序三：计算归一化适应值
%计算归一化适应值的子程序
function fitness=fit(len,maxlen,minlen)
fitness=len;
for i=1:length(len)
    fitness(i,1)=(1-((len(i,1)-minlen)/(maxlen-minlen+0.0001))).^2;
    %fitness(i,1)=1/len(i,1);
end

%选择子程序
function [newpop]=selection(pop,fitvalue,ps)
totalfit=sum(fitvalue); %求适应值之和
fitvalue=fitvalue/totalfit; %单个个体被选择的概率
fitvalue=cumsum(fitvalue); %如 fitvalue=[1 2 3 4]，则 cumsum(fitvalue)=[1 3 6 10] 
[px,py]=size(pop);
ms=sort(rand(px,1)); %从小到大排列
fitin=1;
newin=1;
while newin<=px %蒙特卡洛方法抽样
    if(ms(newin))<fitvalue(fitin)
        newpop(newin,:)=pop(fitin,:);
        newin=newin+1;
    else
        fitin=fitin+1;
    end
end


%交叉算法采用的是由Goldberg和Lingle于1985年提出的PMX(部分匹配交叉)
function [newpop]=crossover(pop,pc) %交叉
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
p=unidrnd(L-W+1);%随机选择交叉范围，从p到p+W
for i=1:W  %交叉
    x=find(a==b(1,p+i-1));
    y=find(b==a(1,p+i-1));
    [a(1,p+i-1),b(1,p+i-1)]=exchange(a(1,p+i-1),b(1,p+i-1));
    [a(1,x),b(1,y)]=exchange(a(1,x),b(1,y));   
end
function [x,y]=exchange(x,y)
temp=x;
x=y;
y=temp;


%变异
function [newpop]=mutation(pop,pm)
[px,py]=size(pop);
newpop=pop;
for i=1:px
    if(rand<pm)
        mpoint=randperm(py);
        newpop(i,mpoint(1):mpoint(2))=fliplr(newpop(i,mpoint(1):mpoint(2)));
    end
end

%计算路径
%该路径长度是一个闭合的路径的长度
function len=myLength(D,p)
[N,NN]=size(D);
len=D(p(1,N),p(1,1));
for i=1:(N-1)
    len=len+D(p(1,i),p(1,i+1));
end

%%绘制路径示意图
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

%中国31城市坐标数据
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