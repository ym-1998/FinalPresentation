%% 随机优化案例
clc;
%% 随机优化模型+有CHP机组+电热锅炉+储能技术---下层包括电力市场和天然气市场
Input;
tic;%将建模的时间包括在内
%% 决策变量
disp("创建上层问题决策变量...");
%CHP机组
Pchp=sdpvar(numCluster,Time);%CHP机组电功率/MW
Hchp=sdpvar(numCluster,Time);%CHP机组热功率/MW
Ichp=binvar(numCluster,Time);%CHP机组状态
SUchp=sdpvar(numCluster,Time);%CHP机组开启成本/$
SDchp=sdpvar(numCluster,Time);%CHP机组关停成本/$
GCchp=sdpvar(numCluster,Time);%CHP机组消耗天然气量/MWh
%电热锅炉机组
Peb=sdpvar(numCluster,Time);%电功率/MWh
Heb=sdpvar(numCluster,Time);%热功率/MWh
%新能源机组
Pwt=sdpvar(numCluster,Time);%风电上网功率，MW
%电储能
Ipch=binvar(numCluster,Time);%储能充电状态
Ipdis=binvar(numCluster,Time);%储能放电状态
Pch=sdpvar(numCluster,Time);%储能充电功率/MW
Pdis=sdpvar(numCluster,Time);%储能放电功率/MW
Ap=sdpvar(numCluster,Time);%储能电量/MWh
%热储能
Ihch=binvar(numCluster,Time);%储能充电状态
Ihdis=binvar(numCluster,Time);%储能放电状态
Hch=sdpvar(numCluster,Time);%储能充电功率/MW
Hdis=sdpvar(numCluster,Time);%储能放电功率/MW
Ah=sdpvar(numCluster,Time);%储能电量/MWh
%气储能
Igch=binvar(numCluster,Time);%储能充电状态
Igdis=binvar(numCluster,Time);%储能放电状态
Gch=sdpvar(numCluster,Time);%储能充电功率/MW
Gdis=sdpvar(numCluster,Time);%储能放电功率/MW
Ag=sdpvar(numCluster,Time);%储能电量/MWh
%购电量
Pmeo=sdpvar(Time,1);%购能量，MWh
%报价参数
Cmeo=sdpvar(Time,1);%报价，$/MWh
%收益各项参数
T1=sdpvar(numCluster,Time);%售能收入，$
T2=sdpvar(numCluster,Time);%储能运维成本，$
Z=sdpvar(numCluster,1);%各场景的收益，$
Obj_MPEC=sdpvar(1,1);%最终的期望收益，$
disp("创建下层电力市场问题决策变量...");
PG=sdpvar(numGenerators,Time);%发电量,MW
Delta=sdpvar(numBuses,Time);%节点相角，rad
disp("创建下层电力市场问题对偶问题决策变量...");
Lanta=sdpvar(numBuses,Time);%功率平衡约束的对偶变量
MiuMax=sdpvar(numGenerators,Time);%发电机最大输出功率对偶变量
MiuMin=sdpvar(numGenerators,Time);%发电机最小输出功率对偶变量
Vmin=sdpvar(numBranches,Time);%线路最小功率对偶变量
Vmax=sdpvar(numBranches,Time);%线路最大功率对偶变量
EpsilonMin=sdpvar(numBuses,Time);%相角最小值对偶变量
EpsilonMax=sdpvar(numBuses,Time);%相角最大值对偶变量
Epsilon=sdpvar(1,Time);%参考相角对偶变量
MiuMeoMax=sdpvar(Time,1);%能源零售商最小中标电量对偶变量
MiuMeoMin=sdpvar(Time,1);%能源零售商最大中标电量对偶变量
disp("创建下层电力市场问题互补松弛决策变量...");
Fai=sdpvar(1,1);%购电成本替代变量
U_MiuMax=binvar(numGenerators,Time);%发电机最大输出功率0-1变量
U_MiuMin=binvar(numGenerators,Time);%发电机最小输出功率0-1变量
U_Vmin=binvar(numBranches,Time);%线路最小功率0-1变量
U_Vmax=binvar(numBranches,Time);%线路最大功率0-1变量
U_EpsilonMin=binvar(numBuses,Time);%相角最小值0-1变量
U_EpsilonMax=binvar(numBuses,Time);%相角最大值0-1变量
U_MiuMeoMax=binvar(Time,1);%能源零售商最小中标电量0-1变量
U_MiuMeoMin=binvar(Time,1);%能源零售商最大中标电量0-1变量
disp("创建下层天然气市场问题决策变量...");
Qw=sdpvar(Nw,Time);%气源产气量
Qlg=sdpvar(Nlg,Time);%主动管道输气量
Qc=sdpvar(Nc,Time);%被动管道输气量
Gmeo=sdpvar(Time,1);%天然气购气量
disp("创建下层天然气市场对偶问题决策变量...");
Faigas=sdpvar(1,1);%购气成本替代变量---new
Cmeogas=sdpvar(Time,1);%报价变量---new
Beta1Min=sdpvar(Nw,Time);%最小产气约束的对偶变量
Beta1Max=sdpvar(Nw,Time);%最大产气约束的对偶变量
Beta2Min=sdpvar(Nlg,Time);%最小主动线路输气约束的对偶变量
Beta2Max=sdpvar(Nlg,Time);%最大主动线路输气约束的对偶变量
Beta3Min=sdpvar(Nc,Time);%最小被动线路输气约束的对偶变量
Beta3Max=sdpvar(Nc,Time);%最大被动线路输气约束的对偶变量
Beta4Min=sdpvar(Time,1);%最小中标气量约束的对偶变量
Beta4Max=sdpvar(Time,1);%最大中标气量约束的对偶变量
Gama=sdpvar(Nn,Time);%天然气节点平衡约束的对偶变量
disp("创建下层天然气市场互补松弛决策变量...");
U_Beta1Min=binvar(Nw,Time);%最小产气约束的对偶变量
U_Beta1Max=binvar(Nw,Time);%最大产气约束的对偶变量
U_Beta2Min=binvar(Nlg,Time);%最小主动线路输气约束的对偶变量
U_Beta2Max=binvar(Nlg,Time);%最大主动线路输气约束的对偶变量
U_Beta3Min=binvar(Nc,Time);%最小被动线路输气约束的对偶变量
U_Beta3Max=binvar(Nc,Time);%最大被动线路输气约束的对偶变量
U_Beta4Min=binvar(Time,1);%最小中标气量约束的对偶变量
U_Beta4Max=binvar(Time,1);%最大中标气量约束的对偶变量
%% 创建目标函数
disp('创建MPEC问题目标函数...');
Cons_MPEC=[];
%各项收益建模
    Cons_MPEC=[Cons_MPEC,T1==PSP*PLScenario+GSP*GLScenario+HSP*HLScenario];%售能收入，$
    Cons_MPEC=[Cons_MPEC,T2==Ches*(Pch+Pdis)+Ctes*(Hch+Hdis)+Cges*(Gch+Gdis)];%HES+TES+GES运维成本，$
    Cons_MPEC=[Cons_MPEC,Faigas==sum(sum(repmat(Cwgas,1,Time).*Qw))+sum(sum(Beta1Max.*repmat(Qwmax,1,Time)))+sum(sum((Beta2Min+Beta2Max).*repmat(Qlgmax,1,Time)))+sum(sum(Beta3Max.*repmat(Qcmax,1,Time)))-sum(sum(Gama.*(Andg*Qdg)))];%购天然气成本，$,不写到一行会报错，极为傻逼
    Cons_MPEC=[Cons_MPEC,Fai==sum(sum(CG_temp.*PG))+sum(sum(MiuMax.*PGmax_temp))-sum(sum(Lanta.*Pd))+pi*(sum(sum(EpsilonMax+EpsilonMin)))+sum(sum(Vmin.*PLmax_temp))+sum(sum(Vmax.*PLmax_temp))];%购电成本，$   
    for s=1:numCluster
    Cons_MPEC=[Cons_MPEC,Z(s)==sum(T1(s,:))-sum(T2(s,:))-Fai-Faigas];%各场景各时段的收益，$---inconsistent，因为没有sum
    end
    Cons_MPEC=[Cons_MPEC,Obj_MPEC==dot(Probability,Z)];%%最终的期望收益，$
%% 创建约束条件
disp('创建上层问题约束条件...');
%CHP机组建模
for s=1:numCluster
    for t=1:Time
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)-PA-(PA-PB)*(Hchp(s,t)-HA)/(HA-HB)<=0];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)-PB-(PB-PC)*(Hchp(s,t)-HB)/(HB-HC)>=-(1-Ichp(s,t))*M];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)-PC-(PC-PD)*(Hchp(s,t)-HC)/(HC-HD)>=-(1-Ichp(s,t))*M];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)>=Pchpmin*Ichp(s,t)];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)<=Pchpmax*Ichp(s,t)];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Hchp(s,t)>=0];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Hchp(s,t)<=Hchpmax*Ichp(s,t)];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,GCchp(s,t)==Pchp(s,t)/ITAchp+SUchp(s,t)+SDchp(s,t)];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,SUchp(s,t)>=0];
        Cons_MPEC=[Cons_MPEC,SDchp(s,t)>=0];
    end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)-Pchp(t-1)<=RUchp];
        Cons_MPEC=[Cons_MPEC,Pchp(t-1)-Pchp(s,t)<=RDchp];
        Cons_MPEC=[Cons_MPEC,SUchp(s,t)>=SUGchp*(Ichp(s,t)-Ichp(t-1))];
        Cons_MPEC=[Cons_MPEC,SDchp(s,t)>=SDGchp*(Ichp(t-1)-Ichp(s,t))];
    end
end
%电热锅炉机组建模
for s=1:numCluster
    for t=1:Time
        Cons_MPEC=[Cons_MPEC,Heb(s,t)==ITAeb*Peb(s,t)];
        Cons_MPEC=[Cons_MPEC,Peb(s,t)>=Pebmin];
        Cons_MPEC=[Cons_MPEC,Peb(s,t)<=Pebmax];
    end
end
%自备新能源机组建模
for s=1:numCluster
    for t=1:Time
           Cons_MPEC=[Cons_MPEC,Pwt(s,t)>=0];%最小风电出力约束
%            Cons_MPEC=[Cons_MPEC,Pwt(s,t)<=WindScenario(s,t)];%最大风电出力约束
           Cons_MPEC=[Cons_MPEC,Pwt(s,t)<=WindPre(t)];%最大风电出力约束
    end
end
%电储能建模
for s=1:numCluster
 for t=1:Time
    Cons_MPEC=[Cons_MPEC,Ipch(s,t)+Ipdis(s,t)<=1];%压缩空气储能系统只能使用充电、放电模式中的一种
    Cons_MPEC=[Cons_MPEC,Pch(s,t)>=Pkmin*Ipch(s,t)];%充电模式功率约束
    Cons_MPEC=[Cons_MPEC,Pch(s,t)<=Pkmax*Ipch(s,t)];
    Cons_MPEC=[Cons_MPEC,Pdis(s,t)>=Pkmin*Ipdis(s,t)];%放电模式功率约束
    Cons_MPEC=[Cons_MPEC,Pdis(s,t)<=Pkmax*Ipdis(s,t)];
    Cons_MPEC=[Cons_MPEC,Ap(s,t)<=Akmax];%电量约束
    Cons_MPEC=[Cons_MPEC,Ap(s,t)>=Akmin];
end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Ap(s,t)==Ap(s,t-1)+ITAkch*Pch(s,t)-Pdis(s,t)/ITAkdis];%相邻时段电量约束
    end
     Cons_MPEC=[Cons_MPEC,Ap(s,1)==0.5*Akmax+ITAkch*Pch(s,1)-Pdis(s,1)/ITAkdis];%初始时段电量约束
     Cons_MPEC=[Cons_MPEC,Ap(s,Time)==0.5*Akmax];%结束时段电量约束
end
%热储能建模
for s=1:numCluster
 for t=1:Time
    Cons_MPEC=[Cons_MPEC,Ihch(s,t)+Ihdis(s,t)<=1];%压缩空气储能系统只能使用充电、放电模式中的一种
    Cons_MPEC=[Cons_MPEC,Hch(s,t)>=0];%充电模式功率约束
    Cons_MPEC=[Cons_MPEC,Hch(s,t)<=Hchmax*Ihch(s,t)];
    Cons_MPEC=[Cons_MPEC,Hdis(s,t)>=0];%放电模式功率约束
    Cons_MPEC=[Cons_MPEC,Hdis(s,t)<=Hdismax*Ihdis(s,t)];
    Cons_MPEC=[Cons_MPEC,Ah(s,t)<=Hmax];%电量约束
    Cons_MPEC=[Cons_MPEC,Ah(s,t)>=Hmin];
end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Ah(s,t)==Ah(s,t-1)+ITAhch*Hch(s,t)-Hdis(s,t)/ITAhdis];%相邻时段电量约束
    end
     Cons_MPEC=[Cons_MPEC,Ah(s,1)==0.5*Hmax+ITAhch*Hch(s,1)-Hdis(s,1)/ITAhdis];%初始时段电量约束
     Cons_MPEC=[Cons_MPEC,Ah(s,Time)==0.5*Hmax];%结束时段电量约束
end
%气储能建模
for s=1:numCluster
 for t=1:Time
    Cons_MPEC=[Cons_MPEC,Igch(s,t)+Igdis(s,t)<=1];%压缩空气储能系统只能使用充电、放电模式中的一种
    Cons_MPEC=[Cons_MPEC,Gch(s,t)>=0];%充电模式功率约束
    Cons_MPEC=[Cons_MPEC,Gch(s,t)<=Gchmax*Igch(s,t)];
    Cons_MPEC=[Cons_MPEC,Gdis(s,t)>=0];%放电模式功率约束
    Cons_MPEC=[Cons_MPEC,Gdis(s,t)<=Gdismax*Igdis(s,t)];
    Cons_MPEC=[Cons_MPEC,Ag(s,t)<=Gmax];%电量约束
    Cons_MPEC=[Cons_MPEC,Ag(s,t)>=Gmin];
end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Ag(s,t)==Ag(s,t-1)+ITAgch*Gch(s,t)-Gdis(s,t)/ITAgdis];%相邻时段电量约束
    end
     Cons_MPEC=[Cons_MPEC,Ag(s,1)==0.5*Gmax+ITAgch*Gch(s,1)-Gdis(s,1)/ITAgdis];%初始时段电量约束
     Cons_MPEC=[Cons_MPEC,Ag(s,Time)==0.5*Gmax];%结束时段电量约束
end
%能量平衡约束建模
for s=1:numCluster
    for t=1:Time
        Cons_MPEC=[Cons_MPEC,Pmeo(t)+Pwt(s,t)+Pchp(s,t)+Pdis(s,t)-Pch(s,t)-Peb(s,t)==PLScenario(s,t)];%电平衡
        Cons_MPEC=[Cons_MPEC,Gmeo(t)+Gdis(s,t)-Gch(s,t)-GCchp(s,t)==GLScenario(s,t)];%天然气平衡
        Cons_MPEC=[Cons_MPEC,Hchp(s,t)+Heb(s,t)+Hdis(s,t)-Hch(s,t)==HLScenario(s,t)];%热平衡
    end
end
%电力市场报价约束
Cons_MPEC=[Cons_MPEC,Cmeo>=Cmeomin];
Cons_MPEC=[Cons_MPEC,Cmeo<=Cmeomax];
%天然气市场报价约束
Cons_MPEC=[Cons_MPEC,Cmeogas>=Cmeogasmin];
Cons_MPEC=[Cons_MPEC,Cmeogas<=Cmeogasmax];
disp("创建下层电力问题约束条件...");
Cons_MPEC=[Cons_MPEC,(PG>=0):'MiuMin'];%最小发电功率约束
Cons_MPEC=[Cons_MPEC,(PG<=PGmax_temp):'MiuMax'];%最大发电功率约束
Cons_MPEC=[Cons_MPEC,(Bf*Delta>=-PLmax_temp):'Vmin'];%最小线路功率约束
Cons_MPEC=[Cons_MPEC,(Bf*Delta<=PLmax_temp):'Vmax'];%最大线路功率约束
Cons_MPEC=[Cons_MPEC,(Pmeo>=Pmeomin):'MiuMeoMin'];%最小中标电量约束
Cons_MPEC=[Cons_MPEC,(Pmeo<=Pmeomax):'MiuMeoMax'];%最大中标电量约束
Cons_MPEC=[Cons_MPEC,(Delta>=-pi*ones(numBuses,Time)):'EpsilonMin'];%相角最小值约束
Cons_MPEC=[Cons_MPEC,(Delta<=pi*ones(numBuses,Time)):'EpsilonMax'];%相角最大值约束
Cons_MPEC=[Cons_MPEC,(Delta(1,:)==0):'Epsilon'];%参考节点相角约束
Cons_MPEC=[Cons_MPEC,(Abg*PG-Ameo*Pmeo'-Bn*Delta==Pd):'Lanta'];%功率平衡约束，对偶变量为节点数X时间数矩阵形式
disp("创建下层问题对偶问题约束条件...");
%对偶变量范围约束
Cons_MPEC=[Cons_MPEC,MiuMax>=0];
Cons_MPEC=[Cons_MPEC,MiuMin>=0];
Cons_MPEC=[Cons_MPEC,MiuMeoMax>=0];
Cons_MPEC=[Cons_MPEC,MiuMeoMin>=0];
Cons_MPEC=[Cons_MPEC,Vmax>=0];
Cons_MPEC=[Cons_MPEC,Vmin>=0];
Cons_MPEC=[Cons_MPEC,EpsilonMax>=0];
Cons_MPEC=[Cons_MPEC,EpsilonMin>=0];
%决策变量系数为0约束
Cons_MPEC=[Cons_MPEC,(-Abg'*Lanta-MiuMin+MiuMax==-CG_temp):'PG'];%Lagrange函数对PG求偏导为0
Cons_MPEC=[Cons_MPEC,(-Cmeo+Lanta'*Ameo-MiuMeoMin+MiuMeoMax==0):'Pmeo'];%Lagrange函数对Pmeo求偏导为0----核心约束
for i=2:numBuses
Cons_MPEC=[Cons_MPEC,(-EpsilonMin(i,:)+EpsilonMax(i,:)+Bn(i,:)*Lanta-Bf_Inv(i,:)*Vmin+Bf_Inv(i,:)*Vmax==0):'Detlait'];%Lagrange函数对参考节点外相角求偏导为0
end
Cons_MPEC=[Cons_MPEC,(-EpsilonMin(1,:)+EpsilonMax(1,:)+Bn(1,:)*Lanta-Bf_Inv(1,:)*Vmin+Bf_Inv(1,:)*Vmax+Epsilon==0):'Detlat'];%Lagrange函数对参考节点相角求偏导为0
disp("创建下层问题强对偶问题互补松弛约束条件...");
Cons_MPEC=[Cons_MPEC,MiuMax<=M*U_MiuMax];
Cons_MPEC=[Cons_MPEC,PGmax_temp-PG<=M*(ones(numGenerators,Time)-U_MiuMax)];
Cons_MPEC=[Cons_MPEC,MiuMin<=M*U_MiuMin];
Cons_MPEC=[Cons_MPEC,PG<=M*(ones(numGenerators,Time)-U_MiuMin)];
Cons_MPEC=[Cons_MPEC,MiuMeoMax<=M*U_MiuMeoMax];
Cons_MPEC=[Cons_MPEC,Pmeomax-Pmeo<=M*(ones(Time,1)-U_MiuMeoMax)];
Cons_MPEC=[Cons_MPEC,MiuMeoMin<=M*U_MiuMeoMin];
Cons_MPEC=[Cons_MPEC,Pmeo-Pmeomin<=M*(ones(Time,1)-U_MiuMeoMin)];
Cons_MPEC=[Cons_MPEC,Vmax<=M*U_Vmax];
Cons_MPEC=[Cons_MPEC,PLmax_temp-Bf*Delta<=M*(ones(numBranches,Time)-U_Vmax)];
Cons_MPEC=[Cons_MPEC,Vmin<=M*U_Vmin];
Cons_MPEC=[Cons_MPEC,Bf*Delta+PLmax_temp<=M*(ones(numBranches,Time)-U_Vmin)];
Cons_MPEC=[Cons_MPEC,EpsilonMax<=M*U_EpsilonMax];
Cons_MPEC=[Cons_MPEC,pi*ones(numBuses,Time)-Delta<=M*(ones(numBuses,Time)-U_EpsilonMax)];
Cons_MPEC=[Cons_MPEC,EpsilonMin<=M*U_EpsilonMin];
Cons_MPEC=[Cons_MPEC,Delta+pi*ones(numBuses,Time)<=M*(ones(numBuses,Time)-U_EpsilonMin)];

disp("创建下层天然气市场原问题约束条件...");
Cons_MPEC=[Cons_MPEC,(Qw>=0):'Beta1Min'];%最小产气约束
Cons_MPEC=[Cons_MPEC,(Qw<=repmat(Qwmax,1,Time)):'Beta1Max'];%最大产气约束
Cons_MPEC=[Cons_MPEC,(Qlg>=-repmat(Qlgmax,1,Time)):'Beta2Min'];%最小主动线路输气约束
Cons_MPEC=[Cons_MPEC,(Qlg<=repmat(Qlgmax,1,Time)):'Beta2Max'];%最大主动线路输气约束
Cons_MPEC=[Cons_MPEC,(Qc>=0):'Beta3Min'];%最小被动线路输气约束
Cons_MPEC=[Cons_MPEC,(Qc<=repmat(Qcmax,1,Time)):'Beta3Max'];%最大动线路输气约束
Cons_MPEC=[Cons_MPEC,(Gmeo>=0):'Beta4Min'];%最小购气约束
Cons_MPEC=[Cons_MPEC,(Gmeo<=Gmeomax):'Beta4Max'];%最大产气约束
Cons_MPEC=[Cons_MPEC,(Anw*Qw-Anlg*Qlg-Anc*Qc-Anm*Gmeo'-Andg*Qdg==0):'Gama'];%节点天然气平衡约束
disp("创建下层天然气市场对偶问题约束条件...");
%对偶变量范围约束
Cons_MPEC=[Cons_MPEC,Beta1Min>=0];
Cons_MPEC=[Cons_MPEC,Beta1Max>=0];
Cons_MPEC=[Cons_MPEC,Beta2Min>=0];
Cons_MPEC=[Cons_MPEC,Beta2Max>=0];
Cons_MPEC=[Cons_MPEC,Beta3Min>=0];
Cons_MPEC=[Cons_MPEC,Beta3Max>=0];
Cons_MPEC=[Cons_MPEC,Beta4Min>=0];
Cons_MPEC=[Cons_MPEC,Beta4Max>=0];
%决策变量系数为0约束
Cons_MPEC=[Cons_MPEC,(repmat(Cwgas,1,Time)-Beta1Min+Beta1Max-Anw'*Gama==0):'Qw'];%Lagrange函数对Qw求偏导为0
Cons_MPEC=[Cons_MPEC,(-Cmeogas-Beta4Min+Beta4Max+Gama'*Anm==0):'Gmer'];%Lagrange函数对Gmer求偏导为0
Cons_MPEC=[Cons_MPEC,(-Beta2Min+Beta2Max+Anlg'*Gama==0):'Qlg'];%Lagrange函数对Qlg求偏导为0
Cons_MPEC=[Cons_MPEC,(-Beta3Min+Beta3Max+Anc'*Gama==0):'Qc'];%Lagrange函数对Qc求偏导为0
disp("创建下层天然气市场互补松弛条件约束...");
Cons_MPEC=[Cons_MPEC,Beta1Min<=M*U_Beta1Min];
Cons_MPEC=[Cons_MPEC,Qw<=M*(ones(Nw,Time)-U_Beta1Min)];
Cons_MPEC=[Cons_MPEC,Beta1Max<=M*U_Beta1Max];
Cons_MPEC=[Cons_MPEC,repmat(Qwmax,1,Time)-Qw<=M*(ones(Nw,Time)-U_Beta1Max)];
Cons_MPEC=[Cons_MPEC,Beta2Min<=M*U_Beta2Min];
Cons_MPEC=[Cons_MPEC,Qlg+repmat(Qlgmax,1,Time)<=M*(ones(Nlg,Time)-U_Beta2Min)];
Cons_MPEC=[Cons_MPEC,Beta2Max<=M*U_Beta2Max];
Cons_MPEC=[Cons_MPEC,repmat(Qlgmax,1,Time)-Qlg<=M*(ones(Nlg,Time)-U_Beta2Max)];
Cons_MPEC=[Cons_MPEC,Beta3Min<=M*U_Beta3Min];
Cons_MPEC=[Cons_MPEC,Qc<=M*(ones(Nc,Time)-U_Beta3Min)];
Cons_MPEC=[Cons_MPEC,Beta3Max<=M*U_Beta3Max];
Cons_MPEC=[Cons_MPEC,repmat(Qcmax,1,Time)-Qc<=M*(ones(Nc,Time)-U_Beta3Max)];
Cons_MPEC=[Cons_MPEC,Beta4Min<=M*U_Beta4Min];
Cons_MPEC=[Cons_MPEC,Gmeo<=M*(ones(Time,1)-U_Beta4Min)];
Cons_MPEC=[Cons_MPEC,Beta4Max<=M*U_Beta4Max];
Cons_MPEC=[Cons_MPEC,repmat(Gmax,Time,1)-Gmeo<=M*(ones(Time,1)-U_Beta4Max)];

%% 参数配置与模型求解
disp('正在求解...');
% tic;
Opt_MPEC=sdpsettings('verbose',1,'debug',1,'solver','cplex','savesolveroutput',1,'savesolverinput',1);
Opt_MPEC.cplex.exportmodel='Opt_MPEC.lp';
MPECModel=optimize(Cons_MPEC,-Obj_MPEC,Opt_MPEC);
if MPECModel.problem==0
    disp('MPEC模型可解');
toc;
tsolve=toc-tic;%求取求解时间
else
disp('MPEC模型不可解');
end
%% 获取计算结果
Res_Case3_MPEC.Pchp=double(Pchp);
Res_Case3_MPEC.Hchp=double(Hchp);
Res_Case3_MPEC.Ichp=double(Ichp);
Res_Case3_MPEC.SUchp=double(SUchp);
Res_Case3_MPEC.SDchp=double(SDchp);
Res_Case3_MPEC.GCchp=double(GCchp);
Res_Case3_MPEC.Peb=double(Peb);
Res_Case3_MPEC.Heb=double(Heb);
Res_Case3_MPEC.Pwt=double(Pwt);
Res_Case3_MPEC.Pch=double(Pch);
Res_Case3_MPEC.Pdis=double(Pdis);
Res_Case3_MPEC.Hch=double(Hch);
Res_Case3_MPEC.Hdis=double(Hdis);
Res_Case3_MPEC.Gch=double(Gch);
Res_Case3_MPEC.Gdis=double(Gdis);
Res_Case3_MPEC.Ap=double(Ap);
Res_Case3_MPEC.Ah=double(Ah);
Res_Case3_MPEC.Ag=double(Ag);
Res_Case3_MPEC.Pmeo=double(Pmeo);

Res_Case3_MPEC.T1=double(T1);
Res_Case3_MPEC.Z=double(Z);
Res_Case3_MPEC.Obj_MPEC=double(Obj_MPEC);
Res_Case3_MPEC.PG=value(PG);
Res_Case3_MPEC.Delta=value(Delta);
Res_Case3_MPEC.Pmeo=value(Pmeo);
Res_Case3_MPEC.Lanta=value(Lanta);
Res_Case3_MPEC.MiuMax=value(MiuMax);
Res_Case3_MPEC.MiuMin=value(MiuMin);
Res_Case3_MPEC.MiuMeoMax=value(MiuMeoMax);
Res_Case3_MPEC.MiuMeoMin=value(MiuMeoMin);
Res_Case3_MPEC.Vmax=value(Vmax);
Res_Case3_MPEC.Vmin=value(Vmin);
Res_Case3_MPEC.EpsilonMax=value(EpsilonMax);
Res_Case3_MPEC.EpsilonMin=value(EpsilonMin);
Res_Case3_MPEC.Epsilon=value(Epsilon);
Res_Case3_MPEC.Cmeo=value(Cmeo);
Res_Case3_MPEC.U_MiuMax=value(U_MiuMax);
Res_Case3_MPEC.U_MiuMin=value(U_MiuMin);
Res_Case3_MPEC.U_MiuMeoMax=value(U_MiuMeoMax);
Res_Case3_MPEC.U_MiuMeoMin=value(U_MiuMeoMin);
Res_Case3_MPEC.U_Vmax=value(U_Vmax);
Res_Case3_MPEC.U_Vmin=value(U_Vmin);
Res_Case3_MPEC.U_EpsilonMax=value(U_EpsilonMax);
Res_Case3_MPEC.U_EpsilonMin=value(U_EpsilonMin);
Res_Case3_MPEC.Fai=value(Fai);
Res_Case3_MPEC.MCP=value(Lanta(5,:))';
Res_Case3_MPEC.Fai=value(Fai);
Res_Case3_MPEC.Faigas=value(Faigas);
Res_Case3_MPEC.Z=value(Z);
Res_Case3_MPEC.ExpectedProfit=value(Obj_MPEC);%期望利润
Res_Case3_MPEC.WorstProfit=min(value(Z));%最差利润



Res_Case3_MPEC.Qw=value(Qw);
Res_Case3_MPEC.Qlg=value(Qlg);
Res_Case3_MPEC.Qc=value(Qc);
Res_Case3_MPEC.Gmeo=value(Gmeo);
Res_Case3_MPEC.Gama=value(Gama);%dual求出对偶变量值有时需要取负号需要自己判断下
Res_Case3_MPEC.Cmeogas=value(Cmeogas);
Res_Case3_MPEC.Beta1Min=value(Beta1Min);Res_Case3_MPEC.Beta1Max=value(Beta1Max);
Res_Case3_MPEC.Beta2Min=value(Beta2Min);Res_Case3_MPEC.Beta2Max=value(Beta2Max);
Res_Case3_MPEC.Beta3Min=value(Beta3Min);Res_Case3_MPEC.Beta3Max=value(Beta3Max);
Res_Case3_MPEC.Beta4Min=value(Beta4Min);Res_Case3_MPEC.Beta4Max=value(Beta4Max);
Res_Case3_MPEC.U_Beta3Min=value(U_Beta3Min);Res_Case3_MPEC.U_Beta3Max=value(U_Beta3Max);


Res_Case3_MPEC.RR=zeros(numCluster,1);
for s=1:numCluster
    Res_Case3_MPEC.RR(s)=1-Res_Case3_MPEC.Z(s)/Z_optimal(s);
end
Res_Case3_MPEC.MRR=max(Res_Case3_MPEC.RR);
