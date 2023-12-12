%% p-robust案例
%% 随机p鲁棒优化程序
clc;

clear all;
Input;%输入数据
Z_optimal=[];%记录每个场景下的最优目标利润值，$
Pmeo_optimal=[];%记录每个场景下的最优购电值
%% 从每个场景寻找最优值
for s=1:numCluster  

disp("创建上层问题决策变量...");
%CHP机组
Pchp=sdpvar(Time,1);%CHP机组电功率/MW
Hchp=sdpvar(Time,1);%CHP机组热功率/MW
Ichp=binvar(Time,1);%CHP机组状态
SUchp=sdpvar(Time,1);%CHP机组开启成本/$
SDchp=sdpvar(Time,1);%CHP机组关停成本/$
GCchp=sdpvar(Time,1);%CHP机组消耗天然气量/MWh
%电热锅炉机组
Peb=sdpvar(Time,1);%电功率/MWh
Heb=sdpvar(Time,1);%热功率/MWh
%新能源机组
Pwt=sdpvar(Time,1);%风电上网功率，MW
%电储能
Ipch=binvar(Time,1);%储能充电状态
Ipdis=binvar(Time,1);%储能放电状态
Pch=sdpvar(Time,1);%储能充电功率/MW
Pdis=sdpvar(Time,1);%储能放电功率/MW
Ap=sdpvar(Time,1);%储能电量/MWh
%热储能
Ihch=binvar(Time,1);%储能充电状态
Ihdis=binvar(Time,1);%储能放电状态
Hch=sdpvar(Time,1);%储能充电功率/MW
Hdis=sdpvar(Time,1);%储能放电功率/MW
Ah=sdpvar(Time,1);%储能电量/MWh
%气储能
Igch=binvar(Time,1);%储能充电状态
Igdis=binvar(Time,1);%储能放电状态
Gch=sdpvar(Time,1);%储能充电功率/MW
Gdis=sdpvar(Time,1);%储能放电功率/MW
Ag=sdpvar(Time,1);%储能电量/MWh
%购能量
Pmeo=sdpvar(Time,1);%购能量，MWh
% Ggm=sdpvar(Time,1);%购买天然气量，MWh
%报价参数
Cmeo=sdpvar(Time,1);%报价，$/MWh
%收益各项参数
T1=sdpvar(Time,1);%售能收入，$
T2=sdpvar(Time,1);%储能运维成本，$
Z=sdpvar(1,1);%各场景的收益，$
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
%     Cons_MPEC=[Cons_MPEC,T1==PSP*PLScenario(s,:)'+GSP*GLScenario(s,:)'+HSP*HLScenario(s,:)'];%售能收入，$
    Cons_MPEC=[Cons_MPEC,T1==PSP*PL+GSP*GL+HSP*HL];%售能收入，$
    Cons_MPEC=[Cons_MPEC,T2==Ches*(Pch+Pdis)+Ctes*(Hch+Hdis)+Cges*(Gch+Gdis)];%HES+TES+GES运维成本，$
    Cons_MPEC=[Cons_MPEC,Faigas==sum(sum(repmat(Cwgas,1,Time).*Qw))+sum(sum(Beta1Max.*repmat(Qwmax,1,Time)))+sum(sum((Beta2Min+Beta2Max).*repmat(Qlgmax,1,Time)))+sum(sum(Beta3Max.*repmat(Qcmax,1,Time)))-sum(sum(Gama.*(Andg*Qdg)))];%购天然气成本，$,不写到一行会报错，极为傻逼
    Cons_MPEC=[Cons_MPEC,Fai==sum(sum(CG_temp.*PG))+sum(sum(MiuMax.*PGmax_temp))-sum(sum(Lanta.*Pd))+pi*(sum(sum(EpsilonMax+EpsilonMin)))+sum(sum(Vmin.*PLmax_temp))+sum(sum(Vmax.*PLmax_temp))];%购电成本，$   
    Cons_MPEC=[Cons_MPEC,Z==sum(T1)-sum(T2)-Fai-Faigas];%各场景各时段的收益，$---inconsistent，因为没有sum
    Cons_MPEC=[Cons_MPEC,Obj_MPEC==Z];%%最终的期望收益，$
%% 创建约束条件
disp('创建上层问题约束条件...');
%CHP机组建模
    for t=1:Time
        Cons_MPEC=[Cons_MPEC,Pchp(t)-PA-(PA-PB)*(Hchp(t)-HA)/(HA-HB)<=0];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Pchp(t)-PB-(PB-PC)*(Hchp(t)-HB)/(HB-HC)>=-(1-Ichp(t))*M];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Pchp(t)-PC-(PC-PD)*(Hchp(t)-HC)/(HC-HD)>=-(1-Ichp(t))*M];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Pchp(t)>=Pchpmin*Ichp(t)];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Pchp(t)<=Pchpmax*Ichp(t)];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Hchp(t)>=0];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,Hchp(t)<=Hchpmax*Ichp(t)];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,GCchp(t)==Pchp(t)/ITAchp+SUchp(t)+SDchp(t)];%CHP运行可行域约束
        Cons_MPEC=[Cons_MPEC,SUchp(t)>=0];
        Cons_MPEC=[Cons_MPEC,SDchp(t)>=0];
    end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Pchp(t)-Pchp(t-1)<=RUchp];
        Cons_MPEC=[Cons_MPEC,Pchp(t-1)-Pchp(t)<=RDchp];
        Cons_MPEC=[Cons_MPEC,SUchp(t)>=SUGchp*(Ichp(t)-Ichp(t-1))];
        Cons_MPEC=[Cons_MPEC,SDchp(t)>=SDGchp*(Ichp(t-1)-Ichp(t))];
    end
%电热锅炉机组建模
    for t=1:Time
        Cons_MPEC=[Cons_MPEC,Heb(t)==ITAeb*Peb(t)];
        Cons_MPEC=[Cons_MPEC,Peb(t)>=Pebmin];
        Cons_MPEC=[Cons_MPEC,Peb(t)<=Pebmax];
    end
%自备新能源机组建模
    for t=1:Time
           Cons_MPEC=[Cons_MPEC,Pwt(t)>=0];%最小风电出力约束
%            Cons_MPEC=[Cons_MPEC,Pwt(t)<=WindScenario(s,t)];%最大风电出力约束
           Cons_MPEC=[Cons_MPEC,Pwt(t)<=WindPre(t)];%最大风电出力约束
    end
%电储能建模
 for t=1:Time
    Cons_MPEC=[Cons_MPEC,Ipch(t)+Ipdis(t)<=1];%压缩空气储能系统只能使用充电、放电模式中的一种
    Cons_MPEC=[Cons_MPEC,Pch(t)>=Pkmin*Ipch(t)];%充电模式功率约束
    Cons_MPEC=[Cons_MPEC,Pch(t)<=Pkmax*Ipch(t)];
    Cons_MPEC=[Cons_MPEC,Pdis(t)>=Pkmin*Ipdis(t)];%放电模式功率约束
    Cons_MPEC=[Cons_MPEC,Pdis(t)<=Pkmax*Ipdis(t)];
    Cons_MPEC=[Cons_MPEC,Ap(t)<=Akmax];%电量约束
    Cons_MPEC=[Cons_MPEC,Ap(t)>=Akmin];
end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Ap(t)==Ap(t-1)+ITAkch*Pch(t)-Pdis(t)/ITAkdis];%相邻时段电量约束
    end
     Cons_MPEC=[Cons_MPEC,Ap(1)==0.5*Akmax+ITAkch*Pch(1)-Pdis(1)/ITAkdis];%初始时段电量约束
     Cons_MPEC=[Cons_MPEC,Ap(Time)==0.5*Akmax];%结束时段电量约束
%热储能建模
 for t=1:Time
    Cons_MPEC=[Cons_MPEC,Ihch(t)+Ihdis(t)<=1];%压缩空气储能系统只能使用充电、放电模式中的一种
    Cons_MPEC=[Cons_MPEC,Hch(t)>=0];%充电模式功率约束
    Cons_MPEC=[Cons_MPEC,Hch(t)<=Hchmax*Ihch(t)];
    Cons_MPEC=[Cons_MPEC,Hdis(t)>=0];%放电模式功率约束
    Cons_MPEC=[Cons_MPEC,Hdis(t)<=Hdismax*Ihdis(t)];
    Cons_MPEC=[Cons_MPEC,Ah(t)<=Hmax];%电量约束
    Cons_MPEC=[Cons_MPEC,Ah(t)>=Hmin];
end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Ah(t)==Ah(t-1)+ITAhch*Hch(t)-Hdis(t)/ITAhdis];%相邻时段电量约束
    end
     Cons_MPEC=[Cons_MPEC,Ah(1)==0.5*Hmax+ITAhch*Hch(1)-Hdis(1)/ITAhdis];%初始时段电量约束
     Cons_MPEC=[Cons_MPEC,Ah(Time)==0.5*Hmax];%结束时段电量约束
%气储能建模
 for t=1:Time
    Cons_MPEC=[Cons_MPEC,Igch(t)+Igdis(t)<=1];%压缩空气储能系统只能使用充电、放电模式中的一种
    Cons_MPEC=[Cons_MPEC,Gch(t)>=0];%充电模式功率约束
    Cons_MPEC=[Cons_MPEC,Gch(t)<=Gchmax*Igch(t)];
    Cons_MPEC=[Cons_MPEC,Gdis(t)>=0];%放电模式功率约束
    Cons_MPEC=[Cons_MPEC,Gdis(t)<=Gdismax*Igdis(t)];
    Cons_MPEC=[Cons_MPEC,Ag(t)<=Gmax];%电量约束
    Cons_MPEC=[Cons_MPEC,Ag(t)>=Gmin];
end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Ag(t)==Ag(t-1)+ITAgch*Gch(t)-Gdis(t)/ITAgdis];%相邻时段电量约束
    end
     Cons_MPEC=[Cons_MPEC,Ag(1)==0.5*Gmax+ITAgch*Gch(1)-Gdis(1)/ITAgdis];%初始时段电量约束
     Cons_MPEC=[Cons_MPEC,Ag(Time)==0.5*Gmax];%结束时段电量约束
%能量平衡约束建模
    for t=1:Time
        Cons_MPEC=[Cons_MPEC,Pmeo(t)+Pwt(t)+Pchp(t)+Pdis(t)-Pch(t)-Peb(t)==PLScenario(s,t)];%电平衡
        Cons_MPEC=[Cons_MPEC,Gmeo(t)+Gdis(t)-Gch(t)-GCchp(t)==GLScenario(s,t)];%天然气平衡
        Cons_MPEC=[Cons_MPEC,Hchp(t)+Heb(t)+Hdis(t)-Hch(t)==HLScenario(s,t)];%热平衡
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
tic;
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
Z_optimal=[Z_optimal;double(Z)] ;%最优值
Pmeo_optimal=[Pmeo_optimal;double(Pmeo)];%最优值
end   


%% 随机p鲁棒优化问题

p=[M;0.038;0.037;0.036;0.035;0.034;0.033];

MRR=[];%最大遗憾度
MRR_Reduction=[];%最大遗憾度减少程度
ExpectedProfit=[];%期望利润,$
WorstProfit=[];%最差利润,$
ExpectedProfit_Scenario=[];%各场景下期望利润,$
Pmeo_temp=[];

for u=1:size(p)
    
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
    Cons_MPEC=[Cons_MPEC,Z>=(1-p(u))*Z_optimal];%随机p鲁棒约束，限制每个场景下的利润不小于最大利润的1-p倍
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
tic;
Opt_MPEC=sdpsettings('verbose',1,'debug',1,'solver','cplex','savesolveroutput',1,'savesolverinput',1);
Opt_MPEC.cplex.exportmodel='Opt_MPEC.lp';
MPECModel=optimize(Cons_MPEC,-Obj_MPEC,Opt_MPEC);
if MPECModel.problem==0
    disp('随机p鲁棒模型可解');
    Zs_optimal=double(Z);%各场景下利润优化值
    RR_temp=zeros(numCluster,1);%各场景的遗憾度
    for s=1:numCluster
    RR_temp(s)=1-Zs_optimal(s)/Z_optimal(s);
    end
    MRR=[MRR;max(RR_temp)];%最大相对遗憾度
    Pmeo_temp=[Pmeo_temp;value(Pmeo)];%记录中标电量
    ExpectedProfit=[ExpectedProfit;double(Obj_MPEC)];%期望利润
    WorstProfit=[WorstProfit;min(Zs_optimal)];%最差利润
    ExpectedProfit_Scenario=[ExpectedProfit_Scenario;Zs_optimal];%各场景下期望利润
toc;
tsolve=toc-tic;%求取求解时间
else
disp('随机p鲁棒模型不可解');
end
end

%% 结果处理
Res_Case_StochasticPRobust.p=p;
Res_Case_StochasticPRobust.Zs_optimal=Zs_optimal;
Res_Case_StochasticPRobust.MRR=MRR;
Res_Case_StochasticPRobust.ExpectedProfit=ExpectedProfit;
Res_Case_StochasticPRobust.WorstProfit=WorstProfit;
Res_Case_StochasticPRobust.ExpectedProfit_Scenario=ExpectedProfit_Scenario;
Res_Case_StochasticPRobust.ScenarioProfit=reshape(ExpectedProfit_Scenario,10,7);