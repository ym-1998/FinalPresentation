%% ����Ż�����
clc;
%% ����Ż�ģ��+��CHP����+���ȹ�¯+���ܼ���---�²���������г�����Ȼ���г�
Input;
tic;%����ģ��ʱ���������
%% ���߱���
disp("�����ϲ�������߱���...");
%CHP����
Pchp=sdpvar(numCluster,Time);%CHP����繦��/MW
Hchp=sdpvar(numCluster,Time);%CHP�����ȹ���/MW
Ichp=binvar(numCluster,Time);%CHP����״̬
SUchp=sdpvar(numCluster,Time);%CHP���鿪���ɱ�/$
SDchp=sdpvar(numCluster,Time);%CHP�����ͣ�ɱ�/$
GCchp=sdpvar(numCluster,Time);%CHP����������Ȼ����/MWh
%���ȹ�¯����
Peb=sdpvar(numCluster,Time);%�繦��/MWh
Heb=sdpvar(numCluster,Time);%�ȹ���/MWh
%����Դ����
Pwt=sdpvar(numCluster,Time);%����������ʣ�MW
%�索��
Ipch=binvar(numCluster,Time);%���ܳ��״̬
Ipdis=binvar(numCluster,Time);%���ܷŵ�״̬
Pch=sdpvar(numCluster,Time);%���ܳ�繦��/MW
Pdis=sdpvar(numCluster,Time);%���ܷŵ繦��/MW
Ap=sdpvar(numCluster,Time);%���ܵ���/MWh
%�ȴ���
Ihch=binvar(numCluster,Time);%���ܳ��״̬
Ihdis=binvar(numCluster,Time);%���ܷŵ�״̬
Hch=sdpvar(numCluster,Time);%���ܳ�繦��/MW
Hdis=sdpvar(numCluster,Time);%���ܷŵ繦��/MW
Ah=sdpvar(numCluster,Time);%���ܵ���/MWh
%������
Igch=binvar(numCluster,Time);%���ܳ��״̬
Igdis=binvar(numCluster,Time);%���ܷŵ�״̬
Gch=sdpvar(numCluster,Time);%���ܳ�繦��/MW
Gdis=sdpvar(numCluster,Time);%���ܷŵ繦��/MW
Ag=sdpvar(numCluster,Time);%���ܵ���/MWh
%������
Pmeo=sdpvar(Time,1);%��������MWh
%���۲���
Cmeo=sdpvar(Time,1);%���ۣ�$/MWh
%����������
T1=sdpvar(numCluster,Time);%�������룬$
T2=sdpvar(numCluster,Time);%������ά�ɱ���$
Z=sdpvar(numCluster,1);%�����������棬$
Obj_MPEC=sdpvar(1,1);%���յ��������棬$
disp("�����²�����г�������߱���...");
PG=sdpvar(numGenerators,Time);%������,MW
Delta=sdpvar(numBuses,Time);%�ڵ���ǣ�rad
disp("�����²�����г������ż������߱���...");
Lanta=sdpvar(numBuses,Time);%����ƽ��Լ���Ķ�ż����
MiuMax=sdpvar(numGenerators,Time);%��������������ʶ�ż����
MiuMin=sdpvar(numGenerators,Time);%�������С������ʶ�ż����
Vmin=sdpvar(numBranches,Time);%��·��С���ʶ�ż����
Vmax=sdpvar(numBranches,Time);%��·����ʶ�ż����
EpsilonMin=sdpvar(numBuses,Time);%�����Сֵ��ż����
EpsilonMax=sdpvar(numBuses,Time);%������ֵ��ż����
Epsilon=sdpvar(1,Time);%�ο���Ƕ�ż����
MiuMeoMax=sdpvar(Time,1);%��Դ��������С�б������ż����
MiuMeoMin=sdpvar(Time,1);%��Դ����������б������ż����
disp("�����²�����г����⻥���ɳھ��߱���...");
Fai=sdpvar(1,1);%����ɱ��������
U_MiuMax=binvar(numGenerators,Time);%���������������0-1����
U_MiuMin=binvar(numGenerators,Time);%�������С�������0-1����
U_Vmin=binvar(numBranches,Time);%��·��С����0-1����
U_Vmax=binvar(numBranches,Time);%��·�����0-1����
U_EpsilonMin=binvar(numBuses,Time);%�����Сֵ0-1����
U_EpsilonMax=binvar(numBuses,Time);%������ֵ0-1����
U_MiuMeoMax=binvar(Time,1);%��Դ��������С�б����0-1����
U_MiuMeoMin=binvar(Time,1);%��Դ����������б����0-1����
disp("�����²���Ȼ���г�������߱���...");
Qw=sdpvar(Nw,Time);%��Դ������
Qlg=sdpvar(Nlg,Time);%�����ܵ�������
Qc=sdpvar(Nc,Time);%�����ܵ�������
Gmeo=sdpvar(Time,1);%��Ȼ��������
disp("�����²���Ȼ���г���ż������߱���...");
Faigas=sdpvar(1,1);%�����ɱ��������---new
Cmeogas=sdpvar(Time,1);%���۱���---new
Beta1Min=sdpvar(Nw,Time);%��С����Լ���Ķ�ż����
Beta1Max=sdpvar(Nw,Time);%������Լ���Ķ�ż����
Beta2Min=sdpvar(Nlg,Time);%��С������·����Լ���Ķ�ż����
Beta2Max=sdpvar(Nlg,Time);%���������·����Լ���Ķ�ż����
Beta3Min=sdpvar(Nc,Time);%��С������·����Լ���Ķ�ż����
Beta3Max=sdpvar(Nc,Time);%��󱻶���·����Լ���Ķ�ż����
Beta4Min=sdpvar(Time,1);%��С�б�����Լ���Ķ�ż����
Beta4Max=sdpvar(Time,1);%����б�����Լ���Ķ�ż����
Gama=sdpvar(Nn,Time);%��Ȼ���ڵ�ƽ��Լ���Ķ�ż����
disp("�����²���Ȼ���г������ɳھ��߱���...");
U_Beta1Min=binvar(Nw,Time);%��С����Լ���Ķ�ż����
U_Beta1Max=binvar(Nw,Time);%������Լ���Ķ�ż����
U_Beta2Min=binvar(Nlg,Time);%��С������·����Լ���Ķ�ż����
U_Beta2Max=binvar(Nlg,Time);%���������·����Լ���Ķ�ż����
U_Beta3Min=binvar(Nc,Time);%��С������·����Լ���Ķ�ż����
U_Beta3Max=binvar(Nc,Time);%��󱻶���·����Լ���Ķ�ż����
U_Beta4Min=binvar(Time,1);%��С�б�����Լ���Ķ�ż����
U_Beta4Max=binvar(Time,1);%����б�����Լ���Ķ�ż����
%% ����Ŀ�꺯��
disp('����MPEC����Ŀ�꺯��...');
Cons_MPEC=[];
%�������潨ģ
    Cons_MPEC=[Cons_MPEC,T1==PSP*PLScenario+GSP*GLScenario+HSP*HLScenario];%�������룬$
    Cons_MPEC=[Cons_MPEC,T2==Ches*(Pch+Pdis)+Ctes*(Hch+Hdis)+Cges*(Gch+Gdis)];%HES+TES+GES��ά�ɱ���$
    Cons_MPEC=[Cons_MPEC,Faigas==sum(sum(repmat(Cwgas,1,Time).*Qw))+sum(sum(Beta1Max.*repmat(Qwmax,1,Time)))+sum(sum((Beta2Min+Beta2Max).*repmat(Qlgmax,1,Time)))+sum(sum(Beta3Max.*repmat(Qcmax,1,Time)))-sum(sum(Gama.*(Andg*Qdg)))];%����Ȼ���ɱ���$,��д��һ�лᱨ����Ϊɵ��
    Cons_MPEC=[Cons_MPEC,Fai==sum(sum(CG_temp.*PG))+sum(sum(MiuMax.*PGmax_temp))-sum(sum(Lanta.*Pd))+pi*(sum(sum(EpsilonMax+EpsilonMin)))+sum(sum(Vmin.*PLmax_temp))+sum(sum(Vmax.*PLmax_temp))];%����ɱ���$   
    for s=1:numCluster
    Cons_MPEC=[Cons_MPEC,Z(s)==sum(T1(s,:))-sum(T2(s,:))-Fai-Faigas];%��������ʱ�ε����棬$---inconsistent����Ϊû��sum
    end
    Cons_MPEC=[Cons_MPEC,Obj_MPEC==dot(Probability,Z)];%%���յ��������棬$
%% ����Լ������
disp('�����ϲ�����Լ������...');
%CHP���齨ģ
for s=1:numCluster
    for t=1:Time
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)-PA-(PA-PB)*(Hchp(s,t)-HA)/(HA-HB)<=0];%CHP���п�����Լ��
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)-PB-(PB-PC)*(Hchp(s,t)-HB)/(HB-HC)>=-(1-Ichp(s,t))*M];%CHP���п�����Լ��
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)-PC-(PC-PD)*(Hchp(s,t)-HC)/(HC-HD)>=-(1-Ichp(s,t))*M];%CHP���п�����Լ��
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)>=Pchpmin*Ichp(s,t)];%CHP���п�����Լ��
        Cons_MPEC=[Cons_MPEC,Pchp(s,t)<=Pchpmax*Ichp(s,t)];%CHP���п�����Լ��
        Cons_MPEC=[Cons_MPEC,Hchp(s,t)>=0];%CHP���п�����Լ��
        Cons_MPEC=[Cons_MPEC,Hchp(s,t)<=Hchpmax*Ichp(s,t)];%CHP���п�����Լ��
        Cons_MPEC=[Cons_MPEC,GCchp(s,t)==Pchp(s,t)/ITAchp+SUchp(s,t)+SDchp(s,t)];%CHP���п�����Լ��
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
%���ȹ�¯���齨ģ
for s=1:numCluster
    for t=1:Time
        Cons_MPEC=[Cons_MPEC,Heb(s,t)==ITAeb*Peb(s,t)];
        Cons_MPEC=[Cons_MPEC,Peb(s,t)>=Pebmin];
        Cons_MPEC=[Cons_MPEC,Peb(s,t)<=Pebmax];
    end
end
%�Ա�����Դ���齨ģ
for s=1:numCluster
    for t=1:Time
           Cons_MPEC=[Cons_MPEC,Pwt(s,t)>=0];%��С������Լ��
%            Cons_MPEC=[Cons_MPEC,Pwt(s,t)<=WindScenario(s,t)];%��������Լ��
           Cons_MPEC=[Cons_MPEC,Pwt(s,t)<=WindPre(t)];%��������Լ��
    end
end
%�索�ܽ�ģ
for s=1:numCluster
 for t=1:Time
    Cons_MPEC=[Cons_MPEC,Ipch(s,t)+Ipdis(s,t)<=1];%ѹ����������ϵͳֻ��ʹ�ó�硢�ŵ�ģʽ�е�һ��
    Cons_MPEC=[Cons_MPEC,Pch(s,t)>=Pkmin*Ipch(s,t)];%���ģʽ����Լ��
    Cons_MPEC=[Cons_MPEC,Pch(s,t)<=Pkmax*Ipch(s,t)];
    Cons_MPEC=[Cons_MPEC,Pdis(s,t)>=Pkmin*Ipdis(s,t)];%�ŵ�ģʽ����Լ��
    Cons_MPEC=[Cons_MPEC,Pdis(s,t)<=Pkmax*Ipdis(s,t)];
    Cons_MPEC=[Cons_MPEC,Ap(s,t)<=Akmax];%����Լ��
    Cons_MPEC=[Cons_MPEC,Ap(s,t)>=Akmin];
end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Ap(s,t)==Ap(s,t-1)+ITAkch*Pch(s,t)-Pdis(s,t)/ITAkdis];%����ʱ�ε���Լ��
    end
     Cons_MPEC=[Cons_MPEC,Ap(s,1)==0.5*Akmax+ITAkch*Pch(s,1)-Pdis(s,1)/ITAkdis];%��ʼʱ�ε���Լ��
     Cons_MPEC=[Cons_MPEC,Ap(s,Time)==0.5*Akmax];%����ʱ�ε���Լ��
end
%�ȴ��ܽ�ģ
for s=1:numCluster
 for t=1:Time
    Cons_MPEC=[Cons_MPEC,Ihch(s,t)+Ihdis(s,t)<=1];%ѹ����������ϵͳֻ��ʹ�ó�硢�ŵ�ģʽ�е�һ��
    Cons_MPEC=[Cons_MPEC,Hch(s,t)>=0];%���ģʽ����Լ��
    Cons_MPEC=[Cons_MPEC,Hch(s,t)<=Hchmax*Ihch(s,t)];
    Cons_MPEC=[Cons_MPEC,Hdis(s,t)>=0];%�ŵ�ģʽ����Լ��
    Cons_MPEC=[Cons_MPEC,Hdis(s,t)<=Hdismax*Ihdis(s,t)];
    Cons_MPEC=[Cons_MPEC,Ah(s,t)<=Hmax];%����Լ��
    Cons_MPEC=[Cons_MPEC,Ah(s,t)>=Hmin];
end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Ah(s,t)==Ah(s,t-1)+ITAhch*Hch(s,t)-Hdis(s,t)/ITAhdis];%����ʱ�ε���Լ��
    end
     Cons_MPEC=[Cons_MPEC,Ah(s,1)==0.5*Hmax+ITAhch*Hch(s,1)-Hdis(s,1)/ITAhdis];%��ʼʱ�ε���Լ��
     Cons_MPEC=[Cons_MPEC,Ah(s,Time)==0.5*Hmax];%����ʱ�ε���Լ��
end
%�����ܽ�ģ
for s=1:numCluster
 for t=1:Time
    Cons_MPEC=[Cons_MPEC,Igch(s,t)+Igdis(s,t)<=1];%ѹ����������ϵͳֻ��ʹ�ó�硢�ŵ�ģʽ�е�һ��
    Cons_MPEC=[Cons_MPEC,Gch(s,t)>=0];%���ģʽ����Լ��
    Cons_MPEC=[Cons_MPEC,Gch(s,t)<=Gchmax*Igch(s,t)];
    Cons_MPEC=[Cons_MPEC,Gdis(s,t)>=0];%�ŵ�ģʽ����Լ��
    Cons_MPEC=[Cons_MPEC,Gdis(s,t)<=Gdismax*Igdis(s,t)];
    Cons_MPEC=[Cons_MPEC,Ag(s,t)<=Gmax];%����Լ��
    Cons_MPEC=[Cons_MPEC,Ag(s,t)>=Gmin];
end
    for t=2:Time
        Cons_MPEC=[Cons_MPEC,Ag(s,t)==Ag(s,t-1)+ITAgch*Gch(s,t)-Gdis(s,t)/ITAgdis];%����ʱ�ε���Լ��
    end
     Cons_MPEC=[Cons_MPEC,Ag(s,1)==0.5*Gmax+ITAgch*Gch(s,1)-Gdis(s,1)/ITAgdis];%��ʼʱ�ε���Լ��
     Cons_MPEC=[Cons_MPEC,Ag(s,Time)==0.5*Gmax];%����ʱ�ε���Լ��
end
%����ƽ��Լ����ģ
for s=1:numCluster
    for t=1:Time
        Cons_MPEC=[Cons_MPEC,Pmeo(t)+Pwt(s,t)+Pchp(s,t)+Pdis(s,t)-Pch(s,t)-Peb(s,t)==PLScenario(s,t)];%��ƽ��
        Cons_MPEC=[Cons_MPEC,Gmeo(t)+Gdis(s,t)-Gch(s,t)-GCchp(s,t)==GLScenario(s,t)];%��Ȼ��ƽ��
        Cons_MPEC=[Cons_MPEC,Hchp(s,t)+Heb(s,t)+Hdis(s,t)-Hch(s,t)==HLScenario(s,t)];%��ƽ��
    end
end
%�����г�����Լ��
Cons_MPEC=[Cons_MPEC,Cmeo>=Cmeomin];
Cons_MPEC=[Cons_MPEC,Cmeo<=Cmeomax];
%��Ȼ���г�����Լ��
Cons_MPEC=[Cons_MPEC,Cmeogas>=Cmeogasmin];
Cons_MPEC=[Cons_MPEC,Cmeogas<=Cmeogasmax];
disp("�����²��������Լ������...");
Cons_MPEC=[Cons_MPEC,(PG>=0):'MiuMin'];%��С���繦��Լ��
Cons_MPEC=[Cons_MPEC,(PG<=PGmax_temp):'MiuMax'];%��󷢵繦��Լ��
Cons_MPEC=[Cons_MPEC,(Bf*Delta>=-PLmax_temp):'Vmin'];%��С��·����Լ��
Cons_MPEC=[Cons_MPEC,(Bf*Delta<=PLmax_temp):'Vmax'];%�����·����Լ��
Cons_MPEC=[Cons_MPEC,(Pmeo>=Pmeomin):'MiuMeoMin'];%��С�б����Լ��
Cons_MPEC=[Cons_MPEC,(Pmeo<=Pmeomax):'MiuMeoMax'];%����б����Լ��
Cons_MPEC=[Cons_MPEC,(Delta>=-pi*ones(numBuses,Time)):'EpsilonMin'];%�����СֵԼ��
Cons_MPEC=[Cons_MPEC,(Delta<=pi*ones(numBuses,Time)):'EpsilonMax'];%������ֵԼ��
Cons_MPEC=[Cons_MPEC,(Delta(1,:)==0):'Epsilon'];%�ο��ڵ����Լ��
Cons_MPEC=[Cons_MPEC,(Abg*PG-Ameo*Pmeo'-Bn*Delta==Pd):'Lanta'];%����ƽ��Լ������ż����Ϊ�ڵ���Xʱ����������ʽ
disp("�����²������ż����Լ������...");
%��ż������ΧԼ��
Cons_MPEC=[Cons_MPEC,MiuMax>=0];
Cons_MPEC=[Cons_MPEC,MiuMin>=0];
Cons_MPEC=[Cons_MPEC,MiuMeoMax>=0];
Cons_MPEC=[Cons_MPEC,MiuMeoMin>=0];
Cons_MPEC=[Cons_MPEC,Vmax>=0];
Cons_MPEC=[Cons_MPEC,Vmin>=0];
Cons_MPEC=[Cons_MPEC,EpsilonMax>=0];
Cons_MPEC=[Cons_MPEC,EpsilonMin>=0];
%���߱���ϵ��Ϊ0Լ��
Cons_MPEC=[Cons_MPEC,(-Abg'*Lanta-MiuMin+MiuMax==-CG_temp):'PG'];%Lagrange������PG��ƫ��Ϊ0
Cons_MPEC=[Cons_MPEC,(-Cmeo+Lanta'*Ameo-MiuMeoMin+MiuMeoMax==0):'Pmeo'];%Lagrange������Pmeo��ƫ��Ϊ0----����Լ��
for i=2:numBuses
Cons_MPEC=[Cons_MPEC,(-EpsilonMin(i,:)+EpsilonMax(i,:)+Bn(i,:)*Lanta-Bf_Inv(i,:)*Vmin+Bf_Inv(i,:)*Vmax==0):'Detlait'];%Lagrange�����Բο��ڵ��������ƫ��Ϊ0
end
Cons_MPEC=[Cons_MPEC,(-EpsilonMin(1,:)+EpsilonMax(1,:)+Bn(1,:)*Lanta-Bf_Inv(1,:)*Vmin+Bf_Inv(1,:)*Vmax+Epsilon==0):'Detlat'];%Lagrange�����Բο��ڵ������ƫ��Ϊ0
disp("�����²�����ǿ��ż���⻥���ɳ�Լ������...");
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

disp("�����²���Ȼ���г�ԭ����Լ������...");
Cons_MPEC=[Cons_MPEC,(Qw>=0):'Beta1Min'];%��С����Լ��
Cons_MPEC=[Cons_MPEC,(Qw<=repmat(Qwmax,1,Time)):'Beta1Max'];%������Լ��
Cons_MPEC=[Cons_MPEC,(Qlg>=-repmat(Qlgmax,1,Time)):'Beta2Min'];%��С������·����Լ��
Cons_MPEC=[Cons_MPEC,(Qlg<=repmat(Qlgmax,1,Time)):'Beta2Max'];%���������·����Լ��
Cons_MPEC=[Cons_MPEC,(Qc>=0):'Beta3Min'];%��С������·����Լ��
Cons_MPEC=[Cons_MPEC,(Qc<=repmat(Qcmax,1,Time)):'Beta3Max'];%�����·����Լ��
Cons_MPEC=[Cons_MPEC,(Gmeo>=0):'Beta4Min'];%��С����Լ��
Cons_MPEC=[Cons_MPEC,(Gmeo<=Gmeomax):'Beta4Max'];%������Լ��
Cons_MPEC=[Cons_MPEC,(Anw*Qw-Anlg*Qlg-Anc*Qc-Anm*Gmeo'-Andg*Qdg==0):'Gama'];%�ڵ���Ȼ��ƽ��Լ��
disp("�����²���Ȼ���г���ż����Լ������...");
%��ż������ΧԼ��
Cons_MPEC=[Cons_MPEC,Beta1Min>=0];
Cons_MPEC=[Cons_MPEC,Beta1Max>=0];
Cons_MPEC=[Cons_MPEC,Beta2Min>=0];
Cons_MPEC=[Cons_MPEC,Beta2Max>=0];
Cons_MPEC=[Cons_MPEC,Beta3Min>=0];
Cons_MPEC=[Cons_MPEC,Beta3Max>=0];
Cons_MPEC=[Cons_MPEC,Beta4Min>=0];
Cons_MPEC=[Cons_MPEC,Beta4Max>=0];
%���߱���ϵ��Ϊ0Լ��
Cons_MPEC=[Cons_MPEC,(repmat(Cwgas,1,Time)-Beta1Min+Beta1Max-Anw'*Gama==0):'Qw'];%Lagrange������Qw��ƫ��Ϊ0
Cons_MPEC=[Cons_MPEC,(-Cmeogas-Beta4Min+Beta4Max+Gama'*Anm==0):'Gmer'];%Lagrange������Gmer��ƫ��Ϊ0
Cons_MPEC=[Cons_MPEC,(-Beta2Min+Beta2Max+Anlg'*Gama==0):'Qlg'];%Lagrange������Qlg��ƫ��Ϊ0
Cons_MPEC=[Cons_MPEC,(-Beta3Min+Beta3Max+Anc'*Gama==0):'Qc'];%Lagrange������Qc��ƫ��Ϊ0
disp("�����²���Ȼ���г������ɳ�����Լ��...");
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

%% ����������ģ�����
disp('�������...');
% tic;
Opt_MPEC=sdpsettings('verbose',1,'debug',1,'solver','cplex','savesolveroutput',1,'savesolverinput',1);
Opt_MPEC.cplex.exportmodel='Opt_MPEC.lp';
MPECModel=optimize(Cons_MPEC,-Obj_MPEC,Opt_MPEC);
if MPECModel.problem==0
    disp('MPECģ�Ϳɽ�');
toc;
tsolve=toc-tic;%��ȡ���ʱ��
else
disp('MPECģ�Ͳ��ɽ�');
end
%% ��ȡ������
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
Res_Case3_MPEC.ExpectedProfit=value(Obj_MPEC);%��������
Res_Case3_MPEC.WorstProfit=min(value(Z));%�������



Res_Case3_MPEC.Qw=value(Qw);
Res_Case3_MPEC.Qlg=value(Qlg);
Res_Case3_MPEC.Qc=value(Qc);
Res_Case3_MPEC.Gmeo=value(Gmeo);
Res_Case3_MPEC.Gama=value(Gama);%dual�����ż����ֵ��ʱ��Ҫȡ������Ҫ�Լ��ж���
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
