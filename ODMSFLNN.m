
clc
clear

run_number =3;%独立运行次数
N = 30000;%样本个数
SA_N = 8;%子带数目
M = 65;
LLL = 40;%输入向量长度
P=2;%三角拓展维数

LL = 2*LLL*P+LLL+1;%拓展后的长度
y = zeros(1,N);
y1=zeros(1,N);
E1=zeros(1,N);
Input = zeros(1,LL);

HK=zeros(SA_N,M);


E2=zeros(1,N);
d_1=zeros(1,N);

%未知系统的冲击响应
%装载分析与综合滤波器组4*64
load filter_bank_8.mat % hk：分析滤波器, fk：综合滤波器，子带数4

for m = 1:M
    HK(:,m)=hk(:,m);
end

for i = 1:run_number
    i
    TFLNN_w=zeros(LL,1);%自适应滤波器权系数
    input1 = unifrnd(-1,1,1,N);%输入一个在（-1，1）之间服从均匀分布的白色输入信号
    input = filter(1,[1,-0.9],input1);
    %非线性系统1
    for ii = 1:N
        if ii==1
            y1(ii)=0.9*input(ii)+0.8*sin(pi*input(ii));
        elseif ii==2
            y1(ii)=0.9*input(ii)+0.8*sin(pi*input(ii))+0.5*cos(pi*input(ii-1));
        elseif ii==3
            y1(ii)=0.9*input(ii)+0.8*sin(pi*input(ii))+0.5*cos(pi*input(ii-1))+0.1*sin(pi*input(ii-2));
        elseif ii==4
            y1(ii)=0.9*input(ii)+0.8*sin(pi*input(ii))+0.5*cos(pi*input(ii-1))+0.1*sin(pi*input(ii-2));
        elseif (5<=ii) && (ii<=15)
            y1(ii)=0.9*input(ii)+0.8*sin(pi*input(ii))+0.5*cos(pi*input(ii-1))+0.1*sin(pi*input(ii-2))+0.08*cos(pi*input(ii-4));
        elseif ii>15
            y1(ii)=0.9*input(ii)+0.8*sin(pi*input(ii))+0.5*cos(pi*input(ii-1))+0.1*sin(pi*input(ii-2))+0.08*cos(pi*input(ii-4))-0.3*cos(pi*input(ii-15));
        end
    end
    
    d_n = awgn(y1,30);
    
    noise = var(d_n-y1);
    for ii = 1:SA_N
        noise_s(ii) = noise * norm(HK(ii,:),2)^2;
    end
    
    flag=0;
    Cw = 100*ones(LL,1)/LL;
    N_e = zeros(1,SA_N);
    V_d = zeros(LL,1);
    V_y = zeros(LL,1);
    N_nois = zeros(SA_N,N);
    v2 = zeros(LL,1);
    
    for p = 3*LLL-1:N
        D_n1 = y1(p:-1:p-M+1);%未加噪声的期望%??????????????????????????要加63？
        D_n = d_n(p:-1:p-M+1);%有噪声的期望
        
        D=zeros(P,2*LLL);%用于存放拓展后的三角元素
        B=zeros(P,LLL);%用于存放拓展后的sin元素
        C=zeros(P,LLL);%用于存放拓展后的cos元素
        S = zeros(LL,M);
        
        for  pp = 0:(M-1)
            %         Input = input(pp:-1:pp-LLL+1);%取出每次需要延迟的分量
            Input = input(p-pp:-1:p-pp-LLL+1);
            
            for j=1:P
                for jj=1:LLL
                    B(j,jj)=sin(j*pi*Input(jj));
                    C(j,jj)=cos(j*pi*Input(jj));
                end
                D(j,:)=[B(j,:) C(j,:)];
            end
            
            for q = 1:LLL
                F(1,q)=D(1,q);
                F(2,q)=D(1,LLL+q);
                F(3,q)=D(2,q);
                F(4,q)=D(2,LLL+q);
            end
            F1 = reshape(F,1,[]);%将拓展后元素按顺序排列
            F2 = [1 Input F1]';
            S(:,pp+1)=F2;%输入矩阵
        end
        
        %主程序
        s_h=S*HK';%子带输入信号
        d_h=HK*D_n';%子带期望信号
        
        if flag == 0 || flag == SA_N
            flag = 1;
            TFLNN_del_w= zeros(LL,1);
            Cd = 0;
            vf = 1-SA_N/(1*LL);
            for kk = 1:SA_N
                
                Y(kk) = TFLNN_w'*s_h(:,kk);
                GR_e(kk,p) = d_h(kk)-Y(kk);
                
                ff = 1-1/(1*LL);
                N_e(kk) = ff*N_e(kk) + (1-ff)*GR_e(kk,p)^2;
                V_d(kk) = ff*V_d(kk) + (1-ff)*d_h(kk)^2;
                V_y(kk) = ff*V_y(kk) + (1-ff)*Y(kk)^2;
                N_nois(kk,p) = V_d(kk)*N_e(kk)/(N_e(kk)+V_y(kk));
                if N_nois(kk,p) <= 0
                    N_nois(kk,p) = -N_nois(kk,p);
                    
                end
                
                g(:,kk) = Cw.*s_h(:,kk) / (sum(s_h(:,kk).*Cw.*s_h(:,kk))+N_nois(kk,p));
                TFLNN_del_w = TFLNN_del_w + g(:,kk)*GR_e(kk,p);
                Cd = Cd + s_h(:,kk).*g(:,kk);
                
                
                
            end
            %%way1
            v2 = 0.985*v2 + 0.015*TFLNN_del_w.^2;
            Cw = Cw - Cd .* Cw + min(norm(TFLNN_del_w)^2/LL,v2);
            
            %%way2
            %             Cw = Cw - Cd .* Cw + norm(TFLNN_del_w)^2/LL;
            %%way3
            %             v2 = 0.9*v2 + 0.1*TFLNN_del_w.^2;
            %             Cw = Cw - Cd .* Cw + v2;
            
            d_1(p)=S(:,1)'*TFLNN_w;
            E1(p)=D_n1(1)-d_1(p);
            TFLNN_w = TFLNN_w + TFLNN_del_w;
        else
            
            flag = flag + 1;
            TFLNN_w = TFLNN_w;
            d_1(p)=S(:,1)'*TFLNN_w;
            E1(p)=D_n1(1)-d_1(p);
        end
    end
    E2 = E2 + E1.*E1;
end

E1_vector1_new = E2/run_number;

figure(3)
plot(10*log10(E1_vector1_new));
hold on
xlabel('Iteration number','fontsize',16)
ylabel('EMSE','fontsize',16)
legend('DMSFLNN-GR')