
clc
clear
% close all;
run_number = 3;%独立运行次数
N = 30000;%样本个数
SA_N = 8;%子带数目
% L = 64;%滤波器的长度
LLL = 40;%输入向量长度
P=2;%三角拓展维数
LL = 2*LLL*P+LLL+1;%拓展后的长度
y = zeros(1,N);
y1=zeros(1,N);
E3=zeros(1,N);
Input = zeros(1,LL);
% S = zeros(LL,64);
HK=zeros(SA_N,LLL);
E3 = zeros(1,N);
u=1;
E2=zeros(1,N);
d_1=zeros(1,N);
E3=zeros(1,N);

%未知系统的冲击响应
%装载分析与综合滤波器组4*64
load filter_bank_8.mat % hk：分析滤波器, fk：综合滤波器，子带数4

for m = 1:LLL
    HK(:,m)=hk(:,m);
end

for i = 1:run_number
    i
    TFLNN_w=zeros(LL,1);%自适应滤波器权系数
    input1 = unifrnd(-1,1,1,N);%输入一个在（-1，1）之间服从均匀分布的白色输入信号
    input = filter(1,[1,-0.1],input1);
    %     input1=randi([-1,1],1,N);

    %******************将输入信号通过unknow nonlinear system
    for ii=1:N
        q(ii) = (2/3)*input(ii) - (3/10)*input(ii)^2;
        if q(ii)>0
            rho=4;
            y1(ii) = 2*(1/(1+exp(-1*rho*q(ii))))-1;
        else
            rho=1/2;
            y1(ii) = 2*(1/(1+exp(-1*rho*q(ii))))-1;
        end
    end

%     for ii=1:N
%         if ii==1
%             y1(ii)=abs(input(ii))^3;
%         else
%             y1(ii)=y1(ii-1)/(1+y1(ii-1)^2)+abs(input(ii))^3;
%         end
%     end


    d_n = awgn(y1,30);
    flag=0;
    for p = LLL:N-(LLL-1)
        D_n1 = y1(p+(LLL-1):-1:p+(LLL-1)-LLL+1);%未加噪声的期望%??????????????????????????要加63？
        D_n = d_n(p+(LLL-1):-1:p+(LLL-1)-LLL+1);%有噪声的期望

        D=zeros(P,2*LLL);%用于存放拓展后的三角元素
        B=zeros(P,LLL);%用于存放拓展后的sin元素
        C=zeros(P,LLL);%用于存放拓展后的cos元素
        S = zeros(LL,LLL);

        for  pp = 0:(LLL-1)
            %         Input = input(pp:-1:pp-LLL+1);%取出每次需要延迟的分量
            Input = input(p+pp:-1:p+pp-LLL+1);

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
            S(:,LLL-pp)=F2;%输入矩阵
        end

        s_h=S*HK';%子带输入信号
        d_h=HK*D_n';%子带期望信号

        if flag == 0 || flag == 8
            flag = 1;
            TFLNN_del_w= zeros(LL,1);
            for kk = 1:SA_N
                Y(kk) = TFLNN_w'*s_h(:,kk);
                e(kk) = d_h(kk)-Y(kk);
                TFLNN_del_w = TFLNN_del_w + (e(kk)*s_h(:,kk))/(s_h(:,kk)'*s_h(:,kk));
            end
            d_1(p)=S(:,1)'*TFLNN_w;
            E3(p)=D_n1(1)-d_1(p);
            TFLNN_w = TFLNN_w + u * TFLNN_del_w;
        else
            flag = flag + 1;
            TFLNN_w = TFLNN_w;
            d_1(p)=S(:,1)'*TFLNN_w;
            E3(p)=D_n1(1)-d_1(p);
        end
    end
    E2 = E2 + E3.*E3;
end

non2_large = E2/run_number;
%性能指标
figure(1)
plot(10*log10(non2_large));
hold off
xlabel('Iteration number','fontsize',16)
ylabel('EMSE','fontsize',16)
legend('DMSFLNN')