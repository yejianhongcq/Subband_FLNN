% 任意参数稳定分布数据产生
function Y=rasd(N,alpha,beta,gama,miu)

% ----- 均匀分布 -----
nuni=rand(1,N);
U=(nuni*pi)-pi/2;
% ----- 指数分布 -----
zhi=rand(1,N);
W=-log(zhi);

if alpha~=1
   X1=S(alpha,beta);
   X2=sin(alpha*(U+B(alpha,beta)))./(cos(U)).^(1/alpha);
   X3=(cos(U-alpha*(U+B(alpha,beta)))./W).^(1/alpha-1);
   X=X1.*X2.*X3;
else
    X=(2/pi)*( (pi/2+beta*U).*tan(U)-beta*log( 0.5*pi*W.*cos(U)./(pi/2+beta*U) ) );
end

if alpha~=1
    Y=gama*X+miu;
else
    Y=gama*X+miu+(2/pi)*beta*gama*log(gama);
end
function y=B(alpha,beta)
y=atan(beta*tan(pi*alpha/2))/alpha;
function y=S(alpha,beta)
y=(1 + beta^2 * (tan(pi*alpha/2)) ^ 2  ) ^   (1/(2*alpha));

