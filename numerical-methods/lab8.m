mu1=1;
mu2=1;
hx=pi*mu1/20;
hy=pi*mu2/20;
ht=(hx^2)/2;
x=0;
y=0;
t=0;
tau=ht;
Tend=ht*10;
alfa=1;
alfa1=alfa;
alfa2=alfa;
alfa3=0;
alfa4=0;
alfa5=0;

U = [];
X = [];
Y = [];
UN = [];

# массивы X и Y
for k=1:1:11
    X(k)=x;
    Y(k)=y;
    x=x+hx;
    y=y+hy;
end
# Нулевая поверхность
for i=1:1:11
    for j=1:1:11
        U(i,j)=(X(i))*(Y(j));
    end
end
U1=[];
U1(:,:)=U(:,:);

# Метод переменных направлений
# Остальные поверхности
k=0;
for t=tau:ht:Tend
    # промежуточная поверхность t=t+ht/2
    for j=2:1:10
        f2=cos(mu2*Y(j))*exp(-(mu1^2+mu2^2)*alfa*t/2);
        f3=0;
        UN(1,j)=f2;
        UN(11,j)=f3;
    end
    for i=1:1:11
        f5=0;
        f4=cos(mu1*X(i))*exp(-(mu1^2+mu2^2)*alfa*t/2);
        UN(i,1)=f4;
        UN(i,11)=f5;
    end
    for j=2:1:10
        y=Y(j);
        for i=2:1:10
            M=zeros(10,10);
            N=zeros(10);
            P=[];
            Q=[];
            x=X(i);
            for k=2:1:9      
                M(k,k-1)= (alfa1-((hx*alfa3)/2)); 
                M(k,k)=((hx^2)-(2*(hx^2))/ht-2*alfa1); 
                M(k,k+1)= (alfa1+((hx*alfa3)/2)); 
                N(k)= (-(2*hx^2)/ht)*U(k,j)-((alfa2*(hx^2))/(hy^2))*(U(k,j+1)-2*U(k,j)+U(k,j-1))-((alfa4*(hx^2))/(2*(hy^2)))*(U(k,j+1)-U(k,j-1))+hx^2*(0);% - вектор правой части системы
            end
            M(1,1)=((hx^2)-(2*(hx^2))/ht-2*alfa1); 
            M(1,2)= (alfa1+((hx*alfa3)/2));
            M(10,9)=(alfa1-((hx*alfa3)/2)); 
            M(10,10)=((hx^2)-(2*(hx^2))/ht-2*alfa1);
            N(1)= (-(2*hx^2)/ht)*U(1,j)-((alfa2*(hx^2))/(hy^2))*(U(1,j+1)-2*U(1,j)+U(1,j-1))-((alfa4*(hx^2))/(2*(hy^2)))*(U(1,j+1)-U(1,j-1))+hx^2*(0)-(alfa1-((hx*alfa3)/2))*U(1,j);
            N(10)= (-(2*hx^2)/ht)*U(10,j)-((alfa2*(hx^2))/(hy^2))*(U(10,j+1)-2*U(10,j)+U(10,j-1))-((alfa4*(hx^2))/(2*(hy^2)))*(U(10,j+1)-U(10,j-1))+hx^2*(0)-(alfa1+((hx*alfa3)/2))*U(11,j);
            
            P(1)= -M(1,2)/M(1,1);
            Q(1)= N(1)/M(1,1);
            for l=2:1:9
                P(l)= -M(l,l+1)/(M(l,l)+M(l,l-1)*P(l-1));
                Q(l)= (N(l)-M(l,l-1)*Q(l-1))/(M(l,l)+M(l,l-1)*P(l-1));
            end
            P(10)=0;
            Q(10)= (N(10)-M(10,9)*Q(9))/(M(10,10)+M(10,9)*P(9));
            UN(10,j)= Q(10);
            for k=2:1:9
                UN(11-k,j)= Q(11-k)+P(11-k)*UN(11-k+1,j);
            end
        end      
    end
    # поверхность для t  
    for j=2:1:10
        f2=cos(mu2*Y(j))*exp(-(mu1^2+mu2^2)*alfa*t/2);
        f3=0;
        U(1,j)=f2;
        U(11,j)=f3;
    end
    for i=1:1:11
        f5=0;
        f4=cos(mu1*X(i))*exp(-(mu1^2+mu2^2)*alfa*t/2);
        U(i,1)=f4;
        U(i,11)=f5;
    end
    for i=2:1:10
        x=X(i);
        for j=2:1:10
            M=zeros(10,10);
            N=zeros(10);
            P=[];
            Q=[];
            y=Y(j);            
            for k=2:1:9        
                M(k,k-1)= (alfa2-((hy*alfa4)/2)); 
                M(k,k)=(-(2*(hy^2))/ht-2*alfa2); 
                M(k,k+1)= (alfa2+((hy*alfa4)/2)); 
                N(k)= ((-(2*hy^2)/ht-(hy^2)))*UN(i,k)-((alfa1*(hy^2))/(hx^2))*(UN(i+1,k)-2*UN(i,k)+UN(i-1,k))-((alfa3*(hy^2))/(2*hx))*(UN(i+1,k)-UN(i-1,k))+hy^2*(0);% - вектор правой части системы
            end
            M(1,1)=(-(2*(hy^2))/ht-2*alfa2); 
            M(1,2)= (alfa2+((hy*alfa4)/2));
            M(10,9)=(alfa2-((hy*alfa4)/2)); 
            M(10,10)=(-(2*(hy^2))/ht-2*alfa2);
            N(1)= ((-(2*hy^2)/ht-(hy^2)))*UN(i,1)-((alfa1*(hy^2))/(hx^2))*(UN(i+1,1)-2*UN(i,1)+UN(i-1,1))-((alfa3*(hy^2))/(2*hx))*(UN(i+1,1)-UN(i-1,1))-hy^2*(0)-(alfa2-((hy*alfa4)/2))*UN(i,1);
            N(10)= ((-(2*hy^2)/ht-(hy^2)))*UN(i,10)-((alfa1*(hy^2))/(hx^2))*(UN(i+1,10)-2*UN(i,k)+UN(i-1,10))-((alfa3*(hy^2))/(2*hx))*(UN(i+1,10)-UN(i-1,10))-hy^2*(0)-(alfa2+((hy*alfa4)/2))*UN(i,11);
            
            P(1)= -M(1,2)/M(1,1);
            Q(1)= N(1)/M(1,1);
            for l=2:1:9
                P(l)= -M(l,l+1)/(M(l,l)+M(l,l-1)*P(l-1));
                Q(l)= (N(l)-M(l,l-1)*Q(l-1))/(M(l,l)+M(l,l-1)*P(l-1));
            end
            P(10)=0;
            Q(10)= (N(10)-M(10,9)*Q(9))/(M(10,10)+M(10,9)*P(9));
            U(i,10)= Q(10);
            for k=2:1:9
                U(i,11-k)= Q(11-k)+P(11-k)*U(i,11-k+1);
            end
        end      
    end
    k=k+1;
end

# метод дробных шагов
# остальные поверхности
k=0;
for t=tau:ht:Tend
    # промежуточная поверхность t=t+ht/2
    for j=2:1:10
        f2=cos(mu2*Y(j))*exp(-(mu1^2+mu2^2)*alfa*t/2);
        f3=0;
        UN1(1,j)=f2;
        UN1(11,j)=f3;
    end
    for i=1:1:11
        f5=0;
        f4=cos(mu1*X(i))*exp(-(mu1^2+mu2^2)*alfa*t/2);
        UN1(i,1)=f4;
        UN1(i,11)=f5;
    end
    for j=2:1:10
        y=Y(j);
        for i=2:1:10
            M=zeros(10,10);
            N=zeros(10);
            P=[];
            Q=[];
            x=X(i);
            for k=2:1:9       
                M(k,k-1)= (alfa1-((hx*alfa3)/2)); 
                M(k,k)=(-2*alfa1-((hx^2)/ht));
                M(k,k+1)= (alfa1+((hx*alfa3)/2)); 
                N(k)= -U1(k,j)*((hx^2)/ht);
            end
            M(1,1)=(-2*alfa1-((hx^2)/ht)); 
            M(1,2)= (alfa1+((hx*alfa3)/2));
            M(10,9)=(alfa1-((hx*alfa3)/2)); 
            M(10,10)=(-2*alfa1-((hx^2)/ht));
            N(1)=-U1(1,j)*((hx^2)/ht)-(alfa1-((hx*alfa3)/2))*U1(1,j);
            N(10)=-U1(10,j)*((hx^2)/ht)-(alfa1+((hx*alfa3)/2))*U1(11,j);
            
            P(1)= -M(1,2)/M(1,1);
            Q(1)= N(1)/M(1,1);
            for l=2:1:9
                P(l)= -M(l,l+1)/(M(l,l)+M(l,l-1)*P(l-1));
                Q(l)= (N(l)-M(l,l-1)*Q(l-1))/(M(l,l)+M(l,l-1)*P(l-1));
            end
            P(10)=0;
            Q(10)= (N(10)-M(10,9)*Q(9))/(M(10,10)+M(10,9)*P(9));
            UN1(10,j)= Q(10);
            for k=2:1:9
                UN1(11-k,j)= Q(11-k)+P(11-k)*UN1(11-k+1,j);
            end
        end      
    end
    # поверхность для t  
    for j=2:1:10
        f2=0;
        f3=cos(mu2*Y(j))*exp(-(mu1^2+mu2^2)*alfa*t);
        U1(1,j)=f2;
        U1(11,j)=f3;
    end
    for i=1:1:11
        f5=0;
        f4=cos(mu1*X(i))*exp(-(mu1^2+mu2^2)*alfa*t);
        U1(i,1)=f4;
        U1(i,11)=f5;
    end
    for i=2:1:10
        x=X(i);
        for j=2:1:10
            M=zeros(10,10);
            N=zeros(10);
            P=[];
            Q=[];
            y=Y(j);            
            for k=2:1:9       
                M(k,k-1)= alfa2; 
                M(k,k)=((-(hy^2))/ht)-2*alfa2; 
                M(k,k+1)= alfa2;
                N(k)= -((hy^2)/ht)*UN1(i,k)+((alfa4*hy)/2)*(UN1(i-1,k)-UN1(i+1,k));
            end
            M(1,1)=((-(hy^2))/ht)-2*alfa2; 
            M(1,2)= alfa2;
            M(10,9)= alfa2; 
            M(10,10)=((-(hy^2))/ht)-2*alfa2;
            N(1)=-((hy^2)/ht)*UN1(i,1)+((alfa4*hy)/2)*(UN1(i-1,1)-UN1(i+1,1))-alfa1*UN1(i,1);
            N(10)=-((hy^2)/ht)*UN1(i,10)+((alfa4*hy)/2)*(UN1(i-1,10)-UN1(i+1,10))-alfa1*UN1(i,11);
            
            P(1)= -M(1,2)/M(1,1);
            Q(1)= N(1)/M(1,1);
            for l=2:1:9
                P(l)= -M(l,l+1)/(M(l,l)+M(l,l-1)*P(l-1));
                Q(l)= (N(l)-M(l,l-1)*Q(l-1))/(M(l,l)+M(l,l-1)*P(l-1));
            end
            P(10)=0;
            Q(10)= (N(10)-M(10,9)*Q(9))/(M(10,10)+M(10,9)*P(9));
            U1(i,10)= Q(10);
            for k=2:1:9
                U1(i,11-k)= Q(11-k)+P(11-k)*U1(i,11-k+1);
            end
        end      
    end
    k=k+1;
end


# аналитическое решения
UA=[];
xa=0;
ya=0;
ta=0;
for i=1:1:11
    for j=1:1:11
        UA(i,j)= cos(mu1*X(i))*cos(mu2*Y(j))*exp(-(mu1^2+mu2^2)*alfa*ta);
    end
    ta=ta+ht;
end

# срез по x для мпн
UPLOTX = [];
for i=1:1:11
    UPLOTX(end+1)=U(i,5);
end
# срез по y для мпн
UPLOTY = [];
for j=1:1:11
    UPLOTY(end+1)=U(5,j);
end

# срез по x для мдш
UPLOTX1 = [];
for i=1:1:11
    UPLOTX1(end+1)=U1(i,5);
end
# срез по y для мдш
UPLOTY1 = [];
for j=1:1:11
    UPLOTY1(end+1)=U1(5,j);
end

# срез по x (аналитическое решение)
UPLOTXA = [];
for i=1:1:11
    UPLOTXA(end+1)=UA(i,5);
end
# серз по y (аналитическое решение)
UPLOTYA = [];
for j=1:1:11
    UPLOTYA(end+1)=UA(5,j);
end

figure('Name','Срезы по y');
subplot(1,2,1);
plot(X,UPLOTY,X,UPLOTYA);
grid on;
xlabel('X');
ylabel('МПН');

subplot(1,2,2);
plot(X,UPLOTY1,X,UPLOTYA);
grid on;
xlabel('X');
ylabel('МДШ');


figure('Name','Срезы по x');
subplot(1,2,1);
plot(Y,UPLOTX,Y,UPLOTXA);
grid on;
xlabel('Y');
ylabel('МПН');

subplot(1,2,2);
plot(Y,UPLOTX1,Y,UPLOTXA);
grid on;
xlabel('Y');
ylabel('МДШ');


figure('Name','3D график МПН');
surf(X,Y,UA);

figure('Name','3D график МДШ');
surf(X,Y,UA);
