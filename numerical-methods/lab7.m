hx=1/10;
hy=1/10;
phi1=1;
phi5=0;
phi3=1;
phi7=0;
phi2=0;
phi6=1;
phi4=0;
phi8=1;
x=0;
y=0;
alfa1=1;
alfa2=0;
alfa3=0;
alfa4=0;
eps=0.001;

U = [];
X = [];
Y = [];
UN = [];

# Нулевая поверхность
for j=2:1:10
    f1=0;
    f3=1-y^2;
    U(1,j)=f1;
    U(11,j)=f3;
    y=y+hy;
end
y=0;
for i=1:1:11
    Y(end+1)=y; 
    X(end+1)=x;
    f2=0;
    f4=x^2-1;
    U(i,1)=f2;
    U(i,11)=f4;
    x=x+hx;
    y=y+hy; 
end
for i=2:1:10
    for j=2:1:10
        U(i,j)=U(1,j)+(U(11,j)-U(1,j))*(i-1)*hx;
    end
end

# Первая поверхность
x=0;
y=hy;
for j=2:1:10
    f1=0;
    f3=1-y^2;
    UN(1,j)=f1;
    UN(11,j)=f3;
    y=y+hy;
end
for i=1:1:11
    f2=0;
    f4=x^2-1;
    UN(i,1)=f2;
    UN(i,11)=f4;
    x=x+hx;
end
for i=2:1:10
    for j=2:1:10
        UN(i,j)= (1/(alfa4*(hx^2)-2*((hy^2+alfa1*(hx^2))/(hy^2))))*(((hx^2)*0)-(U(i+1,j)+U(i-1,j))-(((alfa1*(hx^2))/(hx^2))*(U(i,j+1)+U(i,j-1)))-(((alfa2*hx)/2)*(U(i+1,j)-U(i-1,j)))-(((alfa3*(hx^2))/(2*hy))*(U(i,j+1)-U(i,j-1))));
    end
end

#Следующие поверхности
M = [];
for i=1:1:11
    for j=1:1:11
        M(end+1)=abs(UN(i,j)-U(i,j));
    end
end
m=max(M);

while m > eps
    x=0;
    y=0;
    U=[];
    U(:,:)=UN(:,:);
    for j=2:1:10
        f1=0;
        f3=1-y^2;
        UN(1,j)=f1;
        UN(11,j)=f3;
        y=y+hy;
    end
    for i=1:1:11
        f2=0;
        f4=x^2-1;
        UN(i,1)=f2;
        UN(i,11)=f4;
        x=x+hx;
    end
    for i=2:1:10
        for j=2:1:10
            UN(i,j)= (1/(alfa4*(hx^2)-2*((hy^2+alfa1*(hx^2))/(hy^2))))*(((hx^2)*0)-(U(i+1,j)+U(i-1,j))-(((alfa1*(hx^2))/(hx^2))*(U(i,j+1)+U(i,j-1)))-(((alfa2*hx)/2)*(U(i+1,j)-U(i-1,j)))-(((alfa3*(hx^2))/(2*hy))*(U(i,j+1)-U(i,j-1))));
        end
    end
    M = [];
    for i=1:1:11
        for j=1:1:11
            M(end+1)=abs(UN(i,j)-U(i,j));
        end
    end
    m=max(M);
end

# срез по x
UPLOTX = [];
for i=1:1:11
    UPLOTX(end+1)=U(i,5);
end
# срез по y
UPLOTY = [];
for j=1:1:11
    UPLOTY(end+1)=U(5,j);
end

# аналитическое решения
UA=[];
xa=0;
ya=0;
for i=1:1:11
    for j=1:1:11
        UA(i,j)=(X(i)^2)-(Y(j)^2);
    end
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

figure('Name','Срез по y');
plot(X,UPLOTY,X,UPLOTYA); 
grid on;
xlabel('X');
ylabel('U');

figure('Name','Срез по x');
plot(Y,UPLOTX,Y,UPLOTXA); 
grid on;
xlabel('Y');
ylabel('U');

figure('Name','3D график');
surf(X,Y,U);

# погрешности для срезов по x и y
E1=[];
E2=[];
A=[];
B=[];
XP=[];
YP=[];
yp=0;
xp=0;

for n=1:1:11
    YP(n)=0;
end
for n=12:1:22
    YP(n)=yp;
    yp=yp+hy;
end
for n=1:1:11
    XP(n)=xp;
    xp=xp+hx;
end
for n=12:1:22
    XP(n)=0;
end
for n=1:1:11
    E1(n)=abs(UPLOTYA(n)-UPLOTY(n)); 
end
for n=12:1:22
    E1(n)=0;
end
for n=1:1:11
    E2(n)=0;
end
for n=12:1:22
    E2(n)=abs(UPLOTXA(n-11)-UPLOTX(n-11)); 
end
for k=1:1:11
    A(1,k)=XP(k);
    A(2,k)=YP(k);
    A(3,k)=E1(k);
end
for k=12:1:22
    B(1,k-11)=XP(k);
    B(2,k-11)=YP(k);
    B(3,k-11)=E2(k);        
end

fileID = fopen('absycut.txt','w');
fprintf(fileID,'%7s %7s %12s\n','|       |','|       |','|  absolute |');
fprintf(fileID,'%7s %7s %12s\n','|   x   |','|   y   |','|   error   |');
fprintf(fileID,'%7s %7s %12s\n','|       |','|       |','|   y-cut   |');
fprintf(fileID,'%7s %7s %12s\n','|=======|','|=======|','|===========|');
fprintf(fileID,'|%6.5f| |%6.5f| |%10.9f|\n',A);
fclose(fileID);
type absycut.txt;

fileID = fopen('absxcut.txt','w');
fprintf(fileID,'%7s %7s %12s\n','|       |','|       |','|  absolute |');
fprintf(fileID,'%7s %7s %12s\n','|   x   |','|   y   |','|   error   |');
fprintf(fileID,'%7s %7s %12s\n','|       |','|       |','|   x-cut   |');
fprintf(fileID,'%7s %7s %12s\n','|=======|','|=======|','|===========|');
fprintf(fileID,'|%6.5f| |%6.5f| |%10.9f|\n',B);
fclose(fileID);
type absxcut.txt;




