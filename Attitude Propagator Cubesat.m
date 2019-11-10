%--------------------------------------------------------------------------
%SOLUCIÓN DEL MOVIMIENTO DE ROTACIÓN DE UN CUBESAT EN ÓRBITA 
%LEO ELÍPTICA SOMETIDO A TORQUE GRAVITATORIO CON LA PARAMETRIZACION EN 
%TÉRMINOS DE CUATERNIONES Y VELOCIDADES ANGULARES. 
%PROPAGADOR ORBITAL: TWO BODY

%Autor: Álvaro Fernández Villar                
%--------------------------------------------------------------------------
function propagadorcuaternion
clear;
clc;
global mu;
global A;
global B;
global C;
global n;

%Componentes principales de la matriz de inercia 
A=7.380;
B=7.475;
C=2.155;

%Condiciones iniciales en radianes para los angulos y en rad/s para las velocidades angulares 
wx_inicial=0.0000106;
wy_inicial=0.0000106;
wz_inicial=0.0000106;
roll_inicial=0;
pitch_inicial=0;
yaw_inicial=0;
    
%Parametros orbitales (angulos en radianes)
mu=398600.441; %Parametro gravitacional terrestre (km^3/m^2)
h=54427; %Este valor para el momento cinetico se corresponde con una orbita de semieje mayor a=7078 km. 
e=0.05;
RA=0.6981;
i=1.7104;
w=1.0472;
TA=0;


%Tiempo de simulacion 
T=(2*pi()/mu^2)*(h/(sqrt(1-e^2)))^3; %Periodo orbital
n=2*pi()/T;
np=0.001; %Periodos de simulacion
tsimul=round(np*T); %Tiempo de simulacion

%Elementos orbitales

coe=[h,e,RA,i,w,TA]; %Vector de elementos orbitales
R=6378; %Radio de la Tierra en km

%Calculo del diagrama de Estabilidad para órbita circular
if e==0 
    K1 = ((A) - (B))/(C);
    K2 = ((B) - (C))/(A);
    K3 = -((K1 + K2)/(1 +K1*K2));
 

% Diagrama de Estabilidad
x1=linspace(-1,1,1000);
 
    for i=1:length(x1)
    y1(i)=-x1(i);
    end
 
figure
plot(x1,y1, K1, K2, 'xr','Markersize',10)
axis([-1 1 -1 1])
grid off
hold on
line([-1,1],[0,0])
line([0,0],[-1,1])
hold off
title('Diagrama de Estabilidad')
xlabel('K1')
ylabel('K2')

end

%Obtener el vector de estado a partir de los elementos orbitales
[r,v]=vectorinercial(coe,mu);

%Generacion del vector tiempos
t=0:1:tsimul;

%Propagacion de la orbita
for k=1:tsimul+1
        instante=t(k);
        [r,v]=vectorestadotiempo(r, v, instante);
        coe_instantaneo=coe_vectorestado(r,v,mu);
    for i=1:3
        matrizposicion(k,i)=r(1,i);
    end
    for j=1:6
        matrizcoe(k,j)=coe_instantaneo(1,j);
    end
end


%Conversión angulos de euler al cuaternion segun rotacion 3-2-1

%matriz de cosenos directores
E1=[1 0 0; 0 cos(roll_inicial) sin(roll_inicial); 0 -sin(roll_inicial) cos(roll_inicial)];
E2=[cos(pitch_inicial) 0 -sin(pitch_inicial); 0 1 0; sin(pitch_inicial) 0 cos(pitch_inicial)];
E3=[cos(yaw_inicial) sin(yaw_inicial) 0; -sin(yaw_inicial) cos(yaw_inicial) 0; 0 0 1];
R=E1*E2*E3;
%Obtencion del cuaternion
q = q_dcm(R);

%Vector de condiciones iniciales de integracion
inicial=[wx_inicial wy_inicial wz_inicial q(1) q(2) q(3) q(4)];

%Integracion de las ecuaciones dinamicas de rotacion
for i=1:tsimul
    [t_out,y_out]=ode45(@ecseulergrav,[t(i) t(i)+1],[inicial(1) inicial(2) inicial(3) inicial(4) inicial(5) inicial(6) inicial(7) A B C],odeset(),matrizposicion(i,:));
    [filassolucionparcial,columnassolucionparcial]=size(y_out);
    inicial=[y_out(filassolucionparcial,1) y_out(filassolucionparcial,2) y_out(filassolucionparcial,3) y_out(filassolucionparcial,4) y_out(filassolucionparcial,5) y_out(filassolucionparcial,6) y_out(filassolucionparcial,7)];    
    matrizsolucionparcial=[t_out y_out];
    for tamanosol=1:filassolucionparcial
        for l=1:8
        matrizsolucionfinal(tamanosol,l,i)=matrizsolucionparcial(tamanosol,l,1);
        end
    end
   
end

%transformacion a matrices unicas segun componentes
for i=1:filassolucionparcial
    for j=1:tsimul
        tiempo_matrizsolucionfinal3d(i,1,j)=matrizsolucionfinal(i,1,j);    
        wx_matrizsolucionfinal3d(i,1,j)=matrizsolucionfinal(i,2,j);
        wy_matrizsolucionfinal3d(i,1,j)=matrizsolucionfinal(i,3,j);
        wz_matrizsolucionfinal3d(i,1,j)=matrizsolucionfinal(i,4,j);
        q1_matrizsolucionfinal3d(i,1,j)=matrizsolucionfinal(i,5,j);
        q2_matrizsolucionfinal3d(i,1,j)=matrizsolucionfinal(i,6,j);
        q3_matrizsolucionfinal3d(i,1,j)=matrizsolucionfinal(i,7,j);
        q0_matrizsolucionfinal3d(i,1,j)=matrizsolucionfinal(i,8,j);
    end
end
[filassfinal3d,c]=size(wx_matrizsolucionfinal3d);

%Conversion a matrices de dos dimensiones
for i=1:filassfinal3d
    for k=1:tsimul
        tiempo_matrizsolucionfinal(i,k,1)=tiempo_matrizsolucionfinal3d(i,1,k);
        wx_matrizsolucionfinal(i,k,1)=wx_matrizsolucionfinal3d(i,1,k);
        wy_matrizsolucionfinal(i,k,1)=wy_matrizsolucionfinal3d(i,1,k);
        wz_matrizsolucionfinal(i,k,1)=wz_matrizsolucionfinal3d(i,1,k);
        q1_matrizsolucionfinal(i,k,1)=q1_matrizsolucionfinal3d(i,1,k);
        q2_matrizsolucionfinal(i,k,1)=q2_matrizsolucionfinal3d(i,1,k);
        q3_matrizsolucionfinal(i,k,1)=q3_matrizsolucionfinal3d(i,1,k);
        q0_matrizsolucionfinal(i,k,1)=q0_matrizsolucionfinal3d(i,1,k);
    end
end

%Conversion a vectores
[filasv,cl]=size(tiempo_matrizsolucionfinal);
tiempo_matrizsolucionfinal=tiempo_matrizsolucionfinal(:);
wx_matrizsolucionfinal=wx_matrizsolucionfinal(:);
wy_matrizsolucionfinal=wy_matrizsolucionfinal(:);
wz_matrizsolucionfinal=wz_matrizsolucionfinal(:);
q1_matrizsolucionfinal=q1_matrizsolucionfinal(:);
q2_matrizsolucionfinal=q2_matrizsolucionfinal(:);
q3_matrizsolucionfinal=q3_matrizsolucionfinal(:);
q0_matrizsolucionfinal=q0_matrizsolucionfinal(:);

%Obtencion de los angulos de euler 
[f]=size(q1_matrizsolucionfinal);
for i=1:f
    q=[q1_matrizsolucionfinal(i) q2_matrizsolucionfinal(i) q3_matrizsolucionfinal(i) q0_matrizsolucionfinal(i)];
    Q = dcm_q(q);
    [yaw pitch roll] = dcm_to_ypr(Q);
    roll_matrizsolucionfinal(i)=roll;
    pitch_matrizsolucionfinal(i)=pitch;
    yaw_matrizsolucionfinal(i)=yaw;
end

%Obtencion de la evolucion del modulo del cuaternion
%Modulo del cuaternion
min=3;
suma_componentescuaternion=q0_matrizsolucionfinal.^2+q1_matrizsolucionfinal.^2+q2_matrizsolucionfinal.^2+q3_matrizsolucionfinal.^2;

for i=1:length(suma_componentescuaternion)
    if suma_componentescuaternion(i)<min
        min=suma_componentescuaternion(i);
    end
end


%Resultados
figure 
plot(tiempo_matrizsolucionfinal/T,(roll_matrizsolucionfinal)) %dibuja la posición
title('Alabeo')
ylabel('Angulos (deg)')
xlabel('nº de Periodos')
xlim([0 np])
roll_maximo=max(((abs(roll_matrizsolucionfinal))))

figure 
plot(tiempo_matrizsolucionfinal/T,rad2deg(wx_matrizsolucionfinal)) %dibuja la posición
title('Velocidad angular eje X')
ylabel('(deg/s)')
xlabel('nº de Periodos')
xlim([0 np])
wx_maximo=max((abs(rad2deg(wx_matrizsolucionfinal))))

figure 
plot(tiempo_matrizsolucionfinal/T,(yaw_matrizsolucionfinal)) %dibuja la posición
title('Guiñada')
ylabel('Angulos (deg)')
xlabel('nº de Periodos')
xlim([0 np])
yaw_maximo=max(((abs(yaw_matrizsolucionfinal))))

figure 
plot(tiempo_matrizsolucionfinal/T,rad2deg(wz_matrizsolucionfinal)) %dibuja la posición
title('Velocidad angular eje Z')
ylabel('(deg/s)')
xlabel('nº de Periodos')
xlim([0 np])
wz_maximo=max((abs(rad2deg(wz_matrizsolucionfinal))))

figure 
plot(tiempo_matrizsolucionfinal/T,rad2deg(pitch_matrizsolucionfinal)) %dibuja la posición
title('Cabeceo')
ylabel('Angulos (deg)')
xlabel('nº de Periodos')
xlim([0 np])
pitch_maximo=max(((abs(pitch_matrizsolucionfinal))))



figure 
plot(tiempo_matrizsolucionfinal/T,rad2deg(wy_matrizsolucionfinal)) %dibuja la posición
title('Velocidad angular eje Y')
ylabel('(deg/s)')
xlabel('nº de Periodos')
xlim([0 np])
wy_maximo=max((abs(rad2deg(wy_matrizsolucionfinal))))

figure 
plot(tiempo_matrizsolucionfinal/T,rad2deg(wx_matrizsolucionfinal),'r')
hold on;
plot(tiempo_matrizsolucionfinal/T,rad2deg(wy_matrizsolucionfinal),'m')
plot(tiempo_matrizsolucionfinal/T,rad2deg(wz_matrizsolucionfinal),'b')
title('Componentes velocidad angular')
legend('wx','wy','wz')
ylabel('(deg/s)')
xlabel('nº de Periodos')
xlim([0 np])

figure
plot(tiempo_matrizsolucionfinal/T,log(suma_componentescuaternion));
ylabel('Valor del cuaternion')
xlabel('Periodos')
xlim([0 np])
title('Evolucion logarítmica del cuaternion');
maxima_varmodulocuaternion=max(max(max(suma_componentescuaternion))-1,1-min) %maxima variacion del modulo del cuaternion

if e==0
    
    %Frecuencia de oscilacion estable de cabeceo
 
wf_cabeceo=n*sqrt(3*(A-C)/B)
 
%Frecuencia de oscilacion estable de guiñada
 
wf_guinada=n*sqrt(0.5*((3*((B-C)/A)+((B-A)/C)*((B-C)/A)+1)+sqrt((3*((B-C)/A)+((B-A)/C)*((B-C)/A)+1)^2-4*(4*((B-A)/C)*((B-C)/A)))))
 
%Frecuencia de oscilacion estable de alabeo
 
wf_alabeo=n*sqrt(0.5*((3*((B-C)/A)+((B-A)/C)*((B-C)/A)+1)-sqrt((3*((B-C)/A)+((B-A)/C)*((B-C)/A)+1)^2-4*(4*((B-A)/C)*((B-C)/A)))))
end
end


%Lista de funciones auxiliares
    
function [dw]=ecseulergrav(t,w,matrizposicion)
global A;
global B;
global C;
global mu;
dw=zeros(10,1); 

E=[w(4)^2-w(5)^2-w(6)^2+w(7)^2 2*(w(4)*w(5)+w(6)*w(7)) 2*(w(4)*w(6)-w(5)*w(7))
2*(w(4)*w(5)-w(6)*w(7)) -w(4)^2+w(5)^2-w(6)^2+w(7)^2 2*(w(5)*w(6)+w(4)*w(7))
2*(w(4)*w(6)+w(5)*w(7)) 2*(w(5)*w(6)-w(4)*w(7)) -w(4)^2-w(5)^2+w(6)^2+w(7)^2];

inr=[A 0 0; 0 B 0; 0 0 C];
r=E*matrizposicion';

% Ecuaciones de la dinamica del movimiento

%Ecuaciones de las velocidades angulares
dw(1) = (-(w(2)*w(3))*(C-B)+(C-B)*(3*mu*r(2,1)*r(3,1)/(norm(r)^5)))/A;
dw(2) = (-(w(1)*w(3))*(A-C)+(A-C)*(3*mu*r(1,1)*r(3,1)/(norm(r)^5)))/B;
dw(3) = (-(w(1)*w(2))*(B-A)+(B-A)*(3*mu*r(2,1)*r(1,1)/(norm(r)^5)))/C; 

%Evolucion de los angulos de euler
dw(4) = 0.5*[w(3)*w(5)-w(2)*w(6)+w(1)*w(7)]; 
dw(5) = 0.5*[-w(3)*w(4)+w(1)*w(6)+w(2)*w(7)];
dw(6) = 0.5*[w(2)*w(4)-w(1)*w(5)+w(3)*w(7)];
dw(7) = 0.5*[-w(1)*w(4)-w(2)*w(5)-w(3)*w(6)];

end
 
function s = stumpS(z)

%Evalua las funciones de Stumpff S(z)

    if z > 0
        s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
        s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else
        s = 1/6;
    end

end

function c = stumpC(z)

%Evalua la funcion de Stumff C(z)

    if z > 0
        c = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        c = (cosh(sqrt(-z)) - 1)/(-z);
    else
        c = 1/2;
    end
end

function x = kepler_U(dt, ro, vro, a)

%{
Esta funcion resuleve la ecuacion de Kepler universal para anomalia 
universal a partir del metodo iterativo de Newton-Raphson

  mu   - parametro gravitacional (km^3/s^2)
  x    - anomalia universal (km^0.5)
  dt   - tiempo desde x = 0 (s)
  ro   - posicion radial (km) cuando x = 0
  vro  - velocidad radial (km/s) cuando x = 0
  a    - inverso del semieje mayor (1/km)
  z    - variable auxiliar (z = a*x^2)
  C    - valor de la funcion de Stumpff 
  S    - valor de la funcion de Stumpff 
  n    - numero de iteraciones para la convergencia de la solucion
  nMax - maximo numero de iteraciones permitidas
 
%}

    global mu

    %Establecer la tolerancia y un numero maximo de iteraciones:
    error = 1.e-8;
    nMax = 1000;

    %Valor inicial para x (Chobotov):
    x = sqrt(mu)*abs(a)*dt;
    n = 0;
    ratio = 1;

    %Iteracion hasta convergencia 
    %Error de la tolerancia:

    while abs(ratio) > error && n <= nMax
        n = n + 1;
        C = stumpC(a*x^2);
        S = stumpS(a*x^2);
        F = ro*vro/sqrt(mu)*x^2*C + (1 - a*ro)*x^3*S + ro*x - sqrt(mu)*dt;
        dFdx = ro*vro/sqrt(mu)*x*(1 - a*x^2*S) + (1 - a*ro)*x^2*C + ro;
        ratio = F/dFdx;
        x = x - ratio;
    end

    if n > nMax
        fprintf('\n **No. iterations of Kepler''s equation = %g', n)
        fprintf('\n F/dFdx = %g\n', F/dFdx)
    end
end


function [f, g] = f_y_g(x, t, ro, a)

%{
 
Esta funcion calcula los coeficientes de Lagrange.

  mu - parametro gravitacional (km^3/s^2)
  a  - inverso del semieje mayor (1/km)
  ro - posicion radial para el instante to (km)
  t  - tiempo pasado desde ro (s)
  x  - anomalia universal despues del tiempo t (km^0.5)
  f  - coeficiente de Lagrange f (dimensionless)
  g  - coeficiente de Lagrange g (s)
  
%}

global mu
    z = a*x^2;
    f = 1 - x^2/ro*stumpC(z);
    g = t - 1/sqrt(mu)*x^3*stumpS(z);
end


function [R,V] = vectorestadotiempo(R0, V0, t)

%{
 Esta función calculará el vector de estado [R,V] en el marco de 
 referencia ecuatorial a partir de un vector de estado inicial [r0,v0]
 en el marco de referencia ecuatorial y un intervalo de tiempo 
 transcurrido.
 mu - parametro gravitacional (kmˆ3/sˆ2)
 R0 - vector de posicion inicial (km)
 V0 - vector de velocidad inicial (km/s)
 t - intervalo de tiempo transcurrido (s)
 R - Vector de posicion final (km)
 V - Vector de velocidad final (km/s)

%}


global mu

    %...Magnitudes de R0 y V0:
    r0 = norm(R0);
    v0 = norm(V0);

    %Calculo de la velocidad inicial radial:
    vr0 = dot(R0, V0)/r0;

    %Calculo del inverso del semieje mayor
    alpha = 2/r0 - v0^2/mu;

    %Calculo de la anomalia universal
    x = kepler_U(t, r0, vr0, alpha);

    %Calculo de los coeficientes de Lagrange f and g y 
    %sus derivadas y el vector de posicion final R:
    [f, g] = f_y_g(x, t, r0, alpha);
    R = f*R0 + g*V0;
    r = norm(R);
    [fdot, gdot] = derivadas_f_g(x, r, r0, alpha);

    %Calculo de la velocidad final V :
    V = fdot*R0 + gdot*V0;
end

function coe = coe_vectorestado(R,V,mu)
%{
Esta función calculara el vector (coe) de elementos orbitales a partir
del vector de estado (r,v) 
  mu   - parametro gravitatorio (km^3/s^2)
  R    - vector de posicion en el marco geocentrico ecuatorial (km)
  V    - vector velocidad en el marco geocentrico ecuatorial (km)
  r, v - magnitud de r y v
  vr   - componente radial de la velocidad (km/s)
  H    - vector momento angular (km^2/s)
  h    - magnitud de H (km^2/s)
  i    - inclinacion de la orbita (rad)
  N    - vector linea de nodos (km^2/s)
  n    - magnitud de N
  cp   - producto cruzado de N y R
  RA   - asizension recta del nodo asizendente (rad)
  E    - vector excentricidad
  e    - excentricidad (magnitud de E)
  eps  - un pequeño numero por debajo del cual la excentricidad es
  considerada cero
       
  w    - argumento del perigeo (rad)
  TA   - anomalia verdadera (rad)
  a    - semieje mayor (km)
  pi   - 3.1415926...
  coe  - vector de elementos orbitales [h e RA i w TA a]

%}

    eps = 1.e-10;
    r = norm(R);
    v = norm(V);
    vr = dot(R,V)/r;
    H = cross(R,V);
    h = norm(H);
    incl = acos(H(3)/h);

    N = cross([0 0 1],H);
    n = norm(N);

    if n ~= 0
        RA = acos(N(1)/n);
        if N(2) < 0
            RA = 2*pi - RA;
        end
        else
        RA = 0;
    end
    
    E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
    e = norm(E);

    if n ~= 0
        if e > eps
            w = acos(dot(N,E)/n/e);
                if E(3) < 0
                w = 2*pi - w;
                end
        else
        w = 0;
        end
    else
    w = 0;
    end
    
    if e > eps
        TA = acos(dot(E,R)/e/r);
    if vr < 0
        TA = 2*pi - TA;
    end
    else
        cp = cross(N,R);
        if cp(3) >= 0
    TA = acos(dot(N,R)/n/r);
        else
    TA = 2*pi - acos(dot(N,R)/n/r);
        end
    end
    
a = h^2/mu/(1 - e^2);
coe = [h e RA incl w TA a];

end 



function [yaw pitch roll] = dcm_to_ypr(Q)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
Esta función obtiene la guiñada, el cabeceo y el alabeo a partir de la matriz de cosenos 
directores. 
Q - matriz de cosenos directores
yaw - guiñada(deg)
pitch - cabeceo (deg)
roll - alabeo (deg)
%}

yaw = atan2d_0_360(Q(1,2), Q(1,1));
pitch = asind(-Q(1,3));
roll = atan2d_0_360(Q(2,3), Q(3,3));
end

function t = atan2d_0_360(y,x)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
Esta función calcula el arco tangente de y/x en grados 
para el rango [0,360].
t - ángulo en grados.
%}
if x == 0
if y == 0
t = 0;
elseif y > 0
t = 90;
else
t = 270;
end
elseif x > 0
if y >= 0
t = atand(y/x);
else
t = atand(y/x) + 360;
end
elseif x < 0
if y == 0
t = 180;
else
t = atand(y/x) + 180;
end
end
end

function [r, v] = vectorinercial(coe,mu)

%{
  Esta funcion calcula el vector de estado en el marco geocentrico 
  (r,v) a partir del vector de los elementos orbitales clasicos(coe).
 
  mu   - parammetro gravitatorio terrestre (km^3;s^2)
  coe  - elementos orbitales [h e RA i w TA]
         
             h    = momento angular (km^2/s)
             e    = excentricidad
             RA   = asizension recta del nodo asizendente (rad)
             i    = inclinacion de la orbita (rad)
             w    = argumento del perigeo (rad)
             TA   = anomalia verdadera (rad)
  
  R3_w - matriz de rotacion sobre el eje-z a traves del angulo w
  R1_i - matriz de rotacion sobre el eje-x a traves del angulo i
  R3_W - matriz de rotacion sobre el eje-z a traves del angulo RA
  Q_pX - matriz de transformacion del marco de referencia perifocal al
  marco de referencia ecuatorial
  
  rp   - vector de posicion en el marco perifocal (km)
  vp   - vector de velocidad en el marco perifocal (km/s)
  r    - vector de posicion en el marco geocentrico ecuatorial (km)
  v    - vector velocidad en el marco geocentrico ecuatorial (km/s)
%}

    h = coe(1);
    e = coe(2);
    RA = coe(3);
    incl = coe(4);
    w = coe(5);
    TA = coe(6);
    rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
    vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

    R3_W = [ cos(RA) sin(RA) 0
    -sin(RA) cos(RA) 0
    0 0 1];

    R1_i = [1 0 0
    0 cos(incl) sin(incl)
    0 -sin(incl) cos(incl)];

    R3_w = [ cos(w) sin(w) 0
    -sin(w) cos(w) 0
    0 0 1];

    Q_pX = (R3_w*R1_i*R3_W)';

    r = Q_pX*rp;
    v = Q_pX*vp;

    r = r';
    v = v';
end

function Q = dcm_q(q)
% ~~~~~~~~~~~~~~~~~~~~~~~
%{ 
Esta función calcula la matriz de cosenos directores desde el cuaternion
q - cuaternion (q(4) es la parte escalar)
Q - matriz de cosenos directores
%}

q1 = q(1); q2 = q(2); q3 = q(3); q0 = q(4);
Q = [q1^2-q2^2-q3^2+q0^2, 2*(q1*q2+q3*q0), 2*(q1*q3-q2*q0);
2*(q1*q2-q3*q0), -q1^2+q2^2-q3^2+q0^2, 2*(q2*q3+q1*q0);
2*(q1*q3+q2*q0), 2*(q2*q3-q1*q0), -q1^2-q2^2+q3^2+q0^2 ];
end 


function q = q_dcm(Q)

%{
Esta función calcula el cuaternión a partir de la matriz de 
cosenos directores.
Q - matriz de cosenos directores.
q - cuaternion (donde q(4) es la parte escalar)
%}

K3 = ...
[Q(1,1)-Q(2,2)-Q(3,3), Q(2,1)+Q(1,2), Q(3,1)+Q(1,3), Q(2,3)-Q(3,2);
Q(2,1)+Q(1,2), Q(2,2)-Q(1,1)-Q(3,3), Q(3,2)+Q(2,3), Q(3,1)-Q(1,3);
Q(3,1)+Q(1,3), Q(3,2)+Q(2,3), Q(3,3)-Q(1,1)-Q(2,2), Q(1,2)-Q(2,1);
Q(2,3)-Q(3,2), Q(3,1)-Q(1,3), Q(1,2)-Q(2,1), Q(1,1)+Q(2,2)+Q(3,3)]/3;
[eigvec, eigval] = eig(K3);
[x,i] = max(diag(eigval));
q = eigvec(:,i);
end


function [fdot, gdot] = derivadas_f_g(x, r, ro, a)

%{
Esta funcion calcula las derivadas de las funciones f y g de Lagrange en el
tiempo
%}
    global mu
    z = a*x^2;
    fdot = sqrt(mu)/r/ro*(z*stumpS(z) - 1)*x;
    gdot = 1 - x^2/r*stumpC(z);

end