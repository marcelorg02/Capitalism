# Capitalism
Codes for the simullation of a capitalist system.
clear 
clc 
clf

N = 1000; % Total de personas en la simulacion
M = 100000; % Total de dinero disponible
info = Generar(N,M); % Almacenamiento de la información del sistema, primer renglón indica dinero, segundo empleador , tercero clase, cuarto ingreso anual
%Quinto indica clase, sexto ingreso por propiedad, séptimo ingresos por
%salario
Tiempo = 100; %Tiempo en años de la simulacion
mercado  =8; % Valor de lo que hay en el mercado
salarios = [10 , 90];
frec_clases = [];  %Vector para trackear el tamaño de cada clase a lo largo del timepo
firmtrack= [];
S = [];  %Vector de entropía
bandera =0;
alf = (N/M);
mesesp = zeros(12,N);
mesesw = zeros(12,N);
anosp = zeros(5,N);
anosw = zeros(5,N);
control = 0;
mesest = zeros(12,N);
anost = zeros(5,N);
%Definir muestras 



for i  = 1:12 *Tiempo
        if mod(i,12) ==0 && i ~= Tiempo*12
            info(4,:) = 0;
            info(6,:) = 0;
            info(7,:) =0;
        end
    for j  = 1:N
        a = randi([1 , length(info)]);  %Regla de selección
        info = H1(info,a, salarios);    %Determinar quien lo contrata
        m = rand() * info(1,a);    %Determinar gasto
        info(1,a) = info(1,a) - m; %Quitar dinero del gasto
        mercado = mercado + m; %Agregar gasto al mercado
        [info, mercado] = M1(info,a,mercado);  %Recolectar Ganancias del mercado
        info = F1(info,a,salarios);  %Regla de despidos
        if info(3,a)== 1
            info = W1(info,a,salarios);   %Paga de salarios
        end
         
        %Cáculos necesarios de entropía y bienestar.

        if mod(j,1000) == 1
            bandera = bandera + 1;
            class = frecc(info);
            s = log(class(1,:) + 1);
            if i < 8*12;
            S(j + i -1) = N*log(N) - sum(class(1,:).*s);
            T(j+i -1) = j+i-1;
            end
            t(bandera) = bandera;
            o1 = 1-exp(-alf.*class(2,:));
            O(bandera) = (sum(class(1,:).*o1))/N;
            frec_clases = [frec_clases, EleClase(info).'];  %Linea para calcular frecuencia de los sets en cada año
        end
    end
        if i < 12
            mesesp(i,:) = info(6,:);
            mesesw(i,:) = info(7,:);
            mesest(i,:) = info(4,:);
        end
        if mod(i,20*12) == 0
            control = control+1;
            anosp(control,:) = info(6,:);
            anosw(control,:) = info(7,:);
            anost(control,:) = info(4,:);

        end
end

    


[SC,IPC] = SetClases(info,1);
[Sw,IPW] = SetClases(info,2);
% [Lc,IPLC] = [SetClases(info,2), SetClases(info,3)];
firmas = FirmSize(info);
[firma,masgrande] = max(firmas(2,:));
Mw = DineroClases(info,2);
MC  =DineroClases(info,1);
Employers = VectorEmpleadores(info);


fig1 = figure(1);
tiledlayout(2,2)
nexttile
histogram(frec_clases(1,:),500);
xlim([0, N]);
title("Capitalistas");
nexttile
histogram(frec_clases(2,:),500);
title("Empleados");
xlim([0, N]);
nexttile
histogram(frec_clases(3,:),500);
title("Desempleados");
nexttile
% bar(info(5,:));
% histogram(info(5,:));

fig2 = figure(2);
tiledlayout(2,2)
nexttile
[f,x] = ecdf(info(4,:));
loglog(x,1-f);
xlabel("Ingresos");
ylabel("Pcc")
title("ccdf Ingresos agregado")
nexttile
[f,x] = ecdf(info(6,:));
loglog(x,1-f);
title("ccdf Ingresos desagregados")
xlabel("Ingresos");
ylabel("Pcc")
hold on
[f,x] = ecdf(info(7,:));
loglog(x,1-f);
title("ccdf Ingresos desagregados")
hold off
legend("Ingresos capitalistas","Ingresos salarios","Location","southwest")
nexttile
hold off
histogram(Sw,500);
set(gca, 'xscale','log');
title("Histograma ingresos por sueldos")
xlabel("Ingresos");
ylabel("frecuencia")
nexttile
[f,x] = ecdf(info(6,:));
loglog(x,1-f)
title("Ingresos regimen de propiedad")
xlabel("ingresos")
ylabel("Pcc")



fig3 = figure(3);
tiledlayout(2,2)
nexttile
[f,x] = ecdf(info(1,:));
loglog(x,1-f);
title("distribucion de dinero global")
xlim([0 10^4])
xlabel("Dinero");
ylabel("Pcc")
nexttile 
ylabel("Pcc")
[f,x] = ecdf(Mw);
loglog(x,1-f);
hold on
[f,x] = ecdf(MC);
loglog(x,1-f);
title("distribucion de dinero desagregada")
xlim([0 10^4])
xlabel("Dinero")
ylabel("Pcc")
nexttile 
semilogy(x,1-f);
title("Seccion de distribucion trabajadores")
xlim([0 500])
xlabel("Dinero")
ylabel("Pcc")
nexttile 
[f,x] = ecdf(MC);
loglog(x,1-f);
title("seccion de distribucion capitalistas")
xlabel("Dinero")
ylabel("Pcc")
xlim([200 20000]);
% vec = VectorEmpleadores(info)
% ve = Empleados(info,3)
% frec_clases = EleClase(info);
fig4 = figure(4);
plot(T,S);
xlabel("Meses");
ylabel("Entropía");
title("Evolución de la entropía");

fig5 = figure(5);
t = t./12;
plot(t,O);
xlabel("Años");
ylabel("Bienestar");
title("Evolución del bienestar");


fig6 = figure(6);
[f,x] = ecdf(mesesp(1,:));
loglog(x,1-f);
xlim([0 10^4])
hold on
[f,x] = ecdf(mesesp(3,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(mesesp(6,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(anosp(2,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(anosp(4,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(anosp(5,:));
loglog(x,1-f);
xlim([1 10^4])
title("Ingresos por propiedad");
legend("mes1","mes3","mes6","año 40","año 80","año 100","Location","southwest")
xlabel("Dinero");
ylabel("Pcc");
hold off



hold off
fig7 = figure(7);
[f,x] = ecdf(mesesw(1,:));
loglog(x,1-f);
xlim([0 10^4])
hold on
[f,x] = ecdf(mesesw(3,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(mesesw(6,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(anosw(2,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(anosw(4,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(anosw(5,:));
loglog(x,1-f);
xlim([0 10^4])
title("Ingresos por salarios");
legend("mes1","mes3","mes6","año 40","año 80","año 100","Location","southwest");
xlabel("Dinero");
ylabel("Pcc");
xlim([1 10^4]);
hold off



hold off
fig8 = figure(8);
[f,x] = ecdf(mesest(1,:));
loglog(x,1-f);
xlim([0 10^4])
hold on
[f,x] = ecdf(mesest(3,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(mesest(6,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(anost(2,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(anost(4,:));
loglog(x,1-f);
xlim([0 10^4])
[f,x] = ecdf(anost(5,:));
loglog(x,1-f);
xlim([0 10^4])
title("Ingresos agregados");
legend("mes1","mes3","mes6","año 40","año 80","año 100","Location","southwest")
xlabel("Ingresos")
ylabel("Pcc")
xlim([1 10^4]);
hold off




function inf = Generar(N,M)
    inf = zeros(7,N);
    for i = 1:N
        inf(1,i) = M/N;
        inf(2,i) = 0;
        inf(3,i) = 3;
        inf(4,i) = 0;
        inf(5,i) = Discriminar(inf(1,i));
        inf(6,i) = 0;
        inf(7,i) = 0;
    end
end


function val = truncar(a)
    val = a - mod(a,1);
end

function H = VectorEmpleadores(info)   %Funcion que regresa el conjunto de empleafores, 1.-indice, 2.-dinero, 3.-probabilidad
bandera = 0;

    for i = 1:length(info)
        if info(3,i) == 3 || info(3,i) == 1
            bandera = bandera + 1;
        end
    end
    H= zeros(3,bandera);
    bandera = 0;
    for i = 1:length(info)
        if info(3,i) == 3 || info(3,i) == 1
            bandera = bandera + 1;
            H(1,bandera) = i;
            H(2,bandera) = info(1,i);
            H(3,bandera) = 0;
        end
    end
    H(3,:) = H(2,:)./(sum(H(2,:)));
end

function Vi = VecInt(info) %Funcion de escalamiento de intervalos para transformacion de distribucion uniforme
    H = VectorEmpleadores(info);
    Vi = zeros(1,length(H));
    Vi(1) = H(3,1);
    for i = 2:length(H)
        Vi(i) = Vi(i-1) + H(3,i);
    end
end


function c = Empleador(info) %Funcion que selecciona un empleador aleatorio con distribucion indicada
    n = rand();
    Vi = VecInt(info);
    H = VectorEmpleadores(info);
    for i =  1:length(Vi)
        if n < Vi(1)
            c = H(1,1);
        elseif n < Vi(i) && n > Vi(i-1)
            c = H(1,i);

        end
    end
end


function info = H1(info,a,salarios) %Determinar si se contrata y a quien contrata
        if info(3,a) == 3
            c = Empleador(info);
            if info(1,c) > sum(salarios)/2
             info(2,a) = c;  %Indicar para quien va a trabajar
             info(3,a) = 2;  %Indicar que ahora es trabajador
             info(3,c) = 1;  %Cambiar la clase del empleador
            end
        end
end


function [info,mercado] = M1(info,a,mercado) 
        if info(3,a) ~= 3  %Revisar si el elemento seleciconado no es desempleado
            m = rand() * mercado; %Selecionar una cantidad del mercado a cobrar
            mercado = mercado - m; 
            if info(3,a) == 2
                info(1,info(2,a)) = info(1,info(2,a)) + m; %Transferir dinero al patron
                info(4,info(2,a)) = info(4,info(2,a)) +m; %Sumar a su ingreso anual total
                info(5,info(2,a)) = Discriminar(info(1,info(2,a)));
                info(6,info(2,a)) = info(6,info(2,a)) + m;
            else
                info(1,a)  = info(1,a) + m; %Si el seleccionado es el patron
                info(4,a) = info(4,a) + m;
                info(5,a) = Discriminar(info(1,a));
                info(6,a) = info(6,a) + m;
            end

       end
end

function Wa  = Empleados(info, c)  %Regresa el set de empleados de un empleador "c", cada elemento representa el indice de cada trabajador
Wa  = [];
bandera = 0;
for i = 1:length(info)
    if info(2,i) == c
        bandera  = bandera +1;
        Wa(bandera)  = i;
    end
end
end

function info  = F1(info, a,salarios)  %Funcion para despido de trabajadores
    We  =   Empleados(info,a);
    if info(3,a) == 1
    for i = 1:length(We)
            Wa  =   Empleados(info,a);
        if info(1,a) < length(Wa) * (sum(salarios)/2)
            c = randi([1, length(Wa)]);
            info(2,Wa(c)) = 0;
            info(3,Wa(c)) = 3;
        end
    if isempty(Empleados(info,a))
        info(3,a) = 3;
    end
    end
    end
end

function info = W1(info,a,salarios)   %Función para paga de salarios
       SE = Empleados(info,a);  %SE(i) Representa el indice del empleado en el vector info
       for i = 1:length(SE)
            w = salarios(1) + (rand()*(salarios(2)- salarios(1))); 
            if info(1,a) > w
                info(1,a)  = info(1,a) - w;
                info(5,a) = Discriminar(info(1,a));
                info(1,SE(i)) = info(1,SE(i)) + w; 
                info(4,SE(i)) = info(4,SE(i)) + w;
                info(5,SE(i)) = Discriminar(info(1,SE(i)));
                info(7,SE(i)) = info(7,SE(i)) + w;
            else
                m = rand() * info(1,a);
                info(1,a) = info(1,a) - m;
                info(5,a) = Discriminar(info(1,a));
                info(1,SE(i)) = info(1,SE(i)) + m;
                info(4,SE(i)) = info(4,SE(i)) + m;
                info(5,SE(i)) = Discriminar(info(1,SE(i)));
                info(7,SE(i)) = info(7,SE(i)) + m;
            end
       end
end

function cl = EleClase(info)%Numero de elementos por clase
    cl = [0,0,0];
    for i= 1:length(info)
        for j = 1:3
        if info(3,i) == j
            cl(j) = cl(j) + 1;
        end
        end
    end   
end



function d = Discriminar(m,~) %Función que determina a que clase pertenece cierto elemento.
d = truncar((m)/10) + 1;
if d < 0
    d = 0;
end
end


function [vc,ip] = SetClases(info,p) %Reegresa el ingreso anual de cada elemento de cada clase
    vc =[];
    bandera = 0;
    for i = 1:length(info)
        if info(3,i) == p
            bandera = bandera +1;
            vc(bandera) = info(4,i);
            ip(bandera) = info(6,i);
        end
    end
end

function vc = DineroClases(info,p) %Reegresa el ingreso anual de cada elemento de cada clase
    vc =[];
    bandera = 0;
    for i = 1:length(info)
        if info(3,i) == p
            bandera = bandera +1;
            vc(bandera) = info(1,i);
        end
    end
end

function tf = FirmSize(info)
    VE = VectorEmpleadores(info);
    Ve = [];
    tf = [];
    bandera = 0;
    for i = 1:length(VE)
        Ve = Empleados(info,VE(1,i));
        if isempty(Ve)
        else
        bandera = bandera + 1;
        tf(1,bandera) = VE(1,bandera);
        tf(2,bandera) = length(Ve);
        end
end
end

function fc = frecc(personas) %Función para determinar cuantos miembros pertenecen a cada clase
fc = zeros(2,400);
for i = 1:400
n=0;
p = 0;
for j = 1:length(personas)
    if personas(5,j) == i
        n = n+1;
        p = p + personas(1,j);
    end
    fc(1,i) = n;
    if n ~= 0 
        fc(2,i) = p/n;
    else 
        fc(2,i) = 0;
    end
end
end
end






