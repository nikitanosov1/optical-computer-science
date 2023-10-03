% Исходные данные
n = 1000;
c = 3;
a = -c;
b = c;
beta = 1/10;
m = 1000;
p = -3;
q = 3;
alpha = 1;

% Задаем входной сигнал f(x)=e^(i*betta*x)
h = (b-a)/n;
x = a:h:(b-h/2);
f = exp(1i*beta*x);

% Рисуем входной сигнал
figure;
subplot(2, 1, 1);
plot(x, abs(f));
title('Амплитуда входного сигнала f(x)=e^{i*betta*x}');
subplot(2, 1, 2);
plot(x, angle(f));
title('Фаза входного сигнала f(x)=e^{i*betta*x}');

% Само преобразование
hxi = (q-p)/m;
xi = p:hxi:(q-hxi/2);
[X, XI] = meshgrid(x, xi);
Kernel = exp(-alpha*X.*XI) * laguerre(3, 0, X.*XI);
F = Kernel * f.' * h;

% Построение графика результата преобразования
figure;
subplot(2, 1, 1);
plot(x, abs(F));
title('Амплитуда результата преобразования');
subplot(2, 1, 2);
plot(x, angle(F));
title('Фаза результата преобразования');
