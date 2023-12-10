% Стартовая инициализация переменных
m = 5;
R = 5;
n = 128;
h_r = R/n;
r = 0:h_r:(R-h_r/2);
N = 2*n-1;
M = 2048;


% Получаем значение функции circ
function result = circ(r)
  if r <= 1
    result = 1;
  else
    result = 0;
  end
endfunction


% Функция для построения одномерных графиков
function plot1d(range, fv, title_abs, title_arg, xlabel_s)
  figure;
  plot(range, abs(fv), "linewidth", 2);
  title(title_abs);
  xlabel(xlabel_s);
  ylabel("Амплитуда");

  figure;
  plot(range, angle(fv), "linewidth", 2);
  xlabel(xlabel_s);
  ylabel("Фаза");
  title(title_arg);
endfunction


% Функция для построения двумерных графиков
function plot2d(range1, range2, fv, title_abs, title_arg)
  figure();
  imagesc(range1, range2, abs(fv));
  title(title_abs);
  colormap gray;
  colorbar;

  figure();
  imagesc(range1, range2, arg(fv));
  title(title_arg);
  colormap gray;
  colorbar;
endfunction


% Восстанавливаем изображение функции до двумерного
function restored_func = restore_func(func, n, m, N)
  restored_func = zeros(N);

  for j = 1:N
    for k = 1:N
      alpha = round(sqrt((j-n)^2 + (k-n)^2)) + 1;
      if alpha <= n
        restored_func(j, k) = func(alpha) * exp(1i*m*atan2(k-n, j-n));
      endif
    endfor
  endfor
endfunction


% Численная реализация преобразования Ханкеля
function hankel_result = hankel(func, r, h_r, n, m)
    hankel_result = zeros(1, n);
    for i = 1:n
        hankel_result(i) = sum(func.*besselj(m, 2.*pi.*r.*r(i)).*r.*h_r);
    endfor
    hankel_result = hankel_result.* ((2.*pi)./(1i^m));
endfunction


% Одномерное преобразование Фурье
function F = fourier_bpf(fv, m, n, h_r)
  F = zeros(1, m);

  start_idx = floor((m - n) / 2) + 1;
  end_idx = start_idx + n - 1;
  F(start_idx:end_idx) = fv;
  F = fftshift(F);
  F = fft(F) * h_r;
  F = fftshift(F);

  F = F(length(F) / 2 - n / 2 + 1: length(F) / 2 + n / 2);
endfunction


% Двумерное преобразование Фурье
function F_bpf_2d = fourier_bpf_2d(F_bpf_2d, h_r, N, M)
    for i = 1:N
        F_bpf_2d(i, :) = fourier_bpf(F_bpf_2d(i, :), M, N, h_r);
    endfor
    for i = 1:N
        F_bpf_2d(:, i) = fourier_bpf(F_bpf_2d(:, i), M, N, h_r);
    endfor
endfunction



% Ззначение функции по варианту
result = arrayfun(@circ, r./2) - arrayfun(@circ, 0.95.*r./2);

% Графики входной функции
plot1d(r, result, "Распределение амплитуды входной функции", "Распределение фазы входной функции", "R");

% восстанавливаем изображение входной функции до двумерного
restored_func = restore_func(result, n, m, N);

% Двумерный график восстановленного изображения входной функции
plot2d([-R, R], [-R, R], restored_func, "Амплитуда входной функции", "Фаза входной функции")

% делаем преобразование Ханкеля и замеряем время
tic
hankel_result = hankel(result, r, h_r, n, m);
toc

% строим графики входной функции после преобразования Ханкеля
plot1d(r, hankel_result, "Амплитуда после преобразования Ханкеля", "Фаза после преобразования Ханкеля", "R");

% восстанавливаем изображение входной функции после преобразования Ханкеля
restored_hankel = restore_func(hankel_result, n, m, N);

% строим двумерный график восстановленного изображения входной функции после преобразова-ния Ханкеля
plot2d([-R, R], [-R, R], restored_hankel, "Амплитуда входной функции после преобразования Ханкеля", "Фаза  входной функции после преобразования Ханкеля");

% делаем преобразование Фурье и замеряем время
tic
fourier_result = fourier_bpf_2d(restored_func, h_r, N, M);
toc

% строим двумерный график восстановленного изображения входной функции после преобразова-ния Фурье
plot2d([-R, R], [-R, R], fourier_result, "Амплитуда двумерного БПФ для восстановленной функции", "Фаза двумерного БПФ для восстановленной функции");

