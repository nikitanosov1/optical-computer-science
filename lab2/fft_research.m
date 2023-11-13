clear all
close all

# Одномерное FT

function half_swaped = vector_half_swap(b)
  len_b = length(b); half_swaped = zeros(size(b));

  half_swaped(1 : len_b/2) = b(len_b/2 + 1 : end);
  half_swaped(len_b/2 + 1 : end) = b(1 : len_b/2);
endfunction

function F = FFT(f, N, M)
  zeros_v = zeros(1, M/2 - N/2);
  _f = vector_half_swap([zeros_v, f, zeros_v]);
  F = vector_half_swap(fft(_f))(M/2 - N/2 + 1: M/2 + N/2);
endfunction

N = pow2(9); # pow2(9);
M = pow2(13); # pow2(13);
a = 5;
b = N * N / (4 * a * M);

h_x = (2 * a) / N;
x = -a:h_x:(a - h_x/2);
h_u = 2 * b / N;
u = -b:h_u:(b - h_u/2);
gaussian_beam = exp(-x .* x);

# Строим Гауссовский пучок
figure(1);
subplot(2, 1, 1); plot(x, abs(gaussian_beam)); title("Амплитуда пучка");
subplot(2, 1, 2); plot(x, arg(gaussian_beam)); title("Фаза пучка");

# Получаем FFT от Гауссовского пучка
gaussian_fft = FFT(gaussian_beam, N, M) * h_x;

# Убеждаемся, что gaussian_fft == gaussian_beam, так как e^(-x^2) является
# собственной функцией преобразования Фурье
figure(2);
subplot(2, 1, 1); plot(u, abs(gaussian_fft)); title("Апмлитуда преобразования пучка");
subplot(2, 1, 2); plot(u, arg(gaussian_fft)); title("Фаза преобразования пучка");

# Получаем DFT от Гауссова пучка
[X, U] = meshgrid(x, u);
gaussian_dft = exp(-2 * pi * 1i * X .* U) * gaussian_beam.' * h_x;

figure(3);
subplot(2, 1, 1); plot(u, abs(gaussian_fft), "g", u, abs(gaussian_dft), "r");
title("Амплитуда преобразования пучка: FFT (зеленое), DFT (красное)");
subplot(2, 1, 2); plot(u, arg(gaussian_fft), "g", u, arg(gaussian_dft), "r");
title("Фаза преобразования пучка: FFT (зеленое), DFT (красное)");

f = sinc(x)
F = FFT(f, N, M) * h_x;

figure(4);
subplot(2, 2, 1); plot(x, abs(f)); title("Амплитуда поля");
subplot(2, 2, 2); plot(x, arg(f)); title("Фаза поля");
subplot(2, 2, 3); plot(u, abs(F)); title("Амплитуда преобразования");
subplot(2, 2, 4); plot(u, arg(F)); title("Фаза преобразования");

# F_analytical = ((1 + 16 * pi * pi * u .* u) / 4) .^ (-1);
function y = rect(x)
    y = (abs(x) < 0.5) + (abs(x) == 0.5) * 0.5;
end
# F_analytical = rect(x * pi) / pi
F_analytical = rect(u)

figure(5);
subplot(2, 1, 1); plot(u, abs(F), "g", u, abs(F_analytical), "r");
title("Амплитуда поля: FFT (зеленое), analytical(красное)");
subplot(2, 1, 2); plot(u, arg(F), "g", u, arg(F_analytical), "r");
title("Фаза поля: FFT (зеленое), analytical(красное)");


# Двумерное FT

function F = FFT_2d(f, N, M, h_u, h_v)
  F = zeros(N, N);

  for i = 1:N
    F(i, :) = FFT(f(i, :), N, M) * h_u;
  endfor

  F = F.';
  for j = 1:N
    F(j, :) = FFT(F(j,:), N, M) * h_v;
  endfor
  F = F.';
endfunction

N_2d = pow2(9);
M_2d = pow2(13);
a_2d = 5; b_2d = N_2d * N_2d / (4 * a_2d * M_2d);
h_x_2d = 2 * a_2d / N_2d; x_2d = -a_2d:h_x_2d:(a_2d - h_x_2d / 2);
y_x_2d = 2 * a_2d / N_2d; y_2d = -a_2d:y_x_2d:(a_2d - y_x_2d / 2);
[X_2d, Y_2d] = meshgrid(x_2d, y_2d);
gaussian_beam_2d = exp(-X_2d .* X_2d - Y_2d .* Y_2d);

figure(6);
subplot(2, 1, 1); imagesc(x_2d, y_2d, abs(gaussian_beam_2d)); title("Амплитуда двумерного пучка");
subplot(2, 1, 2); imagesc(x_2d, y_2d, arg(gaussian_beam_2d)); title("Фаза двумерного пучка");

h_u_2d = 2 * b_2d / N_2d; u_2d = -b_2d:h_u_2d:(b_2d - h_u_2d / 2);
h_v_2d = 2 * b_2d / N_2d; v_2d = -b_2d:h_v_2d:(b_2d - h_v_2d / 2);
gaussian_fft_2d = FFT_2d(gaussian_beam_2d, N_2d, M_2d, h_u_2d, h_v_2d);

figure(7);
subplot(2, 1, 1); imagesc(u_2d, v_2d, abs(gaussian_fft_2d));
title("Амплитуда преобразования двумерного пучка");
subplot(2, 1, 2); imagesc(u_2d, v_2d, arg(gaussian_fft_2d));
title("Фаза преобразования двумерного пучка");

f_2d = sinc(X_2d) * sinc(Y_2d);
F_2d = FFT_2d(f_2d, N_2d, M_2d, h_u_2d, h_v_2d);

figure(8);
subplot(2, 1, 1); imagesc(u_2d, v_2d, abs(F_2d));
title("Амплитуда преобразования двумерного поля");
subplot(2, 1, 2); imagesc(u_2d, v_2d, arg(F_2d));
title("Фаза преобразования двумерного поля");


