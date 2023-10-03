function result = laguerre(n, m, x)
result = zeros(size(x));
    for k = 0:n
        result = result + power(-1, k) / factorial(n - k) / factorial (m + k) / factorial(k) * power(x, k);
    end
    result = result * factorial(n + m);
end

