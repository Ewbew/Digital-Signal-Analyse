function y = moving_average_filter(x, M)

y = zeros(size(x)); % Preallocate y for efficiency

for n = 1:length(x)
    if n <= M
        y(n) = sum(x(1:n));
    else
        y(n) = sum(x(n-M+1:n));
    end
end
y=y./M;
end