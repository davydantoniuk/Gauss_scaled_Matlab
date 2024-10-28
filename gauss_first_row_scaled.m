function [x] = gauss_first_row_scaled(A,b)
[m, n] = size(A);
if m == n
    p = 1:n;
    s = max(abs(A'));
    for k = 1:n-1
        [~, d] = max(abs(A(p(k:n),k))./((s(p(k:n)))'));
        d = d+k-1;
        temp = p(k);
        p(k) = p(d);
        p(d) = temp;
        for i = k+1:n
            mnoz = A(p(i), k)/A(p(k), k);
            A(p(i),:) = A(p(i),:) - A(p(k),:) * mnoz;
            b(p(i)) = b(p(i)) - b(p(k)) * mnoz;
        end
    end
    for i = m:-1:1
        suma = 0;
        for j = i + 1:n
            suma = suma + A(p(i), j) * x(j);
        end
        x(i) = (b(p(i)) - suma)/A(p(i), i);
    end
    x = x';
end
end

