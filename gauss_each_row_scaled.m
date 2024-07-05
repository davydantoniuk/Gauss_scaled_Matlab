% Napisz w MATLABIE funkcję rozwiązującą układ równań liniowych
% metodą eliminacji Gaussa ze skalowaniem wierszy w każdym kroku
% iteracyjnym. Zastosuj tę metodę do wybranych przykładów i porównaj
% z klasyczną metodą eliminacji Gaussa.

function [x] = gauss_each_row_scaled(A, b)
[m, n] = size(A);
if m == n
    p = 1:n;
    % Obliczenie początkowych współczynników skalowania dla każdego wiersza
    s = max(abs(A'));
    for k = 1:n-1
        % Ustawienie współczynników skalowania dla przetworzonych wierszy na 1
        s(p(1:k-1)) = 1;
        % Wykonanie częściowego wyboru z skalowaniem
        [~, d] = max(abs(A(p(k:n),k))./((s(p(k:n)))'));
        d = d + k - 1;
        % Zamiana indeksów obrotu
        temp = p(k);
        p(k) = p(d);
        p(d) = temp;
        % Eliminacja Gaussa
        for i = k+1:n
            mnoz = A(p(i), k)/A(p(k), k);
            A(p(i),:) = A(p(i),:) - A(p(k),:) * mnoz;
            b(p(i)) = b(p(i)) - b(p(k)) * mnoz;
        end
        % Aktualizacja współczynników skalowania dla pozostałych wierszy
        for i = k+1:n
            s(p(i)) = max(abs(A(p(i), k:n)));
        end
    end
    % Ustawienie współczynnika skalowania dla ostatniego przetworzonego wiersza na 1
    s(p(1:n-1)) = 1;
    for i = n:-1:1
        suma = 0;
        for j = i+1:n
            suma = suma + A(p(i), j) * x(j);
        end
        x(i) = (b(p(i)) - suma) / A(p(i), i);
    end
    x = x';
end
end

