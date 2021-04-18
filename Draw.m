function Draw(polinom, t, n, left_board, right_board, x_max)
%построение графиков
figure;
ezplot('(cos(x))^2 + x - 1', [left_board right_board]) % рисует реальную функцию, синяя линия
hold on; % остаемся на том же графике
h = left_board:0.001:right_board; % задали малый шаг для построения
f_n = polyval(polinom, h); % для линии полинома
f_t = polyval(polinom, t); % для точек квазиальтернанса (выровненного альтернанса в данном случае)
plot(h, f_n, '-m', t, f_t, 'r*', t, F(t), 'b*'), grid; % рисует полином, сиренивые линии, точки квазиальтернанса - красные звездочки
% синие звездочки - значения функции в точках t
for i = 1:(n+2)
    plot([t(i), t(i)], [f_t(i), F(t(i))], 'Color', 'r', 'LineWidth', 1); % рисуем отклонения на сигму и минус сигму красными отрезками
end
if (x_max >= left_board)
    plot([x_max, x_max], [F(x_max), polyval(polinom, x_max)], 'Color', 'k', 'LineWidth', 1); % рисуем найденное максимальное отклонение
    plot(x_max, F(x_max), 'k*', x_max, polyval(polinom, x_max), 'k*');
axis tight; % подгоняет границы осей под крайние значения функций
legend('исходная функция', 'полином приближения по квазиальтернансу', 'Location', 'NorthWest');
refreshdata
    drawnow
%конец построения графиков
end