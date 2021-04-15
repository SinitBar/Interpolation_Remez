function coefficients = Remez_poly(n, left_board, right_board) % n - степень полинома приближения, границы промежутка, на котором рассматривается функция
t = T(n+2, left_board, right_board); % нашли узлы интерполяции по Чебышеву, отсортированы от меньшего к большему 
number_of_step = 0; % номер шага по алгоритму Ремеза
epsilon = 0.00001; % эпсилон, если модуль разности сигм на последовательных шагах алгоритма меньше него, алгоритм Ремеза завершает работу
f = F(t); % нашли значения в выбранных узлах, функцию оставили прежней

while true % тут начинается выполнение алгоритма Ремеза
    matrix = Matrix(n, t); % вспомогательная матрица для нахождения коэффициентов полинома и сигмы
    %coefficients - первый элемент - сигма, дальше матрица коэффициентов полинома p(n)...p(0)
    coefficients = f/matrix; % искомая строка = строка значений функции умножить на обратную к matrix, тут написано то же самое, что и (f*inv(matrix))
    number_of_step = number_of_step + 1 % выводим номер шага
    if (number_of_step == 1)
        sigma = coefficients(1) % выводим и запоминаем сигму с 1 шага
    else
        Sigma = coefficients(1) % выводим сигму с 2 и других шагов и запоминаем ее
    end
    polinom = coefficients;
    polinom(1) = []; % удаляем из массива 1 элемент - сигму, получаем коэффициенты полинома в нужном порядке
    
    %ищем максимум отклонения и икс, при котором он достигается
    x_max = left_board;
    max_deviation = 0;
    for step = left_board:0.001:right_board
        if (abs(F(step)  -  polyval(polinom, step)) > max_deviation) % если модуль разности f и полинома > макс отклонения
            max_deviation = abs(F(step)  -  polyval(polinom, step));
            x_max = step;
        end
    end
    x_max = x_max
    max_deviation = max_deviation
    % теперь надо найти между какими узлами лежит x_max и изменить массив
    % точек и массив значений
    left_index = 0; % unreal in matlab, индекс узла, лежащего левее икса или индекс последнего элемента, так поймаем крайние иксы
    for i = 1:(n+1)
        if ((x_max > t(i))&&(x_max < t(i+1)))
            left_index = i;
            break;
        end
    end
    was_t = t;
    if (left_index == 0) % тогда икс или левее всех узлов или правее всех
        if (x_max > t(n+2)) % икс правее самого правого узла
            if ((F(t(n+2)) - polyval(polinom, t(n+2)))*(F(x_max) - polyval(polinom, x_max)) > 0) % знаки совпадают
                t(n+2) = x_max; % заменяем последний узел
                f(n+2) = F(x_max); % и меняем значение для него
                %change_last_node = x_max
            else
                t(n+3) = x_max; % добавляем икс в конец последовательности узлов
                f(n+3) = F(x_max);
                t(1) = []; % удаляем первый узел
                f(1) = [];
            end
        else % икс левее самого левого узла
            if ((F(t(1))  -  polyval(polinom, t(1)))*(F(x_max)  -  polyval(polinom, x_max)) > 0) % знаки совпадают
                t(1) = x_max; % заменяем первый узел
                f(1) = F(x_max); % и меняем значение для него
            else
                t = [x_max t]; % добавляем икс в начало последовательности узлов
                f = [(F(x_max)) f];
                t(n+3) = []; % удаляем последний узел
                f(n+3) = [];
            end
        end
    else % икс находится между существующими узлами
        if ((F(t(left_index))  -  polyval(polinom, t(left_index)))*(F(x_max)  -  polyval(polinom, x_max)) > 0) % знаки совпадают
            t(left_index) = x_max;
            f(left_index) = F(x_max);
        else % тогда знаки совпадают с другим соседом, так как знаки значений в узлах чередовались (выровненный альтернанс)
            t(left_index + 1) = x_max;
            f(left_index + 1) = F(x_max);
        end
    end
    
    figure;
    Draw(polinom, was_t, n, left_board, right_board, x_max); % рисует графики
   
    if (number_of_step > 1)
        if (abs(sigma-Sigma) < epsilon) % sigma - с прошлого шага, Sigma - с этого
            delta_sigmas = sigma-Sigma
            break; % завершаем работу алгоритма Ремеза
        else
            sigma = Sigma; % запоминаем текущую сигму как будущую предыдущую
        end
    end
end % конец алгоритма Ремеза

end