# Прямая задача внешней баллистики
# Полёт пули

## Цели работы

Изучить и смоделировать полёт пули в зависимости от заданных начальных условий и встречаемых её на своём пути препятствий.

## Структура проекта

#### Папки

`/build` – автоматически сгенерированная PyCharm-ом папка с единственным файлом `/build/setup.py`.

`/butcher` – папка с текстовыми файлами, содержащими таблицы Бутчера:
* `/butcher/midpoint` – таблица Бутчера для метода средней точки;
* `/butcher/RK5` – таблица Бутчера для метода Рунге-Кутты 4-ого порядка;
* `/butcher/RK5` – таблица Бутчера для метода Рунге-Кутты 5-ого порядка;
* `/butcher/DP5` – таблица Бутчера для метода Дормана-Принца 5-ого порядка;
* `/butcher/DP8` – таблица Бутчера для метода Дормана-Принца 8-ого порядка.

Так же папка `/butcher` содержит изображение `/butcher/butcher_tables.jpg` для более наглядного сравнения таблиц.

![butcher_tables.jpg](/butcher/butcher_tables.jpg)

`/graphs` – папка с некоторыми результами работы программы.

#### Файлы

`/archive.py` – файл с одной из старых версией кода программы,
где ещё не использовался метод библиотеки numpy – np.dot(), который сократил колличество строчек кода во много раз
(файл представлен для того, чтобы любой начинающий математик/программист увидел практическую пользу скалярного произведения).

`/main.py` – файл с текущим кодом программы.

#### Функции `/main.py`

| Название функции                                                    | Аргументы функции                                                                                                                                                                                                                                                                                                                                                           | Что делает и что возвращает                                                                                                                                                                   |
|---------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `def parse(file_name: str)`                                         | `file_name` – название файла таблицы Бутчера.                                                                                                                                                                                                                                                                                                                               | Функция преобразовывает файл с таблицей Бутчера в np.array и возвращает его                                                                                                                   |
| `def cx(mah)`                                                       | `mah` – скорость тела в махах.                                                                                                                                                                                                                                                                                                                                              | cx – вспомогательная функция закона сопротивления воздуха. <br/>Возвращает коэффициент, необходимый `def Force(F_para, S, speed)`                                                             |
| `def Force(F_para, S, speed)`                                       | `F_para` – логическое значение, обозначающее, учитываем ли мы сопротивление воздуха или нет, <br/>`S` – площадь миделева сечения тела, <br/>`speed` – скорость тела (в м/с).                                                                                                                                                                                                | Функция рассчитывает и возвращает модуль силы лобового сопротивления воздуха <br/>(см. Математические выкладки  -> Сопротивление воздуха)                                                     |
| `def air(mass, S, startX, startY, speed, teta, l, butcher, F_para)` | `mass` – масса пули, <br/>`S` – площадь миделева сечения пули, <br/>`startX` – начальный горизонтальное положение, <br/>`startY` – начальная высота, <br/>`speed` – начальная скорость пули, <br/>`teta` – начальный угол, <br/>`l` – горизонтальное расстояние до препятствия, <br/>`butcher` – таблица Бутчера (np.array), <br/>`F_para` – --//--.                        | Функция пошагового расчёта траектории пули в воздухе <br/>(см. Математические выкладки  -> Рассматриваемые системы дифференциальных уравнений). <br/>Возвращает массивы времени и высоты.     |
| `def block(mass, S, startX, startY, speed, teta, ρ, σ, butcher)`    | `mass` – масса пули, <br/>`S` – площадь миделева сечения пули, <br/>`startX` – начальное горизонтальное положение, <br/>`startY` – начальная высота, <br/>`speed` – начальная скорость пули, <br/>`teta` – угол соприкосновения с прептствием, <br/>`ρ` – плотность препятствия, <br/>`σ` – предельное напряжение препятствия, <br/>`butcher` – таблица Бутчера (np.array). | Функция пошагового расчёта траектории пули в препятствии <br/>(см. Математические выкладки  -> Рассматриваемые системы дифференциальных уравнений). <br/>Возвращает массивы времени и высоты. |
| `def dd(i, argX, argY)`                                             | `i` – , <br/>`argX` – , <br/>`argY` – .                                                                                                                                                                                                                                                                                                                                     | Вспомогательная функция                                                                                                                                                                       |
| `def research()`                                                    |                                                                                                                                                                                                                                                                                                                                                                             | Функция задаёт все начальные характеристики тела, системы отсчёта и препятствия, запускает полёт в воздухе до препятствия, запускает встречу с препятствием, отрисовывает графики             |

## Математические выкладки

#### Формулировка задачи

Моделирование траектории движения шарообразной пули, проходящей через заданные препятствия.

|                                | Дано:                                                                                                                                                                                      | 
|--------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Характеристики тела            | m – масса пули; <br/>S – площадь миделева сечения тела; <br/>v₀ – начальная скорость пули;                                                                                                 |
| Характеристики системы отсчёта | y₀ – начальная высота пули; <br/>v₀ – начальная скорость пули; <br/>θ₀ – начальный угол между вектором скорости и плоскостью горизонта; <br/>l – горизонтальное расстояние до препятствия; | 
| Характеристики препятствия     | ρ – плотность препятствия; <br/>σ – предельное напряжение для среды.                                                                                                                       | 

![σ.png](/σ.png)

Задача: построить траекторию движения пули.

#### Сопротивление воздуха

Модуль силы лобового сопротивления в зависимости от скорости движения тела вычисляется по формуле:

F(v) = i * S * p * v^2 / 2 * c_x(v / a),
$$ F(v) = i \cdot S \cdot \cfrac{p \cdot v^2}{2} \cdot c_x(cfrac{v}{a}, $$

где
* S – площадь миделева сечения тела;
* i – коэффициент формы;
* $ c_x $ – функция закона сопротивления воздуха;
* a – скорость звука;
* ρ – плотность атмосферы.

Заметим, что скорость звука в среде a является функцией температуры, а плотность ρ зависит от температуры и давления среды.
Влажность воздуха также влияет на величину скорости звука и плотность среды. <br/>
Форма тела учитывается с помощью коэффициента формы i, который для шара = 0.47. <br/>
Функции закона сопротивления cx (функция 1943 г.) различается для тел разной формы и рассчитывается по следующей формуле:

            0.157,                      при M ≤ 0.73
            0.033∙M + 0.133,            при 0.73 ≤ M < 0.82
            0.161 + 3.9∙(M – 0.823)²,   при 0.82 ≤ M < 0.91
            1.5∙M – 1.176,              при 0.91 ≤ M < 1.00
    c_x =   0.384 – 1.6(M – 1.176)²,    при 1.00 ≤ M < 1.18
            0.384∙sin(1.185/M),         при 1.18 ≤ M < 1.62
            0.29/M + 0.172,             при 1.62 ≤ M < 3.06
            0.316 – 0.016∙M,            при 3.06 ≤ M < 3.53
            0.259,                      при 3.53 ≤ M

#### Встреча с препятствием

В нашей работе так же рассмотрим ситуацию, когда по траектории полёта снаряда встречается препятствие. <br/>
При встрече с преградой пуля обладает кинетической энергией:

E_sum = m * v^2 / 2,

где m – масса пули, кг; v – скорость встречи пули с преградой, м/с.

Кинетическая энергия пули в общем случае расходуется на преодоление сил сопротивления преграды, деформацию и разрушение пули,
нагревание преграды и пули, а также на сообщение запреградной кинетической энергии пули при пробивании тонких преград.
Наибольшая часть энергии расходуется на преодоление сил сопротивления преграды и её пластической деформации.
Около 10% энергии может затрачиваться на упругий удар. Затраты на пластические деформации преграды возрастают с ростом толщины преграды и увеличением ее вязкости. <br/>
Полная энергия, поглощенная средой в процессе пластических деформаций при расширении отверстия рассчитана по формуле Томсона:

E_med_v = π*R^2*h*w * (σ/2 + ρ*(v*R/L));

Дифференциальное уравнение для системы движения будет выглядеть следующим образом:

E_med_v/dh = π*R^2*w * (σ/2 + ρ*(v*R/L));

Для определения угла отклонения при движении в твёрдом теле воспользуемся скоростью пули и скоростью пули по оси х. <br/>
Таким образом, при столкновении с препятствием, изначально необходимо рассчитать начальную кинетическую энергию, а так же начальную кинетическую энергию по оси x:

E_sum_v = m * v^2 / 2

E_sum_u  = m * u^2 / 2

P.S. Строчки `/main.py` 158-159.

#### Рассматриваемые системы дифференциальных уравнений

| Рассмотриваемая среда | Система дифференциальных уравнений                                                                                                                                                                              | , где                                                                                                                                                                                                                                                                                                                                                                      | Код                                 |
|-----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------|
| Воздух                | du/dx = -F(v)/(m*v); <br/>dγ/dx = -g/(u^2); <br/>dy/dx = γ; <br/>dt/dx = 1/u; <br/>v = u * sqrt(1+γ^2)                                                                                                          | x, y – горизонтальная и вертикальная координаты тела соответственно; <br/>t – время; <br/>m – масса тела; <br/>v – скорость движения тела; <br/>u – горизонтальная компонента скорости; <br/>γ – тангенс угла между вектором скорости и плоскостью горизонта; <br/>F(v) – модуль силы лобового сопротивления воздуха; <br/>g – ускорение свободного падения (g = 9,80665). | Строчки `/main.py` 115-119, 124-128 |
| Препятствие (твёрдое) | E_med_v/dx = π*R^2*w * (σ/2 + ρ*(v*R/L)); <br/>E_med_u/dx = π*R^2*w * (σ/2 + ρ*(u*R/L)); <br/> dv/dx = -sqrt(2*E_med_v/m); <br/>du/dx = -sqrt(2*E_med_u/m); <br/>dγ/dx = θ^2; <br/>dy/dx = γ; <br/>dt/dx = 1/u; | --//--; <br/> R – калибр пули, м; <br/>h – толщина пробиваемого слоя преграды, м; <br/>σ – предельное напряжение для среды, Па; <br/>ρ – плотность среды, кг/м3; <br/>v – скорость пули, м/с; <br/>L – длина головной части пули, м.                                                                                                                                       | Строчки `/main.py` 176-182, 188-193 |

## Некоторые результаты работы программы

![/graphs something](/graphs/...)

***

## Список использованных источников

1) Козлитин И.А. – Восстановление входных параметров расчета внешней баллистики тела по результатам траекторных измерений. URL: <br/>
https://www.mathnet.ru/links/f9757b34505280236c2be0aa2743fa93/mm3892.pdf

2) Ефремов А. К. – Аппроксимация закона сопротивления воздуха 1943 г. <br/>
https://cyberleninka.ru/article/n/approksimatsiya-zakona-soprotivleniya-vozduha-1943-g/viewer

    2.1) Точные значения функции закона сопротивления cx "образца 1943 г." <br/>
    https://popgun.ru/viewtopic.php?p=8041476

3) Кривченко А.Л. – Ударно-волновые процессы взаимодействия Высокоскоростных элементов с конденсированными средами. URL: <br/>
http://d21221701.samgtu.ru/sites/d21221701.samgtu.ru/files/aleksenceva_dis.pdf
