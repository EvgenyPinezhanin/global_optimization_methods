# Методы решения задач глобальной оптимизации
Библиотека методов глобальной оптимизации позволяет решать четыре типа задач оптимизации: одномерную задачу глобальной оптимизации, одномерную задачу условной глобальной оптимизации и их многомерные аналоги, многомерную задачу условной глобальной оптимизации, многомерную задачу условной глобальной оптимизации.  
В общем случае задачу глобальной оптимизации с ограничениями можно представить как  

min{φ(y): y ∈ D, g<sub>i</sub>(y) ≤ 0, 1 ≤ i ≤ m},  
D = {y ∈ R<sup>N</sup>: a<sub>j</sub> ≤ y<sub>j</sub> ≤ b<sub>j</sub>, 1 ≤ j ≤ N},
	
где целевая функция φ(y) (далее g<sub>(m+1)</sub>(y)) N-мерная функция, а g<sub>i</sub>(y), 1 ≤ i ≤ m, ограничения.

В библиотеке реализованы три метода: Алгоритм глобального поиска(**АГП**), индексный метод учета ограничений(**ИМГО**) и многомерный обощенный алгоритм глобального поиска(**МОАГП**).  
Важным условием корректной работы каждого метода является **липшицевость** целевой функции и функций-ограничений:  

|g<sub>i</sub>(y<sub>1</sub>) - g<sub>i</sub>(y<sub>2</sub>)| ≤ L<sub>i</sub>‖y<sub>1</sub> - y<sub>2</sub>‖, 1 ≤ i ≤ m + 1.

В папке **samples** расположены простые примеры задач для всех методов, более сложные тестовые функции, а для методов **ИМГО** и **МОАГП** также проведены вычисления **операционной характеристики**(для создания большого объема функций использован генератор **GCGen**).

Папка **GCGen** - реализация глобального генератора задач оптимизации с ограничениями, внутри так же есть реализация функций семейств **Гришагина** и **GKLS**.

Папка **scripts** - скрипты для отрисовки всех графиков, работают с помощью инструмента **Gnuplot**. Домашняя страница проекта **Gnuplot**: [gnuplot.info](http://www.gnuplot.info)
