# раскомментировать при необходимости
# !pip install interval-py
import sympy as sym
import intervalpy as ival
import numpy as np


# функция пример 1
def f(x):
    return x ** 2 - 2


# производная функции 1
def df(x):
    return 2 * x


# функция пример 2
def f1(x):
    return (x + 1) ** 3 / (x ** 2) - 7.1


def f1_symbolic():
    x = sym.Symbol("x")
    f = (x + 1) ** 3 / (x ** 2) - 7.1
    return f, x


# функция пример 3
def f2_symbolic(n=2):
    x = sym.Symbol("x")
    f = 1
    for i in range(n):
        f *= (x - i)
    return f, x


# функция пересечения 2 интервалов
def intersec(a, b):
    if a[1] < b[0] or b[1] < a[0]:
        return None
    else:
        return ival.Interval([max(a[0], b[0]), min(a[1], b[1])])


# функция пересечения множества интервалов
def intersec_mi(mi1, mi2):
    a = []
    ch1 = False
    ch2 = False
    if type(mi1) == ival.Interval:
        ch1 = True
    if type(mi1) == ival.Interval:
        ch2 = True
    if ch1 and ch2:
        intersection = intersec(mi1, mi2)
        if type(intersection) != type(None):
            return intersection
    elif ch1 and not (ch2):
        for ival1 in mi2:
            intersection = intersec(ch1, ival1)
            if type(intersection) != type(None):
                a.append(intersection)
    elif ch2 and not (ch1):
        for ival1 in mi1:
            intersection = intersec(ch2, ival1)
            if type(intersection) != type(None):
                a.append(intersection)
    for ival1 in mi1:
        for ival2 in mi2:
            intersection = intersec(ival1, ival2)
            if type(intersection) != type(None):
                a.append(intersec(ival1, ival2))
    return a


def c_minus(x, df):
    x_right = x[1]
    x_left = x[0]
    df_at_x_right = df[1]
    df_at_x_left = df[0]

    if df_at_x_right <= 0:
        return x_right
    elif df_at_x_left >= 0:
        return x_left
    else:
        return (df_at_x_right * x_left - df_at_x_left * x_right) / (df_at_x_right - df_at_x_left)


def c_plus(x, df):
    x_right = x[1]
    x_left = x[0]
    df_at_x_right = df[1]
    df_at_x_left = df[0]

    if df_at_x_right <= 0:
        return x_left
    elif df_at_x_left >= 0:
        return x_right
    else:
        return (df_at_x_left * x_left - df_at_x_right * x_right) / (df_at_x_left - df_at_x_right)


# функция для бицентрированного оператора
def bicentered_calcul(x, df):
    c_l = c_minus(x, df)
    c_u = c_plus(x, df)
    return c_l, c_u


# оператор Ньютона
def Newton_method(f, df, x, midpoint):
    return midpoint - f(midpoint) / df(x)


# итеративная процедура
def iterative_method(f, df, X_ini, epsilon=1e-6, max_iterations=1000, verbose=False, bicentered=False):
    iterations_count = 0
    roots = []
    df0 = 0
    if verbose:
        print("ini x", X_ini)
    queue = [X_ini]
    while len(queue) > 0:
        if verbose:
            print("Iter", iterations_count)
            print("\tqueue", queue)
        X = queue.pop(0)
        if verbose:
            print("\tX", X)
        if ival.Interval([0, 0]).isIn(f(X)):
            if iterations_count < max_iterations and X.width() > epsilon:
                iterations_count += 1
                if bicentered:
                    c_l, c_h = bicentered_calcul(X, df(X))
                    root_l = Newton_method(f, df, X, c_l)
                    root_h = Newton_method(f, df, X, c_h)
                    root = intersec_mi(root_l, root_h)
                else:
                    midpoint = X.mid()
                    root = Newton_method(f, df, X, midpoint)
                    print(root)
                if verbose:
                    print("\tNewton", root)
                if type(root) != ival.interval.ExtendedInterval and type(root) != list:
                    root = [root]
                for x_newton in root:
                    x_intersec = intersec(x_newton, X)
                    if type(x_intersec) != type(None):
                        if x_intersec[0] == X[0] and x_intersec[1] == X[1]:
                            queue.append(ival.Interval([X[0], midpoint]))
                            queue.append(ival.Interval([midpoint, X[1]]))
                        else:
                            queue.append(x_intersec)

            else:
                if verbose:
                    print("\troot = ", X)
                roots.append(X)
    return roots, iterations_count


# инициализация интервала и вызов итеративной процедуры для функции 1
intervals = ival.Interval([-2, 2])
print("Классический метод Ньютона")
roots, iters = iterative_method(f, df, intervals, bicentered=False)
print("Все корни уравнения:", roots)
print("Количество итераций", iters)

print("Бицентрированный метод Ньютона")
roots, iters = iterative_method(f, df, intervals, bicentered=True)
print("Все корни уравнения:", roots)
print("Количество итераций", iters)


# расчет бицентрических центров
def calculate_bicentered(x, df, coef=1.2):
    c_l, c_h = bicentered_calcul(x, coef * df(x))
    return c_l, c_h


f1_sym, x_sym = f1_symbolic()
df1_sym = f1_sym.diff(x_sym)
df1 = sym.lambdify(x_sym, df1_sym)
# инициализация интервала и вызов итеративной процедуры для функции 2
intervals = ival.Interval([0.2, 7])

roots, iters = iterative_method(f1, df1, intervals, bicentered=False)
print("Все корни уравнения:", roots)
print("Количество итераций", iters)

print("Бицентрированный метод Ньютона")
roots, iters = iterative_method(f1, df1, intervals, bicentered=True)
print("Все корни уравнения:", roots)
print("Количество итераций", iters)
# проверка работы программы для функции 2
for root in roots:
    print(f1(root))

n = 8
f2_sym, x_sym = f2_symbolic(n)
print(f2_sym)
f2 = sym.lambdify(x_sym, f2_sym)
df2_sym = f2_sym.diff(x_sym)
print(df2_sym)
df2 = sym.lambdify(x_sym, df2_sym)
# инициализация интервала и вызов итеративной процедуры для функции 3
intervals = ival.Interval([0, n + 1])
roots, iters = iterative_method(f2, df2, intervals, bicentered=False)
print("Все корни уравнения:", roots)
print("Количество итераций", iters)

print("Бицентрированный метод Ньютона")
roots, iters = iterative_method(f2, df2, intervals, bicentered=True)
print("Все корни уравнения:", roots)
print("Количество итераций", iters)
# проверка работы программы для функции 2
for root in roots:
    print(f2(root))
