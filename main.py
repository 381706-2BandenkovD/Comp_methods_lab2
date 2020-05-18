import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
from math import cos
from tqdm import tqdm
from time import sleep
from datetime import datetime
import PySimpleGUI as sg
import time
plt.rcParams.update({'figure.figsize': '9.5, 6', "figure.facecolor": 'lightgoldenrodyellow', 'axes.edgecolor': 'darkgoldenrod'})

pi = 3.1415926535

def func(x, l, f1, f2):
    return 1/l + f1 * cos((pi*x)/l) + f2 * cos(2*(pi*x)/l)

def bfunc(x, l, b0, b1, b2):
    return b0 + b1 * cos((pi*x)/l) + b2 * cos(2*(pi*x)/l)

def integrate(h, fu): #Функция численного интегрирования
    res = (h/3)*(fu[0] + fu[len(fu) - 1])
    for i in range(1, len(fu) - 1, 2):
        res += (h/3)*(4*fu[i] + 2*fu[i + 1])
    return res

def tridiagAlg(a, b, c, func, count):#Метод прогонки
    A = []
    B = []
    res = [0] * count
    A.append(-c[0]/b[0])
    B.append(func[0]/b[0])
    for i in range(1, count):
        A.append(-c[i] / (a[i] * A[i - 1] + b[i]))
        B.append((func[i] - a[i] * B[i - 1]) / (a[i] * A[i - 1] + b[i]))
    res[count-1] = B[count - 1]
    for i in range(count - 2, -1, -1):
        res[i] = (A[i] * res[i + 1] + B[i])
    return res

def onButtonAddClicked(event):
    global graph_axesm, delta_t, delta_x
    pogreshnostPython = 0.0000001
    k = delta_t / (delta_x * delta_x) + pogreshnostPython
    if (k < 0.25):
        solution(graph_axes)
    else:       
        sg.Popup('Error!', 'Шаг по x и шаг по t не удовлетворяют условию устойчивости:\nτ/(h*h) < 1/4')


def onButtonСreateCliked(event):
    global graph_axes
    alternativeSolution(graph_axes)

def onButtonClearClicked(event):
    global graph_axes
    graph_axes.clear()
    graph_axes.grid()
    plt.draw()

def alternativeSolution(graph_axes):
    global x_val, resA, func_val
    graph_axes.plot(x_val, resA, label = 'Конечная температура Часть(А)', color  ='greenyellow')
    graph_axes.legend()
    plt.draw()

def progressBarSolution():
    global resA, resB, func_val, bfunc_val 
    start_time = datetime.now()
    progressbar = [[sg.ProgressBar(len(func_val), orientation='h', size=(41, 10), key='progressbar')]]
    outputwin = [[sg.Output(size=(60,7))]]

    layout = [
        [sg.Frame('Progress',layout= progressbar)],
        [sg.Frame('Output', layout = outputwin)],
        [sg.Cancel()]
    ]

    window = sg.Window('Custom Progress Meter', layout)
    progress_bar = window['progressbar']

    event, values = window.read(timeout=10)
    for i,item in enumerate(func_val):
        print("φ[", i, "] = ", item, "   b[", i, "] = ", bfunc_val[i] )
        sleep(0.003)
        progress_bar.UpdateBar(i + 1)
    end_time = datetime.now()
    print('Время выполнения: {}'.format(end_time - start_time))
    if event == 'Cancel'  or event is None:
        window.close()


def solution(graph_axes):
    global x_val, resA, resB, func_val, bfunc_val
    func_val = []
    bfunc_val = []
    slices1 = [[]]
    slices2 = [[]]
    count_N = int(_len/delta_x)
    count_T = int(time/delta_t)
    # Вычисление значений функции и заполнение первого слоя сетки
    for i in range(0, count_N):
        func_val.append(func(i*delta_x, _len, f1, f2))
        bfunc_val.append(bfunc(i*delta_x, _len, b0, b1, b2))
        slices1[0].append(func_val[i])
        slices2[0].append(func_val[i])
    # Заполнение матрицы коэффициентов для метода прогонки
    coeff_a = [0.0]
    coeff_b = [1.0]
    coeff_c = [-1.0]
    for i in range(1, count_N - 1):
        coeff_a.append(delta_t / (delta_x * delta_x))
        coeff_b.append(-1 - 2*delta_t / (delta_x * delta_x))
        coeff_c.append(delta_t / (delta_x * delta_x))
    coeff_a.append(-1.0)
    coeff_b.append(1.0)
    coeff_c.append(0.0)
    #Вычисление последующих слоев сетки
    for i in range(1, count_T):
        I = integrate(delta_x, bfunc_val)
        fu = [0]
        fu2 = [0]
        slices1.append([])
        slices2.append([])
        #Вычисляем правую часть системы для прогонки
        for j in range(1, count_N - 1):
            fu.append(-slices1[i - 1][j] * ((bfunc_val[j] - I) * delta_t * delta_t  + 1.0))
            fu2.append(-slices2[i - 1][j] * (bfunc_val[j] * delta_t * delta_t + 1.0))
        fu.append(0)
        fu2.append(0)
        #Метод прогонки для системы из B
        res = tridiagAlg(coeff_a, coeff_b, coeff_c, fu, count_N)
        for j in range(0, count_N):
            slices1[i].append(res[j])
        #Метод прогонки для системы из A
        res2 = tridiagAlg(coeff_a, coeff_b, coeff_c, fu2, count_N)
        for j in range(0, count_N):
            slices2[i].append(res2[j])
    I = integrate(delta_x, slices2[count_T - 1])
    resB = []
    resA = []
    for j in range(0, count_N):
        resA.append(slices2[count_T - 1][j] / I)
        resB.append(slices1[count_T - 1][j])
    x_val = []
    for i in range(0, count_N):
        x_val.append(i * delta_x)
    progressBarSolution()
    graph_axes.plot(x_val, func_val, label = 'Начальная температура', color ='b')
    graph_axes.plot(x_val, slices1[count_T - 1], label = 'Конечная температура', color ='r')
    graph_axes.legend()
    plt.draw()

if __name__ == "__main__":
    global graph_axes, flag, b0, b1, b2, f0, f1, f2, delta_t, delta_x, time, _len
    global button_pb, func_val, resA, resB, x_val, func_val, bfunc_val, slices1, slices2
    global start_time, end_time
    func_val = []
    resA = []
    resB = []
    x_val = []

    func_val = []
    bfunc_val = []
    slices1 = [[]]
    slices2 = [[]]

    fig, graph_axes = plt.subplots()
    graph_axes.grid()
    fig.subplots_adjust(left=0.37, right=0.98, top=0.95, bottom=0.2)
    graph_axes.set_xlabel('X')
    graph_axes.set_ylabel('T')

    def submitTime(text):
        global time
        try:
            time = float(text)
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для параметра 'time' были использованы значения по умолчанию = ", time)

    def submit_Len(text):
        global _len
        try:
            _len = float(text)
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для шага '_len' были использованы значения по умолчанию = ", _len)

    def submitB0(text):
        global b0
        try:
            b0 = float(text)
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для начального 'b0' были использованы значения по умолчанию = ", b0)

    def submitB1(text):
        global b1
        try:
            b1 = float(text)
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для начального 'b1' были использованы значения по умолчанию = ", b1)

    def submitB2(text):
        global b2
        try:
            b2 = float(text)
            return b2
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для колличества точек 'b2' были использованы значения по умолчанию = ", b2)

    def submitF0(text):
        global f0
        try:
            #f0 = float(text)
            f0 = 1 / _len
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для параметра 'f0' были использованы значения по умолчанию = ", f0)

    def submitF1(text):
        global f1
        try:
            f1 = float(text)
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для параметра 'f1' были использованы значения по умолчанию = ", f1)

    def submitF2(text):
        global f2
        try:
            f2 = float(text)
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для параметра 'f2' были использованы значения по умолчанию = ", f2)

    def submitDeltaT(text):
        global delta_t
        try:
            delta_t = float(text)
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для параметра 'delta_t' были использованы значения по умолчанию = ", delta_t)

    def submitDeltaX(text):
        global delta_x
        try:
            delta_x = float(text)
        except ValueError:
            print("Вы пытаетесь ввести не число")
            print("Для параметра 'delta_x' были использованы значения по умолчанию = ", delta_x)

    axes_button_add = plt.axes([0.37, 0.049, 0.29, 0.09])
    button_add = Button(axes_button_add, 'Решение')
    button_add.on_clicked(onButtonAddClicked)

    axes_button_clear = plt.axes([0.05, 0.049, 0.29, 0.09])
    button_clear = Button(axes_button_clear, 'Очистить')
    button_clear.on_clicked(onButtonClearClicked)

    axes_button_create = plt.axes([0.69, 0.049, 0.29, 0.09])
    button_create = Button(axes_button_create, 'Альтернативное решение')
    button_create.on_clicked(onButtonСreateCliked)

    axbox = plt.axes([0.22, 0.88, 0.08, 0.06])
    _len_box = TextBox(axbox, 'Длина стержня =', initial="7")
    _len = 7.
    _len_box.on_submit(submit_Len)
    
    axbox = plt.axes([0.22, 0.78, 0.08, 0.06])
    time_box = TextBox(axbox, 'Время воздействия =', initial= "1")
    time = 1.
    time_box.on_submit(submitTime)
    
    axbox = plt.axes([0.22, 0.68, 0.08, 0.06])
    delta_t_box = TextBox(axbox, 'Шаг по \n времени =', initial="0.01")
    delta_t = 0.01
    delta_t_box.on_submit(submitDeltaT)
    
    axbox = plt.axes([0.05, 0.68, 0.08, 0.06])
    deltax_box = TextBox(axbox, ' Шаг \nпо х =', initial= "0.2")
    delta_x = 0.2
    deltax_box.on_submit(submitDeltaX)
    
    axbox = plt.axes([0.05, 0.55, 0.08, 0.06])
    bo_box = TextBox(axbox, 'b₀=', initial="0")
    b0 = 0.
    bo_box.on_submit(submitB0)
    
    axbox = plt.axes([0.05, 0.45, 0.08, 0.06])
    b1_box = TextBox(axbox, 'b₁=', initial= "-7")
    b1 = -7
    b1_box.on_submit(submitB1)
    
    axbox = plt.axes([0.05, 0.35, 0.08, 0.06])
    b2_box = TextBox(axbox, 'b₂=', initial="0")
    b2 = 0.
    b2_box.on_submit(submitB2)
    
    axbox = plt.axes([0.22, 0.55, 0.08, 0.06])
    f1_box = TextBox(axbox, 'φ₀=', initial= "0")
    f1 = 0.
    f1_box.on_submit(submitF1)
    
    axbox = plt.axes([0.22, 0.45, 0.08, 0.06])
    f2_box = TextBox(axbox, 'φ₁=', initial= "0")
    f2 = 0.
    f2_box.on_submit(submitF2)
    
    plt.show()