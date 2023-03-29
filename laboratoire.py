import matplotlib.pyplot as plt
import numpy as np


def Fx(x):
    if x >= 0:
        return 0.5 * (1 + np.sqrt(1 - np.exp(-1 * x**2 * np.sqrt(np.pi / 8))))
    else:
        return 1 - Fx(-x)

def x(Fx):
    if Fx >= 0.5:
        return np.sqrt(-np.log(-1 * (2 * Fx - 1)**2 + 1)/np.sqrt(np.pi / 8))
    else:
        return -np.sqrt(-np.log(-1 * (2 * Fx - 1)**2 + 1)/np.sqrt(np.pi / 8))

def exercice4(mu = 10, std = 2):
    # Generating random numbers
    N = 10000
    _fx = [Fx(_x) for _x in np.linspace(-5, 5, N)]

    plt.plot(_fx)
    plt.show()

    N = 10000
    vals = np.random.rand(N)

    y = [x(_Fx) for _Fx in vals]
    plt.hist(y, 100)
    plt.show()




def exercice3(mu = 10, std = 2):
    # Generating random numbers
    N  = 10000
    U1 = np.random.rand(N)
    U2 = np.random.rand(N)
    X  = mu + std * np.cos(2 * np.pi * U1) * np.sqrt(-2 * np.log(U2))

    # Displaying histogram and relative frequency plot
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    n = ax1.hist(X, align="left")
    b = n[1]
    n = n[0]
    ax2.plot(b[:-1], n / N, color="red")
    plt.show()

    # Calculating characteristic parameters
    average  = np.average(X)
    stddev   = np.std(X)
    variance = stddev**2
    mse      = (average - mu)**2
    print(average)
    print(stddev)
    print(variance)
    print(mse)

def exercice2(mu = 10, std = 2):
    # Generating random numbers
    N = 10000
    vars = mu + std * np.random.randn(N)

    # Displaying histogram and relative frequency plot
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    n = ax1.hist(vars, align="left")
    b = n[1]
    n = n[0]
    ax2.plot(b[:-1], n / N, color="red")
    plt.show()

    # Calculating characteristic parameters
    average  = np.average(vars)
    stddev   = np.std(vars)
    variance = stddev**2
    mse      = (average - mu)**2
    print(average)
    print(stddev)
    print(variance)
    print(mse)

def exercice1():
    # Generating random numbers
    N = 10000
    bins = 10
    vars = 5 - 10 * np.random.rand(N)

    # Displaying histogram and relative frequency plot
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    n = ax1.hist(vars, bins, align="left")[0]
    ax2.plot([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4], n / N, color="red")
    plt.show()

    # Calculating characteristic parameters
    average  = np.average(vars)
    stddev   = np.std(vars)
    variance = stddev**2
    mse      = average**2
    print(average)
    print(stddev)
    print(variance)
    print(mse)

if __name__ == '__main__':
    exercice4()
