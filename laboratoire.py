import matplotlib.pyplot as plt
import numpy as np

def exercice3(mu, ):
    # Generating random numbers
    N  = 10000
    U1 = np.random.rand(N)
    U2 = np.random.rand(N)
    X  = 1- + 2 * np.cos(2 * np.pi * U1) * np.sqrt(-2 * np.log(U2))

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
    exercice3()
