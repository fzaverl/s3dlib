import matplotlib.pyplot as plt
import numpy as np

# 1. Define function to examine ...............

def f(x,n) :
    y = x**n
    return y

# 2. Setup and map line .......................

x = np.linspace(0, 2, 100)
y1 = f(x,1)
y2 = f(x,2)
y3 = f(x,3)

# 3. Construct figure, add line, and plot .....

plt.xlabel('X')
plt.ylabel('Y = f(X,n)')
plt.title(r'Y = $\ X^n$')

plt.plot(x, y1, label='n = 1')
plt.plot(x, y2, label='n = 2')
plt.plot(x, y3, label='n = 3')
plt.legend()

plt.show()