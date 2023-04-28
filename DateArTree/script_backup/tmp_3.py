
import numpy as np
import matplotlib.pyplot as plt

x, y = np.random.random(size=(2,6))

print(x)
print(y)

for i in range(0, len(x), 2):
    plt.plot(x[i:i+2], y[i:i+2], 'ro-')

plt.show()
