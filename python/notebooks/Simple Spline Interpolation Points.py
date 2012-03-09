# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

### Function for inner ellipse
def f(x):
    return np.sqrt((1.0/64)*(0.25 - x**2))

# <codecell>

X = linspace(0.0,0.5,12)

# <codecell>

Y = -f(X)

# <codecell>

for i,j in zip(X,Y): print (i, j, 0.00)

# <codecell>

### Function for outer ellipse
def g(x):
    return np.sqrt((1.0/64)*(1.0 - x**2))

# <codecell>

A = linspace(0.0,1.0,12)

# <codecell>

B = -g(A)

# <codecell>

for i,j in zip(A,B): print (i, j, 0.00)

# <codecell>


