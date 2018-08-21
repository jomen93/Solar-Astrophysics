import numpy as np 
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.image as mpimg



img=mpimg.imread('campoFinal.png')
imgplot = plt.imshow(img)
plt.show()
