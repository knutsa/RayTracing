import numpy as np
import matplotlib.pyplot as plt


data = np.fromfile("img.bin", dtype="uint8")
shape = np.frombuffer(data[:(3*4)], dtype="int32")
img = np.frombuffer(data[(3*4):], dtype="float64").reshape(*shape)


print("Img shape: ", img.shape)
print("Img max:", img.max())
print("Img min:", img.min())

plt.imshow(img)
plt.savefig("img.png")