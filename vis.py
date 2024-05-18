import numpy as np
import matplotlib.pyplot as plt
import sys
from PIL import Image


inp_fn = "image"
if len(sys.argv) > 1:
    inp_fn = sys.argv[1]

data = np.fromfile(inp_fn + ".bin" , dtype="uint8")
shape = np.frombuffer(data[:(3*4)], dtype="int32")
img = np.frombuffer(data[(3*4):], dtype="float64").reshape(*shape)
counts = None
if img.shape[2] == 4: #Includes count map
    img, counts = img[:, :, :3], img[:, :, 3]

print("Img shape: ", img.shape)
print("Img max:", img.max())
print("Img min:", img.min())

outp_fn = inp_fn
if len(sys.argv) > 2:
    outp_fn = sys.argv[2]
im = Image.fromarray(np.array(256*img).astype(np.uint8))
im.save(outp_fn + ".png")

if counts is not None:
    print("Counts: ", counts.shape, " min: ", counts.min(), " max: ", counts.max())
    vmax=min(counts.mean()*5, counts.max())
    plt.imshow(counts, vmax=vmax)
    # plt.colorbar()
    # plt.show()
    plt.axis("off")
    plt.savefig(outp_fn + "_sample_counts.png", bbox_inches="tight")
    # Image.fromarray(np.repeat((counts/vmax*256)[:,:,None], 3, axis=2).astype(np.uint8)).save(outp_fn + "_sample_counts.png")