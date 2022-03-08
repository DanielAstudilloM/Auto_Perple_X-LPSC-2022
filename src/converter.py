import numpy as np
from PIL import Image

names = np.loadtxt("last_source.txt",dtype=str)

for h in names:
    try:
        im = Image.open(h.replace(".txt",".ps"))
        fig = im.convert("RGBA")
        fig.save(h.replace(".txt",".png"),lossless=True)
    except:
        print(h,"is missing")