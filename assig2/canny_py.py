import cv2
import numpy as np
from matplotlib import pyplot as plt


fd = open('raw_images/leaf.raw', 'rb')
rows = 243
cols = 190
f = np.fromfile(fd, dtype=np.uint8,count=rows*cols)
img = f.reshape((rows, cols)) #notice row, column format
fd.close()

# img = cv2.imread('/home/jhchenjh/4208_test/EE4208ComputerVision/assig2/raw_images/cana.raw')
print(img.shape)
edges = cv2.Canny(img,100,200)

# plt.subplot(121),plt.imshow(img)
# plt.title('Original Image'), plt.xticks([]), plt.yticks([])
# plt.subplot(122),plt.imshow(edges)
# plt.title('Edge Image'), plt.xticks([]), plt.yticks([])

# plt.show()

cv2.imshow('original', img)
cv2.imshow('edge', edges)
cv2.waitKey()
