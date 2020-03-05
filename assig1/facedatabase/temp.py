import os
import shutil

id = 0
for img in os.listdir('new'):
    count = 0
    while count < 16:
        shutil.copy('new/'+img, 'sad/'+str(id)+'.jpg')
        id = id+1
        count += 1
