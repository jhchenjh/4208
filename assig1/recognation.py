import cv2
import numpy as np
import base_faces_lib as bfl
imgNum = 100

#import cv2 Classifiers
face_cascade = cv2.CascadeClassifier('haarcascades/haarcascade_frontalface_default.xml')
video_capture = cv2.VideoCapture(0)
#number = 0

eigenFaces = bfl.calcEigFaces()
eigenVector = bfl.eigenVectorImport()
face_avg = bfl.avgMatrix()
font = cv2.FONT_HERSHEY_SIMPLEX

cnt = 0
while(1):
    # Capture frame-by-frame
	ret, frame = video_capture.read()
	print(cnt)
	gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
	gray = cv2.equalizeHist(gray)
	faces = face_cascade.detectMultiScale(
	    gray,
	    scaleFactor=1.3,
	    minNeighbors=5
	)
	# Draw a rectangle around the faces and eyes
	for (x, y, w, h) in faces:

	    cv2.rectangle(frame,(x,y),(x+w,y+h),(0,255,0),2)
	    #crop the image, resize to 100*100
	    roi_gray = gray[y:y+h, x:x+w]
	    resized = cv2.resize(roi_gray,(100,100))
	    #save image in crop folder,required to create before running code
	    #cv2.imwrite('temp/crop'+str(number)+'.jpg',resized)
	    #number = number +1
	    cv2.waitKey(50)
	    cnt += 1

	# Display the resulting frame
	# cv2.imshow('Video', frame)

	if cnt%8==0:
		cnt = 0
		#calculate nearest eigenFace

		faceVecT = np.resize(resized,(1,10000))
		faceVec = np.resize(faceVecT - face_avg,(10000,1))

		RecEigenFace = np.dot(eigenVector,faceVec)
		RecEigenFaceTrans = RecEigenFace.transpose()
		#print RecEigenFaceTrans.shape

		minDiff = 1000000000000
		for i in range(0,imgNum):
			tempFace = eigenFaces[i]
			#print tempFace.shape

			diffArray = np.absolute(np.array(RecEigenFaceTrans) - np.array(tempFace))
			result = diffArray.astype(int)
			#diff = np.sum(result)
			square = np.square(result)
			SSDdiff = np.sum(square)
			#print diff`
			if SSDdiff<minDiff:
				minDiff = SSDdiff
				nearestFace = i

		if nearestFace in range(80, 84):
			people = 'Jinghuan'
		elif nearestFace in range(85, 89):
			people = 'Jinhai'
		elif nearestFace in range(90, 94):
			people = 'Wenhan'
		elif nearestFace in range(95, 99):
			people = 'Yongyao'
		else:
			people = 'other people'


		print("recognized as ", nearestFace)
		print("minimun difference is ",minDiff)
		#path = 'cap/'+str(nearestFace)+'.jpg'
		#img = cv2.imread(path)
		#cv2.imshow('recognized as: '+str(nearestFace)+'.jpg',img)

		cv2.putText(frame, str(nearestFace)+people, (x+5,y-5), font, 1, (255,255,255) , 2)
		cv2.putText(frame, str(minDiff), (x + 5, y+h - 5), font, 1, (255, 255, 255), 1)
	cv2.imshow('Video', frame)

	if cv2.waitKey(1) & 0xFF == ord('q'):
		break

# When everything is done, release the capture
video_capture.release()
cv2.destroyAllWindows()
