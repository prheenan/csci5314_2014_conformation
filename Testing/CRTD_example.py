import numpy as np
import matplotlib.pyplot as plt

# where the data is located relative to this file
baseFolder = "../output/5x_PBS_25C_/5x_PBS_25C_200ms01 ObjectTrackingData2/Post_Stage3_Analysis/"

CRTD =np.load(baseFolder + "GetModel_CRTD.npy")
timeForCRTD =np.load(baseFolder + "GetModel_Time_For_CRTD.npy")

plt.semilogy(timeForCRTD,CRTD,'ro-')
plt.title('CRTD for P(tau > t), from the .npy files')
plt.ylabel('P(tau > t)')
plt.xlabel("Time (s)")

plt.show()
