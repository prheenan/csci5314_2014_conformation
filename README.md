csci5314_2014_conformation
==========================

CU Boulder // CSCI 5314 Group 2 Project on Protein Conformation Modeling using FRET

TO RUN
=========
(1) Ensure there is at least one '.dat' file in the directory "csci5314_2014_conformation/Data/ObjTr"
    (a) The data can be downloaded from the data section of our google drive
(2) Go into the 'Stage0_Run' directory
(3) Type 'python Main.py'
    (a) You should see something like:
    	Message [Reading: file [../Data/ObjTr/5x_PBS_25C_200ms01001 ObjectTrackingData2.dat]] Received by [Step1::GetListOfObjectData]
	Message [Matching: file [../Data/ObjTr/5x_PBS_25C_200ms01001 ObjectTrackingData2.dat]] Received by [Step1::GetListOfObjectData]
	Message [Processing: [25583] points for [4]] Received by [Step1::GetListOfObjectData]
    (b) If you have problems, you may want to try ensuring you have numpy, matplotlib,  and scipy installed (http://www.scipy.org/)
(4) You should see output information (plots, raw data, etc) in the directory "csci5314_2014_conformation/output/"