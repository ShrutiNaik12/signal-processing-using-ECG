Utilizing the MIT-BIH Arrhythmia Database that is available on www.physionet.org, selected the first10 minutes of an ECG record as described below. Write a computer program to perform the following processing on your selected record. Wrote a computer program (using Matlab) to perform the following processing on your selected record.

1.	Document the signal record that you have used by giving the record number. Plot the ECG data in the record. 
2.	Write a program to determine how many heartbeats are present in the data.
3.	Compute the duration of the time interval between successive heart beats, called RR interval. 
4.	Plot the duration of the intervals obtained in step 3 above as a function of time. 
5.	Obtain and plot the autocorrelation for the time intervals obtained in step 3 above. Allow up to 100 samples shift in either direction. 
6.	Plot the histogram of the RR intervals obtained in step 3.
7.	Plot power density spectrum for the time intervals obtained in step 3 above. 
8.	It is desired to examine the impact of removing high frequency components present in the ECG data on the computation of RR intervals. Design an FIR filter to remove the components with frequencies exceeding 50 Hz. 
9.	Using the filtered signal in part 7, repeat steps 1 and 3 above.  
10.	Compute and plot cross correlation of filtered ECG signal obtained in step 7 with the downloaded signal.
11.	Compare the results in step 8 with the results in steps 1 and 3 above by plotting both ECG raw and filtered on the same plot and, separately,  plot the RR intervals obtain in step 3 and those obtained in step 8 on the same set of axes . 
