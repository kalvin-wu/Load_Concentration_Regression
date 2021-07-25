# Importing modules
import os
import pandas as pd
import numpy as np
from numpy import mean
import matplotlib.pyplot as plt
import statistics as stats
import re
from math import log
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from math import sqrt

# Insert directory for protein detection files below
protein_directory = r'C:\Users\HP\Documents\University\BMEN 501\Final Project\Data for ML\Data for ML\2. Range and caliberation curve\C-Tau'
protein_directory_list = os.listdir(protein_directory)  # List of folders/files in directory

avg_peak_list = []  # List of average currents from each folder
std_err_of_mean = []  # List of standard errors of mean from each folder

# Begin process of reading files from each folder
for i in range(0, len(protein_directory_list)):  # Outer loop to go through each folder of different concentrations
    conc_directory = os.path.join(protein_directory, protein_directory_list[i])  # Entering the folder for one
    # concentration
    conc_directory_list = os.listdir(conc_directory)  # List of files in each concentration folder - these are Excel
    # files
    conc_directory_length = len(conc_directory_list)  # Number of Excel files in each concentration folder
    os.chdir(conc_directory)  # Change program's directory to the concentration folder
    peak_current_list = []  # Empty list for peak currents for all files in the concentration folder

    # Begin reading Excel files in each concentration folder
    for j in range(0, int(1 * conc_directory_length)):  # The multiplying factor adjusts the number of files used
        # in the model. A factor of 1 means all files used.
        data = pd.read_excel(conc_directory_list[j])  # Reading files in the folder. Make sure files are not open.
        df = pd.DataFrame(data, columns=['Applied potential (mV)', 'Current'])  # Naming columns of Excel file
        df = df.dropna(axis=0, how='any')  # Drop any rows (axis=0) that have any blanks.
        df = df.astype(float)  # making sure 'Current' column is a float type
        df = df.nlargest(3, 'Current')  # Narrows dataframe down to rows of 3 largest current values
        peak_current_list.extend(df['Current'].tolist())  # Adding 3 largest current values to list of peak currents
        # for that concentration

    avg_peak = mean(peak_current_list)  # Average of peak currents for ONE concentration
    avg_peak_list.append(avg_peak)  # List of averaged peak currents
    std_err = stats.stdev(peak_current_list) / sqrt(len(peak_current_list))  # Standard error for one concentration
    std_err_of_mean.append(std_err)  # List of standard errors

concentration_labels = protein_directory_list  # Concentration labels taken from folder names
concentrations = []  # Corrected concentrations list

# Start loop to clean concentration labels from folder names
for i in range(0, len(concentration_labels)):
    digit_list = [int(s) for s in re.split('-| ', concentration_labels[i]) if s.isdigit()]  # Extracts the digits from folder
    # name. Make sure the folder names have spaces between digits and letters since is digit only detected digits alone.
    digit = digit_list[0]
    if "pg" in concentration_labels[i]:  # Adjust concentrations to pg for base unit
        pass
    elif "ng" in concentration_labels[i]:
        digit *= 1000  # x1000 if ng
    elif "fg" in concentration_labels[i]:
        digit /= 1000  # divide by 1000 if fg
    elif "ug" in concentration_labels[i]:
        digit *= 1000000  # x1000000 if microgram
    concentrations.append(digit)

# Start obtaining BSA Peak Current

BSA_directory = r'C:\Users\HP\Documents\University\BMEN 501\Final Project\Data for ML\Data for ML\5. Selectivity measurement\Selectivity in cell media\Cell media - Selectivity - C-Tau\1.BSA'
os.chdir(BSA_directory)  # Changing directory to BSA folder
BSA_directory_list = os.listdir(BSA_directory)
BSA_directory_length = len(BSA_directory_list)  # List of BSA files
BSA_current_list = []

# Start reading files
for i in range(0, BSA_directory_length):
    data = pd.read_excel(BSA_directory_list[i])
    BSA_df = pd.DataFrame(data, columns=['Applied potential (mV)', 'Current'])
    BSA_df = BSA_df.dropna(axis=0, how='any')
    BSA_df = BSA_df.astype(float)
    BSA_df = BSA_df.nlargest(3, 'Current')
    BSA_current_list.extend(BSA_df['Current'].tolist())

BSA_avg_peak = mean(BSA_current_list)  # BSA value

# Plotting Figure
BSA_diff_peak = BSA_avg_peak - avg_peak_list  # Take difference between BSA and peak currents
x = concentrations
X_log = [log(i, 10) for i in x]  # Turning concentrations to log scale
Y = BSA_diff_peak
plt.plot(X_log, Y, 'bo')
plt.errorbar(X_log, Y, yerr=std_err_of_mean, fmt='none', ecolor='black', capsize=3)
plt.title("Range and Calibration Curve - C-Tau")
plt.xlabel("Log(10) Concentration (pg)")
plt.ylabel("Peak Current (I0-Ic)(μA)")
#plt.xlim(-0.1, 5.2)

# Linear Regression
X = np.zeros([len(X_log), 1])  # Turn list into arrays for linear regression
for i in range(0, len(X_log)):
    X[i] = X_log[i]
Y = np.array(Y)
model = LinearRegression()
model.fit(X, Y)  # Fitting model with X and Y
slope = model.coef_
slope = slope[0]
intercept = model.intercept_
print("The equation is y=" + str(slope) + "x+" + str(intercept))
LOBF = "y=" + str(round(slope, 2)) + "x+" + str(round(intercept, 2))

# Plot line of best fit
x_fit = np.linspace(min(X_log), max(X_log), 100)
y_fit = slope * x_fit + intercept
plt.plot(x_fit, y_fit, "r", label=LOBF)
plt.legend()
plt.show()

# R2 score
Y_act = Y
X_pred = np.array(X_log)
Y_pred = slope * X_pred + intercept
r2 = r2_score(Y_act, Y_pred)
print('r2=', r2)



