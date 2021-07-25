import os
import pandas as pd
import numpy as np
from numpy import mean
import matplotlib
import matplotlib.pyplot as plt
import statistics as stats
from math import log
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from math import sqrt

# Directory for protein type (Ctau/NEFL)
protein_directory = r'C:\Users\HP\Documents\University\BMEN 501\Final Project\Data for ML\Data for ML\6. Mechanical load\6.1 Torsion\C-Tau'
protein_directory_list = os.listdir(protein_directory)  # list of folders in protein directory (degrees)
protein_directory_len = len(protein_directory_list)  # number of folders in protein folder (5)

avg_peak_list = []  # List of average currents
time_labels = []  # List of times
deg_list = []  # List of degrees

for i in range(0, protein_directory_len):  # outer loop to go through each folder corresponding to degrees
    deg_directory = os.path.join(protein_directory, protein_directory_list[i])  # This will bring us into a degree folder
    deg_directory_list = os.listdir(deg_directory)  # list of time folders
    deg_directory_length = len(deg_directory_list)  # Number of time folders

    for j in range(0, deg_directory_length):  # middle loop to go through each folder corresponding to time
        time_directory = os.path.join(deg_directory, deg_directory_list[j])  # This will bring us into a time folder
        time_directory_list = os.listdir(time_directory)  # list of excel files
        time_directory_length = len(time_directory_list)  # Number of excel files in each time folder
        os.chdir(time_directory)  # enter the directory for each time folder

        peak_current_list = []  # list for peak currents for all files in folder

        # Start file reading loop
        for k in range(0, time_directory_length):
            data = pd.read_excel(time_directory_list[k])  # Reading files in the folder. Make sure files are not open.
            df = pd.DataFrame(data, columns=['Applied potential (mV)', 'Current'])  # Type =dataframe
            df = df.dropna(axis=0,
                           how='any')  # Drop any rows (axis0 =0) that have any blanks. This removes the first blank row and any cells with << error
            df = df.astype(float)  # making sure Current column is a float type, not object
            df = df.nlargest(3, 'Current')  # Narrows dataframe down to rows of 3 largest current values
            peak_current_list.extend(df['Current'].tolist())

        avg_peak = mean(peak_current_list)  # should average N=3*#files
        avg_peak_list.append(avg_peak)  # each input is the average peak for that time
        time_labels.append(deg_directory_list[j])
        deg_list.append(protein_directory_list[i])  # Corresponding list of degrees

deg_labels = []  # corrected degree labels

# Start loop to obtain degrees from file names
for i in range(0, len(deg_list)):
    digit_list = [int(s) for s in deg_list[i].split("d") if
                  s.isdigit()]  # Extracts the digits from the file names in a list. Need to split based on string, so split with d.
    digit = digit_list[0]
    deg_labels.append(digit)  # Corrected labels append

# Final data lists
deg_labels2 = []
avg_peak_list2 = []
time_labels2 = []

for i in range(0, 7):  # Need outer loop 7 times for 7 different time labels
    for j in range(i, 35, 7):  # Every 7th element
        deg_labels2.append(deg_labels[j])  # Every 7th we get our next degree at the same time
        avg_peak_list2.append(avg_peak_list[j])
        time_labels2.append(time_labels[j])

# Making one data table. We have 3 rows, degrees, current, time. Each column is one experiment.
deg_labels2 = np.asarray(deg_labels2)
avg_peak_list2 = np.asarray(avg_peak_list2)
time_labels2 = np.asarray(time_labels2)
table = np.array([deg_labels2, avg_peak_list2, time_labels2]) # Consider using pandas dataframe

# Now we have our data as: deg_labels2=[120, 180, 30, 60, 90] for one time
n_times = len(np.unique(time_labels2))  # Number of unique times
n_intensities = len(np.unique(deg_labels2))  # Number of unique degrees

labels = np.empty([n_times * n_intensities])  # Each label (e.g. 2) classifies each column in the table as the same
# group of times. So label 1 is for all data points where time = 'After 120 min'

for t in range(0, n_times):
    labels[t * n_intensities : (t + 1) * n_intensities] = t + 1

# Plotting Data
plt.figure(figsize=(8, 8))
# Plotting each line
for i in range(0, n_times):
    plot_data = table[:, i * n_intensities: (i+1) * n_intensities]  # Each iteration of the loop plots one time
    # loading = row 1 of the data which are the loading intensities
    # value = row 2 of the data which are the current readings
    loading = plot_data[0, :].astype(float)  # Turn values into floats
    value = plot_data[1, :].astype(float)

    # Sorting the data in ascending loading intensities
    sortmask = np.argsort(loading)  # sortmask is the sorted order of the degrees
    loading = loading[sortmask]
    value = value[sortmask]
    plt.plot(loading, value, marker='o', label=plot_data[2, 0])

plt.title('Regression of Degrees of Torsion vs Average Peak Current')
plt.xlabel('Degrees Torsion')
plt.ylabel('Peak Current (μA)')
plt.legend()
plt.show()

# Now that we have the currents vs loading data, we need the regression equation that relates the current to protein
# concentration.

# Begin linear regression for C-Tau in Cell-Media

# Insert directory for protein detection files below
protein_directory = r'C:\Users\HP\Documents\University\BMEN 501\Final Project\Data for ML\Data for ML\4. Detections - Cell media\C-Tau'
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
    for j in range(0, int(0.5 * conc_directory_length)):  # The multiplying factor adjusts the number of files used
        # in regression. A factor of 1 means all files used.
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
    digit_list = [int(s) for s in concentration_labels[i].split() if s.isdigit()]  # Extracts the digits from folder
    # name
    digit = digit_list[0]
    if "pg" in concentration_labels[i]:  # Adjust concentrations to pg for base unit
        pass
    elif "ng" in concentration_labels[i]:
        digit *= 1000  # x1000 if ng
    elif "fg" in concentration_labels[i]:
        digit /= 1000  # divide by 1000 if fg
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
plt.title("Cell Media - C-Tau")
plt.xlabel("Log(10) Concentration (pg)")
plt.ylabel("Peak Current (I0-Ic)(μA)")
plt.xlim(-0.1, 5.2)

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
x_fit = np.linspace(-1, 6, 100)
y_fit = slope * x_fit + intercept
plt.plot(x_fit, y_fit, "r", label=LOBF)
plt.legend()
plt.show()

# R2 score
Y_act = Y
X_pred = np.array(X_log)
Y_pred = slope * X_pred + intercept
r2 = r2_score(Y_act, Y_pred)

# End of linear regression. Key variables:
# - slope (slope of regression line)
# - intercept (intercept of regression line)


# Plotting Predicted Concentrations. Same thing as plotting the peak currents, but use equation to turn current into
# concentration

# Lists of the predicted concentrations, loading, and times (for bar graph)
pre_conc_list = []
pre_load_list = []
pre_time_list = []

plt.figure(figsize=(8, 8))
# Plotting each line
for i in range(0, n_times):
    plot_data = table[:, i * n_intensities: (i+1) * n_intensities]
    loading = plot_data[0, :].astype(float)
    predicted_conc = plot_data[1, :].astype(float)
    predicted_conc = 10**(((BSA_avg_peak-predicted_conc)-intercept)/slope)  # Convert current to concentration using
    # regression equation

    sortmask = np.argsort(loading)
    loading = loading[sortmask]
    predicted_conc = predicted_conc[sortmask]

    pre_conc_list += list(predicted_conc)
    pre_time_list += list(plot_data[2,:])
    pre_load_list += list(loading)

    plt.plot(loading, predicted_conc, marker='o', label=plot_data[2, 0])
plt.title('Regression of Degrees of Torsion vs Predicted C-Tau Concentration')
plt.xlabel('Degrees Torsion')
plt.ylabel('C-Tau Concentration (pg)')
plt.legend()
plt.show()

# End of plot for predicted concentrations.

# Begin clustered bar chart of results.

# Recall our three lists: pre_conc_list, pre_time_list, and pre_load_list which store the ordered
# data and labels of our concentration predictions.

labels = ['30', '60', '90', '150', '300']
hrs_24 = pre_conc_list[30:35]
hrs_0 = pre_conc_list[25:30]
min_15 = pre_conc_list[5:10]
min_30 = pre_conc_list[10:15]
min_60 = pre_conc_list[15:20]
min_90 = pre_conc_list[20:25]
min_120 = pre_conc_list[0:5]

x = np.arange(len(labels))  # the label locations
width = 0.1  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - 3*width, hrs_24, width, label='before 24 hr')
rects2 = ax.bar(x - 2*width, hrs_0, width, label='0 hr')
rects3 = ax.bar(x - 1*width, min_15, width, label='15 min')
rects4 = ax.bar(x , min_30, width, label='30 min')
rects5 = ax.bar(x + 1*width, min_60, width, label='60 min')
rects6 = ax.bar(x + 2*width, min_90, width, label='90 min')
rects7 = ax.bar(x + 3*width, min_120, width, label='120 min')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('C-Tau Conc. (pg/mL)')
ax.set_xlabel('Applied Torsion (Deg)')
ax.axhline(linewidth = 1, color = 'black')
# ax.set_title('')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

fig.tight_layout()

plt.show()