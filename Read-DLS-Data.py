#! /usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
import glob, os

# The only thing that needs to be changed when running this script is the variable 'path_in'
# The script will go to that directory, read in all text files with the .txt extension
# It will then parse the files, and write out 3 csv files
# These files will then be read in by the script, and it will plot the graphs

path_in = 'C:/Users\isle132/Desktop/For Python/140 nm Silica Data/173 deg/Cell Center'

os.chdir(path_in)  # go to the path specified
# This will loop over all files with the .txt extension in the directory specified
# It will take the prefix for the file, and make a .Diameter.csv, a .correlation.csv and a .residual.csv
# These are extracted from the data in the .txt file
for file in glob.glob('*.txt'):
    exp_data_prefix = os.path.splitext(file)[0]
    exp_data_file = file

    file_in_name_path = path_in+'/'+exp_data_file
    file_diameter_out_name_path = path_in+'/'+exp_data_prefix+'.Diameter.csv'
    file_correlation_out_name_path = path_in+'/'+exp_data_prefix+'.correlation.csv'
    file_residual_out_name_path = path_in+'/'+exp_data_prefix+'.residual.csv'

    with open(file_in_name_path, 'r') as f_in, open(file_diameter_out_name_path, 'w') as f_dia_out, \
            open(file_residual_out_name_path, 'w') as f_res_out, open(file_correlation_out_name_path, 'w') as f_cor_out:
        copy_dia = False
        copy_res = False
        copy_cor = False
        for line in f_in:
            if line.strip() == 'Diameter (nm),Frequency (%),Undersize (%)':
                copy_dia = True
            elif line.strip() == '':
                copy_dia = False
                copy_res = False
                copy_cor = False
            elif line.strip() == 'Delay Time       (µs),Residual':
                copy_res = True
            elif line.strip() == 'Delay Time       (µs),Correlation,Fitting Function':
                copy_cor = True
            if copy_dia:
                f_dia_out.write(line)
            elif copy_res:
                f_res_out.write(line)
            elif copy_cor:
                f_cor_out.write(line)

# Import the data from the .Diameter.csv file
# Plot the diameter vs frequency
dia_files = glob.glob('*.Diameter.csv')

dia_dfs = [pd.read_csv(dia_fp).assign(Data=os.path.splitext(os.path.splitext(dia_fp)[0])[0]) for dia_fp in dia_files]
dia_df = pd.concat(dia_dfs, ignore_index=True)

dia_groups = dia_df.groupby('Data')

fig1, ax1 = plt.subplots()
ax1.margins(0.05)
ax1.set_xscale('log')
ax1.set_xlim(dia_df['Diameter (nm)'].min(), dia_df['Diameter (nm)'].max())
for name, group in dia_groups:
    ax1.plot(group['Diameter (nm)'], group['Frequency (%)'], label=name)
ax1.legend(loc='upper right')

# Import the data from the .residual.csv file
# Plot the Delay time vs Residual
res_files = glob.glob('*.residual.csv')

res_dfs = [pd.read_csv(res_fp, encoding="ISO-8859-1").assign(Data=os.path.splitext(os.path.splitext(res_fp)[0])[0])
           for res_fp in res_files]
res_df = pd.concat(res_dfs, ignore_index=True)

res_groups = res_df.groupby('Data')
fig2, ax2 = plt.subplots()
ax2.margins(0.05)
ax2.set_xscale('log')
ax2.set_xlim(res_df['Delay Time       (µs)'].min(), res_df['Delay Time       (µs)'].max())
for name, group in res_groups:
    ax2.plot(group['Delay Time       (µs)'], group['Residual'], label=name)
ax2.legend(loc='upper right')


# Import the data from the .correlation.csv file
# Plot the delay time vs correlation
cor_files = glob.glob('*.correlation.csv')

cor_dfs = [pd.read_csv(cor_fp, encoding="ISO-8859-1").assign(Data=os.path.splitext(os.path.splitext(cor_fp)[0])[0])
           for cor_fp in cor_files]
cor_df = pd.concat(cor_dfs, ignore_index=True)

cor_groups = cor_df.groupby('Data')
fig3, ax3 = plt.subplots()
ax3.margins(0.05)
ax3.set_xscale('log')
ax3.set_xlim(cor_df['Delay Time       (µs)'].min(), cor_df['Delay Time       (µs)'].max())
for name, group in cor_groups:
    ax3.plot(group['Delay Time       (µs)'], group['Correlation'], label=name)
ax3.legend(loc='upper right')

plt.show()


