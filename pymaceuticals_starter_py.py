#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# * Your analysis here ...

# In[64]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
from scipy.stats import linregress
import numpy as np
# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
combined_data = pd.merge(study_results,mouse_metadata,how='left',on='Mouse ID')
# Display the data table for preview
combined_data.head()


# In[65]:


# Checking the number of mice.
combined_data["Mouse ID"].nunique()


# In[66]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
# Using the duplicated function to comb through our merged data using Mouse ID and Timepoint as identifiers
Duplicated_mice = combined_data[combined_data.duplicated(["Mouse ID","Timepoint"])]
Duplicated_mice


# In[67]:


# Optional: Get all the data for the duplicate mouse ID. 
Duplicated_mice
                              


# In[68]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
combined_data.drop(combined_data[combined_data["Mouse ID"]=="g989"].index, inplace=True)


# In[69]:


# Checking the number of mice in the clean DataFrame.
combined_data["Mouse ID"].nunique()


# ## Summary Statistics

# In[70]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
mean = combined_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].mean()
median = combined_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].median()
variance = combined_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].var()
standard_dv = combined_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].std()
sem = combined_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].sem()

summary_df = pd.DataFrame({"Mean Tumor Volume": mean, "Median Tumor Volume": median, "Tumor Volume Variance": variance, "Tumor Volume Standard Deviation": standard_dv, "Tumor Volume SEM": sem})
summary_df


# In[71]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line


# ## Bar and Pie Charts

# In[72]:


# Generate a bar plot showing the total number of unique mice tested on each drug regimen using pandas.
Drug_Regimen = combined_data.groupby("Drug Regimen")
Drug_Count = pd.DataFrame(Drug_Regimen["Drug Regimen"].count())
Bar_plot = Drug_Count.plot(kind='bar')
plt.title("Drugs tested on Mice")
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Unique Mice Tested")


# In[73]:


# Generate a bar plot showing the total number of unqiue mice tested on each drug regimen using using pyplot.
New_Count = combined_data["Drug Regimen"].value_counts()
x_axis=np.arange(len(New_Count))
tick_locations = [value for value in x_axis]
plt.xticks(tick_locations,New_Count.index.values,rotation=90)
plt.bar(x_axis,New_Count,color="blue",align="center",width=0.5)
plt.title("Drugs Tested on Mice")
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Unique Mice Tested")


# In[74]:


# Generate a pie plot showing the distribution of female versus male mice using pandas
# Creating a tracker for data distribution between male and female:
Gender_Tracker = combined_data["Sex"].value_counts()
# Creating the Pie Chart
Gender_Pie = Gender_Tracker.plot.pie(title="Distribution by Gender",autopct="%1.1f%%")


# In[75]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
Gender_Tracker = combined_data["Sex"].value_counts()
plt.pie(Gender_Tracker,labels=Gender_Tracker.index,autopct="%1.1f%%")
plt.title("Distribution by Gender")


# ## Quartiles, Outliers and Boxplots

# In[76]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
Last_Timepoint = combined_data.groupby("Mouse ID")["Timepoint"].max()
Last_Timepoint_df = pd.DataFrame(Last_Timepoint)
# Merge this group df with the original dataframe to get the tumor volume at the last timepoint
Merged_Timepoint = pd.merge(Last_Timepoint_df,combined_data,on=("Mouse ID","Timepoint"))
Merged_Timepoint.head()


# In[77]:


# Put treatments into a list for for loop (and later for plot labels)
treatment_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_list = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
for drug in treatment_list:
    Drug_df = combined_data.loc[combined_data["Drug Regimen"]==drug]
    Final_Timepoint_df = Drug_df.loc[Drug_df['Timepoint']==Drug_df['Timepoint']]
    Data=Final_Timepoint_df["Tumor Volume (mm3)"]
    tumor_vol_list.append(Data)
    Quartiles = Data.quantile([.25,.5,.75])
    UpperQuartile = Quartiles[0.75]
    LowerQuartile = Quartiles[0.25]
    IQR = UpperQuartile - LowerQuartile
    UpperBound = UpperQuartile + (1.5*IQR)
    LowerBound = LowerQuartile - (1.5*IQR)
    Outliers = (Data.loc[(Final_Timepoint_df["Tumor Volume (mm3)"] >= UpperBound) | (Final_Timepoint_df["Tumor Volume (mm3)"] <= LowerBound)]).count()
    print(f" IQR for {drug} is {IQR}")
    print(f" Upper bound for {drug} is {UpperBound}")
    print(f" Lower bound for {drug} is {LowerBound}")
    print(f" There are {Outliers} outliers for {drug}")
    

    
    


# In[78]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
plt.boxplot(tumor_vol_list,treatment_list)
plt.title("Final Tumor Volume by Drug Regimen")
plt.xticks([1,2,3,4], ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"])
plt.ylabel("Final Tumor Volume (mm3)")


# ## Line and Scatter Plots

# In[79]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin
# Using a random mouse treated with Capomulin
Mice = combined_data.loc[combined_data["Mouse ID"] == "m601"]
# Plotting with tumor volume vs Timepoint
plt.plot(Mice["Timepoint"],Mice["Tumor Volume (mm3)"], marker="o")
plt.title("Mouse M601 Treatment with Capomulin")
plt.xlabel("Time in Days")
plt.ylabel("Tumor Volume (mm3)")


# In[80]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
# Locating within the dataframe with Capomulin as our identifier
Capomulin_df = combined_data.loc[combined_data["Drug Regimen"]=="Capomulin"]
# Averaging the tumor volume of mice treated with Capomulin
Average_Tumor_df = pd.DataFrame(Capomulin_df.groupby("Mouse ID")["Tumor Volume (mm3)"].mean())
# Renaming old Tumor Volume column to become Average Tumor Volume
Renamed_df = Average_Tumor_df.rename(columns={"Tumor Volume (mm3)":"Average Tumor Volume"})
# Merging the two dataframes to get weight and average tumor volume all into one 
Final_Average_df = pd.merge(Capomulin_df,Renamed_df,on="Mouse ID")
Final_Average_df.head()
# Setting values for x and y axis
x=Final_Average_df["Weight (g)"]
y=Final_Average_df["Average Tumor Volume"]
plt.scatter(x,y)
# Labeling and the works
plt.title("Average Tumor Volume by Weight(g)")
plt.xlabel("Weight in Grams")
plt.ylabel("Average Tumor Volume in mm3")



# ## Correlation and Regression

# In[85]:


# Calculate the correlation coefficient and linear regression model 
Correlation = st.pearsonr(x,y)
print(f" The correlation between weight and average tumor volume for Capomulin is {round(Correlation[0],2)}.""")
# Calculating linear regression model
(slope, intercept, rvalue, pvalue, stderr) = linregress(x, y)
Values = x * slope + intercept
Line_Equation = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))
plt.scatter(x,y)
plt.plot(x,Values,color="red")
plt.annotate(Line_Equation,(20,40),color="red")
# for mouse weight and average tumor volume for the Capomulin regimen


# In[ ]:




