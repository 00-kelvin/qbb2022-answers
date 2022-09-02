#!/usr/bin/env python3

# import all necessary packages
import numpy as np 
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm

# create numpy array from formatted data
df = np.genfromtxt("wrangled.txt", dtype=int, encoding=None, names=True)

# plot mutations vs. parent age for maternally and paternally derived mutations
fig, ax = plt.subplots()
ax.scatter(df["Mother_age"], df["mother_count"])
ax.set_ylabel("Count of maternal de novo mutations")
ax.set_xlabel("Maternal age")
plt.savefig("ex2_a.png")

fig, ax = plt.subplots()
ax.scatter(df["Father_age"], df["father_count"])
ax.set_ylabel("Count of paternal de novo mutations")
ax.set_xlabel("Paternal age")
plt.savefig("ex2_b.png")

# test for association between maternal age and maternally inherited mutations
model_mother = smf.ols(formula="mother_count ~ 1 + Mother_age", data=df)
results_mother = model_mother.fit()
# print(results_mother.summary())
# print(results_mother.pvalues)

# test for association between paternal age and paternally inherited mutations
model_father = smf.ols(formula="father_count ~ 1 + Father_age", data=df)
results_father = model_father.fit()
# print(results_father.summary())
# print(results_father.pvalues)

# plot histogram of maternally and paternally derived mutations
fig, ax = plt.subplots()
ax.hist(df["mother_count"], alpha=0.5, label="Maternally inherited")
ax.hist(df["father_count"], alpha=0.5, label="Paternally inherited")
ax.set_xlabel("Number of de novo mutations per proband")
ax.set_ylabel("Frequency")
ax.legend()
plt.savefig("ex2_c.png")

# test for signficant difference between maternally/paternally derived mutations
ttest_result = stats.ttest_ind(df["father_count"], df["mother_count"])
#print(ttest_result)

# predict mutations for a sample proband with father age = 50.5
new_data = df[0]
new_data.fill(0)
new_data["Father_age"] = 50.5
print(model_father.fit().predict(new_data))
