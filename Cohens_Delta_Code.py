
"""
Code written by Lloyd A. Courtenay 
Lloyd A. Courtenay - ladc1995@gmail.com (Universidad de Salamanca [USAL]) #

"""

# load libraries and dependencies

import numpy as np
import pandas as pd
import arviz as az
import pymc3 as pm
import seaborn as sns

print("Package Versions")
print("---------------")
print("Numpy: {}".format(np.__version__))
print("Pandas: {}".format(pd.__version__))
print("ArviZ: {}".format(az.__version__))
print("PyMC3: {}".format(pm.__version__))
print("Seaborn: {}".format(sns.__version__))

az.style.use("arviz-darkgrid")

# load pc scores

data = pd.read_csv("pc_scores.csv", sep = ",")
pc1 = data["PC1"].values
sample = pd.Categorical(data["Labels"],
                       categories = ["O", "C", "sT", "T"]).codes
sample_groups = len(np.unique(sample))
data.head()

# visualise distributions

sns.violinplot(x = "Labels", y = "PC1", data = data)

# fit bayesian model with Gaussian Priors

with pm.Model() as comparing_groups_gaussian:
    μ = pm.Normal("μ", mu = 0, sd = 10, shape = sample_groups)
    σ = pm.HalfNormal("σ", sd = 10, shape = sample_groups)
    y = pm.Normal("y", mu = μ[sample], sd = σ[sample], observed = pc1)
    trace = pm.sample(5000, tune = 1000)
    
# evaluate fit

az.plot_trace(trace)

az.summary(trace)

az.waic(trace)

# sample from posterior

pred = pm.sample_posterior_predictive(trace, 100, comparing_groups_gaussian)

new_data = az.from_pymc3(trace = trace, posterior_predictive = pred)

# calculate Cohen's Delta

means_diff = trace['μ'][:, 0] - trace['μ'][:, 1]
d_cohen_01 = (means_diff / np.sqrt((trace['σ'][:, 0]**2 + trace['σ'][:, 1]**2) / 2)).mean()
means_diff = trace['μ'][:, 0] - trace['μ'][:, 2]
d_cohen_02 = (means_diff / np.sqrt((trace['σ'][:, 0]**2 + trace['σ'][:, 2]**2) / 2)).mean()
means_diff = trace['μ'][:, 0] - trace['μ'][:, 3]
d_cohen_03 = (means_diff / np.sqrt((trace['σ'][:, 0]**2 + trace['σ'][:, 3]**2) / 2)).mean()
means_diff = trace['μ'][:, 1] - trace['μ'][:, 2]
d_cohen_12 = (means_diff / np.sqrt((trace['σ'][:, 1]**2 + trace['σ'][:, 2]**2) / 2)).mean()
means_diff = trace['μ'][:, 1] - trace['μ'][:, 3]
d_cohen_13 = (means_diff / np.sqrt((trace['σ'][:, 1]**2 + trace['σ'][:, 3]**2) / 2)).mean()
means_diff = trace['μ'][:, 2] - trace['μ'][:, 3]
d_cohen_23 = (means_diff / np.sqrt((trace['σ'][:, 2]**2 + trace['σ'][:, 3]**2) / 2)).mean()

print("Cohen's δ")
print(f"S1 vs S2 = {d_cohen_01:.2f}")
print(f"S1 vs S3 = {d_cohen_02:.2f}")
print(f"S1 vs S4 = {d_cohen_03:.2f}")
print(f"S2 vs S3 = {d_cohen_12:.2f}")
print(f"S2 vs S4 = {d_cohen_13:.2f}")
print(f"S3 vs S4 = {d_cohen_23:.2f}")

#

# alternative code ------------------------------------------------------------------

# in case of t-distributed priors the following bayesian model was used
# cohen's delta calculations should also be adjusted to use robust measurements of central tendency
# as well as deviation.

with pm.Model() as comparing_groups_robust:
    μ = pm.Normal("μ", mu = 0, sd = 10, shape = sample_groups)
    σ = pm.HalfNormal("σ", sd = 10, shape = sample_groups)
    ν = pm.Exponential("ν", 1/30, shape = sample_groups)
    y = pm.StudentT("y", mu = μ[sample], sd = σ[sample], nu = ν[sample], observed = pc1)
    trace_robust = pm.sample(10000, tune = 4000, target_accept = 0.95)
    
#