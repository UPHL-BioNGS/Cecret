#!/bin/python

import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("primer_assessment.csv")
df.drop('sample', axis=1, inplace=True)

png = plt.figure()
plt.title('Primer Assessment')
boxplot = df.boxplot(fontsize=5, rot=90, figsize=(15,8), patch_artist=True)
boxplot.set_ylabel('meandepth')
boxplot.set_xlabel('primer name')
png.savefig('primerassessment/primerassessment.png',     format = 'png')
png.savefig('primerassessment/primerassessment_mqc.png', format = 'png')
