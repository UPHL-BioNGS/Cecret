#!/bin/python

import matplotlib.pyplot as plt
import pandas as pd
import os

os.mkdir('primerassessment')

df = pd.read_csv("primer_assessment.csv")
df.drop('sample', axis=1, inplace=True)

boxplot = df.boxplot(fontsize=5, 
                     rot=90, 
                     figsize=(15,8), 
                     grid=False)
boxplot.plot()
plt.title('Primer Assessment')
boxplot.set_ylabel('meandepth')
boxplot.set_xlabel('amplicon name')
boxplot.figure.savefig('primerassessment/primerassessment.png',     bbox_inches='tight')
boxplot.figure.savefig('primerassessment/primerassessment_mqc.png', bbox_inches='tight')

