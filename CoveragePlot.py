import matplotlib.pyplot as plt
import pandas as pd
import glob
import argparse
from matplotlib.offsetbox import AnchoredText
import numpy as np

'''
# How to generate coverage file for 5kb Window using sambamba tool.
sambamba depth window -w 20000 HollysInput.bam -o HollysOutput.txt
'''


def main(args):

  Title = args.Sample #What you want to call the file
  CoverageDF = pd.read_csv(args.CoverageFile,sep='\t') # Path to input file
  CoverageDF['# chrom'] = [str(x) for x in CoverageDF['# chrom'].tolist()]

  # Coverage per gene.
  if len(args.GeneGuide) != 0:

	  Guide = pd.read_csv(args.GeneGuide,
	  					  sep='\t',
						  header=None)

	  UniqueIndexes = Guide.index.tolist()

	  for GIndex in UniqueIndexes:

		  fig = plt.figure(figsize=(20, 10))


		  SubCov  = CoverageDF[(CoverageDF['# chrom'] == str(Guide.loc[GIndex,0]))&
		  						(CoverageDF['chromStart'] >= Guide.loc[GIndex,1])&
								(CoverageDF['chromEnd'] <= Guide.loc[GIndex,2])]

		  x = SubCov.loc[:,'chromEnd']
		  y = SubCov.loc[:,'meanCoverage']
		  plt.title(Guide.loc[GIndex,3])
		  plt.plot(x, y)
		  plt.axhline(y=5, color='red', linestyle='dotted')
		  plt.ylim(0,max(y)*1.1)
		  plt.tight_layout()
		  plt.savefig('{}.png'.format(Guide.loc[GIndex,3]))
		  plt.clf()



  else:
	  UniqueChromosome = CoverageDF['# chrom'].unique().tolist()
	  fig, axs = plt.subplots(len(UniqueChromosome))
	  fig.set_figheight(50)
	  fig.set_figwidth(15)
	  for ChromoIndex in range(len(UniqueChromosome)):
		  x = CoverageDF[CoverageDF['# chrom']==UniqueChromosome[ChromoIndex]].loc[:,'chromEnd']
		  y = CoverageDF[CoverageDF['# chrom']==UniqueChromosome[ChromoIndex]].loc[:,'meanCoverage']
		  axs[ChromoIndex].set_title(UniqueChromosome[ChromoIndex])
		  axs[ChromoIndex].plot(x, y)
		  axs[ChromoIndex].axhline(y=5, color='red', linestyle='dotted')
		  axs[ChromoIndex].set_ylim(0,max(y)*1.1)


  plt.subplots_adjust(hspace = 0.8)
  plt.tight_layout()
  plt.savefig('{}.png'.format(Title))
  plt.clf()




#############################################################
parser = argparse.ArgumentParser()
parser.add_argument('--Sample', type=str, help = 'title of image', required=True)
parser.add_argument('--CoverageFile', type=str, help = 'Output from sambamba e.g. sambamba depth window -w 20000 HollysInput.bam -o HollysOutput.txt', required=True)
parser.add_argument('--GeneGuide', type=str, help = 'Bed file for genes of interest',default='')
parser.set_defaults(func=main)
args = parser.parse_args()
main(args)
