#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json
import matplotlib.pyplot as plt
import numpy
import scipy
from scipy import interpolate
import pandas as pd

def plotXsec(crossSectionData, tag, lineStyles):
  massList = crossSectionData[tag][0]
  crossSectionList = crossSectionData[tag][1]
  crossSectionUncertaintyList = crossSectionData[tag][2]
  label = crossSectionData[tag][3]
  plt.yscale("log")
  style = 'solid'
  color = None
  if tag in lineStyles:
    style = lineStyles[tag][0]
    if len(lineStyles[tag])==2: color = lineStyles[tag][1]
  baseline = plt.plot(massList, crossSectionList, label = label, linestyle = style, color=color)
  # assume symmetric always present
  band = plt.fill_between(massList, numpy.asarray(crossSectionList)-numpy.asarray(crossSectionUncertaintyList), numpy.asarray(crossSectionList)+numpy.asarray(crossSectionUncertaintyList), alpha = 0.2, linewidth=0, facecolor = baseline[0].get_color())

def loadJsonFiles(jsonFilenames, crossSectionData):
  # load data
  for jsonFilename in jsonFilenames:
    # Get input data
    data = json.load(open(jsonFilename))
    df   = pd.DataFrame.from_dict(data["data"], orient = "index")
    # restore mass as column and sort 
    df["mass_GeV"] = df.index.astype(int)
    df = df.sort_values("mass_GeV")
    df.reset_index(inplace = True, drop = True)
    # Set output
    tag = data['process_id']
    crossSectionData[tag] = [[], [], []]
    for iPoint in range(len(df.mass_GeV)):
      mass = df.mass_GeV[iPoint]
      crossSection = df.xsec_pb[iPoint]
      crossSectionUncertainty = df.unc_pb[iPoint]
      crossSectionData[tag][0].append(mass)
      crossSectionData[tag][1].append(crossSection)
      crossSectionData[tag][2].append(crossSectionUncertainty)
    # Set label
    print(tag)
    label = 'no label'
    if tag == 'pp13_glsq': label = '$\\tilde g\\tilde q; m_{\\tilde g} = m_{\\tilde q (u,d,c,s)}$'
    elif tag == 'pp13_sqsq': label = '$\\tilde q\\tilde q^*; m_{\\tilde q} = m_{\\tilde q (u,d,c,s,b)} \\ll m_{\\tilde g, \\tilde t}$'
    elif tag == 'pp13_glgl': label = '$\\tilde g\\tilde g; m_{\\tilde g} \\ll m_{\\tilde q}$'
    elif tag == 'pp13_hino': label = '$\\tilde\\chi\\tilde\\chi; m_{\\tilde {\\chi}^\\pm}=m_{\\tilde {\\chi}^0}$'
    crossSectionData[tag].append(label)

def loadGraphFiles(graphFilenames, crossSectionData):
  for graphFilename in graphFilenames:
    with open(graphFilename) as graphFile:
      tag = graphFile.readline().split(':')[1].rstrip().lstrip()
      label = graphFile.readline().split(':')[1].rstrip().lstrip()
      print(label)
      unitLength = float(graphFile.readline().split(':')[1])
      #print(tag, unitLength)
      massList = []
      crossSectionList = []
      crossSectionUncertaintyList = []
      for line in graphFile:
        mass = int(line.split(':')[0])
        length = float(line.split(':')[1])
        if unitLength == -1:
          crossSection = length
        else:
          crossSection = 10**(length/unitLength)
        massList.append(mass)
        crossSectionList.append(crossSection)
        crossSectionUncertaintyList.append(0)
      crossSectionData[tag] = [massList, crossSectionList, crossSectionUncertaintyList, label]
        #print(mass, length, crossSection)

def combineCrossSections(target_tag, target_label, combine_tags, crossSectionData):
  # Find common mass points
  massList = None
  for tag in combine_tags:
    if massList == None:
      massList = crossSectionData[tag][0]
    else:
      massList = list(set(crossSectionData[tag][0]).intersection(massList))
  massList = sorted(massList)

  # Add common mass points
  crossSectionList = []
  crossSectionUncertaintyList = []
  for target_mass in massList:
    crossSection = 0
    crossSectionUncertainty = 0
    for tag in combine_tags:
      iPoint = crossSectionData[tag][0].index(target_mass)
      crossSection += crossSectionData[tag][1][iPoint]
    crossSectionList.append(crossSection)
    crossSectionUncertaintyList.append(crossSectionUncertainty)

  crossSectionData[target_tag] = [massList, crossSectionList, crossSectionUncertaintyList, target_label]


if __name__ == "__main__":
  # crossSectionData[tag] = [[mass [GeV]], [cross_section [pb]], [cross_section_uncertainty [pb]], label]
  crossSectionData = {}

  jsonFilenames = [
    'json/pp13_gluinosquark_NNLO+NNLL.json', # pp13_glsq: ~g~q; ~g=~q(~u,~d,~c,~s), decoupled:?
    'json/pp13_squark_NNLO+NNLL.json', # pp13_sqsq: ~q~q*; ~q(~u, ~d, ~c, ~s, ~b), decoupled: ~g, ~t
    'json/pp13_gluino_NNLO+NNLL.json', # pp13_glgl: ~g~g; decoupled:~q
    'json/pp13_hino_NLO+NLL.json', #pp13_hino: ~h~h; degenerate(N,C)
  ]
  loadJsonFiles(jsonFilenames, crossSectionData)

  graphFilenames = [
    #'graph/gluinosquark_gleqsq.txt', # pp13_glsq_graph: ~g~q; ~g=~q(~u,~d,~c,~s), decoupled:?
    'graph/squarksquark_gleqsq.txt', # pp13_sqsq_graph: ~g~q; ~g=~q(~u,~d,~c,~s), decoupled:?
    'graph/squarkantisquark_gleqsq.txt', # pp13_sqantisq_graph: ~g~q; ~g=~q(~u,~d,~c,~s), decoupled:?
    'graph/gluinogluino_gleqsq.txt', # pp13_glgl_graph: ~g~q; ~g=~q(~u,~d,~c,~s), decoupled:?
    'graph/squarksquark_nmssm.txt', # pp13_sqsq_nmssm: ~q~q; ~g=1.01*~q(~u,~d,~c,~s) ~~ ~q(~t,~b)
    'graph/mix_nmssm.txt', # pp13_mix_nmssm: ~g~q,~q~q,~g~g; ~g=1.01*~q(~u,~d,~c,~s) ~~ ~q(~t,~b)
  ]
  loadGraphFiles(graphFilenames, crossSectionData)

  combineCrossSections('pp13_mix', '$\\tilde g\\tilde q, \\tilde q\\tilde q, \\tilde g\\tilde g; m_{\\tilde g} = m_{\\tilde q (u,d,c,s)} \\ll m_{\\tilde q (b,t)}$', ['pp13_glsq', 'pp13_sqsq_graph', 'pp13_glgl_graph'], crossSectionData)

  lineStyles = {'pp13_mix': ['dashed', 'Red'], 'pp13_glsq': ['dashed', 'DarkOrange'], 'pp13_sqsq_graph': ['dashed', 'Black'], 'pp13_sqantisq_graph': ['dashed','DarkGreen'], 'pp13_glgl_graph': ['dashed', 'Blue'],
                'pp13_mix_nmssm': ['dashdot', 'Red'], 'pp13_sqsq_nmssm': ['dashdot', 'Black'],
                'pp13_sqsq': ['dotted', 'DarkGreen'],
                'pp13_glgl': ['solid', 'Blue'],
                'pp13_hino': ['dashdot', 'Brown']
               }

  # Plot all data
  plotSequence = [
                  #'pp13_mix', 'pp13_glsq', 'pp13_sqsq_graph', 'pp13_sqantisq_graph', 'pp13_glgl_graph', 
                  'pp13_mix', 'pp13_sqsq_graph', 
                  'pp13_mix_nmssm', 'pp13_sqsq_nmssm',
                  #'pp13_mix',  
                  #'pp13_sqsq', 
                  'pp13_glgl', 
                  #'pp13_hino'
                  ]
  for tag in plotSequence:
    plotXsec(crossSectionData, tag, lineStyles)

  #plt.rc('text', usetex=True)
  #plt.rc('font', size=18)
  plt.rc('legend', fontsize=7.5)
  #plt.rc('text.latex', preamble=r'\usepackage{cmbright}')
  # draw legend and style plot
  plt.xlabel("particle mass [GeV]")
  plt.ylabel("cross section [pb]")
  plt.grid()
  plt.xlim(100, 2500)
  plt.ylim(1e-6, 1e4)
  plt.legend(ncol = 2, framealpha = 1)
  plt.locator_params(axis = "y", base = 100) # for log-scaled axis, it's LogLocator, not MaxNLocator
  plt.savefig('test.pdf')
