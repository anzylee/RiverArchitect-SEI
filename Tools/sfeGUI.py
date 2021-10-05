# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 11:25:32 2021

@author: MM
"""
import os
import pandas as pd
import numpy as np
import tkinter as tk
from datetime import timedelta

## SET INITIAL VALUES FOR ALL BUTTONS AND LISTS ##
def nextButton(parent, options, selection, text, textRow):
    if selection in options:
        parent.destroy()
    else:
        tk.Label(parent,text=text).grid(row=textRow,column=0)

def chooseModel():
    root = tk.Tk()
    root.geometry("300x200")

    tk.Label(root,text="Model Type").grid(row=0,column=0,sticky="NSEW")
    
    modelType = tk.StringVar()
    tk.Radiobutton(root, text="2-D Hydrodynamic", variable=modelType, value='2-d hydrodynamic').grid(sticky="NSEW", row=1, column=0, padx=5)
    tk.Radiobutton(root, text="RCT Discharge", variable=modelType, value='rct discharge').grid(sticky="NSEW", row=2, column=0, padx=5)
    tk.Radiobutton(root, text="RCT Probe", variable=modelType, value='rct probe').grid(sticky="NSEW", row=3, column=0, padx=5)
    
    tk.Button(root, text="Next", command=lambda:nextButton(root, ['2-d hydrodynamic','rct discharge','rct probe'], modelType.get(), 'Choose a model type', 5)).grid(row=4,column=0)
    
    for r in np.arange(5):
        root.grid_rowconfigure(r, weight=1)
        
    root.grid_columnconfigure(0, weight=1)
    
    root.mainloop()
    
    return modelType.get()

def chooseRefSite():
    root = tk.Tk()
    root.geometry("300x200")

    tk.Label(root,text="Streamflow Reference Site").grid(row=0,column=0,sticky="E")
    
    # Import reference site IDs
    dfRef = pd.read_excel('input/All SFE LOI Characteristics with MAF.xlsx', index_col=0, header=0)
    modeledRefs = dfRef[dfRef['Modeled']=='Y']
    refSiteOptions = {}
    for i, r in modeledRefs.iterrows():
        refSiteOptions['Watershed: '+str(i)+', Geomorphic Class: '+r['Geomorphic Class']+', Geologic Setting: '+r['Geologic Setting']+', Mean Annual Flow (cfs): '+'{:0.1f}'.format(r['Mean Annual Flow (cfs)'])] = i
    
    refSite = tk.StringVar()
    tk.OptionMenu(root, refSite, *refSiteOptions.keys()).grid(sticky="NSEW", row=1, column=0, padx=5)
    
    tk.Button(root, text="Next", command=lambda:nextButton(root, refSiteOptions.keys(), refSite.get(), 'Choose a reference streamflow', 3)).grid(row=2,column=0)
    
    for r in np.arange(3):
        root.grid_rowconfigure(r, weight=1)
        
    root.grid_columnconfigure(0, weight=1)
    
    root.mainloop()
    
    return refSiteOptions[refSite.get()]

def chooseCaseSite():
    root = tk.Tk()
    root.geometry("300x200")

    tk.Label(root,text="Topographic Case Site").grid(row=0,column=0,sticky="E")
    
    # Import case site database containing the corresponding LOI for each case site with available data
    dfCase = pd.read_excel('input/case_site_database.xlsx', index_col=0, header=0)
    caseSiteOptions = {}
    for i, r in dfCase.iterrows():
        caseSiteOptions['Site ID: '+str(i)+', Watershed: '+str(r['LOI'])+', Geomorphic Class: '+r['Geomorphic Class']+', Geologic Setting: '+r['Geologic Setting']+', Mean Annual Flow (cfs): '+'{:0.1f}'.format(r['Mean Annual Flow (cfs)'])] = i
    
    caseLOI = tk.StringVar()
    tk.OptionMenu(root, caseLOI, *caseSiteOptions.keys()).grid(sticky="NSEW", row=1, column=0, padx=5)
    
    tk.Button(root, text="Next", command=lambda:nextButton(root, caseSiteOptions.keys(), caseLOI.get(), 'Choose a topographic case site', 3)).grid(row=2,column=0)
    
    for r in np.arange(3):
        root.grid_rowconfigure(r, weight=1)
        
    root.grid_columnconfigure(0, weight=1)
    
    root.mainloop()
    
    return caseSiteOptions[caseLOI.get()], dfCase.loc[caseSiteOptions[caseLOI.get()]]['LOI']

def chooseRCTParadigm():
    root = tk.Tk()
    root.geometry("400x200")

    tk.Label(root,text="RCT Paradigm Model Watershed Site").grid(row=0,column=0,sticky="E")
    
    # Import reference site IDs
    dfRef = pd.read_excel('input/All SFE LOI Characteristics with MAF.xlsx', index_col=0, header=0)
    modeledRefs = dfRef[dfRef['Paradigm']=='Y']
    paradigmOptions = {}
    for i, r in modeledRefs.iterrows():
        paradigmOptions['Watershed: '+str(i)+', Geomorphic Class: '+r['Geomorphic Class']+', Geologic Setting: '+r['Geologic Setting']+', Mean Annual Flow (cfs): '+'{:0.1f}'.format(r['Mean Annual Flow (cfs)'])] = i
    
    paradigmLOI = tk.StringVar()
    tk.OptionMenu(root, paradigmLOI, *paradigmOptions.keys()).grid(sticky="NSEW", row=1, column=0, padx=5)
    
    tk.Button(root, text="Next", command=lambda:nextButton(root, paradigmOptions.keys(), paradigmLOI.get(), 'Choose an RCT paradigm model watershed site', 3)).grid(row=2,column=0)
    
    for r in np.arange(3):
        root.grid_rowconfigure(r, weight=1)
        
    root.grid_columnconfigure(0, weight=1)
    
    root.mainloop()
    
    return paradigmOptions[paradigmLOI.get()]


## Under Construction
def setThresholds(model):
    
    if model.modelType == '2-d hydrodynamic':
        life_stage_period_default = {}
        life_stage_period_default['Chinook Salmon'] = {}
        life_stage_period_default['Chinook Salmon']['fry'] = [[1,1],[9,30], 0.1, np.nan, True, True]
        life_stage_period_default['Chinook Salmon']['juvenile'] = [[1,1],[9,30], 0.1, np.nan, True, True]
        life_stage_period_default['Chinook Salmon']['spawn'] = [[1,1],[9,30], 0.1, np.nan, True, True]
        life_stage_period_default['Rainbow / Steelhead Trout'] = {}
        life_stage_period_default['Rainbow / Steelhead Trout']['fry'] = [[1,1],[9,30], 0.1, np.nan, True, True]
        life_stage_period_default['Rainbow / Steelhead Trout']['juvenile'] = [[6,16],[9,30], 0.1, np.nan, True, True]
        life_stage_period_default['Rainbow / Steelhead Trout']['spawn'] = [[1,1],[9,30], 0.1, np.nan, True, True]

        root = tk.Tk()
        root.geometry("300x200")
        
        tk.Label(root,text="Habitat Ecorisk Threshold").grid(row=0,column=0)
        
        thresholdType = tk.StringVar()
        tk.Radiobutton(root, text="Default", variable=thresholdType, value='default').grid(sticky="NSEW", row=1, column=0)
        tk.Radiobutton(root, text="Custom Entry", variable=thresholdType, value='custom').grid(sticky="NSEW", row=2, column=0)
        
        tk.Button(root, text="Next", command=lambda:nextButton(root, ['default','custom'], thresholdType.get(), 'Choose a threshold type', 4)).grid(row=3,column=0)
        
        for r in np.arange(4):
            root.grid_rowconfigure(r, weight=1)
            
        root.grid_columnconfigure(0, weight=1)
        
        root.mainloop()
        
        if thresholdType.get() == 'default':
            disabled = True
        else:
            disabled = False
            
        life_stage_period = collectLifeStagePeriod(life_stage_period_default, disabled=disabled)
        
    elif (model.modelType == 'rct probe') or (model.modelType == 'rct discharge'):
        
        dfParadigm = pd.read_excel('../Tools/input/RCT_Paradigm_Selection.xlsx', index_col=None, header=0)
        paradigm = dfParadigm[dfParadigm['LOI'] == model.rct_site]
        
        if model.modelType == 'rct discharge':
            t_name = 'Discharge\nThreshold {0:d} (m^3/s)'
            t1 = 'Threshold 1 (cms)'
            t2 = 'Threshold 2 (cms)'
        else:
            t_name = 'RCT Probe\nThreshold {0:d} (m)'
            t1 = 'RCT Probe Threshold 1 (m)'
            t2 = 'RCT Probe Threshold 2 (m)'
        
        life_stage_period_default = {}
        for i, r in paradigm.iterrows():
            try:
                life_stage_period_default[r['Fish name']][r['Life stage']] = [[r['Start Month'],r['Start Day']],[r['End Month'],r['End Day']], r[t1], r[t2], bool(r['Consecutive']), bool(r['Lower Threshold'])]
            except:
                life_stage_period_default[r['Fish name']] = {}
                life_stage_period_default[r['Fish name']][r['Life stage']] = [[r['Start Month'],r['Start Day']],[r['End Month'],r['End Day']], r[t1], r[t2], bool(r['Consecutive']), bool(r['Lower Threshold'])]
        
        life_stage_period = collectLifeStagePeriod(life_stage_period_default, disabled=True, t_name=t_name)
        
    return life_stage_period

def collectLifeStagePeriod(life_stage_period_default, disabled, t_name='Habitat/Bankfull Area\nThreshold {0:d} (0-1)'):
    
    def collectAnswers():
        print('########\nDatabase Report\n########')
        
        answers = {}
        i = 0
        for f in life_stage_period_default:
            answers[f] = {}
            
            for l in life_stage_period_default[f]:
                
                # Set flow thresholds to NaN if provided in the wrong format or not provided
                if root.consecvar[i].get():
                    try:
                        if root.lowervar[i].get():
                            answers[f][l] = [[np.nan, np.nan], [np.nan, np.nan], float(root.t1[i].get()), float(root.t2[i].get()), root.consecvar[i].get(), root.lowervar[i].get()]
                        else:
                            root.t2[i].delete(0,tk.END)
                            root.t2[i].insert(0, "NA")
                            root.t2[i].config(state='disabled')
                            answers[f][l] = [[np.nan, np.nan], [np.nan, np.nan], float(root.t1[i].get()), np.nan, root.consecvar[i].get(), root.lowervar[i].get()]
                        print(f + ' - ' + l + ' has been saved successfully')
                    except:
                        print('Threshold must be entered in number format for ' + f + ' - ' + l + ' if using the consecutive days option. Your database has not been saved correctly and will cause errors.')
                else:
                    try:
                        if root.lowervar[i].get():
                            t1 = float(root.t1[i].get())
                            t2 = float(root.t2[i].get())
                        else:
                            t1 = float(root.t1[i].get())
                            t2 = np.nan
                            root.t2[i].delete(0,tk.END)
                            root.t2[i].insert(0, "NA")
                            root.t2[i].config(state='disabled')
                    except:
                        print('Threshold must be entered in number format for ' + f + ' - ' + l + '. Your database has not been saved correctly and will cause errors.')
                    try:
                        answers[f][l] = [[int(root.startm[i].get()), int(root.startd[i].get())], [int(root.endm[i].get()), int(root.endd[i].get())], t1, t2, root.consecvar[i].get(), root.lowervar[i].get()]
                        print(f + ' - ' + l + ' has been saved successfully')
                    except:
                        print('Start and End Months and Days must be entered for ' + f + ' - ' + l + ' in integer format. Your database has not been saved correctly and will cause errors.')
                
                i += 1
                
        life_stage_period = answers
        
        return life_stage_period
        
    root = tk.Tk()
    
    tk.Label(root,text="Start Month").grid(row=0,column=2,padx=5)
    tk.Label(root,text="Start Day").grid(row=0,column=3,padx=5)
    tk.Label(root,text="End Month").grid(row=0,column=4,padx=5)
    tk.Label(root,text="End Day").grid(row=0,column=5,padx=5)
    tk.Label(root,text="Consecutive").grid(row=0,column=6,padx=5)
    tk.Label(root,text="Threshold").grid(row=0,column=7,padx=5)
    tk.Label(root,text="Lower Threshold").grid(row=0,column=8,padx=5)
    tk.Label(root,text=t_name.format(1)).grid(row=0,column=7,padx=5)
    tk.Label(root,text=t_name.format(2)).grid(row=0,column=9,padx=5)
    
    i = 0
    root.startm = []
    root.startd = []
    root.endm = []
    root.endd = []
    root.consecvar = []
    root.consec = []
    root.t1 = []
    root.t2 = []
    root.lowervar = []
    root.lower = []
    
    for f in life_stage_period_default:
        tk.Label(root,text=f + ':').grid(row=i+1,column=0,padx=5)
        
        for l in life_stage_period_default[f]:
            tk.Label(root,text=l).grid(row=i+1,column=1,padx=5)
            
            root.startm.append(tk.Entry(root))
            root.startd.append(tk.Entry(root))
            root.endm.append(tk.Entry(root))
            root.endd.append(tk.Entry(root))
            root.t1.append(tk.Entry(root))
            root.t2.append(tk.Entry(root))
                            
            root.startm[i].grid(row=i+1,column=2)
            root.startd[i].grid(row=i+1,column=3)
            root.endm[i].grid(row=i+1,column=4)
            root.endd[i].grid(row=i+1,column=5)
            root.t1[i].grid(row=i+1,column=7)
            root.t2[i].grid(row=i+1,column=9,padx=20)
            
            # Whether or not the risk analysis should only count the first set of consecutive days that meet the criteria or any days that meet the criteria
            root.consecvar.append(tk.IntVar())
            root.consec.append(tk.Checkbutton(root, variable=root.consecvar[i], onvalue=True, offvalue=False))
            root.consec[i].grid(row=i+1,column=6)
            
            root.lowervar.append(tk.IntVar())
            root.lower.append(tk.Checkbutton(root, variable=root.lowervar[i], onvalue=True, offvalue=False))
            root.lower[i].grid(row=i+1,column=8)
            
            
            if disabled:
                root.startm[i].insert(0, '{:0.0f}'.format(life_stage_period_default[f][l][0][0]))
                root.startm[i].config(state='disabled')
                root.startd[i].insert(0, '{:0.0f}'.format(life_stage_period_default[f][l][0][1]))
                root.startd[i].config(state='disabled')
                root.endm[i].insert(0, '{:0.0f}'.format(life_stage_period_default[f][l][1][0]))
                root.endm[i].config(state='disabled')
                root.endd[i].insert(0, '{:0.0f}'.format(life_stage_period_default[f][l][1][1]))
                root.endd[i].config(state='disabled')
                root.t1[i].insert(0, '{:0.3f}'.format(life_stage_period_default[f][l][2]))
                root.t1[i].config(state='disabled')
                root.t2[i].insert(0, '{:0.3f}'.format(life_stage_period_default[f][l][3]))
                root.t2[i].config(state='disabled')
                root.consec[i].config(state='disabled')
                root.lower[i].config(state='disabled')
            else:
                root.startm[i].insert(0, "MM")
                root.startd[i].insert(0, "DD")
                root.endm[i].insert(0, "MM")
                root.endd[i].insert(0, "DD")
                root.t1[i].insert(0, "0.00")
                root.t2[i].insert(0, "0.00")
                
            i += 1
    
    if disabled:
        life_stage_period = life_stage_period_default
    else:
        life_stage_period = tk.Button(root,text="Save Data",command = collectAnswers).grid(row=i+1,column=1)
    tk.Button(root, text="Next", command=root.destroy).grid(row=i+2,column=1)
    
    root.mainloop()
    
    return life_stage_period

def figureOptions(fig_Path, fig_name):
    def getEntries():
        fig_Path = root.pathEntry.get()
        fig_name = root.nameEntry.get()
        
    root = tk.Tk()
    root.geometry("300x200")

    tk.Label(root,text="Figure Path").grid(row=1,column=0,padx=5)
    tk.Label(root,text="Figure Name").grid(row=2,column=0,padx=5)
    
    root.pathEntry = tk.Entry(root)
    root.nameEntry = tk.Entry(root)
    
    root.pathEntry.grid(row=1,column=1)
    root.nameEntry.grid(row=2,column=1)
    
    root.pathEntry.insert(0, fig_Path)
    root.nameEntry.insert(0, fig_name)
    
    tk.Button(root, text="Save", command=getEntries).grid(row=3,column=0)
    tk.Button(root, text="Next", command=root.destroy).grid(row=3,column=1)
    
    for r in np.arange(4):
        root.grid_rowconfigure(r, weight=1)
    
    for c in np.arange(2):
        root.grid_columnconfigure(c, weight=1)
        
    root.mainloop()
    
    return fig_Path, fig_name

