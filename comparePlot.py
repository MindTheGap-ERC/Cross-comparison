#!/usr/bin/env python3

############################ plot comparision #################################
# plot multiple simulation outputs on a single figure #
# suports output from rythmite, Matlab and lheureux.f

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import h5py

def plotSpatialRhy(filename, showHeaviside):
    '''
    Plot a depth profile for all solution variables at fixed time.
    For Rhythmite output.

    Parameters
    ----------
    filename : STR
        Name of the data file to plot from, format is 'solution_t_%6d.ascii',
        where the number describes the time at which the depth profile is taken.
        
    showHeaviside : BOOL
        Switch for optionally plotting the Heaviside marking the ADZ. 

    Returns
    -------
    None.

    '''
    
    df = (pd.read_csv(filename, delim_whitespace=True)).shift(axis=1).iloc[:,1:]
    
    Xs = 131.9/0.1 # depth scaling constant
    x = np.array(df.x*Xs)
    
    plt.plot(x,np.array(df.AR),label='AR')
    plt.plot(x,np.array(df.CA),label='CA')
    plt.plot(x,np.array(df.phi),label='phi')
    plt.plot(x,np.array(df.ca),label='Ca')
    plt.plot(x,np.array(df.co),label='CO')
    
    if (showHeaviside):
        plotHeaviside(x/Xs)
    
    # return the loaded df, for plotting residuals separately
    return df
    
    
def plotTemporalRhy(filename):
      '''
      Plot the time series at fixed depth for all solution variables from 
      Rhythmite output 

      Parameters
      ----------
      filename : STR
          Name of the data file to plot from, format is 'solution_x_%6d.ascii',
          where the number describes the grid position at which the time series
          is taken.

      Returns
      -------
      None.

      '''
      df = (pd.read_csv(filename, delim_whitespace=True)).shift(axis=1).iloc[:,1:]
      
      Ts = 131.9/0.1**2 # time scaling constant
      t_plot = np.array(df.x*Ts/1000)
      
      plt.plot(t_plot, np.array(df.AR), label=df.columns[1])
      plt.plot(t_plot, np.array(df.CA), label=df.columns[2])
      plt.plot(t_plot, np.array(df.phi), label=df.columns[5])
      plt.plot(t_plot, np.array(df.ca), label=df.columns[3])
      plt.plot(t_plot, np.array(df.co), label=df.columns[4])
      

def floatMarl(df):
    '''
    take the lheureux.f output as a dataframe and float it
    python apparently doesn't recognise Fortran double format anyomore
    so need to replace the D with E first
    '''
    for i in range(0,len(df.columns)):
        for j in range(0,len(df[df.columns[i]])):
            df[df.columns[i]][j] = float(df[df.columns[i]][j].replace('D','E'))


def plotSpatialFt(filename):
    '''
    Plot a depth profile for all solution variables, form the output of 
    lheureux.f

    Parameters
    ----------
    filename : STR
        File containing the Fortran depth profile.

    Returns
    -------
    None.

    '''
    df = pd.read_csv(filename, delim_whitespace=True)
    floatMarl(df)
    
    Xs = 131.9/0.1 # depth scaling constant
    x = df.x*Xs
    
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    plt.plot(x,df.AR,label='AR_ft', linestyle='--',color=colors[0])
    plt.plot(x,df.CA,label='CA_ft', linestyle='--',color=colors[1])
    plt.plot(x,df.Po,label='Po_ft', linestyle='--',color=colors[2])
    plt.plot(x,df.Ca,label='ca_ft', linestyle='--',color=colors[3])
    plt.plot(x,df.CO,label='co_ft', linestyle='--',color=colors[4])
    
    return df
    
    

def plotTemporalFt(filename,depth_ind):
    '''
    Plot the time series at fixed depth for all solution variables from the 
    output of lheureux.f

    Parameters
    ----------
    filename : STR
        Filename for the data to be plotted.
    depth_ind : INT
        Index describing which depth to plot (1 = 0, 4 = base of domain).

    Returns
    -------
    None.

    '''
    df = pd.read_csv(filename,delim_whitespace=True)
    floatMarl(df)
    
    Ts = 131.9/0.1**2 # time scaling constant
    t_plot = df.t*Ts/1000
    
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    plt.plot(t_plot, df[df.columns[depth_ind]], label=df.columns[depth_ind],\
             linestyle='--',color=colors[0])
    plt.plot(t_plot, df[df.columns[depth_ind+4]], label=df.columns[depth_ind+4],\
             linestyle='--',color=colors[1])
    plt.plot(t_plot, df[df.columns[depth_ind+8]], label=df.columns[depth_ind+8],\
             linestyle='--',color=colors[2])
    plt.plot(t_plot, df[df.columns[depth_ind+12]], label=df.columns[depth_ind+12],\
             linestyle='--',color=colors[3])
    plt.plot(t_plot, df[df.columns[depth_ind+16]], label=df.columns[depth_ind+16],\
             linestyle='--',color=colors[4])


def plotSpatialMAT(filename, showHeaviside, t_ind):
    '''
    Plot a depth profile for all solution variables at fixed time.
    For Matlab output.

    Parameters
    ----------
    filename : STR
        Name of the data file to plot from, format is '/path/to/file/Scenario_integrated.h5',
        which contains the full soln.
        
    showHeaviside : BOOL
        Switch for optionally plotting the Heaviside marking the ADZ. 
        
    t_ind : INT
        Time index at which to plot the depth profile.

    Returns
    -------
    df : Numpy array
        the soln data, in case we want to use it later.

    '''
    # load data in the hdf5 file
    sol = h5py.File(filename,'r')
    
    # convert to a np array
    df = np.array(sol['Solutions'])
    
    # this doesn't store a x points array, calc it from L_x at the number points
    nnx = len(df[0,:,0])
    
    Xs = 131.9/0.1 # depth scaling constant
    x = np.linspace(0, 500, nnx)
    
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    plt.plot(x,df[0,:,t_ind],label='AR_m',linestyle='dotted',color=colors[0])
    plt.plot(x,df[1,:,t_ind],label='CA_m',linestyle='dotted',color=colors[1])
    plt.plot(x,df[4,:,t_ind],label='phi_m',linestyle='dotted',color=colors[2])
    plt.plot(x,df[2,:,t_ind],label='Ca_m',linestyle='dotted',color=colors[3])
    plt.plot(x,df[3,:,t_ind],label='CO_m',linestyle='dotted',color=colors[4])
    
    if (showHeaviside):
        plotHeaviside(x/Xs)
    
    # return the loaded df, for plotting residuals separately
    return df

def plotSpatialPDE(filename, showHeaviside, t_ind):
    '''
    Plot a depth profile for all solution variables at fixed time.
    For py-pde output.

    Parameters
    ----------
    filename : STR
        Name of the data file to plot from, format is '/path/to/file/Scenario_integrated.h5',
        which contains the full soln.
        
    showHeaviside : BOOL
        Switch for optionally plotting the Heaviside marking the ADZ. 
        
    t_ind : INT
        Time index at which to plot the depth profile.

    Returns
    -------
    df : Numpy array
        the soln data, in case we want to use it later.

    '''
    # load data in the hdf5 file
    sol = h5py.File(filename,'r')
    
    # convert to a np array
    df = np.array(sol['data'])
    
    # this doesn't store a x points array, calc it from L_x at the number points
    nnx = len(df[0,0,:])
    
    Xs = 131.9/0.1 # depth scaling constant
    x = np.linspace(0, 500, nnx)
    
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    plt.plot(x,df[t_ind,0,:],label='AR_m',linestyle='dashdot',color=colors[0])
    plt.plot(x,df[t_ind,1,:],label='CA_m',linestyle='dashdot',color=colors[1])
    plt.plot(x,df[t_ind,4,:],label='phi_m',linestyle='dashdot',color=colors[2])
    plt.plot(x,df[t_ind,2,:],label='Ca_m',linestyle='dashdot',color=colors[3])
    plt.plot(x,df[t_ind,3,:],label='CO_m',linestyle='dashdot',color=colors[4])
    
    if (showHeaviside):
        plotHeaviside(x/Xs)
    
    # return the loaded df, for plotting residuals separately
    return df

#################### we need tstep data to plot time series ###################

# def plotTemporalMAT(filename, showHeaviside, x_ind):
#     '''
#     Plot a depth profile for all solution variables at fixed time.
#     For Matlab output.

#     Parameters
#     ----------
#     filename : STR
#         Name of the data file to plot from, format is '/path/to/file/Scenario_integrated.h5',
#         which contains the full soln.
        
#     showHeaviside : BOOL
#         Switch for optionally plotting the Heaviside marking the ADZ. 
        
#     x_ind : INT
#         Space index at which to plot the time series.

#     Returns
#     -------
#     df : Numpy array
#         the soln data, in case we want to use it later.

#     '''
#     # load data in the hdf5 file
#     sol = h5py.File(filename,'r')
    
#     # convert to a np array
#     df = np.array(sol['Solutions'])
    
#     # this doesn't store a x points array, calc it from L_x at the number points
#     nnx = len(df[0,:,0])
    
#     t = ?????
    
#     colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
#     plt.plot(x,df[0,x_ind,:],label='AR_m',linestyle='dotted',color=colors[0])
#     plt.plot(x,df[1,x_ind,:],label='CA_m',linestyle='dotted',color=colors[1])
#     plt.plot(x,df[4,x_ind,:],label='phi_m',linestyle='dotted',color=colors[2])
#     plt.plot(x,df[2,x_ind,:],label='Ca_m',linestyle='dotted',color=colors[3])
#     plt.plot(x,df[3,x_ind,:],label='CO_m',linestyle='dotted',color=colors[4])
    
    
#     # return the loaded df, for plotting residuals separately
#     return df

def plotFig3e():
    '''
    plot the digitized Figure 3e vals from L'Heureux (2018)
    
    This is the 'steady-state' case, phi_0 = 0.6, phi_init = 0.5
    as an addition to a plot from code output
    '''
    
    # use default colour sequence from matplotlib to match data
    # order should be AR, CA, Po, Ca, Co
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    bm = pd.read_csv('fig3e.csv')
    plt.scatter(np.array(bm.ARX), np.array(bm.ARY), label='bm_AR',marker='x',color=colors[0])
    plt.scatter(np.array(bm.CAX), np.array(bm.CAY), label='bm_CA',marker='x',color=colors[1])
    plt.scatter(np.array(bm.PoX), np.array(bm.PoY), label='bm_phi',marker='x',color=colors[2])
    plt.scatter(np.array(bm.CaX), np.array(bm.CaY), label='bm_ca',marker='x',color=colors[3])
    plt.scatter(np.array(bm.CoX), np.array(bm.CoY), label='bm_co',marker='x',color=colors[4])


def plotHeaviside(x):
    '''
    Plot the Heaviside function that is used to define the ADZ.
    An annotation to a plot.

    Parameters
    ----------
    x : ARRAY
        x values from the file that this is plotted with.

    Returns
    -------
    None.

    '''
    # plot the function used to define the ADZ
    
    # for now we hard-code these parameter values
    ADZ_top = 50
    ADZ_bot = 150
    x_scale = 131.9/0.1
    
    smoothK = 500
    
    h = 0.5**2 * ( 1 + np.tanh(smoothK*(x - (ADZ_top/x_scale))))*\
                 ( 1 + np.tanh(smoothK*((ADZ_bot/x_scale) - x)))
    plt.plot(x*x_scale, h, label='Heaviside', color='black', linestyle='--')
    
###############################################################################


##################### do a Spatial comparision plot ###########################
benchmarkComp = False
showHeaviside = False

savedir = ''
savefilename = 'comp_t_1_ft_py_mat.png'

fig = plt.figure(figsize=(12,10))

# plot python output 
rhy = plotSpatialRhy('%srhythmite_solution_t_000001.ascii'%(savedir), showHeaviside)
# plot the Fortran output
ft = plotSpatialFt('%samarlt1'%(savedir))
# plot the Matlab
#mat = plotSpatialMAT('%sMatlab/Scenario_integrated.h5'%(savedir), showHeaviside, 1)
# plot marlpde output
pde = plotSpatialPDE('%spy-pde/LMAHeureuxPorosityDiff.hdf5'%(savedir), showHeaviside, 1)
    
# if benchmark, plot the Fig3e data for comparison
if (benchmarkComp):
    plotFig3e()

plt.legend(loc='lower right')
plt.xlabel('x (cm)')
plt.ylabel('Concentrations')
plt.xlim(0,500)
plt.ylim(0,1)
plt.savefig('%s%s'%(savedir,savefilename))
plt.clf()


##################### also plot residuals py vs. Ft ###########################
resfilename = 'res_t_1_ft_py.png'
x_scale = 131.9/0.1

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
fig = plt.figure(figsize=(12,10))
plt.plot(rhy.x*x_scale, np.abs(rhy.AR - ft.AR.loc[:199]), label='AR',color=colors[0], marker='x')
plt.plot(rhy.x*x_scale, np.abs(rhy.CA - ft.CA.loc[:199]), label='CA',color=colors[1], marker ='x')
plt.plot(rhy.x*x_scale, np.abs(rhy.phi - ft.Po.loc[:199]),label='Po',color=colors[2], marker='x')
plt.plot(rhy.x*x_scale, np.abs(rhy.ca - ft.Ca.loc[:199]),label='ca',color=colors[3], marker='x')
plt.plot(rhy.x*x_scale, np.abs(rhy.co - ft.CO.loc[:199]),label='co',color=colors[4], marker='x')
plt.xlabel('x (cm)')
plt.ylabel('residuals')
plt.xlim(0,500)
plt.ylim(0, 0.01)
plt.legend()
plt.savefig('%s%s'%(savedir,resfilename))
plt.clf()

########################## pypde vs. ft residuals #############################
resfilename = 'res_t_1_ft_pde.png'
fig = plt.figure(figsize=(12,10))
plt.plot(rhy.x*x_scale, np.abs(pde[1,0,:] - ft.AR.loc[:199]), label='AR',color=colors[0], marker='x')
plt.plot(rhy.x*x_scale, np.abs(pde[1,1,:] - ft.CA.loc[:199]), label='CA',color=colors[1], marker ='x')
plt.plot(rhy.x*x_scale, np.abs(pde[1,4,:] - ft.Po.loc[:199]),label='Po',color=colors[2], marker='x')
plt.plot(rhy.x*x_scale, np.abs(pde[1,2,:] - ft.Ca.loc[:199]),label='ca',color=colors[3], marker='x')
plt.plot(rhy.x*x_scale, np.abs(pde[1,3,:] - ft.CO.loc[:199]),label='co',color=colors[4], marker='x')
plt.xlabel('x (cm)')
plt.ylabel('residuals')
plt.xlim(0,500)
plt.ylim(0, 0.01)
plt.legend()

plt.savefig('%s%s'%(savedir,resfilename))
plt.clf()

########################## pypde vs. rhythmite residuals ######################
resfilename = 'res_t_1_py_pde.png'
fig = plt.figure(figsize=(12,10))
plt.plot(rhy.x*x_scale, np.abs(pde[1,0,:] - rhy.AR), label='AR',color=colors[0], marker='x')
plt.plot(rhy.x*x_scale, np.abs(pde[1,1,:] - rhy.CA), label='CA',color=colors[1], marker ='x')
plt.plot(rhy.x*x_scale, np.abs(pde[1,4,:] - rhy.phi),label='Po',color=colors[2], marker='x')
plt.plot(rhy.x*x_scale, np.abs(pde[1,2,:] - rhy.ca),label='ca',color=colors[3], marker='x')
plt.plot(rhy.x*x_scale, np.abs(pde[1,3,:] - rhy.co),label='co',color=colors[4], marker='x')
plt.xlabel('x (cm)')
plt.ylabel('residuals')
plt.xlim(0,500)
plt.ylim(0, 0.01)
plt.legend()
plt.savefig('%s%s'%(savedir,resfilename))
plt.clf()

########################### do a temporal plot ###############################

# Ts = 131.9/0.1**2 # time scaling constant
# tf = 5 # ka
# savefilename = 'comp_x_199_ft_py.png'

# fig = plt.figure(figsize=(12,10))

# # plot python output
# plotTemporalRhy('solution_x_000199.ascii')
# # plot ft output
# # depth ind, should match the position of the python file
# plotTemporalFt('amarlx', 3)

# plt.legend()
# plt.xlim(0,tf)
# plt.ylim(0,2.3)
# plt.xlabel('t (ka)')
# plt.ylabel('Concentrations') 
# plt.savefig('%s'%(savefilename))
# plt.clf()  