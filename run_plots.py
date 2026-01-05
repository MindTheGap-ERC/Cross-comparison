##################### do a Spatial comparision plot ###########################
benchmarkComp = False
showHeaviside = False

savedir = ''
savefilename = 'comp_Tstar_1_py_pde.png'

fig = plt.figure(figsize=(12,10))

# plot python output 
# rhy = plotSpatialRhy('%sdt10/rhytmite/solution_t_000006.ascii'%(savedir), showHeaviside)
# plot the Fortran output
ft = plotSpatialFt('%samarlt1'%(savedir))
# plot the Matlab
mat = plotSpatialMAT('%sdt10/pdepe/Scenario_integrated_10_steps.h5'%(savedir), showHeaviside, 9)
# plot marlpde output
pde = plotSpatialPDE('%sdt10/marlpde/LMAHeureuxPorosityDiff.hdf5'%(savedir), showHeaviside, 1)

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
plt.title('residuals rhythmite vs fortran')
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
plt.title('residuals between pypde and fortran')
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
plt.title('residuals between pypde and rhythmite')
plt.xlim(0,500)
plt.ylim(0, 0.01)
plt.legend()
plt.savefig('%s%s'%(savedir,resfilename))
plt.clf()

########################## pdepe vs. rhythmite residuals #####################

resfilename = 'res_t_1_rhythmite_pdepe.png'
t_index = 9 # index of time step
fig = plt.figure(figsize=(12, 10))
def avg(x):
    '''
    take average of successive vector entries to account for differences in size of grid between implementations
    '''
    avg = 0.5 * (x[1:] + x[:-1])
    return avg

plt.plot(rhy.x * x_scale, np.abs(rhy.AR - avg(mat[0, :, t_index])), label='AR', color=colors[0], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(rhy.CA - avg(mat[1, :, t_index])), label='CA', color=colors[1], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(rhy.phi - avg(mat[4, :, t_index])), label='Po', color=colors[2], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(rhy.ca - avg(mat[2, :, t_index])), label='ca', color=colors[3], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(rhy.co - avg(mat[3, :, t_index])), label='co', color=colors[4], marker = 'x')
plt.xlabel('x (cm)')
plt.ylabel('residuals')
plt.xlim(0, 500)
plt.ylim(0, 0.01)
plt.legend()
plt.title("Residuals between rhythmite and matlab (pdepe)")
plt.savefig('%s%s'%(savedir, resfilename))

########################## pdepe vs. fortran residuals #####################

resfilename = 'res_t_1_fortran_pdepe.png'
t_index = 9 # index of time step
fig = plt.figure(figsize=(12, 10))

plt.plot(rhy.x * x_scale, np.abs(ft.AR.loc[:199] - avg(mat[0, :, t_index])), label='AR', color=colors[0], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(ft.CA.loc[:199] - avg(mat[1, :, t_index])), label='CA', color=colors[1], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(ft.Po.loc[:199] - avg(mat[4, :, t_index])), label='Po', color=colors[2], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(ft.Ca.loc[:199] - avg(mat[2, :, t_index])), label='ca', color=colors[3], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(ft.CO.loc[:199] - avg(mat[3, :, t_index])), label='co', color=colors[4], marker = 'x')
plt.xlabel('x (cm)')
plt.ylabel('residuals')
plt.xlim(0, 500)
plt.ylim(0, 0.01)
plt.legend()
plt.title("Residuals between fortran and matlab (pdepe)")
plt.savefig('%s%s'%(savedir, resfilename))

########################## pdepe vs. pypde residuals #####################

resfilename = 'res_t_1_pypde_pdepe.png'
t_index = 9 # index of time step
fig = plt.figure(figsize=(12, 10))

plt.plot(rhy.x * x_scale, np.abs(pde[1,0,:] - avg(mat[0, :, t_index])), label='AR', color=colors[0], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(pde[1,1,:] - avg(mat[1, :, t_index])), label='CA', color=colors[1], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(pde[1,4,:]- avg(mat[4, :, t_index])), label='Po', color=colors[2], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(pde[1,2,:] - avg(mat[2, :, t_index])), label='ca', color=colors[3], marker = 'x')
plt.plot(rhy.x * x_scale, np.abs(pde[1,3,:] - avg(mat[3, :, t_index])), label='co', color=colors[4], marker = 'x')
plt.xlabel('x (cm)')
plt.ylabel('residuals')
plt.xlim(0, 500)
plt.ylim(0, 0.01)
plt.legend()
plt.title("Residuals between pypde and matlab (pdepe)")
plt.savefig('%s%s'%(savedir, resfilename))

########################### do a temporal plot ###############################

# Ts = 131.9/0.1**2 # time scaling constant
# tf = 5 # ka
# savefilename = 'comp_x_199_ft_py.png'

# fig = plt.figure(figsize=(12,10))

# # plot python output
# plotTemporalRhy('solution_x_000199.ascii')
# # plot ft output
# # depth ind, should match the position of the python file
plotTemporalFt('amarlx', 3)

# plt.legend()
# plt.xlim(0,tf)
# plt.ylim(0,2.3)
# plt.xlabel('t (ka)')
# plt.ylabel('Concentrations') 
# plt.savefig('%s'%(savefilename))
# plt.clf()  

# Plot temporal output of marlpde at some intermediate depth

savedir = ''
savefilename = 'marlpde-temporal.png'

fig = plt.figure(figsize=(12,10))

pde_t = plotTemporalPDE('%sdt10/marlpde/LMAHeureuxPorosityDiff.hdf5'%(savedir), 150)
