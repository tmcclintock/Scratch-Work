"""
Do the fitting for all the samples in the SDSS sample.
"""
import fitsio, sys, os
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
plt.rc("text",usetex=True,fontsize=24)
savepath = "figures/comparison_index%d.png"

#Get the input data
fname = "dr8_run_runpos.fit"
data,header = fitsio.read(fname,header=True)
lam_trues = data['LAMBDA_CHISQ']
z_trues = data['Z_LAMBDA']
lam_arrays = data['LAMBDA_CHISQS']
zname = 'redshift_list.txt'
zs = np.genfromtxt(zname).flatten()

#The comparison functions
def get_lam_model(z,sigmaz,z_true,lam_true):
    return lam_true*np.exp(-0.5*(z_true-z)**2/sigmaz**2)

def total_diff(params,zs,lams):
    z_true,sigmaz,lam_true = params
    if any(params < 0.0): return np.inf
    if z_true > 2.0: return np.inf
    if lam_true > 300: return np.inf
    if sigmaz < 0.005: return np.inf #Avoids numerical issues
    lam_model = get_lam_model(zs,sigmaz,z_true,lam_true)
    #if any(lam_model == 0.0): return np.inf #Avoids numerical issues
    inds = lams > max(lams)/2.0 #Only the top half
    X = np.fabs(lams-lam_model)**2/lam_model**2
    ret_value = sum(X[inds])
    return ret_value

#The plotting functio
def make_comparison(sigmaz,z_true,zs,lam_data,lam_true,see_plots,index=0):
    plt.plot(zs,lam_data)
    plt.scatter(z_true,lam_true)
    plt.scatter(z_trues[index],lam_trues[index],marker="^")
    plt.plot([z_true,z_true],[0,lam_true])
    plt.plot(zs,get_lam_model(zs,sigmaz,z_true,lam_true))
    plt.ylabel("Richness",fontsize=24)
    plt.xlabel("Redshift",fontsize=24)
    plt.title("Cluster #%d"%index)
    plt.gcf().savefig(savepath%index)
    ylim = plt.gca().get_ylim()
    plt.ylim(-10,ylim[1])
    plt.subplots_adjust(bottom=0.15,left=0.15)
    if see_plots:
        plt.show()
    plt.clf()
    return

#Flow control
do_plots = True
see_plots = True
save_outputs = False

N_samples = len(lam_trues)
sigma_z = np.zeros((N_samples))
z_best = np.zeros((N_samples))
lam_best = np.zeros((N_samples))
print N_samples
inds = np.random.randint(0,N_samples,100)
#print inds
#for i in inds:
for i in xrange(0,N_samples):
    lam_data = lam_arrays[i]
    lam_data[lam_data<0.0] = 0.0
    x0 = [z_trues[i],0.03,lam_trues[i]]
    theargs = (zs,lam_data)
#    if do_plots or 0.2<z_trues[i]<0.21:
    if do_plots or 0.3<z_trues[i]<0.31:
        result = minimize(total_diff,x0=x0,args=theargs,method='Nelder-Mead')
        z_best[i],sigma_z[i],lam_best[i] = result['x']
        print "Cluster %d sigmaz = %f"%(i,sigma_z[i])
        print "Creating figure for cluster %d"%i
        make_comparison(sigma_z[i],z_best[i],zs,lam_data,lam_best[i],see_plots,index=i)
    sys.exit()
    continue #end i

if save_outputs:
    np.savetxt("sigma_z_all.txt",sigma_z)
    np.savetxt("z_best_all.txt",z_best)
    np.savetxt("lam_best_all.txt",lam_best)

