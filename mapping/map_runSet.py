import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares
import matplotlib.pyplot as plt

from map_tools import *
from map_cycle import *
from map_run import *

# plot parameters
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12,8]
plt.rcParams["figure.dpi"]=200


class RunSet:
    
    def __init__(self, runrange, cycrange, runtype=None, dtformat='v3', calib=False, trim=False):
        """
        settype = 'rings' or 'zscans'
        rawdttype = 'new' or 'old'
        """
        runs = {}
        for run in runrange:
            runs[run] = Run(run, cycrange, runtype=runtype, dtformat=dtformat, calib=calib, trim=trim)           
        self.runs = runs
        self.runtype = runtype
        self.nr = len(self.runs)
        print("found {} runs of '{}' runtype".format(self.nr, self.runtype))
        
    def saveVars(self):
        
        V = vars(self)
        oldVars = {}
        for var in V:
            oldVars[var] = getattr(self, var)
        self.oldVars = oldVars
        
        runs = self.runs
        for run in runs:
            runs[run].saveVars()
            
    def restoreVars(self):
        
        try:
            oldVars = self.oldVars
        except AttributeError:
            print("no saved variables to restore")
            return 0
        
        for var in oldVars:
            setattr(self, var, oldVars[var])
            
        runs = self.runs
        for run in runs:
            runs[run].restoreVars()
    
    def substractMeanField(self):
        runs = self.runs
        for run in runs:
            runs[run].substractMeanField()
    
    def runsPreview(self, cyc):
        runs = self.runs
        for run in runs:
            print("run {}".format(run))
            try:
                cycle = self.runs[run].cycles[cyc]
                print("scantype :  {}".format(cycle.scantype))
                meanBz = np.mean(cycle.B[2])
                print("mean Bz = {:.2f}".format(meanBz))
            except FileNotFoundError:
                print("file not found")
                continue
                
    def getGlmS(self, lmax, source='rho', norm=True, out=False):
        if self.runtype=='rings':
            runs = self.runs
            Gi = []
            Gi_err = []
            for run in runs:
                print("run {}".format(run))
                gi, gi_err = runs[run].getGlm(lmax=lmax, source=source, norm=norm)
                Gi.append(gi)
                Gi_err.append(gi_err)
            G = np.stack(Gi, axis=2)
            G_err = np.stack(Gi_err, axis=2)
            self.G = G
            self.G_err = G_err
            self.lmax = lmax
            if out==True:
                return self.G, self.G_err
        
        else:
            print("can only extract Glms for 'rings' runtype. Try getG0 for 'zscans' runtype.")
            return 0

    def getG0(self, nfit=9):
        if self.runtype=='zscans':
            runs = self.runs
            nr = self.nr
            G_mean = []
            G_std = []
            for run in runs:
                print("run {}".format(run))
                pars, errs = runs[run].getPars(nfit=nfit)
                g = np.array([ [ pars[i][2][l] for i in range(runs[run].nc) ] for l in range(nfit) ])
                g_mean = np.mean(g, axis=1)
                g_std = np.std(g, axis=1)
                G_mean.append(g_mean)
                G_std.append(g_std)
            self.G_mean = np.transpose(G_mean)
            self.G_std = np.transpose(G_std)
            return self.G_mean, self.G_std
        else:
            print("Better use Glm for 'rings' runtype.")
            return 0
  

    def plotGmS(self, l, lmax=4, dm=0.1, colormap=None):
    
        runs = self.runs
        nr = self.nr
        try:
            Gs = self.G
            Gs_err = self.G_err
            lmax = self.lmax
        except AttributeError:
            Gs, Gs_err = self.getGlmS(lmax=lmax)
        
        fig, ax = plt.subplots(1, 1)
        if colormap==None:
            cmap = plt.get_cmap("tab10")
        else:
            cmap = colormap
        
        M = np.arange(-l-1, l+2)
        moff = dm*(nr//2)
        m1 = 0
        c = 0
        for r in runs:
            G, G_err = runs[r].G, runs[r].G_err
            for m in M:
                if m==0:
                    bar = ax.bar(-moff+m+m1, G[l][lmax+1+m], yerr=G_err[l][lmax+1+m], align='center', alpha=0.5, ecolor='black', capsize=5, color=cmap(c%10), width=dm, label='Run {}'.format(r))
                else:
                    bar = ax.bar(-moff+m+m1, G[l][lmax+1+m], yerr=G_err[l][lmax+1+m], align='center', alpha=0.5, ecolor='black', capsize=5, color=cmap(c%10), width=dm)
            m1 += dm
            c += 1
        ax.grid()
        ax.set_xticks(M)
        ax.set_xlabel(r"$m$", size=15)
        if l==0:
            ax.set_ylabel(r"$G_{{ {}m }}$ (pT)".format(l), size=15)
        elif l==1:
            ax.set_ylabel(r"$G_{{ {}m }}$ (pT/cm)".format(l), size=15)
        elif l>1:
            ax.set_ylabel(r"$G_{{ {}m }}$ (pT/cm$^{{{}}}$)".format(l, l), size=15)
        ax.legend()
        
        return fig, ax
    
    
    ### old version
    def GlmAlternateDegauss_old(self, cycrange, iB=2, iG=1, nfit=9):  

        L6 = []
        L6X = []
        G10_6 = []
        G10_6_std = []
        G10_6X = []
        G10_6X_std = []

        n = len(cycrange)
        runs = self.runs

        for run in runs:
            r = runs[run]
            r.getPars(nfit=nfit)
            G10s = np.array([r.pars[j][iB][iG] for j in range(runs[run].nc)])
            if (run%4)<2:
                print("run {} --> L6".format(run))
                L6.append(run)
                G10_6.append(np.mean(G10s))
                G10_6_std.append(np.std(G10s))
            else:
                print("run {} --> L6X".format(run))
                L6X.append(run)
                G10_6X.append(np.mean(G10s))
                G10_6X_std.append(np.std(G10s))

        L6 = np.array(L6)
        L6X = np.array(L6X)
        G10_6 = np.array(G10_6)
        G10_6X = np.array(G10_6X)
        G10_6_std = np.array(G10_6_std)
        G10_6X_std = np.array(G10_6X_std)
        G6err = G10_6_std/np.sqrt(n)
        G6Xerr = G10_6X_std/np.sqrt(n)

        return L6, G10_6, G10_6_std, G6err, L6X, G10_6X, G10_6X_std, G6Xerr

    
    def GlmAlternateDegauss(self, cycrange, iB=2, iG=1, nfit=9, plot=True):

        X = []
        G0 = []
        G0_std = []
        n = len(cycrange)
        runs = self.runs
        for run in runs:
            r = runs[run]
            print("run {}".format(run))
            r.getPars(nfit=nfit)
            G0s = np.array([r.pars[j][iB][iG] for j in range(runs[run].nc)])
            X.append(run)
            G0.append(np.mean(G0s))
            G0_std.append(np.std(G0s))
        X = np.array(X)
        G0 = np.array(G0)
        G0_std = np.array(G0_std)
        G0_err = G0_std/np.sqrt(n)
        
        if plot==True:
            fig, ax = plt.subplots(1, 1)
            i=0
            for i in range(len(X)):
                if (i%4)<2:
                    if i==0:
                        ax.errorbar(X[i], G0[i], G0_err[i], 0, 'o', capsize=5, ms=6, color='tab:blue', label="after degaussing L6 (20 cycle average)")
                    ax.errorbar(X[i], G0[i], G0_err[i], 0, 'o', capsize=5, ms=6, color='tab:blue')
                else:
                    if i==2:
                        ax.errorbar(X[i], G0[i], G0_err[i], 0, 'o', capsize=5, ms=6, color='tab:orange', label="after degaussing L6X (20 cycle average)")
                    ax.errorbar(X[i], G0[i], G0_err[i], 0, 'o', capsize=5, ms=6, color='tab:orange')
            ax.set_xlabel('Run', size=15)
            ax.set_ylabel(r'$G_{10}$ (pT/cm)', size=15)
            ax.legend()
            ax.set_xticks(np.linspace(X[0], X[-1], X[-1]-X[0]+1))
            ax.set_title(r"$G_{10}$ from $z$-scans under degaussing pattern AABB. $B_0$ down.", size=15)
            ax.grid()
            plt.show()

        return X, G0, G0_std, G0_err
    
    # deprecated
    def getFourierCoefs_old(self, ab, n, T=360):
        runs = self.runs
        nr = self.nr
        s = 0
        for run in runs:
            if s==0:
                C = runs[run].getFourierCoefs(ab, n, T)
                nc = runs[run].nc
                s = 1
            else:
                if runs[run].nc==nc:
                    c = runs[run].getFourierCoefs(ab, n, T)
                    for cyc in runs[run].cycles:
                        C[cyc] += c[cyc]/nr
                #else:
                    #nr += -1
        return C
    
    def avgFourierCoef(self, ab, n, T=360):
        """
        ab = 'a' or 'b', n = coef order
        """
        coefs = {}
        for r in self.runs:
            if ab=='a':
                coefs[r] = self.runs[r].getFourierCoefs('a', n, T)
            elif ab=='b':
                coefs[r] = self.runs[r].getFourierCoefs('b', n, T)
            else:
                print("ab = 'a' or 'b'")
                return None
        coefAvg = {}
        coefStd = {}
        key0 = [*coefs][0]
        if self.runs[key0].runtype==None:
            cycles = self.runs[key0].cycles['rings']
        else:
            cycles = self.runs[key0].cycles 
        for cyc in cycles:
            coefAvg[cyc] = np.mean([ coefs[r][cyc] for r in self.runs ], axis=0)
            coefStd[cyc] = np.std([ coefs[r][cyc] for r in self.runs ], axis=0)
        return coefAvg, coefStd
            
        
    def getRingCoefs(self, ring, icoefmax, T=360):
        """
        ring = tuple (rho, z)
        """
        # initialize fourier coefs for all proj with dynamic run size
        As, Bs = {}, {}
        for ib in range(3):
            As[ib] = {}
            Bs[ib] = {}
            for icoef in range(icoefmax):
                As[ib][icoef] = []
                Bs[ib][icoef] = []
        # fill fourier coefs
        for run in self.runs:
            if self.runs[run].runtype==None:
                cyc = self.runs[run].cycles['rings'][ring]
            else:
                cyc = self.runs[run].cycles[ring]
            for icoef in range(icoefmax):
                a = cyc.getFourierCoef('a', icoef)
                b = cyc.getFourierCoef('b', icoef)
                for ib in range(3):
                    As[ib][icoef].append( a[ib] )
                    Bs[ib][icoef].append( b[ib] )
        return As, Bs