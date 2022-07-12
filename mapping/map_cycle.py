import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares
import matplotlib.pyplot as plt
import os.path

from map_tools import *

# plot parameters
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12,8]
plt.rcParams["figure.dpi"]=200

            
class Cycle:
    
    def __init__(self, run, cyc, path=None, dtformat='v3swap', trim=False):
        """
        dttype = 'v1', 'v2', 'v3', 'v3swap'. Current datatype is 'v3swap'.
        """
        # extract data depending on data type
        if dtformat=='non_binary':
            dat0 = np.load("data/run{}_cyc{}.npy".format(run, cyc), allow_pickle=True)
            self.t = dat0.item()['t']
            self.R = dat0.item()['R']
            self.B = dat0.item()['B']
        else:
        # in the case of binary data use getData() method
            if path==None:
                self.path = "../maps_bin/{}_{}_000_smapper.EDMdat".format(str(run).zfill(6), str(cyc).zfill(6)) 
            else:
                self.path = path
            if os.path.exists(self.path)==False:
                #print("cycle {} not found".format(cyc))
                self.file_exists = False
                return None
            else:
                self.file_exists = True
                self.getData(dtformat=dtformat)
        self.getType()
        self.rho, self.phi, self.z = self.R
        self.Brho, self.Bphi, self.Bz = self.B
        if self.iR!=None:
            self.bounds = Tools.gradScan(self.R[self.iR])
        if trim==True:
            self.trimData()
    
    def saveVars(self):        
        V = vars(self)
        oldVars = {}
        for var in V:
            oldVars[var] = getattr(self, var)
        self.oldVars = oldVars
            
    def restoreVars(self):
        try:
            oldVars = self.oldVars
        except AttributeError:
            print("no saved variables to restore")
            return 0
        for var in oldVars:
            setattr(self, var, oldVars[var])
    
    def trimData(self, clen=10):
        # clen is constant period length
        ilen = len(self.t)
        i=0
        if self.scantype=='rho':
            rho = self.rho
            while rho[i]!=rho[i+clen] and i<(ilen-clen-1):
                i+=1
        elif self.scantype=='phi':
            phi = self.phi
            while (phi[i]!=phi[i+clen] or phi[i]>350) and i<(ilen-clen-1):
                i+=1
        elif self.scantype=='z':
            z = self.z
            while z[i]!=z[i+clen] and i<(ilen-clen-1):
                i+=1
        imax = i-1
        self.t = self.t[:imax]
        self.R = np.array([ self.R[i][:imax] for i in range(3) ])
        self.B = np.array([ self.B[i][:imax] for i in range(3) ])
        self.rho, self.phi, self.z = self.R
        self.Brho, self.Bphi, self.Bz = self.B
        #print("data trimmed, imax = {}".format(imax))
                    
    def getData(self, dtformat='v3', returnall=False, encoder=True):
        # define data type according to the header. there is an extra column ('other') not specified in the header.
        if dtformat=='v1':
            V_to_muT=1
            dt = np.dtype([("timestamp", np.uint64), 
                       ("MotorPosPhi", np.double), ("CalMotorPosPhi", np.double), ("EncPosPhi", np.double), ("CalEncPosPhi", np.double), 
                       ("MotorPosRho", np.double), ("CalMotorPosRho", np.double), ("EncPosRho", np.double), ("CalEncPosRho", np.double),
                       ("MotorPosZ", np.double), ("CalMotorPosZ", np.double), ("EncPosZ", np.double), ("CalEncPosZ", np.double),
                       ("AIN1", np.double), ("AIN2", np.double), ("AIN3", np.double), ("AIN4", np.double), ("other", np.uint64)])
        elif dtformat=='v2':
            V_to_muT=30
            dt = np.dtype([("timestamp", np.uint64), 
                       ("MotorPosPhi", np.double), ("CalMotorPosPhi", np.double), ("EncPosPhi", np.double), ("CalEncPosPhi", np.double), 
                       ("MotorPosRho", np.double), ("CalMotorPosRho", np.double), ("EncPosRho", np.double), ("CalEncPosRho", np.double),
                       ("MotorPosZ", np.double), ("CalMotorPosZ", np.double), ("EncPosZ", np.double), ("CalEncPosZ", np.double),
                       ("AIN1", np.double), ("AIN2", np.double), ("AIN3", np.double), ("AIN4", np.double), ("status", np.uint64), ("other", np.uint64)])
        elif dtformat=='v3' or dtformat=='v3swap':
            V_to_muT=6
            dt = np.dtype([("timestamp", np.uint64), 
                       ("MotorPosRho", np.double), ("CalMotorPosRho", np.double), ("EncPosRho", np.double), ("CalEncPosRho", np.double), 
                       ("MotorPosPhi", np.double), ("CalMotorPosPhi", np.double), ("EncPosPhi", np.double), ("CalEncPosPhi", np.double),
                       ("MotorPosZ", np.double), ("CalMotorPosZ", np.double), ("EncPosZ", np.double), ("CalEncPosZ", np.double),
                       ("0", np.double), ("1", np.double), ("PinPosition", np.double), ("RotationPosition", np.double),
                       ("AIN1", np.double), ("AIN2", np.double), ("AIN3", np.double), ("AIN4", np.double), ("status", np.uint64), ("other", np.uint64)])
            
        V_to_pT = 1e6*V_to_muT
        
        with open(self.path, "rb") as f:
            data = np.fromfile(f, dtype=dt)
            
        t = (np.array(data["timestamp"])-np.array(data["timestamp"])[0])*1e-9 # in s
        R = np.array([data["CalMotorPosRho"]*0.1, data["CalMotorPosPhi"], data["CalMotorPosZ"]*0.1])
        if dtformat=='v2':
            R = np.array([data["CalMotorPosPhi"]*0.1, data["CalMotorPosRho"], data["CalMotorPosZ"]*0.1])
        eR = np.array([data["CalEncPosRho"]*0.1, data["CalEncPosPhi"], data["CalEncPosZ"]*0.1])
        if dtformat=='v3swap':
            B = np.array([data["PinPosition"], data["0"], data["1"]])*V_to_pT # in pT
        else:
            B = np.array([data["AIN3"], data["AIN1"], data["AIN2"]])*V_to_pT # in pT
        Bother = np.array(data["AIN4"])*V_to_pT # in pT
        if dtformat=='v3':
            other = np.array([data["0"], data["1"], data["PinPosition"], data["RotationPosition"], data["status"], data["other"]])
        elif dtformat=='v3swap':
            other = np.array([data["AIN1"], data["AIN2"], data["AIN3"], data["RotationPosition"], data["status"], data["other"]])
        
        if encoder==True:
            self.t, self.R, self.B = t, eR, B
        else:
            self.t, self.R, self.B = t, R, B
            
        if returnall==True:
            return data
        else:
            return t, R, B
        
    def getType(self, minRange=None):
        t = self.t
        rho, phi, z = self.R
        if minRange==None:
            minRhoRange=60
            minPhiRange=350
            minZRange=50
        else:
            minRhoRange=minRange
            minPhiRange=minRange
            minZRange=minRange
        RhoRange = max(rho)-min(rho)
        PhiRange = max(phi)-min(phi)
        ZRange = max(z)-min(z)
        if PhiRange>minPhiRange:
            self.scantype = "phi"
            self.iR = 1
            self.ring = (round(rho[0]), round(z[0]))
        elif ZRange>minZRange:
            self.scantype = "z"
            self.iR = 2
        elif RhoRange>minRhoRange:
            self.scantype = "rho"
            self.iR = 0
        else:
            self.scantype = "static"
            self.iR = None
        return self.scantype

    # old method for determining scan type
    def getType_old(self, sampleLen=100):
        t = self.t
        R = self.R
        B = self.B
        if (len(t)-1)<sampleLen:
            n = len(t)-1
        else:
            n = sampleLen
        i=0 # preliminary test of scan type
        if R[0][i+1]!=R[0][i] and R[1][i+1]==R[1][i] and R[2][i+1]==R[2][i]:
            self.scantype = "rho"
            self.iR = 0
        elif R[0][i+1]==R[0][i] and R[1][i+1]!=R[1][i] and R[2][i+1]==R[2][i]:
            self.scantype = "phi"
            self.iR = 1
            self.ring = (round(R[0][0]), round(R[2][0]))
        elif R[0][i+1]==R[0][i] and R[1][i+1]==R[1][i] and R[2][i+1]!=R[2][i]:
            self.scantype = "z"
            self.iR = 2
        else:
            self.scantype = "static"
            self.iR = None
        return self.scantype
    
    def simplePlot(self, xlabel, ylabel):
        """
        easy to use plotting tool
        xlabel = 't', 'rho', 'phi', 'z'
        ylabel = 'Brho', 'Bphi', 'Bz'
        """
        
        fig, ax = plt.subplots(1, 1)
        
        x = getattr(self, xlabel)
        y = getattr(self, ylabel)
        
        def plotLabel(string):  
            if string=='t':
                label="$t$ (s)"
            elif string=='rho':
                label="$\\rho$ (cm)"
            elif string=='phi':
                label="$\\varphi$ (°)"
            elif string=='z':
                label="$z$ (cm)"
            elif string=='Brho':
                label="$B_\\rho$ (pT)"
            elif string=='Bphi':
                label="$B_\\varphi$ (pT)"
            elif string=='Bz':
                label="$B_z$ (pT)"
            return label
        
        ax.plot(x, y, '-o')
        ax.set_xlabel(r"{}".format(plotLabel(xlabel)), size=15)
        ax.set_ylabel(r"{}".format(plotLabel(ylabel)), size=15)
        ax.grid()
        
        return fig, ax

        
    def cycFit(self, bounds=None, nfit=4, minufit=False, norm=True):
        
        iR = self.iR
        if iR==None:
            print("can only fit rho phi or z scantypes")
            return 0
            
        if bounds==None:
            imin = 0
            imax = len(self.R[iR])-1
        else:
            imin = bounds[0]
            imax = bounds[1]
            
        x = self.R[iR][imin:imax]
        B = np.array([self.B[i][imin:imax] for i in range(3)])
        T = np.max(x)-np.min(x)
        
        if minufit==True:
            par = []
            err = []
            for i in range(3):
                if iR==0:
                    poly_LS = LeastSquares(x, B[i], np.ones_like(x), Functions.fourier)
                elif iR==1:
                    poly_LS = LeastSquares(x, B[i], np.ones_like(x), Functions.fourier)
                elif iR==2:
                    poly_LS = LeastSquares(x, B[i], np.ones_like(x), Functions.polynomial)
                minu = Minuit(poly_LS, np.ones(nfit))
                self.fittable = minu.migrad()
                par.append(np.array(minu.values))
                err.append(np.array(minu.errors))
            self.par = np.array(par)
            self.err = np.array(err)
            
        else:
            if iR==0:
                par = np.zeros((3, nfit+1))
                err = np.zeros((3, nfit+1))
                for i in range(3):
                    par[i], err[i] = Fits.polynomial(B[i], x, N=nfit)
            elif iR==1:
                par = np.zeros((3, 2*nfit+1))
                err = np.zeros((3, 2*nfit+1))
                for i in range(3):
                    par[i], err[i] = Fits.fourier(B[i], x, T=T, N=nfit)
            elif iR==2:
                par = np.zeros((3, nfit+1))
                err = np.zeros((3, nfit+1))
                for i in range(3):
                    par[i], err[i] = Fits.polynomial(B[i], x, N=nfit)
            par = np.array(par)
            err = np.array(err)
            
        if bounds==None:
            self.par = par
            self.err = err
            self.A0 = np.array([ self.par[iB][0] for iB in range(3)])
            self.An = np.array([ self.par[iB][1:nfit+1] for iB in range(3)])
            self.Bn = np.array([ self.par[iB][nfit+1:] for iB in range(3)])
            if iR==0:
                self.nfit = int(len(par[0])-1)
            elif iR==1:
                self.nfit = int((len(par[0])-1)/2)
            elif iR==2:
                self.nfit = int(len(par[0])-1)
            if norm==True and iR==1:
                err = self.errNorm()
            
        return par, err
        
    def getBfit(self, nfit=4, bounds=None, minufit=False):
        
        iR=self.iR
        
        if bounds==None:
            imin = 0
            imax = len(self.R[iR])-1
            try:
                par = self.par
                err = self.err
                nfit = self.nfit
            except AttributeError:
                par, err = self.cycFit(bounds=None, nfit=nfit, minufit=minufit)
                nfit = self.nfit
        else:
            imin = bounds[0]
            imax = bounds[1]
            par, err = self.cycFit(bounds=bounds, nfit=nfit, minufit=minufit)
            
        x = self.R[iR][imin:imax]
        B = np.array([self.B[i][imin:imax] for i in range(3)])
        T = np.max(x)-np.min(x)
        
        if iR==0:
            Bmodel = np.array([Functions.polynomial(x, par[i]) for i in range(3)])
        elif iR==1:
            Bmodel = np.array([Functions.fourier(x, par[i][:nfit+1], par[i][nfit+1:], T=T) for i in range(3)])
        elif iR==2:
            Bmodel = np.array([Functions.polynomial(x, par[i]) for i in range(3)])
            
        if bounds==None:
            self.Bmodel = Bmodel
        
        return x, B, Bmodel
    
    def errNorm(self):
    
        nfit = int((len(self.par[0])-1)/2)
        err = self.err

        if self.iR!=1:
            print("error: only works for ring scans")
            return 0
        else:
            chi2, dof = self.getFitChi2(nscans=0, nfit=nfit)

        self.err = np.array([err[i]*np.sqrt(chi2[i]/dof) for i in range(3)])
        
        return self.err
        
    def parScan(self, iB, l, m=0, N=1, nfit=4):
        """
        usable in the case of multiple scans per cycle. 
        takes l,m of generalized gradient and N number of of scans to average over
        returns average and std of Glm over scans
        """
        
        imin, imax = self.bounds
        ns = len(imin)-1

        G = []
        Gb = []
        for n in range(ns//N):
            par, err = self.cycFit([imin[n*N], imin[(n+1)*N]], nfit)
            if self.iR==2:
                G.append(par[iB][l])
            elif self.iR==1:
                if l==0:
                    if m==1:
                        G.append(par[0][1])
                        Gb.append(-par[1][nfit+1])
                    if m==-1:
                        G.append(par[0][nfit+1])
                        Gb.append(par[1][1])
            else:
                print("only works for z scans m=0, and phi scans l=0")
                break
        
        if self.iR==1:
            return np.array(G), np.array(Gb)
        else:
            return np.array(G)
        
    def getRingG0(self, m, l=0, nfit=6, minufit=False):
        """
        only valid for rho=0, z=0 ring scans
        """
        if self.iR!=1:
            print("ring scans only (iR=1)")
            return 0
        else:
            
            self.cycFit(None, nfit=nfit, minufit=minufit, norm=True)

            if l==0:
                if m==1:
                    G_fromrho = self.par[0][1]
                    Gerr_fromrho = self.errn[0][1]
                    G_fromphi = -self.par[1][nfit+1]
                    Gerr_fromphi = -self.errn[1][nfit+1]
                elif m==-1:
                    G_fromrho = self.par[0][nfit+1]
                    Gerr_fromrho = self.errn[0][nfit+1]
                    G_fromphi = self.par[1][1]
                    Gerr_fromphi = self.errn[1][1]
                else:
                    print("unallowed l, m")
                    return 0
            else:
                print("unallowed l, m")
                return 0
        
        return G_fromrho, Gerr_fromrho, G_fromphi, Gerr_fromphi
    
    def cycPlot(self, iB, iR=None, nscans=0, nfit=4, dofit=True, removeG0=False, minufit=False):
        """
        iR = -1 for t, 0 for rho, 1 for phi, 2 for z
        """
        if iR==None:
            if self.iR!=None:
                iR = self.iR
            else:
                print("none")
                iR = float(input("iR = (-1 for t, 0 for rho, 1 for phi, 2 for z)"))
            
        if nscans==0:
            istart = 0
            istop = len(self.t)-1
        else:
            imin, imax = Tools.gradScan(self.R[iR])
            istart = imin[0]
            istop = imin[nscans]

        t = np.array(self.t[istart:istop])
        if iR>=0:
            x = np.array(self.R[iR][istart:istop])
        B = np.array([self.B[i][istart:istop] for i in range(3)])

        if dofit==True:
            if nscans==0:
                try:
                    par = self.par
                    err = self.err
                    nfit = self.nfit
                except AttributeError:
                    par, err = self.cycFit(None, nfit=nfit, minufit=minufit)
                    nfit = self.nfit
            else:
                par, err = self.cycFit(bounds=[istart, istop], nfit=nfit, minufit=minufit)
            
            if iR==0:
                Y = np.array([Functions.polynomial(x, par[i]) for i in range(3)])
            elif iR==1:
                T = np.max(x)-np.min(x)
                Y = np.array([Functions.fourier(x, par[i][:nfit+1], par[i][nfit+1:], T=T) for i in range(3)])
            elif iR==2:
                Y = np.array([Functions.polynomial(x, par[i]) for i in range(3)])
                
        ################### EDIT PLOTS HERE ######################

        fig, ax = plt.subplots(1, 1)
        
        if iR==-1:
            ax.plot(t, B[iB], color='tab:blue', marker='o', linestyle="", markersize=3)
        else:
            ax.plot(x, B[iB], color='tab:blue', marker='o', linestyle="", markersize=3)
        
        if dofit==True:
            
            if removeG0==False:
                if iR==0:
                    ax.plot(x, Y[iB], color='tab:orange', marker='', linestyle="-", label='polynomial fit (order {})'.format(nfit))
                elif iR==1:
                    ax.plot(x, Y[iB], color='tab:orange', marker='', linestyle="-", label='fourier fit (order {})'.format(nfit))
                elif iR==2:
                    ax.plot(x, Y[iB], color='tab:orange', marker='', linestyle="-", label='polynomial fit (order {})'.format(nfit))
            else:
                if iR==0:
                    y0 = Functions.polynomial(x, par[iB][:1])
                elif iR==1:
                    y0 = Functions.fourier(x, par[iB][:2], par[iB][nfit+1:nfit+2], T=T)
                elif iR==2:
                    y0 = Functions.polynomial(x, par[iB][:1])
                ax.plot(x, B[iB]-y0, color='tab:blue', marker='.', linestyle="-", markersize=8)
                ax.plot(x, Y[iB]-y0, color='tab:blue', marker='.', linestyle="-", markersize=8) 

            if iR==1:
                ax.text(0.2, 0.82, ' a_0 = {:.2e} \n a_1 = {:.2e}'.format(*par[iB][0:2]), color='tab:orange', transform = ax.transAxes, size=12)
                ax.text(0.35, 0.82, ' b_1 = {:.2e}'.format(*par[iB][nfit+1:nfit+2]), color='tab:orange', transform = ax.transAxes, size=12)
            elif iR==2: 
                ax.text(0.06, 0.5, ' $G_{{00}}$ = {:.2e} pT \n $G_{{10}}$ = {:.2e} pT/cm \n $G_{{20}}$ = {:.2e} pT/cm$^2$ \n $G_{{30}}$ = {:.2e} pT/cm$^3$'.format(*par[iB][0:4]), color='tab:orange', transform = ax.transAxes, size=12)
        
        ax.grid()
        if iR==1:
            ax.set_title(r'Ring scan at $\rho =$ {:.2f} cm, $z =$ {:.2f} cm'.format(self.R[0][0], self.R[2][0]), size=15)
        elif iR==2:
            ax.set_title(r'Vertical scan at $\rho =$ {:.2f} cm, $\varphi =$ {:.2f} °'.format(self.R[0][0], self.R[1][0]), size=15)

        if iR==-1:
            ax.set_xlabel(r"$t$ (s)", size=15)
        elif iR==0:
            ax.set_xlabel(r"$\rho$ (cm)", size=15)
        elif iR==1:
            ax.set_xlabel(r"$\varphi$ (°)", size=15)
            ax.set_xlim(-5, 365)
            ax.set_xticks(np.arange(0, 361, 30))
        elif iR==2:
            ax.set_xlabel(r"$z$ (cm)", size=15)

        if iB==0:
            ax.set_ylabel(r"$B_\rho$ (pT)", size=15)
        elif iB==1:
            ax.set_ylabel(r"$B_\varphi$ (pT)", size=15)
        elif iB==2:
            ax.set_ylabel(r"$B_z$ (pT)", size=15)

        ttime = t[-1] 
        ax.text(0.06, 0.16, 'scan time = {:.0f} seconds'.format(ttime), transform = ax.transAxes)
        ax.legend()
        
        return fig, ax
    
    def getFitRes(self, nscans=0, nfit=4):
        
        iR = self.iR
        
        if nscans==0:
            istart = 0
            istop = len(self.t)-1
            try:
                par = self.par
                err = self.err
                nfit = self.nfit
            except AttributeError:
                par, err = self.cycFit(bounds=None, nfit=nfit)  
                nfit = self.nfit
            x, B, Bmodel = self.getBfit(nfit=nfit, bounds=None)
            
        else:
            imin, imax = Tools.gradScan(self.R[iR])
            istart = imin[0]
            istop = imin[nscans]
            par, err = self.cycFit(bounds=[istart, istop], nfit=nfit)
            x, B, Bmodel = self.getBfit(nfit=nfit, bounds=[istart, istop])
        
        residue = np.array([B[i] - Bmodel[i] for i in range(3)])
        
        if nscans==0:
            self.residue = residue
        
        return x, residue
    
    def getFitChi2(self, nscans, nfit):
        
        iR = self.iR
        
        if nscans==0:
            istart = 0
            istop = len(self.t)-1
            x, B, Bmodel = self.getBfit(nfit=nfit, bounds=None)
        else:
            imin, imax = Tools.gradScan(self.R[iR])
            istart = imin[0]
            istop = imin[nscans]
            x, B, Bmodel = self.getBfit(nfit=nfit, bounds=[istart, istop])

        dof = len(x)-1
        chi2 = np.array([Tools.chi2stat(B[i], Bmodel[i], np.ones_like(dof+1)) for i in range(3)])
        
        if nscans==0:
            self.chi2 = chi2
            self.dof = dof
        
        return chi2, dof
    
    def getFitRMS(self, nfit=4):
            
        x, res = self.getFitRes(nscans=0, nfit=nfit)
        RMS = np.array([ np.mean(res[i]**2)**0.5 for i in range(3) ])
        
        return RMS
    
    def getFourierCoef(self, ab, n, T=360):
        # works for a ring scan going +phi then -phi
        imax = Tools.slopeChange(self.phi)[0]
        C = np.zeros(3)
        for i in range(3):
            C[i] += 0.5*fourierCoef(ab, n, self.phi[:imax], self.B[i][:imax], T)
            C[i] += 0.5*fourierCoef(ab, n, -self.phi[imax:], self.B[i][imax:], T)
        return C
        
    
    
