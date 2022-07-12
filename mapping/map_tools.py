import numpy as np


def factorial(n):
    f=1
    for k in range(1, n+1):
        f *= k
    return f

def binomialCoef(n, k):
    return factorial(n) / (factorial(k) * factorial(n-k))

def numericalIntegral(x, y):
    n = min(len(x), len(y))
    dx = x[1:n]-x[:n-1]
    S = sum(y[:n-1]*dx)
    return S

def fourierCoef(ab, n, x, y, T):
    """
    takes: ab (string, 'a' or 'b'), y data (array), n order of coef (int), T period (float)
    returns: Fourier coef (float)
    """        
    if ab=='a':
        C = (2/T)*numericalIntegral(x, y * np.cos(2*np.pi*n*x/T))
    elif ab=='b':
        C = (2/T)*numericalIntegral(x, y * np.sin(2*np.pi*n*x/T))
    if n==0:
        C *= 0.5
    return C

class Geometry:
    
    R=40
    H=12
    Hp=18
    
    ### coefficient <rho B_rho> / (-R**2/4) in double chamber (odd) or single chamber (even)
    a1 = 1
    a2 = Hp/2
    a3 = -R**2/2 + (H**2 + 3*Hp**2)/4  
    a4 = Hp/2 * ( -R**2 + (H**2 + Hp**2)/2 )
    a5 = (5/8) * ( 8*R**4 - 2*R**2*(H**2 + 3*Hp**2)/3 + Hp**2*H**2 + (5*Hp**4 + H**4)/10 )
    
    ### GTB coef and L^n's
    L3 = 32.9**2
    C3 = 4*L3/(R**2 + 2*Hp**2)
    L5 = 32.7**4
    C5 = 48*L5/(15*R**4 + 10*R**2*(3*Hp**2-H**2) - 4*Hp**2*(3*Hp**2+5*H**2))
    
    def volAvg(i, j, k, R=40, H=12, Hp=18, chamber='double'):
        """
        Takes integers i, j, k, chamber geometry (R, H, H'), and chamber type ('single' or 'double')
        Returns < x^i y^j z^k > in given chamber geometry 
        """
        ### treat a = <xy> and b = <z> seperately

        if k%2==1:
            b = 0
        else:
            if chamber=='single':
                b = H**k / (2**k * (k+1))
            elif chamber=='double':
                b = 0
                for n in range(k//2 + 1):
                    #print(n)
                    b += binomialCoef(k, 2*n) * Hp**(k-2*n) * H**(2*n) / (2*n+1)
                b *= 1/(2**k)

        l, m = min(i, j), max(i, j)

        if l%2==1 or m%2==1:
            a = 0
        elif l==0:
            if m==0:
                a = 1
            elif m==2:
                a = R**2/4
            elif m==4:
                a = R**4/8
            elif m==6:
                a = 5*R**6/64
            elif m==8:
                a = 7*R**8/128
            else:
                print("i+j too high")
        elif l==2:
            if m==2:
                a = R**4/24
            elif m==4:
                a = R**6/64
            elif m==6:
                a = R*8/128
            else:
                print("i+j too high")
        elif l==4:
            if m==4:
                a = 3*R**8/640
            else:
                print("i+j too high")

        return a*b

    def bz2(l, m, R=40, H=12, Hp=18):
    
        if l==0:
            if m==0:
                b2 = 1
            else:
                b2 = 0

        elif l==1:
            if m==0:
                b2 = Geometry.volAvg(0, 0, 2)
            elif abs(m)==1:
                b2 = Geometry.volAvg(2, 0, 0)
            else:
                b2 = 0

        elif l==2:
            if m==0:
                b2 = Geometry.volAvg(0, 0, 4) - 2*Geometry.volAvg(2, 0, 2) + 0.5*(Geometry.volAvg(4, 0, 0) + Geometry.volAvg(2, 2, 0))                            
            elif abs(m)==1:
                b2 = 4*Geometry.volAvg(0, 2, 2)
            elif abs(m)==2:
                b2 = 2*Geometry.volAvg(4, 0, 0) - 2*Geometry.volAvg(2, 2, 0)
            else:
                b2 = 0

        return b2

    def dratio1(g1, g3):
        return g1 / (a3*g3)

    def dratio2(g2, g4):
        return g2 / (a4*g4)

    def dratio2(g3, g5):
        return g3 / (a5*g5)
    

class Tools:
    
    def cylindrical(x, y, z):
        rho, phi, z = (x**2 + y**2)**0.5, np.arctan2(y, x), z
        return rho, phi, z
    
    def cylindricalField(phi, Bx, By, Bz, degrees=True):
        if degrees==True:
            phi = phi*np.pi/180
        Brho = Bx*np.cos(phi) + By*np.sin(phi)
        Bphi = -Bx*np.sin(phi) + By*np.cos(phi)
        return Brho, Bphi, Bz

    def cartesian(rho, phi, z, degrees=True):
        if degrees==True:
            phi = np.pi*phi/180
        x, y, z = rho*np.cos(phi), rho*np.sin(phi), z
        return x, y, z
    
    def cartesianField(phi, Brho, Bphi, Bz, degrees=True):
        if degrees==True:
            phi = phi*np.pi/180
        Bx = Brho*np.cos(phi) - Bphi*np.sin(phi)
        By = Brho*np.sin(phi) + Bphi*np.cos(phi)
        return Bx, By, Bz
    
    def gradScan(x):
        imax=[]
        imin=[]
        if (x[1]-x[0])>0:
            imin.append(0)
        elif (x[1]-x[0])<0:
            imax.append(0)
        for i in range(len(x)-2):
            if np.sign(x[i+1]-x[i])!=np.sign(x[i+2]-x[i+1]):
                if (x[i+1]-x[i])>0:
                    imax.append(i+1)
                elif (x[i+1]-x[i])<0:
                    imin.append(i+1)
        return np.array(imin), np.array(imax)
    
    def slopeChange(x):
        I = []
        for i in range(len(x)-2):
            if np.sign(x[i+1]-x[i])!=np.sign(x[i+2]-x[i+1]):
                I.append(i+1)
        return np.array(I)
    
    def chi2stat(ydata, ymodel, yerr):
        c = np.sum(((ydata-ymodel)/yerr)**2)
        return c
    
    def zero_crossing(Y, up=True):
        i=0
        if up==True:
            while Y[i]>0:
                i+=1
        else:
            while Y[i]<0:
                i+=1
        return i
    
        
class Functions:
    
    def polynomial(t, a):
        """
        a = array of polynomial coefficients
        """
        f=0
        for n in range(len(a)):
            f += a[n]*t**n
        return f
    
    def fourier(t, a, b, T=2*np.pi):
        """
        takes: t scalar or array, fourier coefs a n-dim array /!\, b (n-1)-dim array /!\
        returns: fourier series
        """
        f = 0
        for n in range(len(a)):
            f += a[n]*np.cos(n*t*2*np.pi/T) 
            if n>0:
                f += b[n-1]*np.sin(n*t*2*np.pi/T)
        return f
    
    
class Fits:
    
    def polynomial(y0, t, N, sigma=None):
        """
        takes: y0 array (values to be fitted), t array (fit points), N order of the fit, variance of the measure points
        returns: N size array of polynomial coefficients (fit parameters) [a_0, a_1, ..., a_N], with their associated errors
        """

        if np.any(sigma==None):
            sigma = np.ones_like(y0)
        
        u = np.zeros(N+1)
        U = np.zeros((N+1, N+1))
        I = len(t)

        for i in range(I):
            for n in range(N+1):
                u[n] += y0[i]*t[i]**n / sigma[i]**2
                for m in range(N+1):
                    U[n][m] += t[i]**(n+m) / sigma[i]**2

        Uinv = np.linalg.inv(U)

        par = np.dot(Uinv, u)
        err = np.array([Uinv[i][i]**0.5 for i in range(N+1)])

        return par, err
    
    def fourier(y0, t, N, T=2*np.pi, sigma=None):
        """
        takes: y0 array (values to be fitted), t array (fit points), N order of the fit, variance of the measure points
        returns: 2*N+1 size array of fit parameters [a_0, a_1, ..., a_N, b_1, ..., b_N], along with errors on the parameters
        """

        if np.any(sigma==None):
            sigma = np.ones_like(y0)
            
        u = np.zeros(2*N+1)
        U = np.zeros((2*N+1, 2*N+1))
        I = len(t)

        for i in range(I):
            for n in range(2*N+1):

                if n<=N:
                    u[n] += y0[i]*np.cos(n*t[i]*2*np.pi/T) / sigma[i]**2
                    for m in range(2*N+1):
                        if m<=N:
                            U[n][m] += np.cos(m*t[i]*2*np.pi/T)*np.cos(n*t[i]*2*np.pi/T) / sigma[i]**2
                        else:
                            U[n][m] += np.sin((m-N)*t[i]*2*np.pi/T)*np.cos(n*t[i]*2*np.pi/T) / sigma[i]**2

                else:
                    u[n] += y0[i]*np.sin((n-N)*t[i]*2*np.pi/T) / sigma[i]**2
                    for m in range(2*N+1):
                        if m<N:
                            U[n][m] += np.cos(m*t[i]*2*np.pi/T)*np.sin((n-N)*t[i]*2*np.pi/T) / sigma[i]**2
                        else:
                            U[n][m] += np.sin((m-N)*t[i]*2*np.pi/T)*np.sin((n-N)*t[i]*2*np.pi/T) / sigma[i]**2

        Uinv = np.linalg.inv(U)

        par = np.dot(Uinv, u)
        err = np.array([Uinv[i][i]**0.5 for i in range(2*N+1)])

        return par, err
    
    def GlmPoly(iB, a, a_err, ring, m, lmin, lmax):
        """
        takes: iB field proj, a array of fourier coefs of order m to be fitted, ring array of tuples (rho_ring, z_ring), m order of fourier coef and Glms, L order of the fit, sigma std of the measured values
        returns: l-polynomial coefficient of Pi_lm coefficients (fixed m)
        """

        u = np.zeros(lmax+1-lmin)
        U = np.zeros((lmax+1-lmin, lmax+1-lmin))
        I = len(ring)

        for i in range(I):
            rho = ring[i][0]
            z =  ring[i][1]
            for l1 in range(lmin, lmax+1):
                u[l1-lmin] += Harmonic.PiI(iB, rho, z, l1, m) * a[i] / a_err[i]**2
                for l2 in range(lmin, lmax+1):
                    U[l1-lmin][l2-lmin] += Harmonic.PiI(iB, rho, z, l1, m) * Harmonic.PiI(iB, rho, z, l2, m) / a_err[i]**2

        Uinv = np.linalg.inv(U)

        g = np.dot(Uinv, u)
        g_err = np.array([Uinv[i][i]**0.5 for i in range(lmax+1-lmin)])

        return g, g_err

class Dipole:
    
    def field(r, rd, m):
        "r, rd, m shape(3) in cartesian coordinates"
        C = 1.256637062*10**(-1)/(4*np.pi) # pT.cm/nA
        N = np.dot(r-rd,r-rd)**(1/2)
        u = (r-rd)/N 
        Bd = C*(3*np.dot(m,u)*u - m)/N**3
        return Bd
    
    def fields(R, rd, md):
        "R array shape (3, n) in cartesian coordinates"
        RT = np.transpose(R)
        BT = []
        for r in RT:
            BT.append(Dipole.field(r, rd, md))
        B = np.transpose(BT)
        return B
    
    def multiFields(R, Rd, Md):
        "R, Rd, M arrays shape (3, n) in cartesian coordinates"
        RdT, MdT = np.transpose(Rd), np.transpose(Md)
        ndip = len(RdT)
        for i in range(ndip):
            if i==0:
                B = Dipole.fields(R, RdT[i], MdT[i])
            else:
                B += Dipole.fields(R, RdT[i], MdT[i])
        return B 
    
    def matrix(r, rd, inv=False):
        "r and rd shape (3) in cartesian coordinates"
        C = 1.256637062*10**(-1)/(4*np.pi) # pT.cm/nA
        x, y, z = r-rd
        N = np.dot(r-rd,r-rd)**(1/2)
        A = -N**(-3)*np.identity(3) + 3*N**(-5)*np.outer([x, y, z], [x, y, z])
        Ainv = -N**3*np.identity(3) + 1.5*N*np.outer([x, y, z], [x, y, z])
        if inv==False:
            return C*A
        else:
            return (1/C)*Ainv
        
    def multiMatrix(r, Rd):
        "r of shape (3), Rd of shape (3, n) in cartesian coordinates"
        RdT = np.transpose(Rd)
        F = np.zeros((3, 3))
        for i in range(len(RdT)):
            F += Dipole.matrix(r, RdT[i])
        return F
    
    def fit(R, B, rd, sigma=1):
        "R, B of shape (3, n), rd of shape(3) in cartesian coordinates"
        RT = np.transpose(R)
        BT = np.transpose(B)
        F2sum = np.zeros((3,3))
        bFsum = np.zeros(3)
        for i in range(len(RT)):
            F = Dipole.matrix(RT[i], rd)
            F2sum += np.dot(F, F)
            bFsum += np.dot(F, BT[i])
        F2suminv = np.linalg.inv(F2sum)
        mpar = np.dot(F2suminv, bFsum)
        merr = sigma*np.array([F2suminv[i][i]**0.5 for i in range(3)])
        return mpar, merr
    
    def multiFit(R, B, Rd, sigma=1):
        "R, B, Rd of shape (3, n) in cartesian coordinates"
        RT = np.transpose(R)
        BT = np.transpose(B)
        F2sum = np.zeros((3,3))
        bFsum = np.zeros(3)
        for i in range(len(RT)):
            F = Dipole.multiMatrix(RT[i], Rd)
            F2sum += np.dot(F, F)
            bFsum += np.dot(F, BT[i])
        F2suminv = np.linalg.inv(F2sum)
        mpar = np.dot(F2suminv, bFsum)
        merr = sigma*np.array([F2suminv[i][i]**0.5 for i in range(3)])
        return mpar, merr
    
    
class Harmonic:
    
    def PiI(iB, rho, z, l, m):
        if iB==0:
            return Harmonic.PiRho(rho, z, l, m)
        elif iB==1:
            return Harmonic.PiPhi(rho, z, l, m)
        elif iB==2:
            return Harmonic.PiZ(rho, z, l, m)
    
    def PiRho(rho, z, l, m):
        
        if l==0:
            
            if m==-1:
                P = 1
            elif m==0:
                P = 0
            elif m==1:
                P = 1
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
                
        elif l==1:
            
            if m==-2:
                P = rho
            elif m==-1:
                P = z
            elif m==0:
                P = -0.5*rho
            elif m==1:
                P = z
            elif m==2:
                P = rho
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
                
        elif l==2:
            
            if m==-3:
                P = rho**2
            elif m==-2:
                P = 2*rho*z
            elif m==-1:
                P = 0.25*(4*z**2 - 3*rho**2)
            elif m==0:
                P = -rho*z
            elif m==1:
                P = 0.25*(4*z**2 - 3*rho**2)
            elif m==2:
                P = 2*rho*z
            elif m==3:
                P = rho**2
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==3:
            
            if m==-4:
                P = rho**3
            elif m==-3:
                P = 3*rho**2*z
            elif m==-2:
                P = rho*(3*z**2 - rho**2)
            elif m==-1:
                P = 0.25*z*(4*z**2 - 9*rho**2)
            elif m==0:
                P = (3/8)*rho*(rho**2 - 4*z**2)
            elif m==1:
                P = 0.25*z*(4*z**2 - 9*rho**2)
            elif m==2:
                P = rho*(3*z**2 - rho**2)
            elif m==3:
                P = 3*rho**2*z
            elif m==4:
                P = rho**3
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==4:
            
            if m==-5:
                P = rho**4
            elif m==-4:
                P = 4*rho**3*z
            elif m==-3:
                P = 0.25*(24*rho**2*z**2 - 5*rho**4)
            elif m==-2:
                P = 4*(rho*z**3 - rho**3*z)
            elif m==-1:
                P = 0.125*(8*z**4 - 36*rho**2*z**2 + 5*rho**4)
            elif m==0:
                P = 0.5*(3*rho**3*z - 4*rho*z**3)
            elif m==1:
                P = 0.125*(8*z**4 - 36*rho**2*z**2 + 5*rho**4)
            elif m==2:
                P = 4*(rho*z**3 - rho**3*z)
            elif m==3:
                P = 0.25*(24*rho**2*z**2 - 5*rho**4)
            elif m==4:
                P = 4*rho**3*z
            elif m==5:
                P = rho**4
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
        
        elif l==5:
            
            if m==-6:
                P = rho**5
            elif m==-5:
                P = 5*rho**4*z
            elif m==-4:
                P = 0.5*(20*rho**3*z**2 - 3*rho**5)
            elif m==-3:
                P = 1.25*(8*rho**2*z**3 - 5*rho**4*z)
            elif m==-2:
                P = (5/16)*(16*rho*z**4 - 32*rho**3*z**2 + 3*rho**5)
            elif m==-1:
                P = 0.125*(8*z**5 - 60*rho**2*z**3 + 25*rho**4*z)
            elif m==0:
                P = (5/16)*(-8*rho*z**4 + 12*rho**3*z**2 - rho**5)
            elif m==1:
                P = 0.125*(8*z**5 - 60*rho**2*z**3 + 25*rho**4*z)
            elif m==2:
                P = (5/16)*(16*rho*z**4 - 32*rho**3*z**2 + 3*rho**5)
            elif m==3:
                P = 1.25*(8*rho**2*z**3 - 5*rho**4*z)
            elif m==4:
                P = 0.5*(20*rho**3*z**2 - 3*rho**5)
            elif m==5:
                P = 5*rho**4*z
            elif m==6:
                P = rho**5
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
        
        elif l==6:
            
            if m==-7:
                P = rho**6
            elif m==-6:
                P = 6*rho**5*z
            elif m==-5:
                P = 0.25*rho**4*(60*z**2 - 7*rho**2)
            elif m==-4:
                P = rho**3*z*(20*z**2 - 9*rho**2)
            elif m==-3:
                P = (3/16)*rho**2*(80*z**4 - 100*rho**2*z**2 + 7*rho**4)
            elif m==-2:
                P = 0.125*rho*z*(48*z**4 - 160*rho**2*z**2 + 45*rho**4)
            elif m==-1:
                P = (1/64)*(64*z**6 - 720*rho**2*z**2 + 600*rho**4*z**2 - 35*rho**6)
            elif m==0:
                P = (3/8)*rho*(-8*z**5 + 20*rho**2*z**3 - 5*rho**4*z)
            elif m==1:
                P = (1/64)*(64*z**6 - 720*rho**2*z**2 + 600*rho**4*z**2 - 35*rho**6)
            elif m==2:
                P = 0.125*rho*z*(48*z**4 - 160*rho**2*z**2 + 45*rho**4)
            elif m==3:
                P = (3/16)*rho**2*(80*z**4 - 100*rho**2*z**2 + 7*rho**4)
            elif m==4:
                P = rho**3*z*(20*z**2 - 9*rho**2)
            elif m==5:
                P = 0.25*rho**4*(60*z**2 - 7*rho**2)
            elif m==6:
                P = 6*rho**5*z
            elif m==7:
                P = rho**6
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        else:
            print("l too high")
            return None
        
        return P
    
    def PiPhi(rho, z, l, m):
        
        if l==0:
            
            if m==-1:
                P = 1
            elif m==0:
                P = 0
            elif m==1:
                P = -1
        elif l==1:
            
            if m==-2:
                P = rho
            elif m==-1:
                P = z
            elif m==0:
                P = 0
            elif m==1:
                P = -z
            elif m==2:
                P = -rho
                
        elif l==2:
            
            if m==-3:
                P = rho**2
            elif m==-2:
                P = 2*rho*z
            elif m==-1:
                P = 0.25*(4*z**2 - 3*rho**2)
            elif m==0:
                P = 0
            elif m==1:
                P = -0.25*(4*z**2 - 3*rho**2)
            elif m==2:
                P = -2*rho*z
            elif m==3:
                P = -rho**2
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==3:
            
            if m==-4:
                P = rho**3
            elif m==-3:
                P = 3*rho**2*z
            elif m==-2:
                P = 0.5*rho*(6*z**2 - rho**2)
            elif m==-1:
                P = 0.25*z*(4*z**2 - 3*rho**2)
            elif m==0:
                P = 0
            elif m==1:
                P = -0.25*z*(4*z**2 - 3*rho**2)
            elif m==2:
                P = -0.5*rho*(6*z**2 - rho**2)
            elif m==3:
                P = -3*rho**2*z
            elif m==4:
                P = -rho**3
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==4:
            
            if m==-5:
                P = rho**4
            elif m==-4:
                P = 4*rho**3*z
            elif m==-3:
                P = 0.75*(8*rho**2*z**2 - rho**4)
            elif m==-2:
                P = 2*(2*rho*z**3 - rho**3*z)
            elif m==-1:
                P = 0.125*(8*z**4 - 12*rho**2*z**2 + rho**4)
            elif m==0:
                P = 0
            elif m==1:
                P = -0.125*(8*z**4 - 12*rho**2*z**2 + rho**4)
            elif m==2:
                P = -2*(2*rho*z**3 - rho**3*z)
            elif m==3:
                P = -0.75*(8*rho**2*z**2 - rho**4)
            elif m==4:
                P = -4*rho**3*z
            elif m==5:
                P = -rho**4
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
        
        elif l==5:
            
            if m==-6:
                P = rho**5
            elif m==-5:
                P = 5*rho**4*z
            elif m==-4:
                P = rho**3*(10*z**2 - rho**2)
            elif m==-3:
                P = 1.25*(8*rho**2*z**3 - 3*rho**4*z)
            elif m==-2:
                P = (5/16)*(16*rho*z**4 - 16*rho**3*z**2 + rho**5)
            elif m==-1:
                P = 0.125*(8*z**5 - 20*rho**2*z**3 + 5*rho**4*z)
            elif m==0:
                P = 0
            elif m==1:
                P = -0.125*(8*z**5 - 20*rho**2*z**3 + 5*rho**4*z)
            elif m==2:
                P = -(5/16)*(16*rho*z**4 - 16*rho**3*z**2 + rho**5)
            elif m==3:
                P = -1.25*(8*rho**2*z**3 - 3*rho**4*z)
            elif m==4:
                P = -rho**3*(10*z**2 - rho**2)
            elif m==5:
                P = -5*rho**4*z
            elif m==6:
                P = -rho**5
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==6:
            
            if m==-7:
                P = rho**6
            elif m==-6:
                P = 6*rho**5*z
            elif m==-5:
                P = 5*rho**4*(12*z**2 - rho**2)
            elif m==-4:
                P = 2*rho**3*z*(10*z**2 - 3*rho**2)
            elif m==-3:
                P = (3/16)*rho**2*(80*z**4 - 60*rho**2*z**2 + 3*rho**4)
            elif m==-2:
                P = 0.125*rho*z*(48*z**4 - 80*rho**2*z**2 + 15*rho**4)
            elif m==-1:
                P = (1/64)*(64*z**6 - 240*rho**2*z**2 + 120*rho**4*z**2 - 5*rho**6)
            elif m==0:
                P = 0
            elif m==1:
                P = -(1/64)*(64*z**6 - 240*rho**2*z**2 + 120*rho**4*z**2 - 5*rho**6)
            elif m==2:
                P = -0.125*rho*z*(48*z**4 - 80*rho**2*z**2 + 15*rho**4)
            elif m==3:
                P = -(3/16)*rho**2*(80*z**4 - 60*rho**2*z**2 + 3*rho**4)
            elif m==4:
                P = -2*rho**3*z*(10*z**2 - 3*rho**2)
            elif m==5:
                P = -5*rho**4*(12*z**2 - rho**2)
            elif m==6:
                P = -6*rho**5*z
            elif m==7:
                P = -rho**6
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        else:
            print("l too high")
            return None
        
        return P
    
    def PiZ(rho, z, l, m):
        
        if l==0:
            
            if m==-1:
                P = 0
            elif m==0:
                P = 1
            elif m==1:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==1:
            
            if m==-2:
                P = 0
            elif m==-1:
                P = rho
            elif m==0:
                P = z
            elif m==1:
                P = rho
            elif m==2:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
                
        elif l==2:
            
            if m==-3:
                P = 0
            elif m==-2:
                P = rho**2
            elif m==-1:
                P = 2*rho*z
            elif m==0:
                P = -0.5*rho**2 + z**2
            elif m==1:
                P = 2*rho*z
            elif m==2:
                P = rho**2
            elif m==3:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
        
        elif l==3:
            
            if m==-4:
                P = 0
            elif m==-3:
                P = rho**3
            elif m==-2:
                P = 3*rho**2*z**2
            elif m==-1:
                P = rho*(3*z**2 - 0.75*rho**2)
            elif m==0:
                P = 0.5*z*(2*z**2 - 3*rho**2)
            elif m==1:
                P = rho*(3*z**2 - 0.75*rho**2)
            elif m==2:
                P = 3*rho**2*z**2
            elif m==3:
                P = rho**3
            elif m==-4:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==4:
            
            if m==-5:
                P = 0
            elif m==-4:
                P = rho**4
            elif m==-3:
                P = 4*rho**3*z
            elif m==-2:
                P = 6*rho**2*z**2 - rho**4
            elif m==-1:
                P = 4*rho*z**3 - 3*rho**3*z
            elif m==0:
                P = 0.125*(8*z**4 - 24*rho**2*z**2 + 3*rho**4)
            elif m==1:
                P = 4*rho*z**3 - 3*rho**3*z
            elif m==2:
                P = 6*rho**2*z**2 - rho**4
            elif m==3:
                P = 4*rho**3*z
            elif m==4:
                P = rho**4
            elif m==5:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==5:
            
            if m==-6:
                P = 0
            elif m==-5:
                P = rho**5
            elif m==-4:
                P = 5*rho**4*z
            elif m==-3:
                P = 1.25*(8*rho**3*z**2 - rho**5)
            elif m==-2:
                P = 5*(2*rho**2*z**3 - rho**4*z)
            elif m==-1:
                P = (5/8)*(8*rho*z**4 - 12*rho**3*z**2 + rho**5)
            elif m==0:
                P = (1/8)*(8*z**5 - 40*rho**2*z**3 + 15*rho**4*z)
            elif m==1:
                P = (5/8)*(8*rho*z**4 - 12*rho**3*z**2 + rho**5)
            elif m==2:
                P = 5*(2*rho**2*z**3 - rho**4*z)
            elif m==3:
                P = 1.25*(8*rho**3*z**2 - rho**5)
            elif m==4:
                P = 5*rho**4*z
            elif m==5:
                P = rho**5
            elif m==6:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==6:
            
            if m==-7:
                P = 0
            elif m==-6:
                P = rho**6
            elif m==-5:
                P = 6*rho**5*z
            elif m==-4:
                P = 1.5*rho**4*(10*z**2 - rho**2)
            elif m==-3:
                P = 2.5*rho**3*z*(8*z**2 - 3*rho**2)
            elif m==-2:
                P = (15/16)*rho**2*(16*z**4 - 16*rho**2*z**2 + rho**4)
            elif m==-1:
                P = 0.75*rho*z*(8*z**4 - 20*rho**2*z**2 + 5*rho**4)
            elif m==0:
                P = (1/16)*(16*z**6 - 120*rho**2*z**4 + 90*rho**4*z**2 - 5*rho**6)
            elif m==1:
                P = 0.75*rho*z*(8*z**4 - 20*rho**2*z**2 + 5*rho**4)
            elif m==2:
                P = (15/16)*rho**2*(16*z**4 - 16*rho**2*z**2 + rho**4)
            elif m==3:
                P = 2.5*rho**3*z*(8*z**2 - 3*rho**2)
            elif m==4:
                P = 1.5*rho**4*(10*z**2 - rho**2)
            elif m==5:
                P = 6*rho**5*z
            elif m==6:
                P = rho**6
            elif m==7:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        else:
            print("l too high")
            return None
        
        return P
                   