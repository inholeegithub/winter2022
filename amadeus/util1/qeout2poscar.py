import numpy as np
import os
# Written by In-Ho Lee, KRISS, May 3, 2020.
def cart2drt(b,uec):
    vec=np.matmul(b,uec)
    for i in range(3):
        vec[i]=vec[i]-round(vec[i])
        if vec[i] < 0.:
           vec[i]=vec[i]+1.
    return vec
# Written by In-Ho Lee, KRISS, May 3, 2020.
def cart2drtunit(a1,a2,a3,uec):
    aa1=a1
    aa2=a2
    aa3=a3
    dvd=aa1[0]*aa2[1]*aa3[2]-aa1[1]*aa2[0]*aa3[2]-aa1[0]*aa2[2]*aa3[1]
    dvd=dvd+aa1[2]*aa2[0]*aa3[1]+aa1[1]*aa2[2]*aa3[0]-aa1[2]*aa2[1]*aa3[0]
    b=np.zeros((3,3))
    b[0,0]=-aa2[2]*aa3[1]+aa2[1]*aa3[2]
    b[1,0]= aa1[2]*aa3[1]-aa1[1]*aa3[2]
    b[2,0]=-aa1[2]*aa2[1]+aa1[1]*aa2[2]
    b[0,1]= aa2[2]*aa3[0]-aa2[0]*aa3[2]
    b[1,1]=-aa1[2]*aa3[0]+aa1[0]*aa3[2]
    b[2,1]= aa1[2]*aa2[0]-aa1[0]*aa2[2]
    b[0,2]=-aa2[1]*aa3[0]+aa2[0]*aa3[1]
    b[1,2]= aa1[1]*aa3[0]-aa1[0]*aa3[1]
    b[2,2]=-aa1[1]*aa2[0]+aa1[0]*aa2[1]
    b[:,:]=b[:,:]/dvd
    return cart2drt(b,uec)
# Written by In-Ho Lee, KRISS, May 3, 2020.
def dotproduct3(a,b):
    c=a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    return c
# Written by In-Ho Lee, KRISS, May 3, 2020.
def cross3(a,b):
    c=np.zeros(3)
    c[0]=a[1]*b[2]-a[2]*b[1]
    c[1]=a[2]*b[0]-a[0]*b[2]
    c[2]=a[0]*b[1]-a[1]*b[0]
    return c
# Written by In-Ho Lee, KRISS, May 3, 2020.
def meshijk0(dk,ga):
    gatmp=ga/(2.0*np.pi)
    gatmp=ga
    ka=int(round(gatmp/dk))+1
    if ka == 0: 
       ka=1
    dum=gatmp/float(ka)
    if dum >= dk :
       for i in range(15):
           ka=ka+1
           dum=gatmp/float(ka)
           if dum <= dk:
             break 
    return ka
# Written by In-Ho Lee, KRISS, May 3, 2020.
def write_kpoints0(dk,a1,a2,a3):
    iswitch=0
    if dk < 0.0 :
       dk=0.06
    if dk == 0.0 :
       dk=0.015
       iswitch=1
    b1=cross3(a2,a3)
    b2=cross3(a3,a1)
    b3=cross3(a1,a2)
    omega=np.abs(dotproduct3(b1,a1))
    b1[:]=b1[:]*(2.0*np.pi/omega) 
    b2[:]=b2[:]*(2.0*np.pi/omega) 
    b3[:]=b3[:]*(2.0*np.pi/omega)
    ga=np.sqrt(dotproduct3(b1,b1)) 
    gb=np.sqrt(dotproduct3(b2,b2)) 
    gc=np.sqrt(dotproduct3(b3,b3))
    ka=meshijk0(dk,ga)
    kb=meshijk0(dk,gb)
    kc=meshijk0(dk,gc)
    if ka > 15: 
       ka=15
    if kb > 15: 
       kb=15
    if kc > 15: 
       kc=15
    if iswitch == 1:
       ka=ka+1+ka/2
       kb=kb+1+kb/2
       kc=kc+1+kc/2
    return ka,kb,kc
# Written by In-Ho Lee, KRISS, May 3, 2020.
def atomicmass(symbl):
    arr = [None]*110
    symbl0 = [None]*110
    atma=1.0079
    arr[  1]=  1.0079   ; symbl0[  1]='H'
    arr[  2]=  4.0026   ; symbl0[  2]='He'
    arr[  3]=  6.941    ; symbl0[  3]='Li'
    arr[  4]=  9.0122   ; symbl0[  4]='Be'
    arr[  5]=  10.811   ; symbl0[  5]='B'
    arr[  6]=  12.0107  ; symbl0[  6]='C'
    arr[  7]=  14.0067  ; symbl0[  7]='N'
    arr[  8]=  15.9994  ; symbl0[  8]='O'
    arr[  9]=  18.9984  ; symbl0[  9]='F'
    arr[ 10]=  20.1797  ; symbl0[ 10]='Ne'
    arr[ 11]=  22.9897  ; symbl0[ 11]='Na'
    arr[ 12]=  24.305   ; symbl0[ 12]='Mg'
    arr[ 13]=  26.9815  ; symbl0[ 13]='Al'
    arr[ 14]=  28.0855  ; symbl0[ 14]='Si'
    arr[ 15]=  30.9738  ; symbl0[ 15]='P'
    arr[ 16]=  32.065   ; symbl0[ 16]='S'
    arr[ 17]=  35.453   ; symbl0[ 17]='Cl'
    arr[ 18]=  39.948   ; symbl0[ 18]='Ar'
    arr[ 19]=  39.0983  ; symbl0[ 19]='K'
    arr[ 20]=  40.078   ; symbl0[ 20]='Ca'
    arr[ 21]=  44.9559  ; symbl0[ 21]='Sc'
    arr[ 22]=  47.867   ; symbl0[ 22]='Ti'
    arr[ 23]=  50.9415  ; symbl0[ 23]='V'
    arr[ 24]=  51.9961  ; symbl0[ 24]='Cr'
    arr[ 25]=  54.938   ; symbl0[ 25]='Mn'
    arr[ 26]=  55.845   ; symbl0[ 26]='Fe'
    arr[ 27]=  58.9332  ; symbl0[ 27]='Co'
    arr[ 28]=  58.6934  ; symbl0[ 28]='Ni'
    arr[ 29]=  63.546   ; symbl0[ 29]='Cu'
    arr[ 30]=  65.39    ; symbl0[ 30]='Zn'
    arr[ 31]=  69.723   ; symbl0[ 31]='Ga'
    arr[ 32]=  72.64    ; symbl0[ 32]='Ge'
    arr[ 33]=  74.9216  ; symbl0[ 33]='As'
    arr[ 34]=  78.96    ; symbl0[ 34]='Se'
    arr[ 35]=  79.904   ; symbl0[ 35]='Br'
    arr[ 36]=  83.8     ; symbl0[ 36]='Kr'
    arr[ 37]=  85.4678  ; symbl0[ 37]='Rb'
    arr[ 38]=  87.62    ; symbl0[ 38]='Sr'
    arr[ 39]=  88.9059  ; symbl0[ 39]='Y'
    arr[ 40]=  91.224   ; symbl0[ 40]='Zr'
    arr[ 41]=  92.9064  ; symbl0[ 41]='Nb'
    arr[ 42]=  95.94    ; symbl0[ 42]='Mo'
    arr[ 43]=  98.00    ; symbl0[ 43]='Tc'
    arr[ 44]=  101.07   ; symbl0[ 44]='Ru'
    arr[ 45]=  102.9055 ; symbl0[ 45]='Rh'
    arr[ 46]=  106.42   ; symbl0[ 46]='Pd'
    arr[ 47]=  107.8682 ; symbl0[ 47]='Ag'
    arr[ 48]=  112.411  ; symbl0[ 48]='Cd'
    arr[ 49]=  114.818  ; symbl0[ 49]='In'
    arr[ 50]=  118.71   ; symbl0[ 50]='Sn'
    arr[ 51]=  121.76   ; symbl0[ 51]='Sb'
    arr[ 52]=  127.6    ; symbl0[ 52]='Te'
    arr[ 53]=  126.9045 ; symbl0[ 53]='I'
    arr[ 54]=  131.293  ; symbl0[ 54]='Xe'
    arr[ 55]=  132.9055 ; symbl0[ 55]='Cs'
    arr[ 56]=  137.327  ; symbl0[ 56]='Ba'
    arr[ 57]=  138.9055 ; symbl0[ 57]='La'
    arr[ 58]=  140.116  ; symbl0[ 58]='Ce'
    arr[ 59]=  140.9077 ; symbl0[ 59]='Pr'
    arr[ 60]=  144.24   ; symbl0[ 60]='Nd'
    arr[ 61]=  145.00   ; symbl0[ 61]='Pm'
    arr[ 62]=  150.36   ; symbl0[ 62]='Sm'
    arr[ 63]=  151.964  ; symbl0[ 63]='Eu'
    arr[ 64]=  157.25   ; symbl0[ 64]='Gd'
    arr[ 65]=  158.9253 ; symbl0[ 65]='Tb'
    arr[ 66]=  162.5    ; symbl0[ 66]='Dy'
    arr[ 67]=  164.9303 ; symbl0[ 67]='Ho'
    arr[ 68]=  167.259  ; symbl0[ 68]='Er'
    arr[ 69]=  168.9342 ; symbl0[ 69]='Tm'
    arr[ 70]=  173.04   ; symbl0[ 70]='Yb'
    arr[ 71]=  174.967  ; symbl0[ 71]='Lu'
    arr[ 72]=  178.49   ; symbl0[ 72]='Hf'
    arr[ 73]=  180.9479 ; symbl0[ 73]='Ta'
    arr[ 74]=  183.84   ; symbl0[ 74]='W'
    arr[ 75]=  186.207  ; symbl0[ 75]='Re'
    arr[ 76]=  190.23   ; symbl0[ 76]='Os'
    arr[ 77]=  192.217  ; symbl0[ 77]='Ir'
    arr[ 78]=  195.078  ; symbl0[ 78]='Pt'
    arr[ 79]=  196.9665 ; symbl0[ 79]='Au'
    arr[ 80]=  200.59   ; symbl0[ 80]='Hg'
    arr[ 81]=  204.3833 ; symbl0[ 81]='Tl'
    arr[ 82]=  207.2    ; symbl0[ 82]='Pb'
    arr[ 83]=  208.9804 ; symbl0[ 83]='Bi'
    arr[ 84]=  209.     ; symbl0[ 84]='Po'
    arr[ 85]=  210.     ; symbl0[ 85]='At'
    arr[ 86]=  222.     ; symbl0[ 86]='Rn'
    arr[ 87]=  223.     ; symbl0[ 87]='Fr'
    arr[ 88]=  226.     ; symbl0[ 88]='Ra'
    arr[ 89]=  227.     ; symbl0[ 89]='Ac'
    arr[ 90]=  232.0381 ; symbl0[ 90]='Th'
    arr[ 91]=  231.0359 ; symbl0[ 91]='Pa'
    arr[ 92]=  238.0289 ; symbl0[ 92]='U'
    arr[ 93]=  237.     ; symbl0[ 93]='Np'
    arr[ 94]=  244.     ; symbl0[ 94]='Pu'
    arr[ 95]=  243.     ; symbl0[ 95]='Am'
    arr[ 96]=  247.     ; symbl0[ 96]='Cm'
    arr[ 97]=  247.     ; symbl0[ 97]='Bk'
    arr[ 98]=  251.     ; symbl0[ 98]='Cf'
    arr[ 99]=  252.     ; symbl0[ 99]='Es'
    arr[100]=  257.     ; symbl0[100]='Fm'
    arr[101]=  258.     ; symbl0[101]='Md'
    arr[102]=  259.     ; symbl0[102]='No'
    arr[103]=  262.     ; symbl0[103]='Lr'
    arr[104]=  261.     ; symbl0[104]='Rf'
    arr[105]=  262.     ; symbl0[105]='Db'
    arr[106]=  266.     ; symbl0[106]='Sg'
    arr[107]=  264.     ; symbl0[107]='Bh'
    arr[108]=  277.     ; symbl0[108]='Hs'
    arr[109]=  268.     ; symbl0[109]='Mt'
    for i in range(1,110):
        if symbl == symbl0[i] :
           atma=arr[i]
           break
    return atma
#   https://en.wikipedia.org/wiki/Covalent_radius
# Written by In-Ho Lee, KRISS, May 3, 2020.
def covlaentrr(symbol):
    rr=1.
    if symbol == 'H'  : 
       rr=0.31e0
    if symbol == 'He' : 
       rr=0.28e0
    if symbol == 'Li' : 
       rr=1.28e0
    if symbol == 'Be' : 
       rr=0.96e0
    if symbol == 'B'  : 
       rr=0.84e0
    if symbol == 'C'  : 
       rr=0.69e0
    if symbol == 'N'  : 
       rr=0.71e0
    if symbol == 'O'  : 
       rr=0.66e0
    if symbol == 'F'  : 
       rr=0.57e0
    if symbol == 'Ne' : 
       rr=0.58e0
    if symbol == 'Na' : 
       rr=1.66e0
    if symbol == 'Mg' : 
       rr=1.41e0
    if symbol == 'Al' : 
       rr=1.21e0
    if symbol == 'Si' : 
       rr=1.11e0
    if symbol == 'P'  : 
       rr=1.07e0
    if symbol == 'S'  : 
       rr=1.05e0
    if symbol == 'Cl' : 
       rr=1.02e0
    if symbol == 'Ar' : 
       rr=1.06e0
    if symbol == 'K'  : 
       rr=2.03e0
    if symbol == 'Ca' : 
       rr=1.76e0
    if symbol == 'Sc' : 
       rr=1.70e0
    if symbol == 'Ti' : 
       rr=1.60e0
    if symbol == 'V'  : 
       rr=1.53e0
    if symbol == 'Cr' : 
       rr=1.39e0
    if symbol == 'Mn' : 
       rr=1.39e0
    if symbol == 'Fe' : 
       rr=1.32e0
    if symbol == 'Co' : 
       rr=1.26e0
    if symbol == 'Ni' : 
       rr=1.24e0
    if symbol == 'Cu' : 
       rr=1.32e0
    if symbol == 'Zn' : 
       rr=1.22e0
    if symbol == 'Ga' : 
       rr=1.22e0
    if symbol == 'Ge' : 
       rr=1.20e0
    if symbol == 'As' : 
       rr=1.19e0
    if symbol == 'Se' : 
       rr=1.20e0
    if symbol == 'Br' : 
       rr=1.20e0
    if symbol == 'Kr' : 
       rr=1.16e0
    if symbol == 'Rb' : 
       rr=2.20e0
    if symbol == 'Sr' : 
       rr=1.95e0
    if symbol == 'Y'  : 
       rr=1.90e0
    if symbol == 'Zr' : 
       rr=1.75e0
    if symbol == 'Nb' : 
       rr=1.64e0
    if symbol == 'Mo' : 
       rr=1.54e0
    if symbol == 'Tc' : 
       rr=1.47e0
    if symbol == 'Ru' : 
       rr=1.46e0
    if symbol == 'Rh' : 
       rr=1.42e0
    if symbol == 'Pd' : 
       rr=1.39e0
    if symbol == 'Ag' : 
       rr=1.45e0
    if symbol == 'Cd' : 
       rr=1.44e0
    if symbol == 'In' : 
       rr=1.42e0
    if symbol == 'Sn' : 
       rr=1.39e0
    if symbol == 'Sb' : 
       rr=1.39e0
    if symbol == 'Te' : 
       rr=1.38e0
    if symbol == 'I'  : 
       rr=1.39e0
    if symbol == 'Xe' : 
       rr=1.40e0
    if symbol == 'Cs' : 
       rr=2.44e0
    if symbol == 'Ba' : 
       rr=2.15e0
    if symbol == 'La' : 
       rr=2.07e0
    if symbol == 'Lu' : 
       rr=1.87e0
    if symbol == 'Hf' : 
       rr=1.75e0
    if symbol == 'Ta' : 
       rr=1.70e0
    if symbol == 'W'  : 
       rr=1.62e0
    if symbol == 'Re' : 
       rr=1.51e0
    if symbol == 'Os' : 
       rr=1.44e0
    if symbol == 'Ir' : 
       rr=1.41e0
    if symbol == 'Pt' : 
       rr=1.36e0
    if symbol == 'Au' : 
       rr=1.36e0
    if symbol == 'Hg' : 
       rr=1.32e0
    if symbol == 'Tl' : 
       rr=1.45e0
    if symbol == 'Pb' : 
       rr=1.46e0
    if symbol == 'Bi' : 
       rr=1.48e0
    if symbol == 'Po' : 
       rr=1.40e0
    if symbol == 'At' : 
       rr=1.50e0
    if symbol == 'Rn' : 
       rr=1.50e0
    if symbol == 'Fr' : 
       rr=2.60e0
    if symbol == 'Ra' : 
       rr=2.21e0
    if symbol == 'Ce' : 
       rr=2.04e0
    if symbol == 'Pr' : 
       rr=2.03e0
    if symbol == 'Nd' : 
       rr=2.01e0
    if symbol == 'Pm' : 
       rr=1.99e0
    if symbol == 'Sm' : 
       rr=1.98e0
    if symbol == 'Eu' : 
       rr=1.98e0
    if symbol == 'Gd' : 
       rr=1.96e0
    if symbol == 'Tb' : 
       rr=1.94e0
    if symbol == 'Dy' : 
       rr=1.92e0
    if symbol == 'Ho' : 
       rr=1.92e0
    if symbol == 'Er' : 
       rr=1.89e0
    if symbol == 'Tm' : 
       rr=1.90e0
    if symbol == 'Yb' : 
       rr=1.87e0
    if symbol == 'Ac' : 
       rr=2.15e0
    if symbol == 'Th' : 
       rr=2.06e0
    if symbol == 'Pa' : 
       rr=2.00e0
    if symbol == 'U'  : 
       rr=1.96e0
    if symbol == 'Np' : 
       rr=1.90e0
    if symbol == 'Pu' : 
       rr=1.87e0
    if symbol == 'Am' : 
       rr=1.80e0
    if symbol == 'Cm' : 
       rr=1.69e0
    tmp=rr
    return tmp
# Written by In-Ho Lee, KRISS, May 3, 2020.
ipseudo=2
afile=open('vc-relax.out','r')
lalat= False  
lsites= False    
lcart= False 
a1=np.zeros(3)
a2=np.zeros(3)
a3=np.zeros(3)
etot=0.
for line in afile:
    if len(line.split()) == 6:
       if line.split()[1]  == 'total':
          if line.split()[2]  == 'energy':
             if line.split()[5]  == 'Ry':
                etot=float(line.split()[4])
    if len(line.split()) == 5:
       if line.split()[0]  == 'Final':
          if line.split()[1]  == 'enthalpy':
             if line.split()[4]  == 'Ry':
                enth=float(line.split()[3])
    if len(line.split()) == 6:
       if line.split()[0]  == 'number':
          if line.split()[1]  == 'of':
             if line.split()[2]  == 'atomic':
                if line.split()[3]  == 'types':
                   nspecies=int(line.split()[5])
                   symbl=['  ' for j in range(nspecies)]
                   nelements=[ 0 for j in range(nspecies)]
    if len(line.split()) == 5:
       if line.split()[0]  == 'number':
          if line.split()[1]  == 'of':
             if line.split()[2]  == 'atoms/cell':
                natot=int(line.split()[4])
                drt=np.zeros((natot,3))
                xyz=np.zeros((natot,3))
                isymbl=['  ' for j in range(natot)]
    if len(line.split()) == 6:
       if line.split()[0]  == 'lattice':
          if line.split()[1]  == 'parameter':
             if line.split()[2]  == '(alat)':
                alat=float(line.split()[4])
                if line.split()[5]  == 'a.u.':
                   alat=alat*0.529177
    if len(line.split()) == 6:
       if line.split()[0]  == 'site' and line.split()[1] == 'n.' and line.split()[2] == 'atom':
          lsites=True
       if line.split()[4]  == '(alat' and line.split()[5] == 'units)' : 
          lalat=True
    if len(line.split()) == 8:
       if line.split()[0]  == 'crystal' and line.split()[7] == 'alat)' and line.split()[1] == 'axes:':
          lalat=True
    if len(line.split()) == 7:
       if line.split()[0]  == 'a(1)' and line.split()[1] == '=' and line.split()[2] == '(':
          if line.split()[6]  == ')' :
             a1[0]=float(line.split()[3])
             a1[1]=float(line.split()[4])
             a1[2]=float(line.split()[5])
             if lalat :
                a1[:]=a1[:]*alat
    if len(line.split()) == 7:
       if line.split()[0]  == 'a(2)' and line.split()[1] == '=' and line.split()[2] == '(':
          if line.split()[6]  == ')' :
             a2[0]=float(line.split()[3])
             a2[1]=float(line.split()[4])
             a2[2]=float(line.split()[5])
             if lalat :
                a2[:]=a2[:]*alat
    if len(line.split()) == 7:
       if line.split()[0]  == 'a(3)' and line.split()[1] == '=' and line.split()[2] == '(':
          if line.split()[6]  == ')' :
             a3[0]=float(line.split()[3])
             a3[1]=float(line.split()[4])
             a3[2]=float(line.split()[5])
             if lalat :
                a3[:]=a3[:]*alat
                lalat=False
    if len(line.split()) == 8:
       if line.split()[0] ==  'PseudoPot.' :
          if line.split()[1] ==  '#' :
             i=int(line.split()[2])-1
             symbl[i]=line.split()[4]
    if len(line.split()) == 2:
       if line.split()[0] ==  'Cartesian' and line.split()[1] == 'axes': 
          lcart=True
    if len(line.split()) == 10 and lsites :
       if line.split()[9] == ')' and line.split()[2] == 'tau(' :
          if line.split()[5] == '(':
             i=int(line.split()[0])-1
             isymbl[i]=line.split()[1]
             drt[i,0]=float(line.split()[6])
             drt[i,1]=float(line.split()[7])
             drt[i,2]=float(line.split()[8])
             if lalat :
                drt[i,:]=drt[i,:]*alat
             xyz[i,:]=drt[i,:]
             if natot-1 == i:
                lalat=False
    if True:
       if len(line.split()) > 1 and line.split()[0] == 'CELL_PARAMETERS' :
          line=afile.readline()
          a1[0]=float(line.split()[0])
          a1[1]=float(line.split()[1])
          a1[2]=float(line.split()[2])
          line=afile.readline()
          a2[0]=float(line.split()[0])
          a2[1]=float(line.split()[1])
          a2[2]=float(line.split()[2])
          line=afile.readline()
          a3[0]=float(line.split()[0])
          a3[1]=float(line.split()[1])
          a3[2]=float(line.split()[2])
       if len(line.split()) > 1 and line.split()[0] == 'ATOMIC_POSITIONS' :
          for j0 in range(natot):
              line=afile.readline()
              isymbl[j0]=line.split()[0]
              drt[j0,0]=float(line.split()[1])
              drt[j0,1]=float(line.split()[2])
              drt[j0,2]=float(line.split()[3])
afile.close()
if lcart: 
   for k in range(natot):
       drt[k,:]=cart2drtunit(a1,a2,a3,drt[k,:])
       drt[k,0]=drt[k,0]-round(drt[k,0])
       drt[k,1]=drt[k,1]-round(drt[k,1])
       drt[k,2]=drt[k,2]-round(drt[k,2])
       if drt[k,0] < 0.:
          drt[k,0]=drt[k,0]+1.
       if drt[k,1] < 0.:
          drt[k,1]=drt[k,1]+1.
       if drt[k,2] < 0.:
          drt[k,2]=drt[k,2]+1.
for i in range(nspecies):
    j=0
    for k in range(natot):
        if symbl[i] == isymbl[k] :
           j=j+1
    nelements[i]=j
etot=etot*13.6058
afile=open('OUTCAR','w')
astring='  energy  without entropy= '+f'{etot : 24.10f}' +' energy(sigma->0) = '+f'{etot : 24.10f}'+'\n'
afile.write(astring)
afile.close()
lposcar0=False
if os.path.isfile("POSCAR_"):
   lposcar0=True
   afile=open('POSCAR_','r')
   for line in afile:
       iidd=int(line.split()[0])
   afile.close()
if False:
   bfile=open('POSCAR','w')
   if lposcar0:
      bfile.write(str(iidd)+'\n')
   else:
      bfile.write('qe2poscar'+'\n')
   bfile.write('1.' +'\n')
   bfile.write(f"{a1[0] :22.12f}"+   f"{a1[1] :22.12f}"+f"{a1[2] :22.12f}"+'\n')
   bfile.write(f"{a2[0] :22.12f}"+   f"{a2[1] :22.12f}"+f"{a2[2] :22.12f}"+'\n')
   bfile.write(f"{a3[0] :22.12f}"+   f"{a3[1] :22.12f}"+f"{a3[2] :22.12f}"+'\n')
   astring=' '
   bstring=' '
   for i in range(nspecies):
       astring=astring+' '+symbl[i]
       bstring=bstring+' '+str(nelements[i])
   astring=astring+'\n'
   bstring=bstring+'\n'
   bfile.write(astring)
   bfile.write(bstring)
   bfile.write( 'direct' + '\n')
   for i in range(nspecies):
       for k in range(natot):
           if symbl[i] == isymbl[k] :
              bfile.write(f"{drt[k,0] :18.12f}"+f"{drt[k,1] :18.12f}"+f"{drt[k,2] :18.12f}"+'\n')
   bfile.close()
if True:
   bfile=open('CONTCAR','w')
   if lposcar0:
      bfile.write(str(iidd)+'\n')
   else:
      bfile.write('qe2poscar'+'\n')
   bfile.write('1.' +'\n')
   bfile.write(f"{a1[0] :22.12f}"+   f"{a1[1] :22.12f}"+f"{a1[2] :22.12f}"+'\n')
   bfile.write(f"{a2[0] :22.12f}"+   f"{a2[1] :22.12f}"+f"{a2[2] :22.12f}"+'\n')
   bfile.write(f"{a3[0] :22.12f}"+   f"{a3[1] :22.12f}"+f"{a3[2] :22.12f}"+'\n')
   astring=' '
   bstring=' '
   for i in range(nspecies):
       astring=astring+' '+symbl[i]
       bstring=bstring+' '+str(nelements[i])
   astring=astring+'\n'
   bstring=bstring+'\n'
   bfile.write(astring)
   bfile.write(bstring)
   bfile.write( 'direct' + '\n')
   for i in range(nspecies):
       for k in range(natot):
           if symbl[i] == isymbl[k] :
              bfile.write(f"{drt[k,0] :18.12f}"+f"{drt[k,1] :18.12f}"+f"{drt[k,2] :18.12f}"+'\n')
   bfile.close()
if True:
   bfile=open('scf.in_','w')
   bfile.write('&CONTROL'+  '\n')
   bfile.write('prefix= '+ "'rlxallz4'"+  '\n')
   if True:
      bfile.write('wf_collect= .true. '+ '\n')
   bfile.write('calculation= '+ "'scf'"+  '\n') 
   if ipseudo == 1:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/PSEUDOPOTENTIALS_NC/'"+  '\n')
   if ipseudo == 2:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-sr-04_pbe_standard_upf/'"+  '\n')
   if ipseudo == 3:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-fr-04_pbe_standard_upf/'"+  '\n')
   bfile.write('outdir= '+ "'./'"+  '\n')
   bfile.write('tstress= .true.'+  '\n')
   bfile.write('tprnfor= .true.'+  '\n')
   bfile.write('/'+  '\n')
   bfile.write('&SYSTEM'+  '\n')
   bfile.write('nat= '+str(natot)+  '\n')
   bfile.write('ntyp= '+str(nspecies)+  '\n')
   bfile.write('ecutwfc= '+ '100.0'+  '\n')
   bfile.write('occupations= '+ "'smearing'"+  '\n')
   bfile.write('smearing= '+ "'mp'"+  '\n')
   bfile.write('degauss= 0.01'+  '\n')
   bfile.write('ibrav= '+ '0'+  '\n')
   if True:
      bfile.write('use_all_frac= .true.'+ '\n')
      bfile.write('lspinorb= .true.'+ '\n')
      bfile.write('noncolin= .true.'+ '\n')
   bfile.write('/'+  '\n')
   bfile.write('&ELECTRONS'+  '\n')
   bfile.write('conv_thr= 1.0d-8'+  '\n')
   bfile.write('mixing_beta=0.3'+  '\n')
   bfile.write('/'+  '\n')
   bfile.write(''+  '\n')
   bfile.write('ATOMIC_SPECIES'+  '\n')
   if ipseudo == 2 or ipseudo == 3:
      for i in range(nspecies):
          bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.upf'+ '\n')
   if ipseudo == 1:
      for i in range(nspecies):
          if symbl[i] == 'B' or symbl[i] == 'C' or symbl[i] == 'N' or symbl[i] == 'O' :
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-nc.UPF'+ '\n')
          else:    
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-n-nc.UPF'+ '\n')
   bfile.write(''+  '\n')
#  bfile.write('ATOMIC_POSITIONS {crystal}'+  '\n')
   bfile.write('ATOMIC_POSITIONS crystal'+  '\n')
   for i in range(nspecies):
       for k in range(natot):
           if isymbl[k] == symbl[i]:
               bfile.write( isymbl[k]+' '+f"{drt[k,0] :20.12f}"+f"{drt[k,1] :20.12f}"+f"{drt[k,2] :20.12f}"+  '\n')
   bfile.write(''+ '\n')
   bfile.write('CELL_PARAMETERS {angstrom}'+ '\n')
   bfile.write(f"{a1[0] :20.12f}"+   f"{a1[1] :20.12f}"+f"{a1[2] :20.12f}"+ '\n')
   bfile.write(f"{a2[0] :20.12f}"+   f"{a2[1] :20.12f}"+f"{a2[2] :20.12f}"+ '\n')
   bfile.write(f"{a3[0] :20.12f}"+   f"{a3[1] :20.12f}"+f"{a3[2] :20.12f}"+ '\n')
   bfile.write(''+ '\n')
   bfile.write('K_POINTS {automatic}'+ '\n')
   dk=0.3
   dk=0.2
   ka,kb,kc=write_kpoints0(dk,a1,a2,a3)
   bfile.write(f"{ka : 4d}"+f"{kb : 4d}"+f"{kc : 4d}"+f"{0 : 4d} "+f"{0 : 2d} "+f"{0 : 2d} "+ '\n')
   bfile.close()
if True:
   bfile=open('vc-relax.in_','w')
   bfile.write('&CONTROL'+ '\n')
   bfile.write('prefix= '+ "'rlxallz4'"+ '\n')
   if False:
      bfile.write('wf_collect= .true. '+ '\n')
   bfile.write('calculation= '+ "'vc-relax'"+ '\n')
   if ipseudo == 1:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/PSEUDOPOTENTIALS_NC/'"+ '\n')
   if ipseudo == 2:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-sr-04_pbe_standard_upf/'"+  '\n')
   if ipseudo == 3:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-fr-04_pbe_standard_upf/'"+  '\n')
   bfile.write('outdir= '+  "'./'"+ '\n')
   bfile.write('tstress= .true.'+  '\n')
   bfile.write('tprnfor= .true.'+  '\n')
   bfile.write('/'+ '\n')
   bfile.write('&SYSTEM'+ '\n')
   bfile.write('nat= '+str(natot)+ '\n')
   bfile.write('ntyp= '+str(nspecies)+ '\n')
   bfile.write('ecutwfc= '+ '100.0'+ '\n')
   bfile.write('occupations= '+ "'smearing'"+ '\n')
   bfile.write('smearing= '+ "'mp'"+ '\n')
   bfile.write('degauss= 0.01'+ '\n')
   bfile.write('ibrav= '+ '0'+ '\n')
   if False:
      bfile.write('use_all_frac= .true.'+ '\n')
      bfile.write('lspinorb= .true.'+ '\n')
      bfile.write('noncolin= .true.'+ '\n')
   bfile.write('/'+ '\n')
   bfile.write('&ELECTRONS'+ '\n')
   bfile.write('conv_thr= 1.0d-8'+ '\n')
   bfile.write('/'+ '\n')
   bfile.write('&IONS'+ '\n')
   bfile.write('ion_dynamics= '+ "'bfgs'"+ '\n')
   bfile.write('/'+ '\n')
   bfile.write('&CELL'+ '\n')
   bfile.write('cell_dynamics= '+"'bfgs'"+ '\n')
   bfile.write('press= 0.0'+ '\n')
   bfile.write('/'+ '\n')
   bfile.write(''+ '\n')
   bfile.write('ATOMIC_SPECIES'+ '\n')
   if ipseudo ==2 or ipseudo ==3:
      for i in range(nspecies):
          bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.upf'+ '\n')
   if ipseudo ==1:
      for i in range(nspecies):
          if symbl[i] == 'B' or symbl[i] == 'C' or symbl[i] == 'N' or symbl[i] == 'O' :
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-nc.UPF'+ '\n')
          else:    
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-n-nc.UPF'+ '\n')
   bfile.write(''+ '\n')
#  bfile.write('ATOMIC_POSITIONS {crystal}'+ '\n')
   bfile.write('ATOMIC_POSITIONS crystal'+ '\n')
   for i in range(nspecies):
       for k in range(natot):
           if isymbl[k] == symbl[i]:
               bfile.write( isymbl[k]+' '+f"{drt[k,0] :20.12f}"+f"{drt[k,1] :20.12f}"+f"{drt[k,2] :20.12f}"+ '\n')
   bfile.write(''+ '\n')
   bfile.write('CELL_PARAMETERS {angstrom}'+ '\n')
   bfile.write(f"{a1[0] :20.12f}"+   f"{a1[1] :20.12f}"+f"{a1[2] :20.12f}"+ '\n')
   bfile.write(f"{a2[0] :20.12f}"+   f"{a2[1] :20.12f}"+f"{a2[2] :20.12f}"+ '\n')
   bfile.write(f"{a3[0] :20.12f}"+   f"{a3[1] :20.12f}"+f"{a3[2] :20.12f}"+ '\n')
   bfile.write(''+ '\n')
   bfile.write('K_POINTS {automatic}'+ '\n')
   dk=0.3
   dk=0.2
   ka,kb,kc=write_kpoints0(dk,a1,a2,a3)
   bfile.write(f"{ka : 4d}"+f"{kb : 4d}"+f"{kc : 4d}"+f"{0 : 4d} "+f"{0 : 2d} "+f"{0 : 2d} "+ '\n')
   bfile.close()
if True:
   bfile=open('rep.in_','w')
   bfile.write('&CONTROL'+ '\n')
   bfile.write('prefix= '+ "'rlxallz4'"+ '\n')
   if True:
      bfile.write('wf_collect= .true. '+ '\n')
   bfile.write('calculation= '+ "'bands'"+ '\n')
   if ipseudo == 1:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/PSEUDOPOTENTIALS_NC/'"+ '\n')
   if ipseudo == 2:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-sr-04_pbe_standard_upf/'"+ '\n')
   if ipseudo == 3:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-fr-04_pbe_standard_upf/'"+ '\n')
   bfile.write('outdir= '+  "'./'"+ '\n')
   bfile.write('tstress= .true.'+  '\n')
   bfile.write('tprnfor= .true.'+  '\n')
   bfile.write('/'+ '\n')
   bfile.write('&SYSTEM'+ '\n')
   bfile.write('nat= '+str(natot)+ '\n')
   bfile.write('ntyp= '+str(nspecies)+ '\n')
   bfile.write('ecutwfc= '+ '100.0'+ '\n')
   bfile.write('occupations= '+ "'smearing'"+ '\n')
   bfile.write('smearing= '+ "'mp'"+ '\n')
   bfile.write('degauss= 0.01'+ '\n')
   bfile.write('ibrav= '+ '0'+ '\n')
   if True:
      bfile.write('use_all_frac= .true.'+ '\n')
      bfile.write('lspinorb= .true.'+ '\n')
      bfile.write('noncolin= .true.'+ '\n')
   bfile.write('/'+ '\n')
   bfile.write('&ELECTRONS'+ '\n')
   bfile.write('conv_thr= 1.0d-8'+ '\n')
   bfile.write('/'+ '\n')
   bfile.write(''+ '\n')
   bfile.write('ATOMIC_SPECIES'+ '\n')
   if ipseudo == 2 or ipseudo == 3:
      for i in range(nspecies):
          bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.upf'+ '\n')
   if ipseudo == 1:
      for i in range(nspecies):
          if symbl[i] == 'B' or symbl[i] == 'C' or symbl[i] == 'N' or symbl[i] == 'O' :
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-nc.UPF'+ '\n')
          else:    
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-n-nc.UPF'+ '\n')
   bfile.write(''+ '\n')
#  bfile.write('ATOMIC_POSITIONS {crystal}'+ '\n')
   bfile.write('ATOMIC_POSITIONS crystal'+ '\n')
   for i in range(nspecies):
       for k in range(natot):
           if isymbl[k] == symbl[i]:
               bfile.write( isymbl[k]+' '+f"{drt[k,0] :20.12f}"+f"{drt[k,1] :20.12f}"+f"{drt[k,2] :20.12f}"+ '\n')
   bfile.write(''+ '\n')
   bfile.write('CELL_PARAMETERS {angstrom}'+ '\n')
   bfile.write(f"{a1[0] :20.12f}"+   f"{a1[1] :20.12f}"+f"{a1[2] :20.12f}"+ '\n')
   bfile.write(f"{a2[0] :20.12f}"+   f"{a2[1] :20.12f}"+f"{a2[2] :20.12f}"+ '\n')
   bfile.write(f"{a3[0] :20.12f}"+   f"{a3[1] :20.12f}"+f"{a3[2] :20.12f}"+ '\n')
   bfile.write(''+ '\n')
   bfile.write('K_POINTS crystal'+ '\n')
   bfile.write('8'+ '\n')
   bfile.write( '0.0	0.0	0.0	1'+ '\n')
   bfile.write( '0.5	0.0	0.0	1'+ '\n')
   bfile.write( '0.0	0.5	0.0	1'+ '\n')
   bfile.write( '0.0	0.0	0.5	1'+ '\n')
   bfile.write( '0.0	0.5	0.5	1'+ '\n')
   bfile.write( '0.5	0.0	0.5	1'+ '\n')
   bfile.write( '0.5	0.5	0.0	1'+ '\n')
   bfile.write( '0.5	0.5	0.5	1'+ '\n')
   bfile.close()
if True:
   bfile=open('scf.in_tet','w')
   bfile.write('&CONTROL'+  '\n')
   bfile.write('prefix= '+ "'rlxallz4'"+  '\n')
   if False:
      bfile.write('wf_collect= .true. '+ '\n')
   bfile.write('calculation= '+ "'scf'"+  '\n') 
   if ipseudo == 1:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/PSEUDOPOTENTIALS_NC/'"+  '\n')
   if ipseudo == 2:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-sr-04_pbe_standard_upf/'"+  '\n')
   if ipseudo == 3:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-fr-04_pbe_standard_upf/'"+  '\n')
   bfile.write('outdir= '+ "'./'"+  '\n')
   bfile.write('tstress= .true.'+  '\n')
   bfile.write('tprnfor= .true.'+  '\n')
   bfile.write('/'+  '\n')
   bfile.write('&SYSTEM'+  '\n')
   bfile.write('nat= '+str(natot)+  '\n')
   bfile.write('ntyp= '+str(nspecies)+  '\n')
   bfile.write('ecutwfc= '+ '100.0'+  '\n')
#  bfile.write('occupations= '+ "'smearing'"+  '\n')
   bfile.write('occupations ='+ "'tetrahedra_opt'"+ '\n' ) 
#  bfile.write('smearing= '+ "'mp'"+  '\n')
#  bfile.write('degauss= 0.01'+  '\n')
   bfile.write('ibrav= '+ '0'+  '\n')
   if False:
      bfile.write('use_all_frac= .true.'+ '\n')
      bfile.write('lspinorb= .true.'+ '\n')
      bfile.write('noncolin= .true.'+ '\n')
   bfile.write('/'+  '\n')
   bfile.write('&ELECTRONS'+  '\n')
#  bfile.write('conv_thr= 1.0d-8'+  '\n')
#  bfile.write('mixing_beta=0.3'+  '\n')
   bfile.write('/'+  '\n')
   bfile.write(''+  '\n')
   bfile.write('ATOMIC_SPECIES'+  '\n')
   if ipseudo == 2 or ipseudo == 3:
      for i in range(nspecies):
          bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.upf'+ '\n')
   if ipseudo == 1:
      for i in range(nspecies):
          if symbl[i] == 'B' or symbl[i] == 'C' or symbl[i] == 'N' or symbl[i] == 'O' :
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-nc.UPF'+ '\n')
          else:    
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-n-nc.UPF'+ '\n')
   bfile.write(''+  '\n')
#  bfile.write('ATOMIC_POSITIONS {crystal}'+  '\n')
   bfile.write('ATOMIC_POSITIONS crystal'+  '\n')
   for i in range(nspecies):
       for k in range(natot):
           if isymbl[k] == symbl[i]:
               bfile.write( isymbl[k]+' '+f"{drt[k,0] :20.12f}"+f"{drt[k,1] :20.12f}"+f"{drt[k,2] :20.12f}"+  '\n')
   bfile.write(''+ '\n')
   bfile.write('CELL_PARAMETERS {angstrom}'+ '\n')
   bfile.write(f"{a1[0] :20.12f}"+   f"{a1[1] :20.12f}"+f"{a1[2] :20.12f}"+ '\n')
   bfile.write(f"{a2[0] :20.12f}"+   f"{a2[1] :20.12f}"+f"{a2[2] :20.12f}"+ '\n')
   bfile.write(f"{a3[0] :20.12f}"+   f"{a3[1] :20.12f}"+f"{a3[2] :20.12f}"+ '\n')
   bfile.write(''+ '\n')
   bfile.write('K_POINTS {automatic}'+ '\n')
   dk=0.3
   dk=0.2
   ka,kb,kc=write_kpoints0(dk,a1,a2,a3)
   bfile.write(f"{ka : 4d}"+f"{kb : 4d}"+f"{kc : 4d}"+f"{0 : 4d} "+f"{0 : 2d} "+f"{0 : 2d} "+ '\n')
   bfile.close()
if True:
   bfile=open('scf.in_phonopy','w')
   bfile.write('&CONTROL'+  '\n')
   bfile.write('prefix= '+ "'rlxallz4'"+  '\n')
   if False:
      bfile.write('wf_collect= .true. '+ '\n')
   bfile.write('calculation= '+ "'scf'"+  '\n') 
   if ipseudo == 1:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/PSEUDOPOTENTIALS_NC/'"+  '\n')
   if ipseudo == 2:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-sr-04_pbe_standard_upf/'"+  '\n')
   if ipseudo == 3:
      bfile.write('pseudo_dir= '+ "'/home/ihlee/qe/nc-fr-04_pbe_standard_upf/'"+  '\n')
   bfile.write('outdir= '+ "'./'"+  '\n')
   bfile.write('tstress= .true.'+  '\n')
   bfile.write('tprnfor= .true.'+  '\n')
   bfile.write('/'+  '\n')
   bfile.write('&SYSTEM'+  '\n')
   bfile.write('nat= '+str(natot)+  '\n')
   bfile.write('ntyp= '+str(nspecies)+  '\n')
   bfile.write('ecutwfc= '+ '100.0'+  '\n')
   bfile.write('occupations= '+ "'smearing'"+  '\n')
   bfile.write('smearing= '+ "'mp'"+  '\n')
   bfile.write('degauss= 0.01'+  '\n')
   bfile.write('ibrav= '+ '0'+  '\n')
   if False:
      bfile.write('use_all_frac= .true.'+ '\n')
      bfile.write('lspinorb= .true.'+ '\n')
      bfile.write('noncolin= .true.'+ '\n')
   bfile.write('/'+  '\n')
   bfile.write('&ELECTRONS'+  '\n')
   bfile.write('conv_thr= 1.0d-8'+  '\n')
   bfile.write('mixing_beta=0.3'+  '\n')
   bfile.write('/'+  '\n')
   bfile.write(''+  '\n')
   bfile.write('K_POINTS {automatic}'+ '\n')
   dk=0.3
   dk=0.2
   ka,kb,kc=write_kpoints0(dk,a1,a2,a3)
   bfile.write(f"{ka : 4d}"+f"{kb : 4d}"+f"{kc : 4d}"+f"{0 : 4d} "+f"{0 : 2d} "+f"{0 : 2d} "+ '\n')
   bfile.write(''+ '\n')
   bfile.write('ATOMIC_SPECIES'+  '\n')
   if ipseudo == 2 or ipseudo == 3:
      for i in range(nspecies):
          bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.upf'+ '\n')
   if ipseudo == 1:
      for i in range(nspecies):
          if symbl[i] == 'B' or symbl[i] == 'C' or symbl[i] == 'N' or symbl[i] == 'O' :
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-nc.UPF'+ '\n')
          else:    
             bfile.write(symbl[i]+' '+str(atomicmass(symbl[i]))+' '+symbl[i]+'.rel-pbe-n-nc.UPF'+ '\n')
   bfile.write(''+  '\n')
#  bfile.write('ATOMIC_POSITIONS {crystal}'+  '\n')
   bfile.write('ATOMIC_POSITIONS crystal'+  '\n')
   for i in range(nspecies):
       for k in range(natot):
           if isymbl[k] == symbl[i]:
               bfile.write( isymbl[k]+' '+f"{drt[k,0] :20.12f}"+f"{drt[k,1] :20.12f}"+f"{drt[k,2] :20.12f}"+  '\n')
   bfile.write(''+ '\n')
#  bfile.write('CELL_PARAMETERS {angstrom}'+ '\n')
   bfile.write('CELL_PARAMETERS bohr'+ '\n')
   bfile.write(f"{a1[0]/0.529177 :20.12f}"+   f"{a1[1]/0.529177 :20.12f}"+f"{a1[2]/0.529177 :20.12f}"+ '\n')
   bfile.write(f"{a2[0]/0.529177 :20.12f}"+   f"{a2[1]/0.529177 :20.12f}"+f"{a2[2]/0.529177 :20.12f}"+ '\n')
   bfile.write(f"{a3[0]/0.529177 :20.12f}"+   f"{a3[1]/0.529177 :20.12f}"+f"{a3[2]/0.529177 :20.12f}"+ '\n')
   bfile.close()
