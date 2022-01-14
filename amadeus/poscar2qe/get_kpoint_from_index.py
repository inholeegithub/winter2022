# use a small KSPACING = 0.08 in INCAR file
# check the EIGENVAL file and get two integers, iv and ic
# input: EIGENVAL, POSCAR
#        iv, ic
# output: inp_vbm, inp_cbm
iv=54
iv=36
iv=18
ic=iv+1


with open('EIGENVAL','r') as afile:
     k1v=0. ; k2v=0. ; k3v=0. ; wv=0. ; vbm=-1.e20
     for line in afile:
         if len(line.split()) >= 3 :
            if len(line.split()) == 4 :
               k1tmp=float(line.split()[0])
               k2tmp=float(line.split()[1])
               k3tmp=float(line.split()[2])
               wtmp=float(line.split()[3])
            if line.split()[0] == str(iv) and len(line.split()) == 3 :
               if float(line.split()[1]) > vbm:
                  vbm=float(line.split()[1])
                  k1v=k1tmp ; k2v=k2tmp ; k3v=k3tmp ; wv=wtmp
         else :
           continue
print(k1v,k2v,k3v,wv)
print(vbm)
with open('EIGENVAL','r') as afile:
     k1c=0. ; k2c=0. ; k3c=0. ; wc=0. ; cbm=1.e20
     for line in afile:
         if len(line.split()) >= 3 :
            if len(line.split()) == 4 :
               k1tmp=float(line.split()[0])
               k2tmp=float(line.split()[1])
               k3tmp=float(line.split()[2])
               wtmp=float(line.split()[3])
            if line.split()[0] == str(ic) and len(line.split()) == 3 :
               if float(line.split()[1]) < cbm:
                  cbm=float(line.split()[1])
                  k1c=k1tmp ; k2c=k2tmp ; k3c=k3tmp ; wc=wtmp
         else :
           continue
print(k1c,k2c,k3c,wc)
print(cbm)
print(cbm-vbm, ' eV')
jline=0
with open('POSCAR','r') as afile:
     for line in afile:
         jline=jline+1
         if jline == 1: 
            continue
         if jline == 2: 
            scale0=float(line.split()[0])
         if jline == 3: 
            a11=float(line.split()[0])
            a12=float(line.split()[1])
            a13=float(line.split()[2])
         if jline == 4: 
            a21=float(line.split()[0])
            a22=float(line.split()[1])
            a23=float(line.split()[2])
         if jline == 5: 
            a31=float(line.split()[0])
            a32=float(line.split()[1])
            a33=float(line.split()[2])
if scale0 > 0. :
   a11=a11*scale0 ; a12=a12*scale0 ; a13=a13*scale0
   a21=a21*scale0 ; a22=a22*scale0 ; a23=a23*scale0
   a31=a31*scale0 ; a32=a32*scale0 ; a33=a33*scale0
if scale0 < 0. :
   vtest=(a12*a23-a13*a22)*a31+(a13*a21-a11*a23)*a32+(a11*a22-a12*a21)*a33
   if vtest < 0. :
      vtest=-vtest
   vtest=(-scale0)/vtest 
   vtest=vtest**(1.e0/3.e0)
   a11=a11*vtest ; a12=a12*vtest ; a13=a13*vtest
   a21=a21*vtest ; a22=a22*vtest ; a23=a23*vtest
   a31=a31*vtest ; a32=a32*vtest ; a33=a33*vtest
print('')
print('cp inp_vbm   inp')
print('/home/ihlee/csa_vasp/emc-1.50/fortran/emc_gen')
print('mpirun -np 8 /TGM/Apps/VASP/bin/5.4.4/O3/NORMAL/vasp.5.4.4.pl2.O3.NORMAL.ncl.x > stdout.log')
print('/home/ihlee/csa_vasp/emc-1.50/fortran/emc_calc')
print(k1v,k2v,k3v)
print('0.01')
print(iv)
print('V')
print(a11,a12,a13)
print(a21,a22,a23)
print(a31,a32,a33)
print('')
print('cp inp_cbm   inp')
print('/home/ihlee/csa_vasp/emc-1.50/fortran/emc_gen')
print('mpirun -np 8 /TGM/Apps/VASP/bin/5.4.4/O3/NORMAL/vasp.5.4.4.pl2.O3.NORMAL.ncl.x > stdout.log')
print('/home/ihlee/csa_vasp/emc-1.50/fortran/emc_calc')
print(k1c,k2c,k3c)
print('0.01')
print(ic)
print('V')
print(a11,a12,a13)
print(a21,a22,a23)
print(a31,a32,a33)


# Write-Overwrites 
file1 = open("inp_vbm", "w")  # write mode 
astring=str(k1v)+' '+str(k2v)+' '+str(k3v)+' \n'
file1.write(astring)
astring='0.01'+' \n'
file1.write(astring)
astring=str(iv)+' \n'
file1.write(astring)
astring='V'+' \n'
file1.write(astring)
astring=str(a11)+' '+str(a12)+' '+str(a13)+' \n'
file1.write(astring)
astring=str(a21)+' '+str(a22)+' '+str(a23)+' \n'
file1.write(astring)
astring=str(a31)+' '+str(a32)+' '+str(a33)+' \n'
file1.write(astring)
file1.close() 
# Write-Overwrites 
file1 = open("inp_cbm", "w")  # write mode 
astring=str(k1c)+' '+str(k2c)+' '+str(k3c)+' \n'
file1.write(astring)
astring='0.01'+' \n'
file1.write(astring)
astring=str(ic)+' \n'
file1.write(astring)
astring='V'+' \n'
file1.write(astring)
astring=str(a11)+' '+str(a12)+' '+str(a13)+' \n'
file1.write(astring)
astring=str(a21)+' '+str(a22)+' '+str(a23)+' \n'
file1.write(astring)
astring=str(a31)+' '+str(a32)+' '+str(a33)+' \n'
file1.write(astring)
file1.close() 
