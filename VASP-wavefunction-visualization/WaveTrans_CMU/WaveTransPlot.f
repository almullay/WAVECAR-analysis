************************* WaveTransPlot ************************************
*
*   input values from the WAVECAR file produced by VASP, and
*   output wavefunction vs. z at specified x,y, spin value, k value, 
*   and band, to WAVEFCN.txt
*                                                                        
*   version 1.0 - July 3, 2012 - R. M. Feenstra and M. Widom
*   version 1.1 - July 13, 2012 - split into WAVECARin subroutine
*                                 and calling program

      implicit real*8 (a-h, o-z)                                        

*   dimensions below must be sufficiently large to agree with values
*   in WAVECAR file; WAVECARin routine issues error message(s) otherwise

      parameter(nbdim=120,npdim=15000,nwdim=11,nsdim=1)

      complex*8 coeff(npdim,nbdim,nwdim,nsdim)
      complex*16 cener(nbdim,nwdim,nsdim),csum
      dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),sumkg(3),
     &occ(nbdim,nwdim,nsdim),wk(3,nwdim,nsdim),nplane(nwdim,nsdim),
     &igall(3,npdim,nwdim,nsdim)

      pi=4.*atan(1.)
      
*   input

      call WAVECARin(nsdim,nwdim,nbdim,npdim,nspin,nwk,nband,
     &nplane,igall,coeff,cener,occ,wk,a1,a2,a3,b1,b2,b3,Vcell)                                                                        
 
*   specified location of wavefunction plot, and number of points

      ispin=1
      iwk=1
      iband=1
      x=0.
      y=0.
      nz=256
      zmax=a1(3)+a2(3)+a3(3)
      delz=zmax/float(nz)

*   construct wavefunction, and output

      open(unit=12,file='WAVEFCN.txt')
      do 550 iz=1,nz+1
         z=(iz-1)*delz
         csum=cmplx(0.,0.)
         do 500 iplane=1,nplane(iwk,ispin)
            ig1=igall(1,iplane,iwk,ispin)
            ig2=igall(2,iplane,iwk,ispin)
            ig3=igall(3,iplane,iwk,ispin)
            do 470 j=1,3
               sumkg(j)=(wk(1,iwk,ispin)+ig1)*b1(j)+
     &                  (wk(2,iwk,ispin)+ig2)*b2(j)+
     &                  (wk(3,iwk,ispin)+ig3)*b3(j)
470         continue
            csum=csum+coeff(iplane,iband,iwk,ispin)*cdexp(cmplx(0.,1.)*
     &            (sumkg(1)*x+sumkg(2)*y+sumkg(3)*z))
500      continue
         csum=csum/dsqrt(Vcell)
         write(12,*) sngl(z),sngl(real(csum)),sngl(dimag(csum))
550   continue
      write(6,*) ' '
      write(6,*) 'output WAVEFCN.txt'

*   exit

      write(6,*) 'press ENTER to continue'
      read(5,*)
      stop
      end