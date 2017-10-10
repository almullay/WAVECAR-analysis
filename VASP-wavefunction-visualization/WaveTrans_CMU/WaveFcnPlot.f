************************* WaveFcnPlot ************************************
*
*   input values from the GCOEFF file produced by WaveTrans, and
*   outputs wavefunction vs. z at specified x,y, spin value, k value, 
*   and band, to WAVEFCN.txt
*
*   version 1.0 - July 3, 2012 - R. M. Feenstra and M. Widom
                                                                        
      implicit real*8 (a-h, o-z)                                        
                                                                        
*   dimensions below must be sufficiently large to agree with values
*   in GCOEFF file; program issues error message(s) otherwise

      parameter(nbdim=120,npdim=15000,nwdim=11,nsdim=1)

      complex*8 coeff(npdim,nbdim,nwdim,nsdim)
      complex*16 cener(nbdim,nwdim,nsdim),csum
      dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),vtmp(3),sumkg(3),
     &occ(nbdim,nwdim,nsdim),wk(3,nwdim,nsdim),nplane(nwdim,nsdim),
     &igall(3,npdim,nwdim,nsdim)
      pi=4.*atan(1.)
      
*   input
      
      open(unit=11,file='GCOEFF.txt')
      read(11,*) nspin
      read(11,*) nwk
      read(11,*) nband
      write(6,*) 'no. spin values =',nspin
      write(6,*) 'no. k points =',nwk
      write(6,*) 'no. bands =',nband
      if (nspin.gt.nsdim) then
         write(6,*) '*** error - nsdim < nspin'
         write(6,*) 'increase nsdim in PARAMETER statement and ',
     &              'recompile program'
         write(6,*) 'press ENTER to continue'
         read(5,*)
         stop
      end if
      if (nwk.gt.nwdim) then
         write(6,*) '*** error - nwdim < nwk'
         write(6,*) 'increase nwdim in PARAMETER statement and ',
     &              'recompile program'
         write(6,*) 'press ENTER to continue'
         read(5,*)
         stop
      end if
      if (nband.gt.nbdim) then
         write(6,*) '*** error - nbdim < nband'
         write(6,*) 'increase nbdim in PARAMETER statement and ',
     &              'recompile program'
         write(6,*) 'press ENTER to continue'
         read(5,*)
         stop
      end if
      read(11,*) (a1(j),j=1,3)
      read(11,*) (a2(j),j=1,3)
      read(11,*) (a3(j),j=1,3)
      write(6,*) 'real space lattice vectors:'
      write(6,*) 'a1 =',(sngl(a1(j)),j=1,3)
      write(6,*) 'a2 =',(sngl(a2(j)),j=1,3)
      write(6,*) 'a3 =',(sngl(a3(j)),j=1,3)
      write(6,*) ' '
      call vcross(vtmp,a2,a3)
      Vcell=a1(1)*vtmp(1)+a1(2)*vtmp(2)+a1(3)*vtmp(3)
      write(6,*) 'volume unit cell =',sngl(Vcell)
      write(6,*) ' '
      read(11,*) (b1(j),j=1,3)
      read(11,*) (b2(j),j=1,3)
      read(11,*) (b3(j),j=1,3)
      write(6,*) 'reciprocal lattice vectors:'
      write(6,*) 'b1 =',(sngl(b1(j)),j=1,3)
      write(6,*) 'b2 =',(sngl(b2(j)),j=1,3)
      write(6,*) 'b3 =',(sngl(b3(j)),j=1,3)
      write(6,*) ' '
      write(6,*) 'reading plane wave coefficients'

      do 300 ispin=1,nspin
         do 200 iwk=1,nwk
            write(6,*) 'k point #',iwk
            read(11,*) (wk(j,iwk,ispin),j=1,3)
            do 190 iband=1,nband
               read(11,*) itmp,nplane(iwk,ispin)
               read(11,*) cener(iband,iwk,ispin),occ(iband,iwk,ispin)
               do 180 iplane=1,nplane(iwk,ispin)
                  read(11,*) (igall(j,iplane,iwk,ispin),j=1,3),
     &                        coeff(iplane,iband,iwk,ispin)
180            continue               
190         continue
200      continue
300   continue

*   specified (x,y,z) location of wavefunction plot, number of points, and
*   the spin, wavevector, and band values (ispin,iwk and iband)
*   (in this sample program, these are all specified here in the code, whereas
*   in an actual application they would typically be part of an input file)

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
     &          (wk(2,iwk,ispin)+ig2)*b2(j)+(wk(3,iwk,ispin)+ig3)*b3(j)
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
*
*   routine for computing vector cross-product
*
      subroutine vcross(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(3),b(3),c(3)
      
      a(1)=b(2)*c(3)-b(3)*c(2)
      a(2)=-(b(1)*c(3)-b(3)*c(1))
      a(3)=b(1)*c(2)-b(2)*c(1)
      return
      end
      