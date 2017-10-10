************************* WAVECARin ************************************
*
*   inputs the WAVECAR file in binary format from VASP, and outputs 
*   an array of G values and corresponding plane wave coefficients, 
*   together with energy eigenvalues and occupations
*
*   version 1.0 - July 3, 2012 - R. M. Feenstra and M. Widom
*   version 1.1 - July 13, 2012 - split into WAVECARin subroutine
*                                 and calling programs
*   version 1.2 - Sept 30, 2012 - changed estimator for max. no. of
*                                 plane waves
*   version 1.3 - Sept 17, 2014 - updated 'c' value
*
*   inputs:
*   nsdim = dimension for number of spin values
*   nwdim = dimension for number of wavevector (k) values
*   nbdim = dimension for number of bands
*   npdim = dimension for number of plane waves
*
*   dimensions nsdim,nwdim,nbdim,npdim are set in calling program, and
*   must be sufficiently large to agree with values in WAVECAR file; 
*   routine issues error message(s) otherwise
*
*   outputs:
*   nspin = actual number of spin values
*   nwk = actual number of wavevector (k) values
*   nband = actual number of bands
*   nplane = array of numbers of plane waves (one for each k and spin)
*   igall = array of G indices (see below for definition)
*   coeff = array of plane wave coefficients
*   cener = array of energy eigenvalues
*   occ = array of occupation of states
*   wk = array of wavevector (k) values
*   a1,a2,a3 = real-space primitive lattice vectors
*   b1,b2,b3 = reciprocal-space lattice vectors
*   Vcell = volume of real-space primitive unit cell
*   
*   the x,y,z components of each G value are given in terms of the
*   igall values and the components of the recip. lattice vectors
*   according to:
*   ig1*b1_x + ig2*b2_x + ig3*b3_x,
*   ig1*b1_y + ig2*b2_y + ig3*b3_y, and
*   ig1*b1_z + ig2*b2_z + ig3*b3_z, respectively,
*   with
*   ig1=igall(1,iplane,iwk,ispin),
*   ig2=igall(2,iplane,iwk,ispin), and
*   ig3=igall(3,iplane,iwk,ispin),
*   where iplane=1,2,...,nplane(iwk,ispin) is an index incrementing
*   the plane waves for specific k and spin values
*
*   note that the energy eigenvalues are complex, as provided in the
*   WAVECAR file, but the imaginary part is zero (at least for all cases
*   investigated thus far)
*     

      subroutine WAVECARin(nsdim,nwdim,nbdim,npdim,nspin,nwk,nband,
     &nplane,igall,coeff,cener,occ,wk,a1,a2,a3,b1,b2,b3,Vcell)                                                                        

      implicit real*8 (a-h, o-z)                                        

      complex*8 coeff(npdim,nbdim,nwdim,nsdim)
      complex*16 cener(nbdim,nwdim,nsdim)
      character*1 chartmp
      dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),vtmp(3),sumkg(3),
     &occ(nbdim,nwdim,nsdim),wk(3,nwdim,nsdim),nplane(nwdim,nsdim),
     &igall(3,npdim,nwdim,nsdim)
     
*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
*   adjusted in final decimal places to agree with VASP value; program
*   checks for discrepancy of any results between this and VASP values)

      data c/0.262465831d0/ 
*      data c/0.26246582250210965422d0/ 
      pi=4.*atan(1.)
      
*   input (first test for presence of non-empty WAVECAR file)
      
      open(unit=10,file='WAVECAR',iostat=iost,status='old')
      if (iost.ne.0) then
         write(6,*) '*** error - WAVECAR file nonexistent'
         write(6,*) 'press ENTER to continue'
         read(5,*)
         stop
      end if
      read(10,10,end=20) chartmp
10    format(1a1)      
      go to 30
20       write(6,*) '*** error - WAVECAR file empty'
         write(6,*) 'press ENTER to continue'
         read(5,*)
         stop
30    continue      
      close(unit=10)
      
      nrecl=16
      open(unit=10,file='WAVECAR',access='direct',recl=nrecl,
     &iostat=iost,status='old')
      read(unit=10,rec=1) xnrecl,xnspin
      close(unit=10)
      nrecl=nint(xnrecl)
      nspin=nint(xnspin)
      write(6,*) 'record length  =',nrecl
      write(6,*) 'no. spin values =',nspin
      open(unit=10,file='WAVECAR',access='direct',recl=nrecl,
     &iostat=iost,status='old')
      if (iost.ne.0) write(6,*) 'open error - iostat =',iost            
      read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3),
     &(a3(j),j=1,3)
      nwk=nint(xnwk)
      nband=nint(xnband)
      write(6,*) 'no. k points =',nwk
      write(6,*) 'no. bands =',nband
      write(6,*) 'max. energy =',sngl(ecut)
      if (nspin.gt.nsdim) then
         write(6,*) '*** error - nsdim < nspin'
         write(6,*) 'increase nsdim in calling program and ',
     &              'recompile program'
         write(6,*) 'press ENTER to continue'
         read(5,*)
         stop
      end if
      if (nwk.gt.nwdim) then
         write(6,*) '*** error - nwdim < nwk'
         write(6,*) 'increase nwdim in calling program and ',
     &              'recompile program'
         write(6,*) 'press ENTER to continue'
         read(5,*)
         stop
      end if
      if (nband.gt.nbdim) then
         write(6,*) '*** error - nbdim < nband'
         write(6,*) 'increase nbdim in calling program and ',
     &              'recompile program'
         write(6,*) 'press ENTER to continue'
         read(5,*)
         stop
      end if
      write(6,*) 'real space lattice vectors:'
      write(6,*) 'a1 =',(sngl(a1(j)),j=1,3)
      write(6,*) 'a2 =',(sngl(a2(j)),j=1,3)
      write(6,*) 'a3 =',(sngl(a3(j)),j=1,3)
      write(6,*) ' '

      irec=2
      do 220 ispin=1,nspin
         do 200 iwk=1,nwk
            irec=irec+1
            read(unit=10,rec=irec) xnplane,(wk(i,iwk,ispin),i=1,3),
     &      (cener(iband,iwk,ispin),occ(iband,iwk,ispin),iband=1,nband)
            nplane(iwk,ispin)=nint(xnplane)
            write(6,*) 'k point #',iwk,'  input no. of plane waves =',
     &                  nplane(iwk,ispin)
            write(6,*) 'k value =',(sngl(wk(j,iwk,ispin)),j=1,3)
            if (nplane(iwk,ispin).gt.npdim) then
               write(6,*) '*** error - npdim < nplane'
               write(6,*) 'increase npdim in PARAMETER statement and ',
     &                    'recompile program'
               write(6,*) 'press ENTER to continue'
               read(5,*)
               stop
            end if
            do 190 iband=1,nband
               irec=irec+1
               read(unit=10,rec=irec) (coeff(iplane,iband,iwk,ispin),
     &            iplane=1,nplane(iwk,ispin))
190         continue
200      continue
220   continue

*   compute corresponding G values

      write(6,*) ' '
      call vcross(vtmp,a2,a3)
      Vcell=a1(1)*vtmp(1)+a1(2)*vtmp(2)+a1(3)*vtmp(3)
      write(6,*) 'volume unit cell =',sngl(Vcell)
      call vcross(b1,a2,a3)
      call vcross(b2,a3,a1)
      call vcross(b3,a1,a2)
      do 250 j=1,3
         b1(j)=2.*pi*b1(j)/Vcell
         b2(j)=2.*pi*b2(j)/Vcell
         b3(j)=2.*pi*b3(j)/Vcell
250   continue
      write(6,*) 'reciprocal lattice vectors:'
      write(6,*) 'b1 =',(sngl(b1(j)),j=1,3)
      write(6,*) 'b2 =',(sngl(b2(j)),j=1,3)
      write(6,*) 'b3 =',(sngl(b3(j)),j=1,3)
      write(6,*) 'reciprocal lattice vector magnitudes:'
      b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
      b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
      b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)
      write(6,*) sngl(b1mag),sngl(b2mag),sngl(b3mag)
      write(6,*) ' '
      
*   estimate max. values for recip. lattice vectors indices, and
*   the max. no. plane waves (the latter is not actually used,
*   but keep it here for reference with f90 version)

      phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
      call vcross(vtmp,b1,b2)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
      nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
      nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
      nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
      npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
      
      phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
      call vcross(vtmp,b1,b3)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
      phi123=abs(asin(sinphi123))
      nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
      nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
      nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
      npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      
      phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
      call vcross(vtmp,b2,b3)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
      phi123=abs(asin(sinphi123))
      nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
      nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
      nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
      npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

      nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
      nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
      nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
      npmax=min0(npmaxA,npmaxB,npmaxC)
      
      write(6,*) 'max. no. G values; 1,2,3 =',nb1max,nb2max,nb3max
      write(6,*) 'estimated max. no. plane waves =',npmax
      write(6,*) ' '
      
      do 400 ispin=1,nspin
         do 350 iwk=1,nwk
            ncnt=0
            do 300 ig3=0,2*nb3max
               ig3p=ig3
               if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
               do 290 ig2=0,2*nb2max
                  ig2p=ig2
                  if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
                  do 280 ig1=0,2*nb1max
                     ig1p=ig1
                     if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
                     do 270 j=1,3
                        sumkg(j)=(wk(1,iwk,ispin)+ig1p)*b1(j)+
     &                           (wk(2,iwk,ispin)+ig2p)*b2(j)+
     &                           (wk(3,iwk,ispin)+ig3p)*b3(j)
270                  continue
                     gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
                     etot=gtot**2/c
                     if (etot.lt.ecut) then
                        ncnt=ncnt+1
                        igall(1,ncnt,iwk,ispin)=ig1p
                        igall(2,ncnt,iwk,ispin)=ig2p
                        igall(3,ncnt,iwk,ispin)=ig3p
                     end if
280               continue
290            continue
300         continue
            write(6,*) 'k point #',iwk,'   computed no. plane waves =',
     &         ncnt
            if (ncnt.ne.nplane(iwk,ispin)) then
               write(6,*) '*** error - computed no. plane waves <>',
     &                    ' input no.'
               write(6,*) 'press ENTER to exit'
               read(5,*)
               stop
            end if
350      continue 
400   continue

*   exit

      return
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
      
      