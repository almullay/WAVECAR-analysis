************************* WaveTrans ************************************
*
*   call WAVECARin to input values from the WAVECAR file produced by 
*   VASP, and outputs GCOEFF file with G values and corresponding plane
*   wave coefficients, together with energy eigenvalues and occupations
*
*   version 1.0 - July 3, 2012 - R. M. Feenstra and M. Widom
*   version 1.1 - July 13, 2012 - split into WAVECARin subroutine
*                                 and calling program
*
*   format of GCOEFF file:
*     no. spin values (integer)
*     no. wavevector values (integer)
*     no. bands (integer)
*     real space lattice vector (a1) x,y,z coefficients (real)
*     real space lattice vector (a2) x,y,z coefficients (real)
*     real space lattice vector (a3) x,y,z coefficients (real)
*     recip. lattice vector 1 (b1) x,y,z coefficients (real)
*     recip. lattice vector 2 (b2) x,y,z coefficients (real)
*     recip. lattice vector 3 (b3) x,y,z coefficients (real)
*     loop through spin values
*        loop through wavevector values:
*           wavevector x,y,z components (real)
*           energy (complex), occupation (real)
*           loop through bands:
*              index of band (integer), no. plane waves (integer)
*              loop through plane waves:
*                 ig1,ig2,ig3 values (integer) and coefficient (complex)
*              end loop through plane waves
*           end loop through bands
*        end loop through wavevector values
*     end loop through spin values
*                                                                        
*   the x,y,z components of each G value are given in terms of the
*   ig1,ig2,ig3 values and the components of the recip. lattice vectors
*   according to:
*   ig1*b1_x + ig2*b2_x + ig3*b3_x,
*   ig1*b1_y + ig2*b2_y + ig3*b3_y, and
*   ig1*b1_z + ig2*b2_z + ig3*b3_z, respectively
*
*   note that the energy eigenvalues are complex, as provided in the
*   WAVECAR file, but the imaginary part is zero (at least for all cases
*   investigated thus far)
*     
      implicit real*8 (a-h, o-z)                                        

*   dimensions below must be sufficiently large to agree with values
*   in WAVECAR file; WAVECARin routine issues error message(s) otherwise

      parameter(nbdim=208,npdim=400000,nwdim=1,nsdim=1)

      complex*8 coeff(npdim,nbdim,nwdim,nsdim)
      complex*16 cener(nbdim,nwdim,nsdim)
      dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),vtmp(3),sumkg(3),
     &occ(nbdim,nwdim,nsdim),wk(3,nwdim,nsdim),nplane(nwdim,nsdim),
     &igall(3,npdim,nwdim,nsdim)
     
*   read WAVECAR file, using routine WAVECARin

      call WAVECARin(nsdim,nwdim,nbdim,npdim,nspin,nwk,nband,
     &nplane,igall,coeff,cener,occ,wk,a1,a2,a3,b1,b2,b3,Vcell)                                                                        

*   output G values and coefficients

      write(6,*) ' '
      write(6,*) 'writing file GCOEFF.txt'
      open(unit=11,file='GCOEFF.txt')
      write(11,*) nspin
      write(11,*) nwk
      write(11,*) nband
      write(11,*) (sngl(a1(j)),j=1,3)
      write(11,*) (sngl(a2(j)),j=1,3)
      write(11,*) (sngl(a3(j)),j=1,3)
      write(11,*) (sngl(b1(j)),j=1,3)
      write(11,*) (sngl(b2(j)),j=1,3)
      write(11,*) (sngl(b3(j)),j=1,3)
      do 700 ispin=1,nspin
         do 600 iwk=1,nwk
            write(6,*) 'k point #',iwk
            write(11,*) (sngl(wk(j,iwk,ispin)),j=1,3)
            do 590 iband=1,3
               write(11,*) iband,nplane(iwk,ispin)
               write(11,560) cener(iband,iwk,ispin),occ(iband,iwk,ispin)
560            format('( ',g14.6,' , ',g14.6,' ) ',g14.6)            
               do 580 iplane=1,nplane(iwk,ispin)
                  write(11,570) (igall(j,iplane,iwk,ispin),j=1,3),
     &                        coeff(iplane,iband,iwk,ispin)
570            format(3i6,'  ( ',g14.6,' , ',g14.6,' )')     
580            continue               
590         continue
600      continue
700   continue

*   exit

      write(6,*) 'press ENTER to continue'
      read(5,*)
      stop
      end
      
      