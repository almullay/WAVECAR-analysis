!!$************************* WaveTransPlot*********************************
!!$
!!$   input the WAVECAR file in binary format from VASP, and write
!!$   selected real space wavefunction in a3 direction to standard output
!!$   Plane wave coefficients are written to GCOEFF.txt
!!$
!!$   Compile with gfortran or ifort. Flag "-assume byterecl" is required
!!$   for ifort.
!!$
!!$   version 2.0 - July 11, 2012 - R. M. Feenstra and M. Widom
!!$   version 2.1 - Sept 30, 2012 - changed estimator for max. no. of
!!$                                 plane waves
!!$   version 1.3 - Sept 17, 2014 - updated 'c' value
!!$
!!$   options are -f filename -s spin -k k-point -b band -x coord -y coord
!!$   defaults are -f WAVECAR -s 1 -k 1 -b 1 -x 0 -y 0
!!$   coordinates are direct coordinates

implicit real*8 (a-h, o-z)
complex*8, allocatable :: coeff(:)
complex*16, allocatable :: cener(:)
real*8, allocatable :: occ(:)
integer, allocatable :: igall(:,:)
dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3),vtmp(3)
dimension wk(3),xyz(3),wkpg(3),ig(3)
complex*16 csum
integer spin,kpoint,band
character*75 filename

!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value; program
!!$*   checks for discrepancy of any results between this and VASP values)

data c/0.262465831d0/ 
!!$*   data c/0.26246582250210965422d0/ 
pi=4.*atan(1.)

!!$ parse arguments
call parse(filename,spin,kpoint,band,x,y)
xyz(1)=x
xyz(2)=y

!!$*   input

nrecl=24
open(unit=10,file=filename,access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(6,*) 'open error - iostat =',iost            
read(unit=10,rec=1) xnrecl,xnspin,xnprec
close(unit=10)
nrecl=nint(xnrecl)
nspin=nint(xnspin)
nprec=nint(xnprec)
if(nprec.eq.45210) then
   write(0,*) '*** error - WAVECAR_double requires complex*16'
   stop
endif
open(unit=10,file=filename,access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(6,*) 'open error - iostat =',iost
open(unit=11,file='GCOEFF.txt')
read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
     (a3(j),j=1,3)
nwk=nint(xnwk)
nband=nint(xnband)
if (kpoint.gt.nwk) then
   write(0,*) '*** error - selected k=',kpoint,' > max k=',nwk
   stop
endif
if (band.gt.nband) then
   write(0,*) '*** error - selected band=',band,' > max band=',nband
   stop
endif
allocate(occ(nband))
allocate(cener(nband))

!!$*   compute reciprocal properties

call vcross(a2xa3,a2,a3)
Vcell=dot_product(a1,a2xa3)
a3mag=dsqrt(dot_product(a3,a3))
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
   b1=2.*pi*b1/Vcell
   b2=2.*pi*b2/Vcell
   b3=2.*pi*b3/Vcell
write(11,*) (sngl(b1(j)),j=1,3)
write(11,*) (sngl(b2(j)),j=1,3)
write(11,*) (sngl(b3(j)),j=1,3)
b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

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

allocate (igall(3,npmax))
allocate (coeff(npmax))

!!$*   Find desired wavefunction
irec=3+(kpoint-1)*(nband+1)+(spin-1)*nwk*(nband+1)
read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
     (cener(iband),occ(iband),iband=1,nband)
nplane=nint(xnplane)
write(11,*) spin,kpoint,band
write(11,*) (sngl(wk(j)),j=1,3)
write(11,*) cener(band),occ(band)
      
!!$*   Calculate plane waves
ncnt=0
do ig3=0,2*nb3max
   ig3p=ig3
   if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
   do ig2=0,2*nb2max
      ig2p=ig2
      if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
      do ig1=0,2*nb1max
         ig1p=ig1
         if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
         do j=1,3
            sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                 (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
         enddo
         gtot=sqrt(dot_product(sumkg,sumkg))
         etot=gtot**2/c
         if (etot.lt.ecut) then
            ncnt=ncnt+1
            igall(1,ncnt)=ig1p
            igall(2,ncnt)=ig2p
            igall(3,ncnt)=ig3p
         end if
      enddo
   enddo
enddo
if (ncnt.ne.nplane) then
   write(0,*) '*** error - computed ncnt=',ncnt, &
        ' != input nplane=',nplane
   if(ncnt.eq.(2*nplane-1)) write(0,*) '*** suspect Gamma-only WAVECAR'
   stop
endif

irec=irec+band
read(unit=10,rec=irec) (coeff(iplane), iplane=1,nplane)
!!$*   output G value and coefficients
do iplane=1,nplane
   write(11,570) (igall(j,iplane),j=1,3), &
        coeff(iplane)
570 format(3i6,'  ( ',g14.6,' , ',g14.6,' )')     
enddo

do iz=0,2*nb3max
   z=dble(iz)/dble(1+2*nb3max)
   xyz(3)=z
   csum=cmplx(0.,0.)
   do iplane=1,nplane
      ig=igall(:,iplane)
      wkpg=wk+ig
      csum=csum+coeff(iplane)* &
           cdexp(2.*pi*cmplx(0.,1.)*dot_product(wkpg,xyz))
   enddo
   csum=csum/dsqrt(Vcell)
!!$ output z*a3mag for units of Angstroms
   write(6,*) sngl(z),sngl(real(csum)),sngl(dimag(csum))
enddo

deallocate(igall)
deallocate(coeff)
stop
end program

!!$*   routine for computing vector cross-product
subroutine vcross(a,b,c)
  implicit real*8(a-h,o-z)
  dimension a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross      
      
!!$   parse command line arguments
subroutine parse(filename,spin,kpoint,band,x,y)
character*75 filename
integer spin,band,kpoint
real*8 x,y
character*20 option,value
integer iarg,narg,ia
iarg=iargc()
nargs=iarg/2
filename="WAVECAR"
spin = 1
kpoint = 1
band = 1
x = 0.
y = 0.
if(iarg.ne.2*nargs) then
   call help
endif
do ia=1,nargs
   call getarg(2*ia-1,option)
   call getarg(2*ia,value)
   if(option == "-f") then
      read(value,*) filename
   else if(option == "-s") then
      read(value,*) spin
   else if(option == "-k") then
      read(value,*) kpoint
   else if(option == "-b") then
      read(value,*) band
   else if(option == "-x") then
      read(value,*) x
   else if(option == "-y") then
      read(value,*) y
   else if(option =="-h") then
      call help
   else
      call help
   endif
enddo
return
end subroutine parse

subroutine help
  write(6,*) "syntax: WaveTransPlot -f file -s spin -k k-point -b band &
       -x coord -y coord"
  write(6,*) "defaults: -f WAVECAR -s 1 -k 1 -b 1 -x 0.0 -y 0.0"
  write(6,*) "inputs: x and y are direct coordinates on axes a1 and a2"
  write(6,*) "output: wavefunction psi(x,y,z) with z direct coordinate &
       on a3 axis"
  stop
end subroutine help
