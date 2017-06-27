!This program performs a Monte Carlo simulation for the 2D Ising gauge model on a square lattice
!The temperature is included in the coupling strength
!J=K/(k_BT) where K is the bare coupling strength and
!H=B/(k_BT) where B is the magnetic field
program Ising
include 'globIsgaugeMC'
integer i,isite,ix,iy
integer num_sites
character (len=90) :: filename
character (len=5) :: size_str

size=28
write(size_str,'(I5)') size
write (filename, '(A,A,A)' )'configuration',trim(adjustl(size_str)),'.csv'
open(101,file=filename)
!set temperature points
do i=1,tpoints
 temperature(i)=temperature_min+dble(i-1)*temperature_steps
 J=1.0d0/temperature(i)
 call Initialise
 call MonteCarlo
 en_per_site(i)=energy_av
 C_v(i)=Energy_square_av-Energy_av**2
 mag_per_site(i)=mag_av
 xi(i)=J*(mag_square_av-mag_av**2)
end do

close(101)

!call output

end program Ising

subroutine Initialise
include 'globIsgaugeMC'
integer i,ix,iy

H=0.0d0
MCSteps=50000
MCSteps_eq=1000

!weight of single spin flip
s=0.75d0

do ix = 1, Size
 do iy = 1, Size      
  gauges(ix, iy) = 1
  spins_x(ix,iy) = 1
  spins_y(ix,iy) = 1
 end do
end do
Magnetisation = 2*Size**2
mag_av=0.0d0
mag_square_av=0.0d0
Energy = -1.0d0*Size*Size*J
Energy_av=0.0d0
Energy_square_av=0.0d0
do i=1,3
 MCWeights(i) = dexp(-dble(2*(i-1))*J)
end do

end subroutine Initialise

subroutine flipx(ix,iy)
include 'globIsgaugeMC'
integer ix, iy, costfactor

!flip x bond
costfactor=spins_x(ix+1,iy+1)*spins_y(mod(ix+1,size)+1,iy+1)&
   *spins_x(ix+1,mod(iy+1,size)+1)*spins_y(ix+1,iy+1)&
   +spins_x(ix+1,iy+1)*spins_y(mod(ix+1,size)+1,mod(iy-1+size,size)+1)&
   *spins_x(ix+1,mod(iy-1+size,size)+1)*spins_y(ix+1,mod(iy-1+size,size)+1)

randnumber=grnd()
if (costfactor>0) then
 if (randnumber<s*MCWeights(costfactor+1)) then
  spins_x(ix+1,iy+1)=-spins_x(ix+1,iy+1)
  magnetisation=magnetisation+2*spins_x(ix+1,iy+1)
  energy=energy+dble(2*costfactor)*J
 end if
else
 if (randnumber<s) then
  spins_x(ix+1,iy+1)=-spins_x(ix+1,iy+1)
  magnetisation=magnetisation+2*spins_x(ix+1,iy+1)
  energy=energy+dble(2*costfactor)*J
 end if
end if

end subroutine flipx

subroutine flipy(ix,iy)
include 'globIsgaugeMC'
integer ix, iy, costfactor

!flip y bond
costfactor=spins_y(ix+1,iy+1)*spins_x(ix+1,iy+1)&
   *spins_y(mod(ix+1,size)+1,iy+1)*spins_x(ix+1,mod(iy+1,size)+1)&
   +spins_y(ix+1,iy+1)*spins_x(mod(ix-1+size,size)+1,iy+1)&
   *spins_y(mod(ix-1+size,size)+1,iy+1)*spins_x(mod(ix-1+size,size)+1,mod(iy+1,size)+1)

randnumber=grnd()
if (costfactor>0) then
 if (randnumber<s*MCWeights(costfactor+1)) then
  spins_y(ix+1,iy+1)=-spins_y(ix+1,iy+1)
  magnetisation=magnetisation+2*spins_y(ix+1,iy+1)
  energy=energy+dble(2*costfactor)*J
 end if
else
 if (randnumber<s) then
  spins_y(ix+1,iy+1)=-spins_y(ix+1,iy+1)
  magnetisation=magnetisation+2*spins_y(ix+1,iy+1)
  energy=energy+dble(2*costfactor)*J
 end if
end if

end subroutine flipy

subroutine gaugeupdate(ix,iy)
include 'globIsgaugeMC'
integer ix,iy

randnumber=grnd()
if (randnumber<(1.0d0-s)) then
 spins_x(ix+1,iy+1)=-spins_x(ix+1,iy+1)
 spins_y(ix+1,iy+1)=-spins_y(ix+1,iy+1)
 spins_x(mod(ix-1+size,size)+1,iy+1)=-spins_x(mod(ix-1+size,size)+1,iy+1)
 spins_y(ix+1,mod(iy-1+size,size)+1)=-spins_y(ix+1,mod(iy-1+size,size)+1)
end if

end subroutine gaugeupdate

subroutine MonteCarlo
include 'globIsgaugeMC'

integer mciter, ix,iy,costfactor,curspin
integer num_bonds

!equilibrate
do mciter=1,mcsteps_eq
 do ix=0,size-1
  do iy=0,size-1
   !flip x bond
   call flipx(ix,iy)
   !flip y bond
   call flipy(ix,iy)
   !gaugeupdate
   call gaugeupdate(ix,iy)
  end do
 end do
end do

num_bonds=2*size*size

!MC steps after equilibration
do mciter=1,mcsteps
 do ix=0,size-1
  do iy=0,size-1
   !flip x bond
   call flipx(ix,iy)
   !flip y bond
   call flipy(ix,iy)
   !guageupdate
   call gaugeupdate(ix,iy)
  end do
 end do
 mag_av=mag_av+dble(magnetisation)
 energy_av=energy_av+energy
 Energy_square_av=Energy_square_av+(energy/dble(2*size**2))**2
 mag_square_av=mag_square_av+(dble(magnetisation)/dble(2*size**2))**2

 !record configurations
 if(mod(mciter,10)==0) then
  write(101,"(f20.10,<num_bonds>(',', I10))")1.0d0/J,(spins_x(ix,:size),ix=1,size),(spins_y(ix,:size),ix=1,size)
 end if
end do

!average energy persite
energy_av=energy_av/dble(num_bonds*mcsteps)
Energy_square_av=Energy_square_av/dble(mcsteps)

!average magnetisation
mag_av=mag_av/dble(num_bonds*mcsteps)
mag_square_av=mag_square_av/dble(mcsteps)

end subroutine MonteCarlo

subroutine output
include 'globIsgaugeMC'
integer i
character (len=90) :: filename

write(filename, '(A,I2,A)' )'en',size,'.csv'
open(102,file=filename)
do i=1,tpoints
 write(102,"(9(f20.10, ',', f20.10))")temperature(i),en_per_site(i)*temperature(i)
end do
close(102)

write(filename, '(A,I2,A)' )'Cv',size,'.csv'
open(103,file=filename)
do i=1,tpoints
 write(103,"(9(f20.10,',', f20.10))")temperature(i),C_v(i)
end do
close(103)

write(filename, '(A,I2,A)' )'m',size,'.csv'
open(104,file=filename)
do i=1,tpoints
 write(104,"(f20.10, (',', f20.10))")temperature(i),mag_per_site(i)
end do
close(104)

write(filename, '(A,I2,A)' )'xi',size,'.csv'
open(105,file=filename)
do i=1,tpoints
 write(105,"(9(f20.10,',', f20.10))")temperature(i),xi(i)
end do
close(105)

end subroutine output

