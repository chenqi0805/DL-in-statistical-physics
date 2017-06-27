!This program performs a Monte Carlo simulation for the 2D Ising model on a square lattice
!The temperature is included in the coupling strength
!J=K/(k_BT) where K is the bare coupling strength and
!H=B/(k_BT) where B is the magnetic field
program Ising
include 'globIsMC'
integer i

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

call output

end program Ising

subroutine Initialise
include 'globIsMC'
integer i,ix,iy

H=0.0d0
Size=20
MCSteps=100000
MCSteps_eq=1000
do ix = 1, Size
 do iy = 1, Size      
  Spins(ix, iy) = 1
 end do
end do
Magnetisation = Size**2
mag_av=0.0d0
mag_square_av=0.0d0
Energy = -2.0d0*Size*Size*J
Energy_av=0.0d0
Energy_square_av=0.0d0
do i=1,5
 MCWeights(i) = dexp(-dble(2*(i-1))*J)
end do

end subroutine Initialise

subroutine MonteCarlo
include 'globIsMC'

integer mciter, ix,iy,costfactor,curspin

!equilibrate
do mciter=1,mcsteps_eq
 do ix=0,size-1
  do iy=0,size-1
   curspin=Spins(ix+1,iy+1)
   costfactor=Curspin*(Spins(mod(ix+1,size)+1,iy+1)&
   +Spins(mod(ix-1+size,size)+1,iy+1)&
   +Spins(ix+1,mod(iy+1,size)+1)&
   +Spins(ix+1,mod(iy-1+size,size)+1))
   if (costfactor>0) then
    randnumber=grnd()
    if (randnumber<MCWeights(costfactor+1)) then
     spins(ix+1,iy+1)=-spins(ix+1,iy+1)
     magnetisation=magnetisation+2*spins(ix+1,iy+1)
     energy=energy+dble(2*costfactor)*J
    end if
   else
    spins(ix+1,iy+1)=-spins(ix+1,iy+1)
    magnetisation=magnetisation+2*spins(ix+1,iy+1)
    energy=energy+dble(2*costfactor)*J
   end if
  end do
 end do
end do
!MC steps after equilibration
do mciter=1,mcsteps
 do ix=0,size-1
  do iy=0,size-1
   curspin=Spins(ix+1,iy+1)
   costfactor=Curspin*(Spins(mod(ix+1,size)+1,iy+1)&
                      +Spins(mod(ix-1+size,size)+1,iy+1)&
                      +Spins(ix+1,mod(iy+1,size)+1)&
                      +Spins(ix+1,mod(iy-1+size,size)+1))
   if (costfactor>0) then
    randnumber=grnd()
    if (randnumber<MCWeights(costfactor+1)) then
     spins(ix+1,iy+1)=-spins(ix+1,iy+1)
     magnetisation=magnetisation+2*spins(ix+1,iy+1)
     energy=energy+dble(2*costfactor)*J
    end if
   else
    spins(ix+1,iy+1)=-spins(ix+1,iy+1)
    magnetisation=magnetisation+2*spins(ix+1,iy+1)
    energy=energy+dble(2*costfactor)*J
   end if
  end do
 end do
 mag_av=mag_av+dble(magnetisation)
 energy_av=energy_av+energy
 Energy_square_av=Energy_square_av+(energy/dble(size**2))**2
 mag_square_av=mag_square_av+(dble(magnetisation)/dble(size**2))**2
end do

!average energy persite
energy_av=energy_av/dble(size*size*mcsteps)
Energy_square_av=Energy_square_av/dble(mcsteps)
!average magnetisation
mag_av=mag_av/dble(size*size*mcsteps)
mag_square_av=mag_square_av/dble(mcsteps)

end subroutine MonteCarlo

subroutine output
include 'globIsMC'
integer i

open(102,file="en.dat")
do i=1,tpoints
 write(102,"(9(1x, f20.10, 5x))")temperature(i),en_per_site(i)*temperature(i)
end do
close(102)

open(103,file="Cv.dat")
do i=1,tpoints
 write(103,"(9(1x, f20.10, 5x))")temperature(i),C_v(i)
end do
close(103)

open(104,file="m.dat")
do i=1,tpoints
 write(104,"(9(1x, f20.10, 5x))")temperature(i),mag_per_site(i)
end do
close(104)

open(105,file="xi.dat")
do i=1,tpoints
 write(105,"(9(1x, f20.10, 5x))")temperature(i),xi(i)
end do
close(105)

end subroutine output

