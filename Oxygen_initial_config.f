      FUNCTION  random(iseed)
      ! When first call, iseed must be a large positive integer.
      ! iseed will be changed when exit and be used for next calling.
      ! The range of the generated random number is between 1 and -1
      implicit none
      integer, intent(inout) :: iseed
      real :: random
      iseed = mod(8121*iseed+28411, 134456) ! 0 =< iseed < 134456
      random = real(iseed)/134456. ! 0 < random < 1
      end FUNCTION random

      PROGRAM main
      IMPLICIT NONE

      INTEGER I,J,K,L,M,N

      ! Initiate constants
      REAL*8 NA
      PARAMETER (NA=6.02214076D23)
      REAL*8 PI
      PARAMETER (PI=4.0*atan(1.0))
      REAL*8 BOLTZ
      PARAMETER (BOLTZ=1.38064852D-23) ! in J/K

      ! Set initial parameters for our system
      INTEGER natom ! the number of atoms in the simulation
      PARAMETER (natom=500) !

      REAL*8 massAr ! mass of the particle
      PARAMETER (massAr=0.016/NA) ! in kg from kg/mol rahman64

      REAL*8 density
      PARAMETER (density=22.664) ! in nm^3

      REAL*8 Tref !the reference temprerature for the thermostat
      PARAMETER (Tref=320.00000) ! in K

      REAL*8 Length ! Length of the box
      PARAMETER (Length=(natom/density)**(1.0/3.0)) !

      INTEGER nDOF
      PARAMETER (nDOF=5)

      REAL vrms ! target average final_velocity
      PARAMETER (vrms = sqrt(3.0*BOLTZ*Tref/massAr)/1000) ! velocity in nm/ps

      REAL*8 sigma ! sigma value for calculating LJ potential
      PARAMETER (sigma=0.3006) ! sigma in nm, from rahman64

      REAL*8 dist_part ! distance between particles in unit cell
      PARAMETER (dist_part=Length/10.0)!sigma/sqrt(2.0))

      REAL*8 dist_unit ! distance between unit cells
      PARAMETER (dist_unit=Length/5.0)

      REAL*8 bond_length
      PARAMETER (bond_length = 0.1208) !nm
      INTEGER repeats
      PARAMETER (repeats = 5)

      ! Create unit cell
      REAL*8 unit_cell(4,3) ! sum of the total force on each particle


      REAL phi,theta,costheta ! for generating random vectors


      REAL scaled ! the value that the random vectors are scaled by

      REAL vx,vy,vz !velocity vectors in the x,y,z directions
      REAL rand_mag ! the magnitude of the x,y,z velocity vectors

      CHARACTER*5 resname,atomname1,atomname2
      PARAMETER (resname='   O',atomname1='  Ox1',atomname2='  Ox2')

      INTEGER icount,rescount
      REAL x,y,z,dx1,dy1,dz1
      REAL dx,dy,dz,total_length
      REAL random
      INTEGER iseed

      print*,Length
      print*,natom/density
      print*,vrms
      ! Create unit cell
      DO I=1,4
         DO J=1,3
           unit_cell(I,J) = 0.0
         END DO
      END DO
      unit_cell(1,2) = dist_part
      unit_cell(1,3) = dist_part
      unit_cell(3,1) = dist_part
      unit_cell(3,3) = dist_part
      unit_cell(4,1) = dist_part
      unit_cell(4,2) = dist_part


      iseed=10

      OPEN(91,FILE="Oxygen_500_initial_coordinates.gro") !open output file

      WRITE(91,*)'A box of liquid O2'
      WRITE(91,*)natom*2

      OPEN(92,FILE="Oxygen_500_initial_velocity.gro") !open output file
      WRITE(92,*)'A box of liquid O2'
      WRITE(92,*)natom*2

      icount=1      ! stores which atom we are on
      rescount=1
      ! Loop over all repeats of the unit cell
      DO I=1,repeats
         dx=REAL(I-1)*dist_unit
         DO J=1,repeats
            dy=REAL(J-1)*dist_unit
            DO K=1,repeats
              dz=REAL(K-1)*dist_unit

              ! Now that we know the displacements of the unit cell, move each atom
              DO M=1,4
                x=dx+unit_cell(M,1)
                y=dy+unit_cell(M,2)
                z=dz+unit_cell(M,3)

                ! Get random velocity
                phi = random(iseed)*PI*2.0
                costheta = -1.0 + random(iseed)*2.0

                theta = acos( costheta )
                vx = sin( theta) * cos( phi )
                vy = sin( theta) * sin( phi )
                vz = cos( theta )

                ! calculate the magnitude of the random x,y,z vectors
                rand_mag = sqrt((vx**2)+(vy**2)+(vz**2))

                ! calculate how much the random vector needs to be scaled
                scaled=vrms/rand_mag

                ! Scale the random vector
                vx=vx*vrms ! scale the x vector
                vy=vy*vrms ! scale the y vector
                vz=vz*vrms ! scale the z vector


                WRITE(91,31)rescount,resname,atomname1,icount,x,y,z,vx,
     :vy,vz
                WRITE(92,31)rescount,resname,atomname1,icount,10.0*vx,
     :10.0*vy,10.0*vz,vx,vy,vz
                icount=icount+1

                ! Add second molecule
        
                x = x+bond_length/sqrt(3.0)
                y = y+bond_length/sqrt(3.0)
                z = z+bond_length/sqrt(3.0)
                ! Get random velocity
                phi = random(iseed)*PI*2.0
                costheta = -1.0 + random(iseed)*2.0

                theta = acos( costheta )
                vx = sin( theta) * cos( phi )
                vy = sin( theta) * sin( phi )
                vz = cos( theta )

                ! calculate the magnitude of the random x,y,z vectors
                rand_mag = sqrt((vx**2)+(vy**2)+(vz**2))

                ! calculate how much the random vector needs to be scaled
                scaled=vrms/rand_mag

                ! Scale the random vector
                vx=vx*vrms ! scale the x vector
                vy=vy*vrms ! scale the y vector
                vz=vz*vrms ! scale the z vector

                WRITE(91,31)rescount,resname,atomname2,icount,x,y,z,vx,
     :vy,vz
                WRITE(92,31)rescount,resname,atomname2,icount,10.0*vx,
     :10.0*vy,10.0*vz,vx,vy,vz
                icount=icount+1
                rescount = rescount+1

              END DO
          END DO
        END DO
      END DO

  51  WRITE(91,*)Length,Length,Length
      WRITE(92,*)Length,Length,Length

  31  FORMAT(i5,2a5,i5,3f8.4,3f8.4)

      CLOSE(91)
      CLOSE(92)

      END
