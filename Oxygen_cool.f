      !-------------------------------------------------------------
      ! function kinetic_energy: Calculate kinetic energy based on velocities
      !
      ! Input:
      !   natoms: number of atoms in the system
      !   vel1: (natoms,3) velocities of all atoms
      !   massO: Mass of 1 Ar atom
      !
      ! Output:
      !   ke: Total kinetic energy in the system
      !-------------------------------------------------------------
      function kinetic_energy(nMol,nAtomPer,velold,vel1,mass) result(ke)
        INTEGER, intent(in):: nMol        ! number of molecule
        INTEGER, intent(in):: nAtomPer        ! number of atoms per mol
        REAL*8, intent(in) :: velold(nMol,nAtomPer,3)  ! velocity array
        REAL*8, intent(in) :: vel1(nMol,nAtomPer,3)  ! velocity array

        REAL*8, intent(in):: mass           ! mass of 1 atom
        REAL*8             :: ke              ! output
        INTEGER I,J,K
        REAL*8 vel_mid
        ke=0.0

        ! Loop over all atoms
        DO I=1,nMol

          DO J=1,nAtomPer
            ! Loop over x,y,z components of velocity
            !
            DO K=1,3
                vel_mid=0.5*(velold(I,J,K)+vel1(I,J,K))
                ke=ke+vel_mid**2    ! sum vel^2
              END DO
            END DO
          END DO

        ke = 0.5*ke*mass          ! calculate kinetic energy

      end function


      !-------------------------------------------------------------
      ! function run_SHAKE: "Shake" bonds to fix bond lengths
      !
      ! Input:
      !   natoms:
      !   vel1: (natoms,3) velocities of all atoms
      !   massO: Mass of 1 Ar atom
      !
      ! Output:
      !   ke: Total kinetic energy in the system
      !-------------------------------------------------------------
      function run_SHAKE(nMol,nAPer,pos,pos_old,vel,
     :mass,tbl_SQ,dt) result(mass_term)
        INTEGER, intent(in):: nMol        ! number of molecule
        INTEGER, intent(in):: nAPer        ! number of atoms per mol
        REAL*8:: pos(nMol,nAPer,3)  ! position array
        REAL*8:: pos_old(nMol,nAPer,3)  ! last position array
        REAL*8:: vel(nMol,nAPer,3)  ! velocity array
        REAL*8, intent(in):: mass           ! mass of 1 atom
        REAL*8, intent(in)::tbl_SQ
        REAL*8, intent(in)::dt
        INTEGER I,JA,JB,K
        REAL*8 thisLength_SQ,dot_prod,g
        REAL*8 bond_vector(3), old_bond_vector(3)
        REAL*8 mass_term
        LOGICAL flag

        mass_term=((1.0/mass)+(1.0/mass))

        ! Loop over all molecules
        DO I=1,nMol
          dev_count=0
          flag=.TRUE.
          DO WHILE(flag) ! shake
            flag=.FALSE.

            ! Loop over all atom pairs in molecule
            DO JA=1,nAPer-1
               DO JB=JA+1,nAPer

                  thisLength_SQ=0.0
                  dot_prod=0.0

                  DO K=1,3
                     bond_vector(K)=pos(I,JA,K)-pos(I,JB,K)
                     thisLength_SQ=thisLength_SQ+bond_vector(K)**2

                     old_bond_vector(K)=pos_old(I,JA,K)-pos_old(I,JB,K)
                     dot_prod=dot_prod+bond_vector(K)*old_bond_vector(K)
                  END DO

                  g=tbl_SQ-thisLength_SQ

                  g=g/(2.0*mass_term*dot_prod)

                  thisLength_SQ=0.0
                  DO K=1,3
                     correction=g*old_bond_vector(K)/mass

                     pos(I,JA,K)=pos(I,JA,K)+correction
                     pos(I,JB,K)=pos(I,JB,K)-correction
                     bond_vector(K)=pos(I,JA,K)-pos(I,JB,K)
                     thisLength_SQ=thisLength_SQ+bond_vector(K)**2
                  END DO

                END DO ! JB loop
              END DO ! JA loop
            END DO ! do while loop
            ! recalculating velocities according to changes in shake
            ! positions

            DO J=1,nAPer

              DO K=1,3
                vel(I,J,K)=(pos(I,J,K)-pos_old(I,J,K))/dt
              END DO

            END DO  !J loop

          END DO ! I loop

        end function

      !-------------------------------------------------------------
      ! Argon_simulation
      !
      ! Simulate a set of argon atoms
      ! Based on code from Stephanie Kearing
      !
      ! Input:
      !   Argon_864_initial_coordinates.gro : file with atomic coordinates and
      !       velocities for to initiate simulation
      !
      ! Ouput:
      !   A bunch of files
      !-------------------------------------------------------------


      PROGRAM Oxygen_simulation

      ! Based on code from Stephanie Kearing
      IMPLICIT NONE

      INTEGER I,J,K,L,M,N,time_loop

      ! Set to NVE or NVT
      INTEGER is_NVT ! 0 if no thermostat (i.e. NVE); 1 for NVT


      ! Initiate constants
      REAL*8 NA
      PARAMETER (NA=6.02214076D23)

      REAL*8 BOLTZ
      PARAMETER (BOLTZ=1.38064852D-23) ! in J/K

      ! Declare functions
      REAL*8 kinetic_energy
      REAL*8 run_SHAKE
      ! Set constant parameters for our system
      INTEGER nMol ! the number of atoms in the simulation
      PARAMETER (nMol=500) ! From rahman64

      INTEGER nAtomPer ! the number of atoms per molecule
      PARAMETER (nAtomPer=2) ! From rahman64

      REAL*8 eps ! epsilon value for calculating LJ potential
      PARAMETER (eps=120.0*BOLTZ) ! epsilon in J, from rahman64

      REAL*8 sigma ! sigma value for calculating LJ potential
      PARAMETER (sigma=0.3006000E-9) ! sigma in m, from nm rahman64
      REAL*8 sigma2 ! sigma value for calculating LJ potential
      PARAMETER (sigma2=sigma**2) ! sigma in m, from nm rahman64

      REAL*8 massO ! mass of the particle
      PARAMETER (massO=0.016/NA) ! in kg from kg/mol rahman64

      REAL*8 massO2 ! mass of the particle
      PARAMETER (massO2=massO*2.0) ! in kg from kg/mol rahman64

      INTEGER nstep ! number of steps in the simulation
      PARAMETER (nstep=50000)!1000000

      INTEGER nsave ! frequency to save data
      PARAMETER (nsave=10)

      INTEGER comShiftTime              ! frequency to shift com of velocity
      PARAMETER (comShiftTime=110)

      REAL*8 DT ! time step
      PARAMETER (DT=4.0E-15) !  unit: s

      REAL*8 cutoff,cutoffSQ

      REAL*8 Tref1 !the reference temprerature for the thermostat
      PARAMETER (Tref1=320.00000) ! in K from Rahman64
      REAL*8 Tref2 !the reference temprerature for the thermostat
      PARAMETER (Tref2=77.00000) ! in K from Rahman64
      REAL*8 Temp_step
      PARAMETER (Temp_step=(Tref1-Tref2)/nstep)
      REAL*8 Tref

      REAL*8 Length ! Length of the box

      ! Create arrays for coordinates and velocities
      REAL*8 pos(nMol,nAtomPer,3),vel(nMol,nAtomPer,3) ! 2 atoms each, 3 for x,y,z directions
      REAL*8 vel_old(nMol,nAtomPer,3),vel_half(nMol,nAtomPer,3)
      REAL*8 pos_old(nMol,nAtomPer,3)

      ! calculating VLJ and FLJ in x,y,z directions
      REAL*8 rij(3) ! the vector between Ar I and J (I-J)
      REAL*8 dist_ij ! total distance between Ar I and J
      REAL*8 V_ij ! The LJ potential between Ar I and J
      REAL*8 V_tot ! The total potential energy in the system
      REAL*8 F_ij ! the magnitude of force between Ar I and J
      REAL*8 sr_12 ! (sigma/r)^12
      REAL*8 sr_6 ! (sigma/r)^6
      REAL*8 half_box ! half the length of the box
      REAL*8 eps4
      PARAMETER (eps4 = 4.0*eps)
      REAL*8 eps24
      PARAMETER (eps24=24.0*eps)

      INTEGER nDOF            ! Atomic DOF
      PARAMETER (nDOF=5)
      REAL*8 KE_Temp          ! Conversion between KE and Temp
      PARAMETER (KE_Temp = 2.0/(REAL(nDOF)*BOLTZ*(nMol)))

      ! Variables to attempt leap frog generation of new r,v,a
      REAL*8 force(nMol,nAtomPer,3) ! the force on each particle
      REAL*8 this_force(3)  ! force on this particle
      REAL*8 accel(nMol,nAtomPer,3) ! acceleration
      REAL*8 is_NVT_scale,KE,vel_2 ! for KE calculation


      ! Shake algorithm
      REAL*8 delr(3),delr_old(3),g,correction
      REAL*8 bond_length,new_delr
      PARAMETER (bond_length=0.1208E-9) ! m from Perng01
      REAL*8 blSQ
      PARAMETER (blSQ=bond_length**2)
      REAL*8 mass_term
      PARAMETER (mass_term=((1.0/massO)+(1.0/massO)))
      REAL*8 dot_prod ! the dot product of old and new distance vectors
      REAL*8 delr_magSQ,delr_old_mag
      REAL*8 max_X,max_Y,max_Z

      ! Variables for g(r)
      REAL*8 bin_ends(100)
      REAL*8 g_of_r(100)

      REAL*8 time ! the step multiplied by DT to give the time in sec
      REAL*8 temp ! the temperature of the system to be calc from KE

      ! preventing the flying ice cube
      REAL*8 vel_sum(3),vel_scale(3)

      ! reading in/values from input file
      INTEGER idum ! loop/array variables / dummy integer
      CHARACTER*5 resname,atomname1,atomname2

      ! File for storing energy data
      OPEN(92,FILE="Energies_O2_cool.txt")
      OPEN(91,FILE="O2_positions_cool.gro") ! open output traj file
      OPEN(94,FILE="Temp_cool.txt") ! open output temperature
      OPEN(95,FILE="Velocity_O2_cool.txt")
      OPEN(96,FILE="Velocity_vector_O2_cool.txt")
      OPEN(97,FILE="G_r_all_distances_O2_cool.txt")
      OPEN(98,FILE="Oxygen_cool.gro")

      !Set to NVT or NVE
      is_NVT=1
      Tref = Tref1

      cutoff=100*sigma
      cutoffSQ=cutoff**2 ! cutoffs for energy equations

      !----------------------
      ! Read in intial file
      !----------------------
      OPEN(11,FILE="Oxygen_NVT_320.gro") ! open input file
      READ(11,*)
      READ(11,*)
      DO I=1,nMol
          READ(11,31)idum,resname,atomname1,idum,
     :(pos(I,1,K),K=1,3),(vel(I,1,K),K=1,3)
          READ(11,31)idum,resname,atomname2,idum,
     :(pos(I,2,K),K=1,3),(vel(I,2,K),K=1,3)

      END DO
      READ(11,*)Length,Length,Length
      CLOSE(11)

      !--------------------------------------------
      ! unit conversion for pos/vel/box read in from the gro file; to SI
      ! Also calculate average velocity and initialize some variables
      !--------------------------------------------
      vel_sum(1) = 0
      vel_sum(2) = 0
      vel_sum(3) = 0
      Length=Length*1.0E-9 ! from nm to m

      ! Loop over all atoms
      DO I=1,nMol
        DO J=1,nAtomPer
          DO K=1,3                       ! Loop over x,y,z components
              pos(I,J,K)=pos(I,J,K)*1.0E-9    ! from nm to m
              vel(I,J,K)=vel(I,J,K)*1.0E2     ! from nm/ps to m/s
              accel(I,J,K) = 0.0            ! initialize acceleration to 0
              force(I,J,K) = 0.0            ! initialize forces to 0
              vel_sum(K)=vel_sum(K)+vel(I,J,K)  ! Calculate total velocity of system
              pos_old(I,J,K) = pos(I,J,K)
            END DO
         END DO
      END DO
      half_box = Length/2.0D0
      vel_scale(1) = vel_sum(1)/(nMol*nAtomPer)
      vel_scale(2) = vel_sum(2)/(nMol*nAtomPer)
      vel_scale(3) = vel_sum(3)/(nMol*nAtomPer)
      print*,Length


      ! Calculate kinetic energy and temperature
      KE = kinetic_energy(nMol, nAtomPer, vel,vel,massO)
      temp = (KE*KE_Temp)
      print*,temp
      is_NVT_scale = sqrt(Tref/temp)  ! Scale for velocity

      ! Loop over all atoms and scale velocity
      DO I=1,nMol
        DO J=1,nAtomPer
          DO K=1,3
            vel(I,J,K)=vel(I,J,K)*is_NVT_scale
          end do
        end do
      end do

      ! shift velocity to keep whole system from moving in a direction
      vel_sum(1) = 0
      vel_sum(2) = 0
      vel_sum(3) = 0
      ! Loop over all atoms
      DO I=1,nMol
        DO J=1,nAtomPer
          vel_2 = 0     ! new velocity of this atom
          DO K=1,3                           ! Loop over x,y,z components
            vel(I,J,K) = vel(I,J,K)-vel_scale(K) ! shift velocity
            vel_sum(K)=vel_sum(K)+vel(I,J,K)   ! Calculate new total velocity of system
            vel_2 = vel_2 + vel(I,J,K)*vel(I,J,K)
          END DO
        END DO
       END DO


       ! Initialize g_of_r data
       dist_ij = sigma*0.05D0
       DO I=1,100
         g_of_r(i) = 0.0
         bin_ends(i) = (0.05D0*sigma*I)**2
         dist_ij = dist_ij + sigma*0.05D0
       END DO

       ! Calculate kinetic energy and temperature
       KE = kinetic_energy(nMol, nAtomPer, vel,vel,massO)
       temp = (KE*KE_Temp)
       print*,temp

      !--------------------------------------------
      ! Progress system in time
      !--------------------------------------------
      time = 0.0
      DO time_loop=1,nstep
        time = time+dt                  ! Calculate current time

        !--------------------------------------------
        ! Calculate forces and potential energy
        !--------------------------------------------

        V_tot = 0.0                     ! Initialize total velocity to 0

        ! loop over all pairs of atoms
        DO I=1,nMol-1
          DO J=i+1,nMol
            DO M=1,2                  ! Atoms in molecule 1
              DO N=1,2                ! Atoms in molecule 2

                dist_ij = 0.0             ! reset dist_ij to 0


                ! Loop over x,y,z components
                DO K=1,3
                  ! get distance vector between atoms
                  rij(K)= pos(I,M,K)-pos(J,N,K)
                  rij(K)= rij(K)-length*ANINT(rij(K)/Length)

                  dist_ij = dist_ij + rij(K)*rij(K)
                END DO ! end K: calcualte dist_ij

                ! If within a certain distance, calcualte force
                if (dist_ij <= cutoffSQ) then

                  sr_6 = (sigma2/dist_ij)**3
                  sr_12 = sr_6**2

                  ! Calculate potential energy and add to total of system
                  V_ij=eps4*(sr_12-sr_6)
                  V_tot=V_tot+V_ij

                  ! Calculate magnitude of force
                  F_ij=-eps24*(-2.0*sr_12+sr_6)/dist_ij

                  ! Find unit vector of rij and use to add force vector to i and j
                  DO K=1,3
                    this_force(K) = F_ij*rij(K)
                    force(i,M,K) = force(i,M,K) + this_force(K)
                    force(j,N,K) = force(j,N,K) - this_force(K)
                  END DO ! end K: force vectors
                END if  ! end if dist_ij < cutoffSQ
              end do ! end N
            end do  ! end M
          END do    ! end J
        END DO      ! end I


        !--------------------------------------------
        ! Use the force to update the positions
        ! This is done using leap frog algorithm
        ! Remember that vel is really vel at t-0.5dt
        !--------------------------------------------
        max_X = 0.0
        max_Y = 0.0
        max_Z = 0.0
        DO I=1,nMol
          DO J=1,nAtomPer
           DO K=1,3
             pos_old(I,J,K) = pos(I,J,K)
             accel(I,J,K) =force(I,J,K)/massO                  ! Acceleration based on force
             vel_old(I,J,K) = vel(I,J,K)                        ! Store old velocity
             vel(I,J,K)=vel(I,J,K) + DT*accel(I,J,K)              ! Get new velocity (at t+0.tdt)
             pos(I,J,K) = pos(I,J,K) + DT*vel(I,J,K)              ! Get new position (at t+dt)
             vel_half(I,J,K) = (vel(I,J,K) + vel_old(I,J,K))*0.5  ! Get velocity at (t+dt) for KE
             force(I,J,K) = 0.0                               ! reset force array for next time loop
           END DO
          end do

        END DO

        ! Run SHAKE
        temp = run_SHAKE(nMol,nAtomPer,pos,pos_old,vel,massO,blSQ,dt)

        ! Calculate kinetic energy and temperature
        KE = kinetic_energy(nMol,nAtomPer,vel_old,vel,massO)
        temp = (KE*KE_Temp)
        WRITE(94,*)temp

        !--------------------------------------------
        ! If we're in NVT, scale velocity to keep set temperature
        !--------------------------------------------
        IF (is_NVT == 1) THEN

            is_NVT_scale = sqrt(Tref/temp)  ! Scale for velocity
            ! Modify target temp
            Tref = Tref - Temp_step

            ! Loop over all atoms and scale velocity
            DO I=1,nMol
              DO J=1,nAtomPer
                DO K=1,3
                  vel(I,J,K)=vel(I,J,K)*is_NVT_scale
                end do
              end do
            end do

          !---------------------------------------------------------------------
          ! Every 100 time points, shift com of veloctiy to be 0
          !---------------------------------------------------------------------
          IF(MOD(time_loop,comShiftTime).EQ.0) THEN
            vel_sum(1) = 0.0
            vel_sum(2) = 0.0
            vel_sum(3) = 0.0

            ! Loop over all atoms and find center of mass of velocties
            DO I=1,nMol
              DO J=1,nAtomPer
               DO K=1,3                           ! Loop over x,y,z components
                  vel_sum(K)=vel_sum(K)+vel(I,J,K)  ! Calculate total velocity of system
               END DO
             END DO
            END DO

            ! Calculate how much system needs to shift to keep from flying ice cube
            vel_sum(1) = vel_sum(1)/(nAtomPer*nMol)
            vel_sum(2) = vel_sum(2)/(nAtomPer*nMol)
            vel_sum(3) = vel_sum(3)/(nAtomPer*nMol)

            ! shift velocity to keep whole system from moving in a direction
            DO I=1,nMol                         ! Loop over all atoms
              DO J=1,nAtomPer
               DO K=1,3                           ! Loop over x,y,z components
                 vel(I,J,K) = vel(I,J,K)-vel_sum(K) ! shift velocity
               END DO
              END DO
            END DO
          END IF

        end if

        !--------------------------------------------
        ! Store system data if we're in a save time
        !--------------------------------------------
        IF(MOD(time_loop,nsave).EQ.0) THEN
           resname='   O'
           IF(is_NVT <2) THEN
             WRITE(91,*)'After step ',time_loop
             WRITE(91,*)(nMol*nAtomPer)

             DO I=1,nMol
                WRITE(91,31)I,resname,atomname1,I,
     :(pos(I,1,K)*1.0E9,K=1,3),
     :(vel(I,1,K)*1.0E-4,K=1,3)
                 WRITE(91,31)I,resname,atomname2,I,
     :(pos(I,2,K)*1.0E9,K=1,3),
     :(vel(I,2,K)*1.0E-4,K=1,3)
              END DO

            WRITE(91,*)Length*1.0E9,Length*1.0E9,Length*1.0E9
            END IF
        END IF
        WRITE(92,*)time,V_tot,KE,V_tot+KE,temp
        
      END DO

      ! Write g_of_r data
      DO I=1,100
        WRITE(97,*)bin_ends(I),g_of_r(I)
      END DO

      WRITE(95,*)'After step ',time_loop
      WRITE(95,*)(nMol*nAtomPer)

      DO I=1,nMol
         WRITE(95,31)I,resname,atomname1,I,
     :(vel(I,1,K)*1.0E-3,K=1,3),
     :(vel(I,1,K)*1.0E-4,K=1,3)
          WRITE(95,31)I,resname,atomname2,I,
     :(vel(I,2,K)*1.0E-3,K=1,3),
     :(vel(I,2,K)*1.0E-4,K=1,3)
       END DO

       WRITE(91,*)Length*1.0E9,Length*1.0E9,Length*1.0E9

      ! Write final coordinates
      WRITE(98,*)'A box of liquid O2'
      WRITE(98,*)(nMol*nAtomPer)
      DO I=1,nMol
        WRITE(98,31)I,resname,atomname1,I,
     :(pos(I,1,K)*1.0E9,K=1,3),
     :(vel(I,1,K)*1.0E-2,K=1,3)
        WRITE(98,31)I,resname,atomname2,I,
     :(pos(I,2,K)*1.0E9,K=1,3),
     :(vel(I,2,K)*1.0E-2,K=1,3)
      END DO
      WRITE(98,*)Length*1.0E9,Length*1.0E9,Length*1.0E9

  31  FORMAT(i5,2a5,i5,3f8.4,3f8.4) ! format of input file
      print*,Length
      CLOSE(92)
      CLOSE(91)
      CLOSE(94)
      CLOSE(95)
      CLOSE(96)
      CLOSE(97)
      END
