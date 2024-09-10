!******* Subroutine added by Adolfo (magnetic) **********
SUBROUTINE colloid_compute_magnetic_torque(this, time_current, magnetic_B_amp, &
     coll_torque, stat_info)
        !----------------------------------------------------
        ! Compute the magnetic torque due to the misalignment between the magnetic moment
        ! and the magnetic field.
        !     T = m x B
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), INTENT(IN)                    :: time_current
        REAL(MK), INTENT(IN)                    :: magnetic_B_amp
        REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: coll_torque
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------

        INTEGER                                 :: stat_info_sub
        
        REAL(MK), DIMENSION(3)          :: B
        REAL(MK) :: alpha
        INTEGER :: num, dim
        INTEGER :: I
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        dim   = this%num_dim
        num   = this%num_colloid

        !-- Alternating field --
        alpha = this%cc_magnet_rot_freq * time_current        
        !-- Following is counterclockwise, initially along z-axis --
        B(1) = -magnetic_B_amp * sin(alpha)
        B(2) = 0.0_MK
        B(3) = magnetic_B_amp * cos(alpha)

        !---- T = m x B -----
        DO I = 1, num
           coll_torque(1, I) = coll_torque(1, I) +  this%magnet_mom(2,I) * B(3) - &
                this%magnet_mom(3,I) * B(2)
           coll_torque(2, I) = coll_torque(2, I) +  this%magnet_mom(3,I) * B(1) - &
                this%magnet_mom(1,I) * B(3)
           coll_torque(3, I) = coll_torque(3, I) +  this%magnet_mom(1,I) * B(2) - &
                this%magnet_mom(2,I) * B(1)

        ENDDO


9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_compute_magnetic_torque
      !********************************************************

!****** Subroutine changed by Adolfo (magnetic). The old one is after this one. ****

      SUBROUTINE colloid_compute_magnetic_moment(this, dt, stat_info)
        !----------------------------------------------------
        ! Subroutine to calculate the evolution of the magnetic moment as
        !      dm/dt = Omega x m
        !-----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), INTENT(IN)                    :: dt        
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------

        INTEGER                                 :: stat_info_sub
        REAL(MK), DIMENSION(3) :: dm_dt
        REAL(MK) :: current_magnitude
        INTEGER :: num, dim
        INTEGER :: I

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        dim   = this%num_dim
        num   = this%num_colloid


        !-----  dm/dt = Omega x m --------
        DO I = 1, num
           dm_dt(1) = this%omega(2,I,1) * this%magnet_mom(3,I) - &
                this%omega(3,I,1) * this%magnet_mom(2,I)
           dm_dt(2) = this%omega(3,I,1) * this%magnet_mom(1,I) - &
                this%omega(1,I,1) * this%magnet_mom(3,I)
           dm_dt(3) = this%omega(1,I,1) * this%magnet_mom(2,I) - &
                this%omega(2,I,1) * this%magnet_mom(1,I)

           this%magnet_mom(:,I) = this%magnet_mom(:,I) + dm_dt(:) * dt

           !-- The magnetic moment is renormalized to avoid numerical errors --
           current_magnitude = SQRT(DOT_PRODUCT(this%magnet_mom(:,I), this%magnet_mom(:,I)))
           this%magnet_mom(:,I) = this%magnet_mom(:,I) / current_magnitude * this%magnet_mom_magnitude(I)
        ENDDO

9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_compute_magnetic_moment
      !********************************************************

      SUBROUTINE colloid_compute_magnetism_moment_old(this, dt, stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_magnetism_moment
        !----------------------------------------------------
        !
        ! Purpose     : Compute the magnetism moment
        !               using the rotating frequency of 
        !               the magnetic field.
        !               
        !
        ! Reference   :
        !
        ! Remark      : 
        !
        ! Revision    : V0.1  26.06.2013, original version.
        !
        !----------------------------------------------------
        ! This code is  based on the original MCF code  developed by Xin Bian.
        ! The  current version  has  been developed  in collaboration  between
        ! Marco Ellero,  leader of the  CFD Modelling and Simulation  group at
        ! BCAM (Basque Center  for Applied Mathematics) in  Bilbao, Spain, and
        ! Adolfo Vazquez-Quesada from  the Department of Fundamental Physics
        ! at UNED, in Madrid, Spain.
        !
        ! Developers:
        !     Xin Bian.
        !     Adolfo Vazquez-Quesada.
        !
        ! Contact: a.vazquez-quesada@fisfun.uned.es
        !          mellero@bcamath.org
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !
        ! this           : an object of Particles Class.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), INTENT(IN)                    :: dt
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------

        INTEGER                                 :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        CALL colloid_compute_magnetism_rotation_vector(this,dt,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, __FILE__, __LINE__, &
           "computing magnetism rotation vector failed!"
           stat_info = -1
           GOTO 9999
        END IF

        CALL colloid_compute_magnetism_rotation_matrix(this,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, __FILE__, __LINE__, &
                "computing magnetism rotation matrix failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%cc_magnet_mom(1:3) = &
             MATMUL(this%cc_magnet_rot_matrix(1:3,1:3),this%cc_magnet_mom(1:3))
        
        !PRINT *, "mom: ", this%cc_magnet_mom(1:3)

        CALL colloid_compute_magnetism_accumulation_matrix(this,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, __FILE__, __LINE__, &
           "computing magnetism accumulation matrix failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        CALL colloid_compute_magnetism_accumulation_vector(this,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, __FILE__, __LINE__, &
                "computing magnetism accumulation vector failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_compute_magnetism_moment_old
      
      
