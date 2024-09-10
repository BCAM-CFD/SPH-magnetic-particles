!************ Subroutine changed by Adolfo (magnetic). The subroutine used
!             for the rotating magnetic field is after this one. We have also added
!             the dipole-dipole torque. The original
!             subroutine programmed by Xin Bian is at the end of the file *******************
      SUBROUTINE colloid_compute_magnetism_cc(this,&
           e_ij, r, sid_ip, sid_jp, F_ij, T_ij, T_ji, stat_info)
        !***** Check Landecker 1999, and careful with some signs ******
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_magnetism_cc
        !----------------------------------------------------
        !
        ! Purpose     : 
        !
        ! Routines    :
        !
        ! References  : 
        !
        ! Remarks     : V0.1 18.06 2013, original version,
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
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:),INTENT(IN)       :: e_ij
        REAL(MK), INTENT(INOUT)                 :: r
        INTEGER, INTENT(IN)                     :: sid_ip
        INTEGER, INTENT(IN)                     :: sid_jp
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: F_ij
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: T_ij
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: T_ji
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: dim
        REAL(MK)                        :: a, F0, T0
        REAL(MK)                        :: c1
        REAL(MK), DIMENSION(3)          :: mom_i, mom_j
        REAL(MK)                        :: mi_dot_e, mj_dot_e, m_dot_m
        REAL(MK), DIMENSION(3)          :: mi_vec_e, mj_vec_e
        REAL(MK), DIMENSION(3)          :: mi_vec_mj
       
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        
        dim = this%num_dim
        
        a   = this%radius(1,sid_ip)
        F0 = this%cc_magnet_F0
        T0 = F0 / 3.0_MK 

!        !-- We consider the same moment for both particles --
!        mom(1:dim) = this%cc_magnet_mom(1:dim) 

        !-- We consider that each particle has a different magnetic moment --
        mom_i(1:dim) = this%magnet_mom(1:dim, sid_ip)
        mom_j(1:dim) = this%magnet_mom(1:dim, sid_jp)
                
        !----------------------------------------------------
        ! Calculate the magnetic force 
        !----------------------------------------------------

        ! **** Note that in Landecker 1999, the direction of r is r_{beta,alpha} here (or r_ji)
        mi_dot_e = DOT_PRODUCT(mom_i(1:dim), e_ij(1:dim))
        mj_dot_e = DOT_PRODUCT(mom_j(1:dim), e_ij(1:dim))
        m_dot_m   = DOT_PRODUCT(mom_i(1:dim), mom_j(1:dim))
        !-- mi_vec_e is m_i x e_ij
        mi_vec_e(1) = mom_i(2) * e_ij(3) - mom_i(3) * e_ij(2)
        mi_vec_e(2) = mom_i(3) * e_ij(1) - mom_i(1) * e_ij(3)
        mi_vec_e(3) = mom_i(1) * e_ij(2) - mom_i(2) * e_ij(1)
        !-- mj_vec_e is m_j x e_ij
        mj_vec_e(1) = mom_j(2) * e_ij(3) - mom_j(3) * e_ij(2)
        mj_vec_e(2) = mom_j(3) * e_ij(1) - mom_j(1) * e_ij(3)
        mj_vec_e(3) = mom_j(1) * e_ij(2) - mom_j(2) * e_ij(1)
        !-- mi_vec_mj is m_i x m_j
        mi_vec_mj(1) = mom_i(2) * mom_j(3) - mom_i(3) * mom_j(2)
        mi_vec_mj(2) = mom_i(3) * mom_j(1) - mom_i(1) * mom_j(3)
        mi_vec_mj(3) = mom_i(1) * mom_j(2) - mom_i(2) * mom_j(1)
        
        !**** Force on i by j ***** F0 = 3 * mu0 / (4 * pi)
        F_ij(1:dim) = -F0 / r**4.0_MK * (mi_dot_e * mom_j(1:dim) + &
             mj_dot_e * mom_i(1:dim) - &
             (5.0_MK * mj_dot_e * mi_dot_e - m_dot_m)*e_ij(1:dim))

        !**** Torque *****
        ! T0 = mu0 / (4 * pi)
        ! Torque on j by i
        T_ij(1:dim) = T0 / r**3.0_MK * (3.0_MK * mi_dot_e * mj_vec_e(1:dim) + &
             mi_vec_mj(1:dim))
        ! Torque on i by j (Note that mi x mj = -mj x mi)
        T_ji(1:dim) = T0 / r**3.0_MK * (3.0_MK * mj_dot_e * mi_vec_e(1:dim) - &
             mi_vec_mj(1:dim))
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_compute_magnetism_cc
!***********************************************************

!************ Subroutine changed by Adolfo (magnetic). The original is
!       after this one *******************
!************ This subroutine was used for the paper with a rotating magnetic field ************
      SUBROUTINE colloid_compute_magnetism_cc_rotating(this,&
           e_ij,r,sid_ip,F_ij,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_magnetism_cc
        !----------------------------------------------------
        !
        ! Purpose     : 
        !
        ! Routines    :
        !
        ! References  : 
        !
        ! Remarks     : V0.1 18.06 2013, original version,
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
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:),INTENT(IN)       :: e_ij
        REAL(MK), INTENT(INOUT)                 :: r
        INTEGER, INTENT(IN)                     :: sid_ip
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: F_ij
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: dim
        REAL(MK)                        :: a, F0
        REAL(MK)                        :: c1
        REAL(MK), DIMENSION(3)          :: mom
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        
        dim = this%num_dim
        
        a   = this%radius(1,sid_ip)
        F0 = this%cc_magnet_F0
        !-- We consider the same moment for both particles --
        mom(1:dim) = this%cc_magnet_mom(1:dim) 

        F_ij(1:dim) = 0.0_MK
        
        
        !----------------------------------------------------
        ! Calculate the magnetic force 
        !----------------------------------------------------

        r          = r / a
           
        c1 = DOT_PRODUCT(mom(1:dim), e_ij(1:dim))
        
        F_ij(1:dim) = F0 / r**4.0_MK * (2.0_MK * c1 * mom(1:dim) - &
             (5.0_MK * c1 * c1 - 1.0_MK) * e_ij(1:dim))
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_compute_magnetism_cc_rotating
!********************************************************************

      
      SUBROUTINE colloid_compute_magnetism_cc_old(this,&
           x_ip,x_jp,sid_ip,sid_jp,F_ij,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_magnetism_cc
        !----------------------------------------------------
        !
        ! Purpose     : 
        !
        ! Routines    :
        !
        ! References  : Sing et al. PNAS 2009, Supporting info
        !               [S8].
        !
        ! Remarks     : V0.1 18.06 2013, original version,
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
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:),INTENT(IN)       :: x_ip
        REAL(MK), DIMENSION(:),INTENT(IN)       :: x_jp
        INTEGER, INTENT(IN)                     :: sid_ip
        INTEGER, INTENT(IN)                     :: sid_jp
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: F_ij
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: dim,num
        REAL(MK)                        :: a, aa, r, F0, cut_off, cut_on, s
        REAL(MK)                        :: vc, mu, cof
        REAL(MK)                        :: c1, c2, c3, c4
        REAL(MK), DIMENSION(3)          :: r12
        REAL(MK), DIMENSION(3)          :: mom1,mom2
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        
        dim = this%num_dim
        num = this%num_colloid
        
        a   = this%radius(1,sid_ip)
        aa  = 2.0_MK*a
        !mu  = this%cc_magnet_mu
        !vc  = 4.0_MK*mcf_pi*a**3*this%cc_magnet_f/3.0_MK
        !cof =  vc * this%cc_magnet_chi/ mu
        !mom1(1:dim) = cof * this%cc_magnet_B(1:dim)
        F0 = this%cc_magnet_F0
        mom1(1:dim) = this%cc_magnet_mom(1:dim)
        mom2(1:dim) = mom1(1:dim)
        cut_off     = this%cc_magnet_cut_off
        cut_on      = this%cc_magnet_cut_on
        
        F_ij(1:dim) = 0.0_MK
        
        
        !----------------------------------------------------
        ! Calculate the magnetic force 
        !----------------------------------------------------

        r12(1:dim) = x_ip(1:dim) - x_jp(1:dim)
        r = SQRT(DOT_PRODUCT(r12(1:dim), r12(1:dim)))
        
        s = r - aa
        
        IF ( s <= cut_off ) THEN
           
           IF ( s <= cut_on ) THEN
              
              r12(1:dim) = r12(1:dim) * (aa+cut_on)/ r
              r = aa + cut_on
              
           END IF

           r12(1:dim) = r12(1:dim)/a
           r          = r/a
           
           c1 = DOT_PRODUCT(mom1(1:dim), r12(1:dim))/r
           c2 = DOT_PRODUCT(mom2(1:dim), r12(1:dim))/r
           c3 =  DOT_PRODUCT(mom1(1:dim), mom2(1:dim))
           
           F_ij(1:dim) = F0 / r**4 * ( c1 * mom2(1:dim) + c2*mom1(1:dim) &
                - (5.0_MK*c1*c2-c3)* r12(1:dim) / r )
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_compute_magnetism_cc_old
      
      
      
