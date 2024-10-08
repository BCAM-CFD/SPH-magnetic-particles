      SUBROUTINE colloid_init_accumulation_matrix(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_init_accumulation_matrix
        !----------------------------------------------------
        !
        ! Purpose     : Using current accumulative rotaiton 
        !               vector to compute
        !               accumulative rotation matrix for colloid.
        !
        ! Referecen   : Chen et al. 
        !               Physics of Fluids, 18, 103605, 2006.
        !
        ! Revision    : V.01  2.12.2011
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
        
        TYPE(Colloid), INTENT(OUT)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: i,dim,stat_info_sub
        REAL(MK), DIMENSION(3,3)        :: rot_matrix
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        dim           = this%num_dim        
        
        DO i = 1, this%num_colloid
           
           rot_matrix(:,:) = 0.0_MK

           CALL tool_rotation_matrix(this%tool,&
                dim,this%acc_vector(1:3,i),this%acc_vector(4,i),&
                rot_matrix(1:3,1:3),stat_info_sub)
           
           this%acc_matrix(:,:,i) = rot_matrix(:,:)
            
        END DO ! i = 1, num_colloid
        
           
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_init_accumulation_matrix
      
