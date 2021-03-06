module modGraviton
  use ModParameters
  implicit none
  private
  integer, parameter  :: dp = KIND(1d0)
  real(dp), private, parameter :: tol = 0.00000010_dp
  
  logical :: goodhel_gg(1:16), goodhel_qq(1:16)
  integer :: ncalls_gg = 0, ncalls_qq = 0
  integer, private, parameter :: maxcalls = 100
  
  !----- notation for subroutines
  public :: EvalAmp_gg_G_VV,EvalAmp_qqb_G_VV !! EvalAmp_gg_G_VV_old
  
contains
  
  
  !----- a subroutinefor gg -> G -> ZZ/WW
  !----- all outgoing convention and the following momentum assignment
  !-----  0 -> g(p1) + g(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
  subroutine EvalAmp_gg_G_VV(p,MY_IDUP,sum)                         ! modify p -> pp
    use ModMisc
    implicit none
    real(dp), intent(out) ::  sum
    real(dp), intent(in) :: p(4,6)
    integer, intent(in) :: MY_IDUP(6:9)
    real(dp) :: s, pin(4,4)
    complex(dp) :: A(2), sp(4,4), propG, propZ1, propZ2
    integer :: i1,i2,i3,i4,ihel
    real(dp) :: aL1,aR1,aL2,aR2
    real(dp) :: gZ_sq
    real(dp) :: prefactor, Lambda_inv
    
    if( ncalls_gg .eq. 0 ) then
       do ihel=1,16
          goodhel_gg(ihel)=.false.
       enddo
    endif
    if( ncalls_gg .le. maxcalls ) then
       ncalls_gg = ncalls_gg + 1
    endif
    
    gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2
    
    !---- the 1/Lambda coupling
    Lambda_inv = 1.0d0/Lambda
    
    !---- full prefactor; 8 is  the color factor
    prefactor = 8d0*(Lambda_inv**2)**2*gZ_sq**2
    
    
    if( DecayMode1.le.3 ) then!  Z decay
       if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
          aL1=aL_lep
          aR1=aR_lep
       elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
          aL1=aL_neu
          aR1=aR_neu
       elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
          aL1=aL_QUp
          aR1=aR_QUp
       elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
          aL1=aL_QDn
          aR1=aR_QDn
       else
          aL1=0d0
          aR1=0d0
       endif
       prefactor = prefactor *(one/two*M_V*Ga_V)**2
       s = 2d0 * scr(p(:,3),p(:,4))
       propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
       s = 2d0 * scr(p(:,5),p(:,6))
       propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
    elseif( DecayMode1.ge.4 .and. DecayMode1.le.6 ) then !  W decay
       aL1 = bL
       aR1 = bR
       prefactor = prefactor *(one/two*M_V*Ga_V)**2
       s = 2d0 * scr(p(:,3),p(:,4))
       propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
       s = 2d0 * scr(p(:,5),p(:,6))
       propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
    elseif( DecayMode1.eq.7 ) then !  photon "decay"
       aL1=1d0
       aR1=1d0
       propZ1 = 1d0
       propZ2 = 1d0
       prefactor = prefactor/gZ_sq**2! cancel the overall z coupling
    else
       aL1=0d0
       aR1=0d0            
    endif
    
    if( DecayMode2.le.3 ) then!  Z decay
       if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) .or. abs(MY_IDUP(8)).eq.abs(TaM_) ) then
          aL2=aL_lep
          aR2=aR_lep
       elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
          aL2=aL_neu
          aR2=aR_neu
       elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
          aL2=aL_QUp
          aR2=aR_QUp
       elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
          aL2=aL_QDn
          aR2=aR_QDn
       else
          aL2=0d0
          aR2=0d0
       endif
    elseif( DecayMode2.ge.4 .and. DecayMode2.le.6 ) then !  W decay
       aL2 = bL
       aR2 = bR
    elseif( DecayMode2.eq.7 ) then !  photon "decay"
       aL2=1d0
       aR2=1d0  
    else
       aL2=0d0
       aR2=0d0
    endif
    
    ! print *, 'DEBUG:', M_V, GA_V, DecayMode1, DecayMode2
    
    s  = 2d0 * scr(p(:,1),p(:,2))
    propG = one/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)
    
    sum = zero
    
    ! @TODO implement helicity selection ...
    ! @TODO implement only use non-zero helicity combinations
    
    ihel = 0

    do i1=1,2
       do i2 = 1,2
          do i3 = 1,2
             do i4 = 1,2
                ihel = ihel + 1
                if( ncalls_gg .ge. maxcalls .and. .not. goodhel_gg(ihel) ) then
                   cycle
                endif

                ! if( ihel > 8 ) then
                !    exit !!! amplitudes are mirrored in helicity combinations -- dschoute 05/09/2012
                ! endif

                if( SelectedHelComb .ge. 1 ) then
                   if( ihel .ne. SelectedHelComb ) then
                      cycle
                   endif
                endif
                
                pin(1,:) = p(:,1)
                pin(2,:) = p(:,2)
                sp(1,:) = pol_mless2(dcmplx(p(:,1)),-3+2*i1,'in')  ! gluon
                sp(2,:) = pol_mless2(dcmplx(p(:,2)),-3+2*i2,'in')  ! gluon
                
                !-------- -1 == left, 1 == right
                if( DecayMode1.ne.7 ) then 
                   pin(3,:) = p(:,3)+p(:,4)
                   pin(4,:) = p(:,5)+p(:,6)
                   sp(3,:) = pol_dk2mom(dcmplx(p(:,3)),dcmplx(p(:,4)),-3+2*i3)  !e-,e+
                   sp(4,:) = pol_dk2mom(dcmplx(p(:,5)),dcmplx(p(:,6)),-3+2*i4)  !mu-,mu+ / Q,Qbar
                elseif( DecayMode1.eq.7 ) then 
                   pin(3,:) = p(:,3)
                   pin(4,:) = p(:,5)
                   sp(3,:) = pol_mless2(dcmplx(pin(3,:)),-3+2*i3,'out')  ! photon
                   sp(4,:) = pol_mless2(dcmplx(pin(4,:)),-3+2*i4,'out')  ! photon
                endif
                
                call ggGZZampl(pin,sp,A(1))

                if( ncalls_gg .lt. maxcalls .and. abs(A(1)) .gt. tol .and. .not. goodhel_gg(ihel) ) then
                   goodhel_gg(ihel) = .true.
                   print *,'added helicity (gg->G)', ihel,-3+2*i1,-3+2*i2,-3+2*i3,-3+2*i4,abs(A(1))
                endif
                
                if (i3.eq.1) then
                   A(1) = aL1*A(1)
                elseif(i3.eq.2) then
                   A(1) = aR1*A(1)
                endif
                if (i4.eq.1) then
                   A(1) = aL2*A(1)
                elseif(i4.eq.2) then
                   A(1) = aR2*A(1)
                endif
                
                sum = sum + abs(propG*propZ1*propZ2*A(1))**2
                
             enddo
          enddo
       enddo
    enddo
    
    sum = sum*prefactor

    sum = sum
    
  end subroutine EvalAmp_gg_G_VV
  
  
  
  
  !----- a subroutine for q qbar -> G -> Z -> lept + Z --> 2 lepts
  !----- all outgoing convention and the following momentum assignment
  !-----  0 -> bq(p1) + q(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
  subroutine EvalAmp_qqb_G_VV(p,MY_IDUP,sum)
    implicit none
    real(dp), intent(out) ::  sum
    real(dp), intent(in) :: p(4,6)
    integer, intent(in) :: MY_IDUP(6:9)
    real(dp) :: s, pin(4,4)
    complex(dp) :: A(2), sp(4,4), propG, propZ1, propZ2
    integer :: i1,i2,i3,i4,ihel
    real(dp) :: aL1,aR1,aL2,aR2
    real(dp) :: gZ_sq
    real(dp) :: prefactor, Lambda_inv
    
    if( ncalls_qq .eq. 0 ) then
       do ihel=1,16
          goodhel_qq(ihel)=.false.
       enddo
    endif
    if( ncalls_qq .lt. maxcalls ) then
       ncalls_qq = ncalls_qq + 1
    endif

    !---- electroweak couplings
    !       aL = -one + two*sitW**2
    !       aR = aL+one
    gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2
    
    !---- the 1/Lambda coupling
    Lambda_inv = 1.0_dp/Lambda
    
    !---- full prefactor; 3 is  the color factor
    prefactor = 3d0*(Lambda_inv**2)**2*gZ_sq**2
    
    
    if( DecayMode1.le.3 ) then!  Z decay
       if( abs(MY_IDUP(6)).eq.abs(ElM_) .or. abs(MY_IDUP(6)).eq.abs(MuM_) .or. abs(MY_IDUP(6)).eq.abs(TaM_) ) then
          aL1=aL_lep
          aR1=aR_lep
       elseif( abs(MY_IDUP(6)).eq.abs(NuE_) .or. abs(MY_IDUP(6)).eq.abs(NuM_) .or. abs(MY_IDUP(6)).eq.abs(NuT_) ) then
          aL1=aL_neu
          aR1=aR_neu
       elseif( abs(MY_IDUP(6)).eq.abs(Up_) .or. abs(MY_IDUP(6)).eq.abs(Chm_) ) then
          aL1=aL_QUp
          aR1=aR_QUp
       elseif( abs(MY_IDUP(6)).eq.abs(Dn_) .or. abs(MY_IDUP(6)).eq.abs(Str_) .or. abs(MY_IDUP(6)).eq.abs(Bot_) ) then
          aL1=aL_QDn
          aR1=aR_QDn
       else
          aL1=0d0
          aR1=0d0
       endif
       prefactor = prefactor *(one/two*M_V*Ga_V)**2
       s = 2d0 * scr(p(:,3),p(:,4))
       propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
       s = 2d0 * scr(p(:,5),p(:,6))
       propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
    elseif( DecayMode1.ge.4 .and. DecayMode1.le.6 ) then !  W decay
       aL1 = bL
       aR1 = bR
       prefactor = prefactor *(one/two*M_V*Ga_V)**2
       s = 2d0 * scr(p(:,3),p(:,4))
       propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
       s = 2d0 * scr(p(:,5),p(:,6))
       propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
    elseif( DecayMode1.eq.7 ) then !  photon decay
       aL1=1d0
       aR1=1d0
       propZ1 = 1d0
       propZ2 = 1d0
       prefactor = prefactor/gZ_sq**2! cancel the overall z coupling
    else
       aL1=0d0
       aR1=0d0            
    endif
    
    if( DecayMode2.le.3 ) then!  Z decay
       if( abs(MY_IDUP(8)).eq.abs(ElM_) .or. abs(MY_IDUP(8)).eq.abs(MuM_) .or. abs(MY_IDUP(8)).eq.abs(TaM_) ) then
          aL2=aL_lep
          aR2=aR_lep
       elseif( abs(MY_IDUP(8)).eq.abs(NuE_) .or. abs(MY_IDUP(8)).eq.abs(NuM_) .or. abs(MY_IDUP(8)).eq.abs(NuT_) ) then
          aL2=aL_neu
          aR2=aR_neu
       elseif( abs(MY_IDUP(8)).eq.abs(Up_) .or. abs(MY_IDUP(8)).eq.abs(Chm_) ) then
          aL2=aL_QUp
          aR2=aR_QUp
       elseif( abs(MY_IDUP(8)).eq.abs(Dn_) .or. abs(MY_IDUP(8)).eq.abs(Str_) .or. abs(MY_IDUP(8)).eq.abs(Bot_) ) then
          aL2=aL_QDn
          aR2=aR_QDn
       else
          aL2=0d0
          aR2=0d0
       endif
    elseif( DecayMode2.ge.4 .and. DecayMode2.le.6 ) then !  W decay
       aL2 = bL
       aR2 = bR
    elseif( DecayMode2.eq.7 ) then !  photon decay
       aL2=1d0
       aR2=1d0  
    else
       aL2=0d0
       aR2=0d0  
    endif
    
    s  = two*scr(p(:,1),p(:,2))
    propG = one/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)
    
    ihel = 0

    sum = zero
    do i1=1,2
       do i3 = 1,2
          do i4 = 1,2
             ihel = ihel + 1
             if( ncalls_qq .ge. maxcalls .and. .not. goodhel_qq(ihel) ) then
                cycle
             endif
             
             ! if( ihel > 4 ) then
             !    exit !!! amplitudes are mirrored in helicity combinations -- dschoute 05/09/2012
             ! endif

             if( SelectedHelComb .ge. 1 ) then
                if( ihel .ne. SelectedHelComb ) then
                   cycle
                endif
             endif
             
             pin(1,:) = p(:,1)
             pin(2,:) = p(:,2)
             sp(1,:) = pol_dk2mom(dcmplx(p(:,2)),dcmplx(p(:,1)),-3+2*i1)  !qbq
             sp(2,:) = sp(1,:)  !-- the same, isn't really needed but for uniform bookeeping
             
             !-------- -1 == left, 1 == right
             if( DecayMode1.ne.7 ) then 
                pin(3,:) = p(:,3)+p(:,4)
                pin(4,:) = p(:,5)+p(:,6)
                sp(3,:) = pol_dk2mom(dcmplx(p(:,3)),dcmplx(p(:,4)),-3+2*i3)  !e-,e+
                sp(4,:) = pol_dk2mom(dcmplx(p(:,5)),dcmplx(p(:,6)),-3+2*i4)  !mu-,mu+ / Q,Qbar
             elseif( DecayMode1.eq.7 ) then 
                pin(3,:) = p(:,3)
                pin(4,:) = p(:,5)
                sp(3,:) = pol_mless2(dcmplx(pin(3,:)),-3+2*i3,'out')  ! photon
                sp(4,:) = pol_mless2(dcmplx(pin(4,:)),-3+2*i4,'out')  ! photon
             endif
             
             call qqGZZampl(pin,sp,A(1))

             if( ncalls_qq .lt. maxcalls .and. abs(A(1)) .gt. tol .and. .not. goodhel_qq(ihel) ) then
                goodhel_qq(ihel) = .true.
                print *,'added helicity (qq->G)', ihel,-3+2*i1,-(-3+2*i1),-3+2*i3,-3+2*i4,abs(A(1))
             endif
             
             if (i3.eq.1) then
                A(1) = aL1*A(1)
             elseif(i3.eq.2) then
                A(1) = aR1*A(1)
             endif
             if (i4.eq.1) then
                A(1) = aL2*A(1)
             elseif(i4.eq.2) then
                A(1) = aR2*A(1)
             endif
             
             sum = sum + abs(propG*propZ1*propZ2*A(1))**2
          enddo
       enddo
    enddo
    
    sum = sum*prefactor
    
  end subroutine EvalAmp_qqb_G_VV
  
  
  
  
  
  subroutine qqGZZampl(p,sp,res)
    implicit none
    real(dp), intent(in) :: p(4,4)
    complex(dp), intent(in) :: sp(4,4)
    complex(dp), intent(out) :: res
    complex(dp) :: e1_e2, e1_e3, e1_e4
    complex(dp) :: e2_e3, e2_e4
    complex(dp) :: e3_e4
    complex(dp) :: q_q,q3_q3,q4_q4
    complex(dp) :: q1_q2,q1_q3,q1_q4
    complex(dp) :: q2_q3,q2_q4
    complex(dp) :: q3_q4
    complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
    complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4
    complex(dp) :: e3_q4,e4_q3
    complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
    complex(dp) :: e1(4),e2(4),e3(4),e4(4)
    complex(dp) :: yyy1,yyy2,yyy3,yyy4,yyy5,yyy6, yyy7
    real(dp) :: q34,MG,MZ3,MZ4
    real(dp) :: rr
    
    q1 = dcmplx(p(1,:),0d0)
    q2 = dcmplx(p(2,:),0d0)
    q3 = dcmplx(p(3,:),0d0)
    q4 = dcmplx(p(4,:),0d0)
    
    
    e1 = sp(1,:)
    e2 = sp(2,:)
    e3 = sp(3,:)
    e4 = sp(4,:)
    
    q = -q1-q2
    
    q_q =sc(q,q)
    q3_q3 = sc(q3,q3)
    q4_q4 = sc(q4,q4)
    
    
    q1_q2 = sc(q1,q2)
    q1_q3 = sc(q1,q3)
    q1_q4 = sc(q1,q4)
    q2_q3 = sc(q2,q3)
    q2_q4 = sc(q2,q4)
    q3_q4 = sc(q3,q4)
    
    e1_e2 = sc(e1,e2)
    e1_e3 = sc(e1,e3)
    e1_e4 = sc(e1,e4)
    
    e2_e3 = sc(e2,e3)
    e2_e4 = sc(e2,e4)
    
    e3_e4 = sc(e3,e4)
    
    
    q1_e3 = sc(q1,e3)
    q1_e4 = sc(q1,e4)
    q2_e3 = sc(q2,e3)
    q2_e4 = sc(q2,e4)
    e1_q3 = sc(e1,q3)
    e1_q4 = sc(e1,q4)
    e2_q3 = sc(e2,q3)
    e2_q4 = sc(e2,q4)
    e3_q4 = sc(e3,q4)
    e4_q3 = sc(e4,q3)
    
    MZ3=dsqrt(cdabs(q3_q3))
    MZ4=dsqrt(cdabs(q4_q4))
    if( use_dynamic_MG ) then
       MG = dsqrt(cdabs(q_q))          
    else
       MG = M_Reso
    endif
    
    q34 = (MG**2-MZ3**2-MZ4**2)/2d0
    
    
    if( DecayMode1.le.3 .or. ((DecayMode1.ge.4) .and. (DecayMode1.le.6)) ) then! decay into Z's or W's
       if (generate_bis) then
          rr = q34/Lambda**2
          yyy1 = q34*(b1   + b2*rr*(one + two*M_V**2/q34+ M_V**4/q34**2)  + b5*M_V**2/q34)
          yyy2 = -b1/two + b3*rr*(1d0-M_V**2/q34) + two*b4*rr+b7*rr*M_V**2/q34
          yyy3 = (-b2/two - b3- two*b4)*rr/q34
          yyy4 = -b1 - b2*rr -(b2+b3+b6)*rr*M_V**2/q34
          yyy5 = two*b8*rr*MG**2/q34
          yyy6 = b9
          yyy7 = b10*rr*MG**2/q34
       else
          yyy1 = (q34)*c1
          yyy2 = c2
          yyy3 = c3
          yyy4 = c4
          yyy5 = c5
          yyy6 = c6
          yyy7 = c7
       endif
    elseif( DecayMode1.eq.7 ) then! decay into photons
       if (generate_bis) then
          rr = q34/Lambda**2
          yyy1 = q34*(b1   + b2*rr*(one + two*M_V**2/q34+ M_V**4/q34**2)  + b5*M_V**2/q34)
          yyy2 = -b1/two + b3*rr*(1d0-M_V**2/q34) + two*b4*rr+b7*rr*M_V**2/q34
          yyy3 = (-b2/two - b3- two*b4)*rr/q34
          yyy4 = -b1 - b2*rr -(b2+b3+b6)*rr*M_V**2/q34
          yyy5 = two*b8*rr*MG**2/q34
          yyy6 = 0d0
          yyy7 = b10*rr*MG**2/q34
       else
          yyy1 = q34*c1
          yyy2 = c2
          yyy3 = c3
          yyy4 = c4
          yyy5 = c5
          yyy6 = c6
          yyy7 = c7
       endif
    endif
    
    res = czero
    res =                                                               &
         q1_e3*e1_e4*yyy1 - q1_e3*e1_q3*e4_q3*yyy4 + q1_e4*e1_e3*yyy1 +       &
         q1_e4*e1_q3*e3_q4*yyy4 - q1_q3*e1_e3*e4_q3*yyy4 + q1_q3*e1_e4*       &
         e3_q4*yyy4 + 4*q1_q3*e1_q3*e3_e4*yyy2 + 4*q1_q3*e1_q3*e3_q4*         &
         e4_q3*yyy3 + 1./2.*e1_e3*e4_q3*yyy1 - 1./4.*e1_e3*e4_q3*MG**2*   &
         yyy4 + 1./2.*e1_e4*e3_q4*yyy1 + 1./4.*e1_e4*e3_q4*MG**2*yyy4 +   &
         e1_q3*e3_e4*MG**2*yyy2 + e1_q3*e3_q4*e4_q3*MG**2*yyy3
    
    res = res +                                                        &
         1./2.*et1(e3,e4,q1,q)*e1_q3*yyy6 - 1./2.*et1(e3,e4,q2,q)*e1_q3*  &
         yyy6 + 1./4.*et1(e3,e4,e1,q)*MG**2*yyy6 + et1(e3,e4,e1,q)*q1_q3* &
         yyy6 + 4*et1(e3,e4,q3,q4)*q1_q3*e1_q3*MG**(-2)*yyy5 + et1(e3,e4, &
         q3,q4)*e1_q3*yyy5
    
    res = res +                                                               &
         1./2.*et1(q1,e3,q,q3)*e1_q3*e4_q3*MG**(-2)*yyy7 - 1./2.*et1(q1,        &
         e3,q,q4)*e1_q3*e4_q3*MG**(-2)*yyy7 + 1./2.*et1(q1,e4,q,q3)*e1_q3      &
         *e3_q4*MG**(-2)*yyy7 - 1./2.*et1(q1,e4,q,q4)*e1_q3*e3_q4*             &
         MG**(-2)*yyy7 - 1./2.*et1(q2,e3,q,q3)*e1_q3*e4_q3*MG**(-2)*yyy7   &
         + 1./2.*et1(q2,e3,q,q4)*e1_q3*e4_q3*MG**(-2)*yyy7 - 1./2.*et1(       &
         q2,e4,q,q3)*e1_q3*e3_q4*MG**(-2)*yyy7 + 1./2.*et1(q2,e4,q,q4)*        &
         e1_q3*e3_q4*MG**(-2)*yyy7 + et1(e1,e3,q,q3)*q1_q3*e4_q3*MG**(-2)  &
         *yyy7 + 1./4.*et1(e1,e3,q,q3)*e4_q3*yyy7 - et1(e1,e3,q,q4)*q1_q3          &
         *e4_q3*MG**(-2)*yyy7 - 1./4.*et1(e1,e3,q,q4)*e4_q3*yyy7 + et1(e1      &
         ,e4,q,q3)*q1_q3*e3_q4*MG**(-2)*yyy7 + 1./4.*et1(e1,e4,q,q3)*          &
         e3_q4*yyy7 - et1(e1,e4,q,q4)*q1_q3*e3_q4*MG**(-2)*yyy7 - 1./4.*       &
         et1(e1,e4,q,q4)*e3_q4*yyy7
    
    
    
  end subroutine qqGZZampl
  
  
  
  
  
  subroutine ggGZZampl(p,sp,res)
    implicit none
    real(dp), intent(in) :: p(4,4)
    complex(dp), intent(in) :: sp(4,4)
    complex(dp), intent(out) :: res
    complex(dp) :: e1_e2, e1_e3, e1_e4
    complex(dp) :: e2_e3, e2_e4
    complex(dp) :: e3_e4
    complex(dp) :: q_q
    complex(dp) :: q1_q2,q1_q3,q1_q4,q3_q3,q4_q4
    complex(dp) :: q2_q3,q2_q4
    complex(dp) :: q3_q4
    complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
    complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4
    complex(dp) :: e3_q4,e4_q3
    complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
    complex(dp) :: e1(4),e2(4),e3(4),e4(4)
    complex(dp) :: xxx1,xxx2,xxx3,yyy1,yyy2,yyy3,yyy4,yyy5,yyy6
    complex(dp) :: yyy7
    real(dp) :: q34,MZ3,MZ4,MG
    logical :: new
    real(dp) :: rr_gam, rr
    
    
    new = .true.
    
    q1 = dcmplx(p(1,:),0d0)
    q2 = dcmplx(p(2,:),0d0)
    q3 = dcmplx(p(3,:),0d0)
    q4 = dcmplx(p(4,:),0d0)
    
    
    e1 = sp(1,:)
    e2 = sp(2,:)
    e3 = sp(3,:)
    e4 = sp(4,:)
    
    q = -q1-q2
    q_q = sc(q,q)
    q3_q3 = sc(q3,q3)
    q4_q4 = sc(q4,q4)
    
    
    q1_q2 = sc(q1,q2)
    q1_q3 = sc(q1,q3)
    q1_q4 = sc(q1,q4)
    q2_q3 = sc(q2,q3)
    q2_q4 = sc(q2,q4)
    q3_q4 = sc(q3,q4)
    
    e1_e2 = sc(e1,e2)
    e1_e3 = sc(e1,e3)
    e1_e4 = sc(e1,e4)
    e2_e3 = sc(e2,e3)
    e2_e4 = sc(e2,e4)
    e3_e4 = sc(e3,e4)
    
    q1_e3 = sc(q1,e3)
    q1_e4 = sc(q1,e4)
    q2_e3 = sc(q2,e3)
    q2_e4 = sc(q2,e4)
    e1_q3 = sc(e1,q3)
    e1_q4 = sc(e1,q4)
    e2_q3 = sc(e2,q3)
    e2_q4 = sc(e2,q4)
    e3_q4 = sc(e3,q4)
    e4_q3 = sc(e4,q3)
    
    MZ3=dsqrt(cdabs(q3_q3))
    MZ4=dsqrt(cdabs(q4_q4))
    if( use_dynamic_MG ) then
       MG = dsqrt(cdabs(q_q))
    else
       MG = M_Reso
    endif
    
    !---- define couplings
    q34 = (MG**2-MZ3**2-MZ4**2)/2d0!  = s = pV1.pV2    = (q_q-MZ3^2-MZ4^2)/2
    rr_gam = q_q/two/Lambda**2! kappa for IS
    xxx1 = (a1 + a2*rr_gam)
    xxx2 = -a1/two + (a3+two*a4)*rr_gam
    xxx3 = 4d0*a5*rr_gam
    
    if( DecayMode1.le.3 .or. ((DecayMode1.ge.4) .and. (DecayMode1.le.6)) ) then! decay into Z's or W's
       if (generate_bis) then
          rr = q34/Lambda**2! kappa for FS
          yyy1 = q34*(b1 + b2*rr*(one + two*M_V**2/q34+ M_V**4/q34**2)   + b5*M_V**2/q34)
          yyy2 = -b1/two + b3*rr*(one-M_V**2/q34) + two*b4*rr+b7*rr*M_V**2/q34
          yyy3 = (-b2/two - b3- two*b4)*rr/q34
          yyy4 = -b1 - b2*rr -(b2+b3+b6)*rr*M_V**2/q34
          yyy5 = two*b8*rr*MG**2/q34
          yyy6 = b9
          yyy7 = b10*rr*MG**2/q34
       else
          yyy1 = q34*c1
          yyy2 = c2
          yyy3 = c3
          yyy4 = c4
          yyy5 = c5
          yyy6 = c6
          yyy7 = c7
       endif
    elseif( DecayMode1.eq.7 ) then! decay into photons
       if (generate_bis) then
          rr = q34/Lambda**2! kappa for FS
          yyy1 = q34*(b1 + b2*rr*(one + two*M_V**2/q34+ M_V**4/q34**2)   + b5*M_V**2/q34)
          yyy2 = -b1/two + b3*rr*(one-M_V**2/q34) + two*b4*rr+b7*rr*M_V**2/q34
          yyy3 = (-b2/two - b3- two*b4)*rr/q34
          yyy4 = -b1 - b2*rr -(b2+b3+b6)*rr*M_V**2/q34
          yyy5 = two*b8*rr*MG**2/q34
          yyy6 = 0d0
          yyy7 = b10*rr*MG**2/q34
       else
          yyy1 = q34*c1
          yyy2 = c2
          yyy3 = c3
          yyy4 = c4
          yyy5 = c5
          yyy6 = c6
          yyy7 = c7
       endif
    endif
    
    
    
    res = czero
    
    
    if (new) then
       
       
       res = &
            + 8.D0*q1_e3*q1_e4*e1_e2*yyy1*xxx2 - 8.D0*q1_e3*q1_q3*e1_e2* &
            e4_q3*yyy4*xxx2 + 4.D0*q1_e3*e1_e2*e4_q3*yyy1*xxx2 + 2.D0* &
            q1_e3*e1_e2*e4_q3*MZ4**2*yyy4*xxx2 - 2.D0*q1_e3*e1_e2*e4_q3* &
            MZ3**2*yyy4*xxx2 - 2.D0*q1_e3*e1_e2*e4_q3*MG**2*yyy4*xxx2 + 8.d0 &
            *q1_e4*q1_q3*e1_e2*e3_q4*yyy4*xxx2 + 4.D0*q1_e4*e1_e2*e3_q4 &
            *yyy1*xxx2 - 2.D0*q1_e4*e1_e2*e3_q4*MZ4**2*yyy4*xxx2 + 2.D0* &
            q1_e4*e1_e2*e3_q4*MZ3**2*yyy4*xxx2 + 2.D0*q1_e4*e1_e2*e3_q4* &
            MG**2*yyy4*xxx2 - 8.D0*q1_q3*e1_e2*e3_e4*MZ4**2*yyy2*xxx2 + 8.d0 &
            *q1_q3*e1_e2*e3_e4*MZ3**2*yyy2*xxx2 + 8.D0*q1_q3*e1_e2* &
            e3_e4*MG**2*yyy2*xxx2 - 8.D0*q1_q3*e1_e2*e3_q4*e4_q3*MZ4**2* &
            yyy3*xxx2 + 8.D0*q1_q3*e1_e2*e3_q4*e4_q3*MZ3**2*yyy3*xxx2 + 8.d0 &
            *q1_q3*e1_e2*e3_q4*e4_q3*MG**2*yyy3*xxx2 + 16.D0*q1_q3**2* &
            e1_e2*e3_e4*yyy2*xxx2 + 16.D0*q1_q3**2*e1_e2*e3_q4*e4_q3*yyy3 &
            *xxx2 + 2.D0/3.D0*e1_e2*e3_e4*MZ4**4*yyy2*xxx2 + 1.D0/3.D0* &
            e1_e2*e3_e4*MZ4**4*yyy2*xxx1
       res = res - 4.D0/3.D0*e1_e2*e3_e4*MZ3**2*MZ4**2*yyy2*xxx2 - 2.D0/ &
            3.D0*e1_e2*e3_e4*MZ3**2*MZ4**2*yyy2*xxx1 + 2.D0/3.D0*e1_e2* &
            e3_e4*MZ3**4*yyy2*xxx2 + 1.D0/3.D0*e1_e2*e3_e4*MZ3**4*yyy2* &
            xxx1 + 2.D0/3.D0*e1_e2*e3_e4*MG**2*yyy1*xxx2 - 2.D0/3.D0* &
            e1_e2*e3_e4*MG**2*yyy1*xxx1 - 4.D0/3.D0*e1_e2*e3_e4*MG**2* &
            MZ4**2*yyy2*xxx2 - 2.D0/3.D0*e1_e2*e3_e4*MG**2*MZ4**2*yyy2* &
            xxx1 + 8.D0/3.D0*e1_e2*e3_e4*MG**2*MZ3**2*yyy2*xxx2 - 2.D0/3.D0 &
            *e1_e2*e3_e4*MG**2*MZ3**2*yyy2*xxx1 + 2.D0/3.D0*e1_e2*e3_e4* &
            MG**4*yyy2*xxx2 + 1.D0/3.D0*e1_e2*e3_e4*MG**4*yyy2*xxx1 + 4.D0 &
            /3.D0*e1_e2*e3_q4*e4_q3*yyy1*xxx2 + 2.D0/3.D0*e1_e2*e3_q4* &
            e4_q3*yyy1*xxx1 + 2.D0/3.D0*e1_e2*e3_q4*e4_q3*MZ4**4*yyy3* &
            xxx2 + 1.D0/3.D0*e1_e2*e3_q4*e4_q3*MZ4**4*yyy3*xxx1 - 4.D0/3.D0 &
            *e1_e2*e3_q4*e4_q3*MZ3**2*MZ4**2*yyy3*xxx2 - 2.D0/3.D0*e1_e2 &
            *e3_q4*e4_q3*MZ3**2*MZ4**2*yyy3*xxx1 + 2.D0/3.D0*e1_e2*e3_q4* &
            e4_q3*MZ3**4*yyy3*xxx2 + 1.D0/3.D0*e1_e2*e3_q4*e4_q3*MZ3**4* &
            yyy3*xxx1
       res = res + 2.D0/3.D0*e1_e2*e3_q4*e4_q3*MG**2*yyy4*xxx2 - 2.D0/3.D0 &
            *e1_e2*e3_q4*e4_q3*MG**2*yyy4*xxx1 - 4.D0/3.D0*e1_e2*e3_q4* &
            e4_q3*MG**2*MZ4**2*yyy3*xxx2 - 2.D0/3.D0*e1_e2*e3_q4*e4_q3* &
            MG**2*MZ4**2*yyy3*xxx1 + 8.D0/3.D0*e1_e2*e3_q4*e4_q3*MG**2* &
            MZ3**2*yyy3*xxx2 - 2.D0/3.D0*e1_e2*e3_q4*e4_q3*MG**2*MZ3**2* &
            yyy3*xxx1 + 2.D0/3.D0*e1_e2*e3_q4*e4_q3*MG**4*yyy3*xxx2 + 1.D0 &
            /3.D0*e1_e2*e3_q4*e4_q3*MG**4*yyy3*xxx1 + e1_e3*e2_e4*MG**2* &
            yyy1*xxx1 - e1_e3*e2_q3*e4_q3*MG**2*yyy4*xxx1 + e1_e4*e2_e3* &
            MG**2*yyy1*xxx1 + e1_e4*e2_q3*e3_q4*MG**2*yyy4*xxx1 - e1_q3* &
            e2_e3*e4_q3*MG**2*yyy4*xxx1 + e1_q3*e2_e4*e3_q4*MG**2*yyy4* &
            xxx1 + 4.D0*e1_q3*e2_q3*e3_e4*MG**2*yyy2*xxx1 + 4.D0*e1_q3* &
            e2_q3*e3_q4*e4_q3*MG**2*yyy3*xxx1 + 8.D0*et1(q1,q2,e1,e2)* &
            q1_e3*q1_e4*MG**(-2)*yyy1*xxx3 - 8.D0*et1(q1,q2,e1,e2)*q1_e3* &
            q1_q3*e4_q3*MG**(-2)*yyy4*xxx3 + 4.D0*et1(q1,q2,e1,e2)*q1_e3* &
            e4_q3*MG**(-2)*yyy1*xxx3 + 2.D0*et1(q1,q2,e1,e2)*q1_e3*e4_q3* &
            MG**(-2)*MZ4**2*yyy4*xxx3
       res = res - 2.D0*et1(q1,q2,e1,e2)*q1_e3*e4_q3*MG**(-2)*MZ3**2* &
            yyy4*xxx3 - 2.D0*et1(q1,q2,e1,e2)*q1_e3*e4_q3*yyy4*xxx3 + 8.D0* &
            et1(q1,q2,e1,e2)*q1_e4*q1_q3*e3_q4*MG**(-2)*yyy4*xxx3 + 4.D0* &
            et1(q1,q2,e1,e2)*q1_e4*e3_q4*MG**(-2)*yyy1*xxx3 - 2.D0*et1(q1 &
            ,q2,e1,e2)*q1_e4*e3_q4*MG**(-2)*MZ4**2*yyy4*xxx3 + 2.D0*et1( &
            q1,q2,e1,e2)*q1_e4*e3_q4*MG**(-2)*MZ3**2*yyy4*xxx3 + 2.D0* &
            et1(q1,q2,e1,e2)*q1_e4*e3_q4*yyy4*xxx3 - 8.D0*et1(q1,q2,e1,e2 &
            )*q1_q3*e3_e4*MG**(-2)*MZ4**2*yyy2*xxx3 + 8.D0*et1(q1,q2,e1, &
            e2)*q1_q3*e3_e4*MG**(-2)*MZ3**2*yyy2*xxx3 + 8.D0*et1(q1,q2,e1 &
            ,e2)*q1_q3*e3_e4*yyy2*xxx3 - 8.D0*et1(q1,q2,e1,e2)*q1_q3* &
            e3_q4*e4_q3*MG**(-2)*MZ4**2*yyy3*xxx3 + 8.D0*et1(q1,q2,e1,e2) &
            *q1_q3*e3_q4*e4_q3*MG**(-2)*MZ3**2*yyy3*xxx3 + 8.D0*et1(q1,q2 &
            ,e1,e2)*q1_q3*e3_q4*e4_q3*yyy3*xxx3 + 16.D0*et1(q1,q2,e1,e2)* &
            q1_q3**2*e3_e4*MG**(-2)*yyy2*xxx3 + 16.D0*et1(q1,q2,e1,e2)* &
            q1_q3**2*e3_q4*e4_q3*MG**(-2)*yyy3*xxx3 + 2.D0/3.D0*et1(q1,q2 &
            ,e1,e2)*e3_e4*MG**(-2)*MZ4**4*yyy2*xxx3
       res = res - 4.D0/3.D0*et1(q1,q2,e1,e2)*e3_e4*MG**(-2)*MZ3**2* &
            MZ4**2*yyy2*xxx3 + 2.D0/3.D0*et1(q1,q2,e1,e2)*e3_e4*MG**(-2)* &
            MZ3**4*yyy2*xxx3 + 2.D0/3.D0*et1(q1,q2,e1,e2)*e3_e4*yyy1*xxx3 &
            - 4.D0/3.D0*et1(q1,q2,e1,e2)*e3_e4*MZ4**2*yyy2*xxx3 + 8.D0/3.D0 &
            *et1(q1,q2,e1,e2)*e3_e4*MZ3**2*yyy2*xxx3 + 2.D0/3.D0*et1(q1 &
            ,q2,e1,e2)*e3_e4*MG**2*yyy2*xxx3 + 4.D0/3.D0*et1(q1,q2,e1,e2) &
            *e3_q4*e4_q3*MG**(-2)*yyy1*xxx3 + 2.D0/3.D0*et1(q1,q2,e1,e2)* &
            e3_q4*e4_q3*MG**(-2)*MZ4**4*yyy3*xxx3 - 4.D0/3.D0*et1(q1,q2, &
            e1,e2)*e3_q4*e4_q3*MG**(-2)*MZ3**2*MZ4**2*yyy3*xxx3 + 2.D0/3.D0 &
            *et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**(-2)*MZ3**4*yyy3*xxx3 + 2.D0 &
            /3.D0*et1(q1,q2,e1,e2)*e3_q4*e4_q3*yyy4*xxx3 - 4.D0/3.D0* &
            et1(q1,q2,e1,e2)*e3_q4*e4_q3*MZ4**2*yyy3*xxx3 + 8.D0/3.D0* &
            et1(q1,q2,e1,e2)*e3_q4*e4_q3*MZ3**2*yyy3*xxx3 + 2.D0/3.D0* &
            et1(q1,q2,e1,e2)*e3_q4*e4_q3*MG**2*yyy3*xxx3 - et1(q1,q2,e1, &
            e2)*et1(q1,q,e3,e4)*MG**(-2)*MZ4**2*yyy6*xxx3 + et1(q1,q2,e1, &
            e2)*et1(q1,q,e3,e4)*MG**(-2)*MZ3**2*yyy6*xxx3
       res = res + et1(q1,q2,e1,e2)*et1(q1,q,e3,e4)*yyy6*xxx3 + 4.D0* &
            et1(q1,q2,e1,e2)*et1(q1,q,e3,e4)*q1_q3*MG**(-2)*yyy6*xxx3 - 4.0d0 &
            *et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)*q1_q3*e4_q3*MG**(-4)*yyy7 &
            *xxx3 + et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)*e4_q3*MG**(-4)* &
            MZ4**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)*e4_q3* &
            MG**(-4)*MZ3**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e3,q3)* &
            e4_q3*MG**(-2)*yyy7*xxx3 + 4.D0*et1(q1,q2,e1,e2)*et1(q1,q,e3, &
            q4)*q1_q3*e4_q3*MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1, &
            q,e3,q4)*e4_q3*MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)* &
            et1(q1,q,e3,q4)*e4_q3*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2, &
            e1,e2)*et1(q1,q,e3,q4)*e4_q3*MG**(-2)*yyy7*xxx3 - 4.D0*et1(q1 &
            ,q2,e1,e2)*et1(q1,q,e4,q3)*q1_q3*e3_q4*MG**(-4)*yyy7*xxx3 + &
            et1(q1,q2,e1,e2)*et1(q1,q,e4,q3)*e3_q4*MG**(-4)*MZ4**2*yyy7* &
            xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e4,q3)*e3_q4*MG**(-4)*MZ3**2 &
            *yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e4,q3)*e3_q4*MG**(-2)* &
            yyy7*xxx3
       res = res + 4.D0*et1(q1,q2,e1,e2)*et1(q1,q,e4,q4)*q1_q3*e3_q4* &
            MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q1,q,e4,q4)*e3_q4* &
            MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q1,q,e4,q4)* &
            e3_q4*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q1,q, &
            e4,q4)*e3_q4*MG**(-2)*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1(q2,q, &
            e3,e4)*MG**(-2)*MZ4**2*yyy6*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q, &
            e3,e4)*MG**(-2)*MZ3**2*yyy6*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q, &
            e3,e4)*yyy6*xxx3 - 4.D0*et1(q1,q2,e1,e2)*et1(q2,q,e3,e4)* &
            q1_q3*MG**(-2)*yyy6*xxx3 + 4.D0*et1(q1,q2,e1,e2)*et1(q2,q,e3, &
            q3)*q1_q3*e4_q3*MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2, &
            q,e3,q3)*e4_q3*MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)* &
            et1(q2,q,e3,q3)*e4_q3*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2, &
            e1,e2)*et1(q2,q,e3,q3)*e4_q3*MG**(-2)*yyy7*xxx3 - 4.D0*et1(q1 &
            ,q2,e1,e2)*et1(q2,q,e3,q4)*q1_q3*e4_q3*MG**(-4)*yyy7*xxx3 + &
            et1(q1,q2,e1,e2)*et1(q2,q,e3,q4)*e4_q3*MG**(-4)*MZ4**2*yyy7* &
            xxx3
       res = res - et1(q1,q2,e1,e2)*et1(q2,q,e3,q4)*e4_q3*MG**(-4)* &
            MZ3**2*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e3,q4)*e4_q3* &
            MG**(-2)*yyy7*xxx3 + 4.D0*et1(q1,q2,e1,e2)*et1(q2,q,e4,q3)* &
            q1_q3*e3_q4*MG**(-4)*yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e4 &
            ,q3)*e3_q4*MG**(-4)*MZ4**2*yyy7*xxx3 + et1(q1,q2,e1,e2)*et1( &
            q2,q,e4,q3)*e3_q4*MG**(-4)*MZ3**2*yyy7*xxx3 + et1(q1,q2,e1,e2 &
            )*et1(q2,q,e4,q3)*e3_q4*MG**(-2)*yyy7*xxx3 - 4.D0*et1(q1,q2, &
            e1,e2)*et1(q2,q,e4,q4)*q1_q3*e3_q4*MG**(-4)*yyy7*xxx3 + et1( &
            q1,q2,e1,e2)*et1(q2,q,e4,q4)*e3_q4*MG**(-4)*MZ4**2*yyy7*xxx3 &
            - et1(q1,q2,e1,e2)*et1(q2,q,e4,q4)*e3_q4*MG**(-4)*MZ3**2* &
            yyy7*xxx3 - et1(q1,q2,e1,e2)*et1(q2,q,e4,q4)*e3_q4*MG**(-2)* &
            yyy7*xxx3 - 1.D0/3.D0*et1(q1,q2,e1,e2)*et1(q,e3,e4,q3)*yyy6* &
            xxx3 + 1.D0/3.D0*et1(q1,q2,e1,e2)*et1(q,e3,e4,q4)*yyy6*xxx3 &
            + 2.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*MG**(-4)* &
            MZ4**4*yyy5*xxx3 - 4.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4 &
            )*MG**(-4)*MZ3**2*MZ4**2*yyy5*xxx3
       res = res + 2.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*MG**(-4)* &
            MZ3**4*yyy5*xxx3 - 4.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)* &
            MG**(-2)*MZ4**2*yyy5*xxx3 + 8.D0/3.D0*et1(q1,q2,e1,e2)*et1(e3 &
            ,e4,q3,q4)*MG**(-2)*MZ3**2*yyy5*xxx3 + 2.D0/3.D0*et1(q1,q2,e1 &
            ,e2)*et1(e3,e4,q3,q4)*yyy5*xxx3 - 8.D0*et1(q1,q2,e1,e2)*et1( &
            e3,e4,q3,q4)*q1_q3*MG**(-4)*MZ4**2*yyy5*xxx3 + 8.D0*et1(q1,q2 &
            ,e1,e2)*et1(e3,e4,q3,q4)*q1_q3*MG**(-4)*MZ3**2*yyy5*xxx3 + 8.D0 &
            *et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*q1_q3*MG**(-2)*yyy5*xxx3 &
            + 16.D0*et1(q1,q2,e1,e2)*et1(e3,e4,q3,q4)*q1_q3**2*MG**(-4)* &
            yyy5*xxx3 + 4.D0*et1(q1,q,e3,e4)*q1_q3*e1_e2*yyy6*xxx2 - et1( &
            q1,q,e3,e4)*e1_e2*MZ4**2*yyy6*xxx2 + et1(q1,q,e3,e4)*e1_e2* &
            MZ3**2*yyy6*xxx2 + et1(q1,q,e3,e4)*e1_e2*MG**2*yyy6*xxx2 - 4.D0 &
            *et1(q1,q,e3,q3)*q1_q3*e1_e2*e4_q3*MG**(-2)*yyy7*xxx2 + et1( &
            q1,q,e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q1,q, &
            e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q1,q,e3,q3 &
            )*e1_e2*e4_q3*yyy7*xxx2
       res = res + 4.D0*et1(q1,q,e3,q4)*q1_q3*e1_e2*e4_q3*MG**(-2)*yyy7* &
            xxx2 - et1(q1,q,e3,q4)*e1_e2*e4_q3*MG**(-2)*MZ4**2*yyy7*xxx2 + &
            et1(q1,q,e3,q4)*e1_e2*e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 + et1( &
            q1,q,e3,q4)*e1_e2*e4_q3*yyy7*xxx2 - 4.D0*et1(q1,q,e4,q3)* &
            q1_q3*e1_e2*e3_q4*MG**(-2)*yyy7*xxx2 + et1(q1,q,e4,q3)*e1_e2* &
            e3_q4*MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q1,q,e4,q3)*e1_e2*e3_q4 &
            *MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q1,q,e4,q3)*e1_e2*e3_q4*yyy7 &
            *xxx2 + 4.D0*et1(q1,q,e4,q4)*q1_q3*e1_e2*e3_q4*MG**(-2)*yyy7* &
            xxx2 - et1(q1,q,e4,q4)*e1_e2*e3_q4*MG**(-2)*MZ4**2*yyy7*xxx2 &
            + et1(q1,q,e4,q4)*e1_e2*e3_q4*MG**(-2)*MZ3**2*yyy7*xxx2 + &
            et1(q1,q,e4,q4)*e1_e2*e3_q4*yyy7*xxx2 - 4.D0*et1(q2,q,e3,e4)* &
            q1_q3*e1_e2*yyy6*xxx2 + et1(q2,q,e3,e4)*e1_e2*MZ4**2*yyy6* &
            xxx2 - et1(q2,q,e3,e4)*e1_e2*MZ3**2*yyy6*xxx2 - et1(q2,q,e3, &
            e4)*e1_e2*MG**2*yyy6*xxx2 + 4.D0*et1(q2,q,e3,q3)*q1_q3*e1_e2* &
            e4_q3*MG**(-2)*yyy7*xxx2 - et1(q2,q,e3,q3)*e1_e2*e4_q3* &
            MG**(-2)*MZ4**2*yyy7*xxx2
       res = res + et1(q2,q,e3,q3)*e1_e2*e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 &
            + et1(q2,q,e3,q3)*e1_e2*e4_q3*yyy7*xxx2 - 4.D0*et1(q2,q,e3, &
            q4)*q1_q3*e1_e2*e4_q3*MG**(-2)*yyy7*xxx2 + et1(q2,q,e3,q4)* &
            e1_e2*e4_q3*MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q2,q,e3,q4)*e1_e2 &
            *e4_q3*MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q2,q,e3,q4)*e1_e2* &
            e4_q3*yyy7*xxx2 + 4.D0*et1(q2,q,e4,q3)*q1_q3*e1_e2*e3_q4* &
            MG**(-2)*yyy7*xxx2 - et1(q2,q,e4,q3)*e1_e2*e3_q4*MG**(-2)* &
            MZ4**2*yyy7*xxx2 + et1(q2,q,e4,q3)*e1_e2*e3_q4*MG**(-2)* &
            MZ3**2*yyy7*xxx2 + et1(q2,q,e4,q3)*e1_e2*e3_q4*yyy7*xxx2 - 4.D0 &
            *et1(q2,q,e4,q4)*q1_q3*e1_e2*e3_q4*MG**(-2)*yyy7*xxx2 + et1( &
            q2,q,e4,q4)*e1_e2*e3_q4*MG**(-2)*MZ4**2*yyy7*xxx2 - et1(q2,q, &
            e4,q4)*e1_e2*e3_q4*MG**(-2)*MZ3**2*yyy7*xxx2 - et1(q2,q,e4,q4 &
            )*e1_e2*e3_q4*yyy7*xxx2 - et1(q,e1,e3,e4)*e2_q3*MG**2*yyy6* &
            xxx1 + et1(q,e1,e3,q3)*e2_q3*e4_q3*yyy7*xxx1 - et1(q,e1,e3,q4 &
            )*e2_q3*e4_q3*yyy7*xxx1 + et1(q,e1,e4,q3)*e2_q3*e3_q4*yyy7* &
            xxx1
       res = res - et1(q,e1,e4,q4)*e2_q3*e3_q4*yyy7*xxx1 - et1(q,e2,e3, &
            e4)*e1_q3*MG**2*yyy6*xxx1 + et1(q,e2,e3,q3)*e1_q3*e4_q3*yyy7* &
            xxx1 - et1(q,e2,e3,q4)*e1_q3*e4_q3*yyy7*xxx1 + et1(q,e2,e4,q3 &
            )*e1_q3*e3_q4*yyy7*xxx1 - et1(q,e2,e4,q4)*e1_q3*e3_q4*yyy7* &
            xxx1 - 1.D0/3.D0*et1(q,e3,e4,q3)*e1_e2*MG**2*yyy6*xxx2 + 1.D0/ &
            3.D0*et1(q,e3,e4,q3)*e1_e2*MG**2*yyy6*xxx1 + 1.D0/3.D0*et1(q, &
            e3,e4,q4)*e1_e2*MG**2*yyy6*xxx2 - 1.D0/3.D0*et1(q,e3,e4,q4)* &
            e1_e2*MG**2*yyy6*xxx1 - 8.D0*et1(e3,e4,q3,q4)*q1_q3*e1_e2* &
            MG**(-2)*MZ4**2*yyy5*xxx2 + 8.D0*et1(e3,e4,q3,q4)*q1_q3*e1_e2 &
            *MG**(-2)*MZ3**2*yyy5*xxx2 + 8.D0*et1(e3,e4,q3,q4)*q1_q3* &
            e1_e2*yyy5*xxx2 + 16.D0*et1(e3,e4,q3,q4)*q1_q3**2*e1_e2* &
            MG**(-2)*yyy5*xxx2 + 2.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2* &
            MG**(-2)*MZ4**4*yyy5*xxx2 + 1.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2* &
            MG**(-2)*MZ4**4*yyy5*xxx1 - 4.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2* &
            MG**(-2)*MZ3**2*MZ4**2*yyy5*xxx2 - 2.D0/3.D0*et1(e3,e4,q3,q4) &
            *e1_e2*MG**(-2)*MZ3**2*MZ4**2*yyy5*xxx1
       res = res + 2.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2*MG**(-2)*MZ3**4*yyy5 &
            *xxx2 + 1.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2*MG**(-2)*MZ3**4*yyy5* &
            xxx1 - 4.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2*MZ4**2*yyy5*xxx2 - 2.D0 &
            /3.D0*et1(e3,e4,q3,q4)*e1_e2*MZ4**2*yyy5*xxx1 + 8.D0/3.D0* &
            et1(e3,e4,q3,q4)*e1_e2*MZ3**2*yyy5*xxx2 - 2.D0/3.D0*et1(e3,e4 &
            ,q3,q4)*e1_e2*MZ3**2*yyy5*xxx1 + 2.D0/3.D0*et1(e3,e4,q3,q4)* &
            e1_e2*MG**2*yyy5*xxx2 + 1.D0/3.D0*et1(e3,e4,q3,q4)*e1_e2* &
            MG**2*yyy5*xxx1 + 4.D0*et1(e3,e4,q3,q4)*e1_q3*e2_q3*yyy5*xxx1
       
       
       
    else
       
       
       
       res = q1_e3*q2_e4*e1_e2 * ( M_Reso**2 )
       res = res + q1_e3*q2_q4*e1_e2*e4_q3 * (  - 2.0 )
       res = res + q1_e4*q2_e3*e1_e2 * ( M_Reso**2 )
       res = res + q1_e4*q2_q3*e1_e2*e3_q4 * (  - 2.0 )
       res = res + q1_q3*q2_e4*e1_e2*e3_q4 * (  - 2.0 )
       res = res + q1_q3*q2_q4*e1_e2*e3_e4 * ( 2.0 )
       res = res + q1_q3*e1_e2*e3_e4 * ( 1./2.*M_Reso**2 )
       res = res + q1_q3*e1_e2*e3_q4*e4_q3 * (  - 1. )
       res = res + q1_q4*q2_e3*e1_e2*e4_q3 * (  - 2. )
       res = res + q1_q4*q2_q3*e1_e2*e3_e4 * ( 2 )
       res = res + q1_q4*e1_e2*e3_e4 * ( 1./2.*M_Reso**2 )
       res = res + q1_q4*e1_e2*e3_q4*e4_q3 * (  - 1. )
       res = res + q2_q3*e1_e2*e3_e4 * ( 1./2.*M_Reso**2 )
       res = res + q2_q3*e1_e2*e3_q4*e4_q3 * (  - 1 )
       res = res + q2_q4*e1_e2*e3_e4 * ( 1./2.*M_Reso**2 )
       res = res + q2_q4*e1_e2*e3_q4*e4_q3 * (  - 1. )
       res = res + e1_e2*e3_e4 * ( M_Reso**2*M_V**2 - 1./2.*M_Reso**4 )
       res = res + e1_e2*e3_q4*e4_q3 * ( M_Reso**2 )
       res = res + e1_e3*e2_e4 * ( 1./2.*M_Reso**4 )
       res = res + e1_e3*e2_q4*e4_q3 * (  - M_Reso**2 )
       res = res + e1_e4*e2_e3 * ( 1./2.*M_Reso**4 )
       res = res + e1_e4*e2_q3*e3_q4 * (  - M_Reso**2 )
       res = res + e1_q3*e2_e4*e3_q4 * (  - M_Reso**2 )
       res = res + e1_q3*e2_q4*e3_e4 * ( M_Reso**2 )
       res = res + e1_q4*e2_e3*e4_q3 * (  - M_Reso**2 )
       res = res + e1_q4*e2_q3*e3_e4 * ( M_Reso**2 )
       
    endif
    
    
  end subroutine ggGZZampl
  
  
  
  
  
  
  
  
  
  
  subroutine qqGZZampl_old(p,sp,res)
    implicit none
    real(dp), intent(in) :: p(4,4)
    complex(dp), intent(in) :: sp(4,4)
    complex(dp), intent(out) :: res
    complex(dp) :: e1_e2, e1_e3, e1_e4
    complex(dp) :: e2_e3, e2_e4
    complex(dp) :: e3_e4
    complex(dp) :: q_q
    complex(dp) :: q1_q2,q1_q3,q1_q4
    complex(dp) :: q2_q3,q2_q4
    complex(dp) :: q3_q4
    complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
    complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4
    complex(dp) :: e3_q4,e4_q3
    complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
    complex(dp) :: e1(4),e2(4),e3(4),e4(4)
    
    q1 = dcmplx(p(1,:),0d0)
    q2 = dcmplx(p(2,:),0d0)
    q3 = dcmplx(p(3,:),0d0)
    q4 = dcmplx(p(4,:),0d0)
    
    
    e1 = sp(1,:)
    e2 = sp(2,:)
    e3 = sp(3,:)
    e4 = sp(4,:)
    
    q = -q1-q2
    
    q_q =sc(q,q)
    
    
    q1_q2 = sc(q1,q2)
    q1_q3 = sc(q1,q3)
    q1_q4 = sc(q1,q4)
    q2_q3 = sc(q2,q3)
    q2_q4 = sc(q2,q4)
    q3_q4 = sc(q3,q4)
    
    e1_e2 = sc(e1,e2)
    e1_e3 = sc(e1,e3)
    e1_e4 = sc(e1,e4)
    
    e2_e3 = sc(e2,e3)
    e2_e4 = sc(e2,e4)
    
    e3_e4 = sc(e3,e4)
    
    
    q1_e3 = sc(q1,e3)
    q1_e4 = sc(q1,e4)
    q2_e3 = sc(q2,e3)
    q2_e4 = sc(q2,e4)
    e1_q3 = sc(e1,q3)
    e1_q4 = sc(e1,q4)
    e2_q3 = sc(e2,q3)
    e2_q4 = sc(e2,q4)
    e3_q4 = sc(e3,q4)
    e4_q3 = sc(e4,q3)
    
    res = czero
    
    
    res =                                                            &
         + 1./4.*q1_e3*q_q*e1_e4 - 1./2.*q1_e3*e1_q4*e4_q3 + 1./4.*q1_e4  &
         *q_q*e1_e3 - 1./2.*q1_e4*e1_q3*e3_q4 - 1./2.*q1_q3*e1_e4*      &
         e3_q4 + 1./2.*q1_q3*e1_q4*e3_e4 - 1./2.*q1_q4*e1_e3*e4_q3 + 1. &
         /2.*q1_q4*e1_q3*e3_e4 - 1./4.*q2_e3*q_q*e1_e4 + 1./2.*q2_e3*   &
         e1_q4*e4_q3 - 1./4.*q2_e4*q_q*e1_e3 + 1./2.*q2_e4*e1_q3*e3_q4  &
         + 1./2.*q2_q3*e1_e4*e3_q4 - 1./2.*q2_q3*e1_q4*e3_e4 + 1./2.*  &
         q2_q4*e1_e3*e4_q3 - 1./2.*q2_q4*e1_q3*e3_e4
    
    
  end subroutine qqGZZampl_old
  
  
  
  
  
  
  !----- a subroutine for q qbar -> G -> Z -> lept + Z --> 2 lepts
  !----- all outgoing convention and the following momentum assignment
  !-----  0 -> bq(p1) + q(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
  subroutine EvalAmp_qqb_G_VV_old(p,sum)
    implicit none
    real(dp), intent(out) ::  sum
    real(dp), intent(in) :: p(4,6)
    real(dp) :: s, pin(4,4)
    complex(dp) :: A(2), sp(4,4), propG, propZ1, propZ2
    integer :: i1,i2,i3,i4
    real(dp) :: aL,aR
    real(dp) :: gZ_sq
    real(dp) :: prefactor, Lambda_inv
    
    !---- electroweak couplings
    
    aL = -one + two*sitW**2
    aR = aL+one
    
    gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2
    
    !---- the 1/Lambda coupling
    
    
    Lambda_inv = 1.0_dp/Lambda
    
    !---- full prefactor; 3 is  the color factor
    
    prefactor = 3d0*(Lambda_inv**2)**2*(one/two*M_V*Ga_V)**2*gZ_sq**2
    
    sum = zero
    
    s  = two*scr(p(:,1),p(:,2))
    propG = one/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)
    
    s = two*scr(p(:,3),p(:,4))
    propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
    
    s = two*scr(p(:,5),p(:,6))
    propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
    
    do i1=1,2
       do i3 = 1,2
          do i4 = 1,2
             
             
             pin(1,:) = p(:,1)
             pin(2,:) = p(:,2)
             
             sp(1,:) = pol_dk2mom(dcmplx(p(:,2)),dcmplx(p(:,1)),-3+2*i1)  !qbq
             sp(2,:) = sp(1,:)  !-- the same, isn't really needed but for uniform
             ! bookeeping
             
             pin(3,:) = p(:,3) + p(:,4)
             pin(4,:) = p(:,5) + p(:,6)
             
             !-------- -1 == left, 1 == right
             
             sp(3,:) = pol_dk2mom(dcmplx(p(:,3)),dcmplx(p(:,4)),-3+2*i3)  !e-,e+
             sp(4,:) = pol_dk2mom(dcmplx(p(:,5)),dcmplx(p(:,6)),-3+2*i4)  !mu-.mu+
             
             call qqGZZampl_old(pin,sp,A(1))
             
             if (i3.eq.1) then
                A(1) = aL*A(1)
             elseif(i3.eq.2) then
                A(1) = aR*A(1)
             endif
             
             if (i4.eq.1) then
                A(1) = aL*A(1)
             elseif(i4.eq.2) then
                A(1) = aR*A(1)
             endif
             
             
             
             sum = sum + abs(propG*propZ1*propZ2*A(1))**2
             
             
          enddo
       enddo
    enddo
    
    
    sum = sum*prefactor
    
  end subroutine EvalAmp_qqb_G_VV_old
  
  
  
  
  
  !----- a subroutinefor gg -> G -> Z -> lept + Z --> 2 lepts
  !----- all outgoing convention and the following momentum assignment
  !-----  0 -> g(p1) + g(p2) + e-(p3) + e+(p4) +mu-(p5) +mu+(p6)
  subroutine EvalAmp_gg_G_VV_old(p,sum)
    use ModMisc
    implicit none
    real(dp), intent(out) ::  sum
    real(dp), intent(in) :: p(4,6)
    real(dp) :: s, pin(4,4)
    complex(dp) :: A(2), sp(4,4), propG, propZ1, propZ2
    integer :: i1,i2,i3,i4
    real(dp) :: aL,aR
    real(dp) :: gZ_sq
    real(dp) :: prefactor, Lambda_inv
    
    
    !---- electroweak couplings
    aL = -one + two*sitW**2
    aR = aL+one
    
    
    gZ_sq = 4.0_dp*pi*alpha_QED/4.0_dp/(one-sitW**2)/sitW**2
    
    !---- the 1/Lambda coupling
    Lambda_inv = 1.0d0/Lambda
    
    !---- full prefactor; 8 is  the color factor
    prefactor = 8d0*(Lambda_inv**2)**2*(one/two*M_V*Ga_V)**2*gZ_sq**2
    
    
    sum = zero
    
    s  = 2d0 * scr(p(:,1),p(:,2))
    propG = one/dcmplx(s - M_Reso**2,M_Reso*Ga_Reso)
    
    s = 2d0 * scr(p(:,3),p(:,4))
    propZ1 = s/dcmplx(s - M_V**2,M_V*Ga_V)
    
    s = 2d0 * scr(p(:,5),p(:,6))
    propZ2 = s/dcmplx(s - M_V**2,M_V*Ga_V)
    
    do i1=1,2
       do i2 = 1,2
          do i3 = 1,2
             do i4 = 1,2
                
                
                pin(1,:) = p(:,1)
                pin(2,:) = p(:,2)
                
                sp(1,:) = pol_mless2(dcmplx(p(:,1)),-3+2*i1,'in')  ! gluon
                sp(2,:) = pol_mless2(dcmplx(p(:,2)),-3+2*i2,'in')  ! gluon
                
                pin(3,:) = p(:,3) + p(:,4)
                pin(4,:) = p(:,5)+p(:,6)
                
                !-------- -1 == left, 1 == right
                sp(3,:) = pol_dk2mom(dcmplx(p(:,3)),dcmplx(p(:,4)),-3+2*i3)  !e-,e+
                sp(4,:) = pol_dk2mom(dcmplx(p(:,5)),dcmplx(p(:,6)),-3+2*i4)  !mu-.mu+
                
                
                call ggGZZampl_old(pin,sp,A(1))
                
                if (i3.eq.1) then
                   A(1) = aL*A(1)
                elseif(i3.eq.2) then
                   A(1) = aR*A(1)
                endif
                
                if (i4.eq.1) then
                   A(1) = aL*A(1)
                elseif(i4.eq.2) then
                   A(1) = aR*A(1)
                endif
                
                
                sum = sum + abs(propG*propZ1*propZ2*A(1))**2
                
                
             enddo
          enddo
       enddo
    enddo
    
    sum = sum*prefactor
    
  end subroutine EvalAmp_gg_G_VV_old
  
  
  
  
  
  subroutine ggGZZampl_old(p,sp,res)
    implicit none
    real(dp), intent(in) :: p(4,4)
    complex(dp), intent(in) :: sp(4,4)
    complex(dp), intent(out) :: res
    complex(dp) :: e1_e2, e1_e3, e1_e4
    complex(dp) :: e2_e3, e2_e4
    complex(dp) :: e3_e4
    complex(dp) :: q_q
    complex(dp) :: q1_q2,q1_q3,q1_q4
    complex(dp) :: q2_q3,q2_q4
    complex(dp) :: q3_q4
    complex(dp) :: q1_e3,q1_e4,q2_e3,q2_e4
    complex(dp) :: e1_q3,e1_q4,e2_q3,e2_q4
    complex(dp) :: e3_q4,e4_q3
    complex(dp) :: q1(4),q2(4),q3(4),q4(4),q(4)
    complex(dp) :: e1(4),e2(4),e3(4),e4(4)
    
    q1 = dcmplx(p(1,:),0d0)
    q2 = dcmplx(p(2,:),0d0)
    q3 = dcmplx(p(3,:),0d0)
    q4 = dcmplx(p(4,:),0d0)
    
    
    e1 = sp(1,:)
    e2 = sp(2,:)
    e3 = sp(3,:)
    e4 = sp(4,:)
    
    q = -q1-q2
    
    q_q =sc(q,q)
    
    
    q1_q2 = sc(q1,q2)
    q1_q3 = sc(q1,q3)
    q1_q4 = sc(q1,q4)
    q2_q3 = sc(q2,q3)
    q2_q4 = sc(q2,q4)
    q3_q4 = sc(q3,q4)
    
    e1_e2 = sc(e1,e2)
    e1_e3 = sc(e1,e3)
    e1_e4 = sc(e1,e4)
    
    e2_e3 = sc(e2,e3)
    e2_e4 = sc(e2,e4)
    
    e3_e4 = sc(e3,e4)
    
    
    q1_e3 = sc(q1,e3)
    q1_e4 = sc(q1,e4)
    q2_e3 = sc(q2,e3)
    q2_e4 = sc(q2,e4)
    e1_q3 = sc(e1,q3)
    e1_q4 = sc(e1,q4)
    e2_q3 = sc(e2,q3)
    e2_q4 = sc(e2,q4)
    e3_q4 = sc(e3,q4)
    e4_q3 = sc(e4,q3)
    
    res = czero
    
    
    res =  &
         + M_Reso**(-2) * (  - q1_q3*q_q*e1_e2*e3_q4*e4_q3 + 1./2.*q1_q3*    &
         q_q**2*e1_e2*e3_e4 - q1_q4*q_q*e1_e2*e3_q4*e4_q3 + 1./2.*       &
         q1_q4*q_q**2*e1_e2*e3_e4 - q2_q3*q_q*e1_e2*e3_q4*e4_q3 + 1./2.  &
         *q2_q3*q_q**2*e1_e2*e3_e4 - q2_q4*q_q*e1_e2*e3_q4*e4_q3 + 1./   &
         2.*q2_q4*q_q**2*e1_e2*e3_e4 - q_q**2*e1_e2*e3_q4*e4_q3 + 1./2.  &
         *q_q**3*e1_e2*e3_e4 )
    
    res = res + q1_e3*q2_e4*q_q*e1_e2 - 2*q1_e3*q2_q4*e1_e2*e4_q3   &
         + q1_e4*q2_e3*q_q*e1_e2 - 2*q1_e4*q2_q3*e1_e2*e3_q4 - 2*       &
         q1_q3*q2_e4*e1_e2*e3_q4 + 2*q1_q3*q2_q4*e1_e2*e3_e4 - 2*q1_q4  &
         *q2_e3*e1_e2*e4_q3 + 2*q1_q4*q2_q3*e1_e2*e3_e4 + q_q*e1_e2*    &
         e3_e4*M_V**2 + 2*q_q*e1_e2*e3_q4*e4_q3 - q_q*e1_e3*e2_q4*e4_q3  &
         - q_q*e1_e4*e2_q3*e3_q4 - q_q*e1_q3*e2_e4*e3_q4 + q_q*e1_q3*   &
         e2_q4*e3_e4 - q_q*e1_q4*e2_e3*e4_q3 + q_q*e1_q4*e2_q3*e3_e4    &
         - q_q**2*e1_e2*e3_e4 + 1./2.*q_q**2*e1_e3*e2_e4 + 1./2.*       &
         q_q**2*e1_e4*e2_e3
    
  end subroutine ggGZZampl_old
  
  
  
  
  
  
  
  
  double complex function et1(e1,e2,e3,e4)
    implicit none
    complex(dp), intent(in) :: e1(4), e2(4), e3(4), e4(4)
    
    et1 =  e1(1)*e2(2)*e3(3)*e4(4)-e1(1)*e2(2)*e3(4)*e4(3) &
         -e1(1)*e2(3)*e3(2)*e4(4)+e1(1)*e2(3)*e3(4)*e4(2) &
         +e1(1)*e2(4)*e3(2)*e4(3)-e1(1)*e2(4)*e3(3)*e4(2) &
         -e1(2)*e2(1)*e3(3)*e4(4)+e1(2)*e2(1)*e3(4)*e4(3) &
         +e1(2)*e2(3)*e3(1)*e4(4)-e1(2)*e2(3)*e3(4)*e4(1) &
         -e1(2)*e2(4)*e3(1)*e4(3)+e1(2)*e2(4)*e3(3)*e4(1) &
         +e1(3)*e2(1)*e3(2)*e4(4)-e1(3)*e2(1)*e3(4)*e4(2) &
         -e1(3)*e2(2)*e3(1)*e4(4)+e1(3)*e2(2)*e3(4)*e4(1) &
         +e1(3)*e2(4)*e3(1)*e4(2)-e1(3)*e2(4)*e3(2)*e4(1) &
         -e1(4)*e2(1)*e3(2)*e4(3)+e1(4)*e2(1)*e3(3)*e4(2) &
         +e1(4)*e2(2)*e3(1)*e4(3)-e1(4)*e2(2)*e3(3)*e4(1) &
         -e1(4)*e2(3)*e3(1)*e4(2)+e1(4)*e2(3)*e3(2)*e4(1)
    
    return
  end function et1
  
  
  
  
  double complex function sc(q1,q2)
    complex(dp), intent(in) :: q1(4)
    complex(dp), intent(in) :: q2(4)
    
    sc = q1(1)*q2(1) - q1(2)*q2(2)-q1(3)*q2(3) -q1(4)*q2(4)
    
  end function sc
  
  double precision function scr(p1,p2)
    real(dp), intent(in) :: p1(4),p2(4)
    
    scr = p1(1)*p2(1) - p1(2)*p2(2)-p1(3)*p2(3) -p1(4)*p2(4)
    
  end function scr
  
  !---- THESE ARE POLARIZATION ROUTINES
  
  ! -- massless vector polarization subroutine
  function pol_mless(p,i,outgoing)
    complex(dp), intent(in)    :: p(4)
    integer, intent(in)          :: i
    logical, intent(in),optional :: outgoing
    ! -------------------------------
    integer :: pol
    real(dp) :: p0,px,py,pz
    real(dp) :: pv,ct,st,cphi,sphi
    complex(dp) :: pol_mless(4)
    
    !^^^IFmp
    !    p0=(p(1)+conjg(p(1)))/two
    !    px=(p(2)+conjg(p(2)))/two
    !    py=(p(3)+conjg(p(3)))/two
    !    pz=(p(4)+conjg(p(4)))/two
    !^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    !^^^END
    
    
    pv=sqrt(abs(p0**2))
    ct=pz/pv
    st=sqrt(abs(1.0_dp-ct**2))
    
    if (st < tol) then
       cphi=1.0_dp
       sphi=0.0_dp
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif
    
    
    ! -- distinguish between positive and negative energies
    if ( p0 > 0.0_dp) then
       pol=i
    else
       pol=-i
    endif
    
    ! -- take complex conjugate for outgoing
    if (present(outgoing)) then
       if (outgoing) pol = -pol
    endif
    
    pol_mless(1)=czero
    pol_mless(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
    pol_mless(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
    pol_mless(4)=-st/sqrt2
    
  end function pol_mless
  
  
  function pol_mless2(p,i,out)
    integer, intent(in)       :: i
    complex(dp), intent(in) :: p(4)
    character(len=*), intent(in):: out
    complex(dp)             :: pol_mless2(4)
    ! -------------------------------------
    
    if (out == 'out') then
       pol_mless2 = pol_mless(p,i,outgoing=.true.)
    else
       pol_mless2 = pol_mless(p,i,outgoing=.false.)
    endif
  end function pol_mless2
  
  
  function pol_dk2mom(plepton,antilepton,i,outgoing)
    integer, intent(in) :: i
    integer :: j
    complex(dp), intent(in) :: plepton(:),antilepton(:)
    logical, intent(in),optional :: outgoing
    complex(dp) :: pol_dk2mom(4),Ub(4),V(4),q(4),qsq
    
    
    q=plepton+antilepton
    qsq=q(1)**2-q(2)**2-q(3)**2-q(4)**2
    
    Ub(:)=ubar0(plepton,i)
    V(:)=v0(antilepton,-i)
    
    !---Now return in Kirill's notation  1=E,2=px,3=py,4=pz
    !   This is an expression for (-i)/qsq* (-i) Ub(+/-)) Gamma^\mu V(-/+)
    pol_dk2mom(1)=-(Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3))
    pol_dk2mom(2)=-(-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3))
    pol_dk2mom(3)=-ci*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
    pol_dk2mom(4)=-(Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3))
    
    do j=1,4
       pol_dk2mom(j)=pol_dk2mom(j)/qsq
    enddo
    
    ! -- do nothing in this case
    if (present(outgoing)) then
       !if (outgoing) pol_dk2mom = conjg(pol_dk2mom)
    endif
  end function pol_dk2mom
  
  
  !     ubar spinor, massless
  
  function ubar0(p,i)
    complex(dp), intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------
    complex(dp) :: ubar0(4)
    complex(dp) :: fc, fc2
    real(dp)    :: p0,px,py,pz
    
    !^^^IFmp
    !    p0=(p(1)+conjg(p(1)))/two
    !    px=(p(2)+conjg(p(2)))/two
    !    py=(p(3)+conjg(p(3)))/two
    !    pz=(p(4)+conjg(p(4)))/two
    !^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    !^^^END
    
    fc2 = p0 + pz
    fc=sqrt(fc2)
    
    if (abs(fc2).gt. tol) then
       
       if (i.eq.1) then
          ubar0(1)=czero
          ubar0(2)=czero
          ubar0(3)=fc
          ubar0(4)=(px-ci*py)/fc
       elseif (i.eq.-1) then
          ubar0(1)=(px+ci*py)/fc
          ubar0(2)=-fc
          ubar0(3)=czero
          ubar0(4)=czero
       else
          stop 'ubar0: i out of range'
       endif
       
    else
       if (i.eq.1) then
          ubar0(1) = czero
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = sqrt(cone*two*p0)
       elseif (i.eq.-1) then
          ubar0(1) = sqrt(cone*(two*p0))
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = czero
       else
          stop 'ubar0: i out of range'
       endif
    endif
    
    
  end function ubar0
  
  
  
  ! -- v0  spinor, massless
  function v0(p,i)
    complex(dp), intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------
    complex(dp) :: v0(4)
    complex(dp) :: fc2, fc
    real(dp)    :: p0,px,py,pz
    
    !^^^IFmp
    !    p0=(p(1)+conjg(p(1)))/two
    !    px=(p(2)+conjg(p(2)))/two
    !    py=(p(3)+conjg(p(3)))/two
    !    pz=(p(4)+conjg(p(4)))/two
    !^^^ELSE
    p0=real(p(1),dp)
    px=real(p(2),dp)
    py=real(p(3),dp)
    pz=real(p(4),dp)
    !^^^END
    
    fc2 = p0 + pz
    fc=sqrt(fc2)
    
    if (abs(fc2).gt. tol) then
       
       if (i.eq.1) then
          v0(1)=czero
          v0(2)=czero
          v0(3)=(px-ci*py)/fc
          v0(4)=-fc
       elseif (i.eq.-1) then
          v0(1)=fc
          v0(2)=(px+ci*py)/fc
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range'
       endif
       
    else
       
       if (i.eq.1) then
          v0(1)=czero
          v0(2)=czero
          v0(3)=sqrt(cone*two*p0)
          v0(4)=czero
       elseif (i.eq.-1) then
          v0(1)=czero
          v0(2)=sqrt(cone*two*p0)
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range'
       endif
       
    endif
    
  end function v0
  
  
  
  
end module modGraviton


