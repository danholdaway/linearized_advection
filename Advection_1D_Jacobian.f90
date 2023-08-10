! ================================================================
!                     MAIN PROGRAM
! ================================================================

PROGRAM advection

IMPLICIT NONE

!Parameters 
!----------
INTEGER,          PARAMETER :: N  = 64                    !Grid size
DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384 !Pi
DOUBLE PRECISION, PARAMETER :: dx = 1.0/N                 !Grid spacing
DOUBLE PRECISION, PARAMETER :: C  = 0.1                   !Courant number
DOUBLE PRECISION, PARAMETER :: dt = C*dx                  !Time step
INTEGER,          PARAMETER :: N1 = ceiling((1.0/(dt)))   !Number of time steps
LOGICAL,          PARAMETER :: COMPUTE_J = .true.

INTEGER,          PARAMETER :: NSTOP = 292

! Initial Conditions
!-------------------
!INTEGER, PARAMETER :: traceval=1 ! Sine wave
INTEGER, PARAMETER :: traceval=2 ! Step function
!INTEGER, PARAMETER :: traceval=3 ! One Cloud
!INTEGER, PARAMETER :: traceval=4 ! Many clouds

!Locals
!------
INTEGER ::  tt 
!INTEGER, PARAMETER :: FVCHOICE = 1
INTEGER, PARAMETER :: FVCHOICE = 2
DOUBLE PRECISION :: x(N), twodt

!Grid
!----
INTEGER :: fx(N), f1x(N), f2x(N), f3x(N), f4x(N), g1x(N),  g2x(N), g3x(N), g4x(N)

!Schemes - free
!--------------
DOUBLE PRECISION :: q_free(N)
DOUBLE PRECISION :: u_free
!Linear schemes
DOUBLE PRECISION :: q1st_free(N),   q1stnew_free(N), q1sthalf_free(N)              !1st Order Upwind Scheme (requires u>1)
DOUBLE PRECISION :: q2nd_free(N),   q2ndnew_free(N), q2ndhalf_free(N)              !2nd Order Lax-Wendroff Scheme
DOUBLE PRECISION :: q3rd_free(N),   q3rdnew_free(N), q3rdhalf_free(N)              !3rd Order Upwind Scheme (requires u>1)
!PPM Schemes
DOUBLE PRECISION :: qppm_unl_free(N),  qppm_unlnew_free(N),  qppm_unlhalf_free(N)  !Unlimited PPM
DOUBLE PRECISION :: qppm_lin_free(N),  qppmnew_lin_free(N),  qppmhalf_lin_free(N)  !PPM with Lin RMC etc
DOUBLE PRECISION :: qppm_cwl_free(N),  qppm_cwlnew_free(N),  qppm_cwlhalf_free(N)  !PPM Colella and Woodward JCP 1984
DOUBLE PRECISION :: qppm_csl_free(N),  qppm_cslnew_free(N),  qppm_cslhalf_free(N)  !PPM Colella and Sekora JCP 2008
!FD with limiters
DOUBLE PRECISION :: q3rdul_free(N), q3rdulnew_free(N), q3rdulhalf_free(N)          !3rd Order Scheme with Universal Limiter (requires u>1)
DOUBLE PRECISION :: q3rdrl_free(N), q3rdrlnew_free(N), q3rdrlhalf_free(N)          !3rd Order Scheme with Quasi-Monotonic Relaxed Limiter (requires u>1)
!Semi-lagrangian schemes
DOUBLE PRECISION :: qsl_free(N),      qslnew_free(N)                               !Simple Semi-Lagrangian Scheme
DOUBLE PRECISION :: qslice_free(N),   qslicenew_free(N)                            !SLICE Scheme
DOUBLE PRECISION :: qslicebs_free(N), qslicebsnew_free(N)                          !SLICE Scheme (with BS limiter)
!Spectral
DOUBLE PRECISION :: qspec_free(N), qspecold_free(N)                                !Spectral Transform Method with Leapfrog Timestep

!Schemes - replay
!----------------
DOUBLE PRECISION :: q_repl(N)
DOUBLE PRECISION :: u_repl
!Linear schemes
DOUBLE PRECISION :: q1st_repl(N), q1stnew_repl(N), q1sthalf_repl(N)
DOUBLE PRECISION :: q2nd_repl(N), q2ndnew_repl(N), q2ndhalf_repl(N)
DOUBLE PRECISION :: q3rd_repl(N), q3rdnew_repl(N), q3rdhalf_repl(N)
!PPM Schemes
DOUBLE PRECISION :: qppm_unl_repl(N), qppm_unlnew_repl(N), qppm_unlhalf_repl(N)
DOUBLE PRECISION :: qppm_lin_repl(N), qppmnew_lin_repl(N), qppmhalf_lin_repl(N)
DOUBLE PRECISION :: qppm_cwl_repl(N), qppm_cwlnew_repl(N), qppm_cwlhalf_repl(N)
DOUBLE PRECISION :: qppm_csl_repl(N), qppm_cslnew_repl(N), qppm_cslhalf_repl(N)
!FD with limiters
DOUBLE PRECISION :: q3rdul_repl(N), q3rdulnew_repl(N), q3rdulhalf_repl(N)
DOUBLE PRECISION :: q3rdrl_repl(N), q3rdrlnew_repl(N), q3rdrlhalf_repl(N)
!Semi-lagrangian schemes
DOUBLE PRECISION :: qsl_repl(N),      qslnew_repl(N)
DOUBLE PRECISION :: qslice_repl(N),   qslicenew_repl(N)
DOUBLE PRECISION :: qslicebs_repl(N), qslicebsnew_repl(N) 
!Spectral
DOUBLE PRECISION :: qspec_repl(N), qspecold_repl(N)

!Schemes tangent linear
!----------------------
DOUBLE PRECISION :: q_traj(N)
DOUBLE PRECISION :: q_tlm(N)
DOUBLE PRECISION :: u_traj
DOUBLE PRECISION :: u_tlm
DOUBLE PRECISION :: qppm_lin_traj(N), qppmnew_lin_traj(N), qppmhalf_lin_traj(N)
DOUBLE PRECISION :: qppm_lin_tlm(N), qppmnew_lin_tlm(N), qppmhalf_lin_tlm(N), qppm_lin_tlm_new(N)
DOUBLE PRECISION :: qppm_unl_traj(N), qppm_unlnew_traj(N), qppm_unlhalf_traj(N)
DOUBLE PRECISION :: qppm_unl_tlm(N), qppm_unlnew_tlm(N), qppm_unlhalf_tlm(N), qppm_unl_tlm_new(N)
DOUBLE PRECISION :: qppm_cwl_traj(N), qppm_cwlnew_traj(N), qppm_cwlhalf_traj(N)
DOUBLE PRECISION :: qppm_cwl_tlm(N), qppm_cwlnew_tlm(N), qppm_cwlhalf_tlm(N), qppm_cwl_tlm_new(N)
DOUBLE PRECISION :: qppm_csl_traj(N), qppm_cslnew_traj(N), qppm_cslhalf_traj(N)
DOUBLE PRECISION :: qppm_csl_tlm(N), qppm_cslnew_tlm(N), qppm_cslhalf_tlm(N), qppm_csl_tlm_new(N)
DOUBLE PRECISION :: q1st_traj(N), q1stnew_traj(N), q1sthalf_traj(N)
DOUBLE PRECISION :: q1st_tlm(N), q1stnew_tlm(N), q1sthalf_tlm(N), q1st_tlm_new(N)
DOUBLE PRECISION :: q2nd_traj(N), q2ndnew_traj(N), q2ndhalf_traj(N)
DOUBLE PRECISION :: q2nd_tlm(N), q2ndnew_tlm(N), q2ndhalf_tlm(N), q2nd_tlm_new(N)
DOUBLE PRECISION :: q3rd_traj(N), q3rdnew_traj(N), q3rdhalf_traj(N)
DOUBLE PRECISION :: q3rd_tlm(N), q3rdnew_tlm(N), q3rdhalf_tlm(N), q3rd_tlm_new(N)
DOUBLE PRECISION :: q3rdul_traj(N), q3rdulnew_traj(N), q3rdulhalf_traj(N)
DOUBLE PRECISION :: q3rdul_tlm(N), q3rdulnew_tlm(N), q3rdulhalf_tlm(N), q3rdul_tlm_new(N)
DOUBLE PRECISION :: q3rdrl_traj(N), q3rdrlnew_traj(N), q3rdrlhalf_traj(N)
DOUBLE PRECISION :: q3rdrl_tlm(N), q3rdrlnew_tlm(N), q3rdrlhalf_tlm(N), q3rdrl_tlm_new(N)
DOUBLE PRECISION :: qsl_traj(N),      qslnew_traj(N)
DOUBLE PRECISION :: qsl_tlm(N),      qslnew_tlm(N), qsl_tlm_new(N)
DOUBLE PRECISION :: qslice_traj(N),   qslicenew_traj(N)
DOUBLE PRECISION :: qslice_tlm(N),   qslicenew_tlm(N), qslice_tlm_new(N)
DOUBLE PRECISION :: qslicebs_traj(N), qslicebsnew_traj(N)
DOUBLE PRECISION :: qslicebs_tlm(N), qslicebsnew_tlm(N), qslicebs_tlm_new(N)
DOUBLE PRECISION :: qspec_traj(N), qspecold_traj(N)
DOUBLE PRECISION :: qspec_tlm(N), qspecold_tlm(N), qspec_tlm_new(N)

!TLM test
!--------
DOUBLE PRECISION :: nlm_pert_traj_qppm(N)
DOUBLE PRECISION :: nlm_pert_traj_qppm_unl(N)
DOUBLE PRECISION :: nlm_pert_traj_qppm_cwl(N)
DOUBLE PRECISION :: nlm_pert_traj_qppm_csl(N)
DOUBLE PRECISION :: nlm_pert_traj_q1st(N)
DOUBLE PRECISION :: nlm_pert_traj_q2nd(N)
DOUBLE PRECISION :: nlm_pert_traj_q3rd(N)
DOUBLE PRECISION :: nlm_pert_traj_q3rdul(N)
DOUBLE PRECISION :: nlm_pert_traj_q3rdrl(N)
DOUBLE PRECISION :: nlm_pert_traj_qsl(N)
DOUBLE PRECISION :: nlm_pert_traj_qslice(N)
DOUBLE PRECISION :: nlm_pert_traj_qslicebs(N)
DOUBLE PRECISION :: nlm_pert_traj_qspec(N)

!Jacobian
!--------
INTEGER :: i,j,k
DOUBLE PRECISION :: JAC(N,N), JACDT(N,N)
DOUBLE PRECISION :: EYE(N,N)
DOUBLE PRECISION :: qtraj(N), qnewtraj(N), qhalftraj(N)
DOUBLE PRECISION :: qpert(N), qnewpert(N), qhalfpert(N)
DOUBLE PRECISION :: utraj, upert
DOUBLE PRECISION :: maxeval_3rd(N1)
DOUBLE PRECISION :: maxeval_ppm_unl(N1)
DOUBLE PRECISION :: maxeval_ppm_lin(N1), maxeval_ppm_cwl(N1), maxeval_ppm_csl(N1)
DOUBLE PRECISION :: maxeval_3rdul(N1), maxeval_3rdrl(N1)
DOUBLE PRECISION :: maxeval_slice(N1), maxeval_slicebs(N1)

!Eigenvalue computation
!----------------------
integer, parameter :: LDA = N
REAL(8) :: WR(N), WI(N)
INTEGER, PARAMETER :: LDVL = N !Left vectors
REAL(8) :: VL(LDVL,N)
INTEGER, PARAMETER :: LDVR = N !Right vectors
REAL(8) :: VR(LDVR,N)
INTEGER :: INFO, LWORK
INTEGER, PARAMETER :: LWMAX = 10000
REAL(8) :: WORK(LWMAX)
EXTERNAL :: DGEEV
REAL*8 :: maxeval

!This program runs tangent linear test on a series of advection routines.
!Two nonlinear integrations give the nonlinear perturbation trajectory.
!One tangnet linear run recoveres the nonlinear perturbation trajectory using the tangent linear model.

 print*,'PROGRAM START'

 !Setup grid
 !----------
 call gridsetup(x,dx,fx,f1x,f2x,f3x,f4x,g1x,g2x,g3x,g4x,N)

 !Initial Conditions - free
 !-------------------------
 u_free = 1.0
 call initialconditions(q_free,x,dx,pi,fx,N,traceval)

 qppm_lin_free(fx)=q_free(fx)
 qppm_unl_free(fx)=q_free(fx)
 qppm_cwl_free(fx)=q_free(fx)
 qppm_csl_free(fx)=q_free(fx)
 q1st_free(fx)=q_free(fx)
 q2nd_free(fx)=q_free(fx)
 q3rd_free(fx)=q_free(fx)
 q3rdul_free(fx)=q_free(fx)
 q3rdrl_free(fx)=q_free(fx)
 qsl_free(fx)=q_free(fx)
 qslice_free(fx)=q_free(fx)
 qslicebs_free(fx)=q_free(fx)
 qspec_free(fx)=q_free(fx)
 qspecold_free(fx)=q_free(fx)

 ! Initial Conditions - replay
 !----------------------------
 u_repl = 1.0
 call initialconditions(q_repl,x,dx,pi,fx,N,traceval)

 !Add a constant to the whole profile (testing only)
 !q_repl = q_repl+1e-4

 !Multiply by the whole profile to make a positive perturbation
 q_repl = q_repl*(1+(1-0.99999))+1e-6

 !Make a nonlinear perturbation by adding between two points
 !q_repl(25:39) = q_repl(25:39)+1e-4

 !Perturb one location
 !q_repl(10) = q_repl(10)+1e-3

 qppm_lin_repl(fx)=q_repl(fx)
 qppm_unl_repl(fx)=q_repl(fx)
 qppm_cwl_repl(fx)=q_repl(fx)
 qppm_csl_repl(fx)=q_repl(fx)
 q1st_repl(fx)=q_repl(fx)
 q2nd_repl(fx)=q_repl(fx)
 q3rd_repl(fx)=q_repl(fx)
 q3rdul_repl(fx)=q_repl(fx)
 q3rdrl_repl(fx)=q_repl(fx)
 qsl_repl(fx)=q_repl(fx)
 qslice_repl(fx)=q_repl(fx)
 qslicebs_repl(fx)=q_repl(fx)
 qspec_repl(fx)=q_repl(fx)
 qspecold_repl(fx)=q_repl(fx)

 !Initial Conditions - tlm
 !-------------------------
 u_traj = u_free
 q_traj = q_free
 u_tlm = u_repl - u_free
 q_tlm = q_repl - q_free

 qppm_lin_traj(fx)=q_free(fx)
 qppm_lin_tlm(fx)=qppm_lin_repl(fx)-qppm_lin_free(fx)

 qppm_unl_traj(fx)=q_free(fx)
 qppm_unl_tlm(fx)=qppm_unl_repl(fx)-qppm_unl_free(fx)

 qppm_cwl_traj(fx)=q_free(fx)
 qppm_cwl_tlm(fx)=qppm_cwl_repl(fx)-qppm_cwl_free(fx)

 qppm_csl_traj(fx)=q_free(fx)
 qppm_csl_tlm(fx)=qppm_csl_repl(fx)-qppm_csl_free(fx)

 q1st_traj(fx)=q_free(fx)
 q1st_tlm (fx)=q1st_repl(fx)-q1st_free(fx)

 q2nd_traj(fx)=q_free(fx)
 q2nd_tlm (fx)=q2nd_repl(fx)-q2nd_free(fx)

 q3rd_traj(fx)=q_free(fx)
 q3rd_tlm (fx)=q3rd_repl(fx)-q3rd_free(fx)

 q3rdul_traj(fx)=q_free(fx)
 q3rdul_tlm (fx)=q3rdul_repl(fx)-q3rdul_free(fx)

 q3rdrl_traj(fx)=q_free(fx)
 q3rdrl_tlm (fx)=q3rdrl_repl(fx)-q3rdrl_free(fx)

 qsl_traj(fx)=q_free(fx)
 qsl_tlm(fx) =q_repl(fx)-q_free(fx)

 qslice_traj(fx)=q_free(fx)
 qslice_tlm(fx) =q_repl(fx)-q_free(fx)

 qslicebs_traj(fx)=q_free(fx)
 qslicebs_tlm(fx) =q_repl(fx)-q_free(fx)

 qspec_traj(fx)=q_free(fx)
 qspec_tlm(fx) =q_repl(fx)-q_free(fx)
 qspecold_traj(fx)=q_free(fx)
 qspecold_tlm(fx) =q_repl(fx)-q_free(fx)

 !Zero out the eigenvalues
 maxeval_ppm_unl = 0.0
 maxeval_ppm_lin = 0.0
 maxeval_ppm_cwl = 0.0
 maxeval_ppm_csl = 0.0
 maxeval_3rd = 0.0
 maxeval_3rdul = 0.0
 maxeval_3rdrl = 0.0
 maxeval_slice = 0.0
 maxeval_slicebs = 0.0

 !Identity matrix
 EYE = 0.0
 DO i = 1,N
    EYE(i,i) = 1.0
 ENDDO

 !Run Advection
 !-------------
 DO tt=1, NSTOP

    !Compute time step for spectral method
    if (tt .eq. 1) then
       twodt=dt
    else
       twodt=2.0*dt
    endif

    !CALL THE (NON)LINEAR ADVECTION ROUTINES FOR THE FREE RUN
    !--------------------------------------------------------
    !PPM without Limiter
    call unlimitedppm1d(qppm_unlhalf_free,qppm_unl_free,u_free,dt,dx,N)
    qppm_unlnew_free(fx) = qppm_unl_free(fx) - (u_free*dt/dx)*(qppm_unlhalf_free(f1x)-qppm_unlhalf_free(fx))
    qppm_unl_free=qppm_unlnew_free

    !PPM with Lin limiter
    call ppmlin1d(qppmhalf_lin_free,qppm_lin_free,u_free,dt,dx,N)
    qppmnew_lin_free(fx) = qppm_lin_free(fx) - (u_free*dt/dx)*(qppmhalf_lin_free(f1x)-qppmhalf_lin_free(fx))
    qppm_lin_free=qppmnew_lin_free

    !PPM with Colella-Woodward (JCP 1984) limiter
    call ppmoriginal1d(qppm_cwlhalf_free,qppm_cwl_free,u_free,dt,dx,N)
    qppm_cwlnew_free(fx) = qppm_cwl_free(fx) - (u_free*dt/dx)*(qppm_cwlhalf_free(f1x)-qppm_cwlhalf_free(fx))
    qppm_cwl_free=qppm_cwlnew_free

    !PPM with Colella-Sekora (JCP 2008) limiter
    call ppmcs1d(qppm_cslhalf_free,qppm_csl_free,u_free,dt,dx,N)
    qppm_cslnew_free(fx) = qppm_csl_free(fx) - (u_free*dt/dx)*(qppm_cslhalf_free(f1x)-qppm_cslhalf_free(fx))
    qppm_csl_free=qppm_cslnew_free
    
    !FD/FV unlimited 3rd order upwind
    call getqhalf(q3rdhalf_free,q3rd_free,3,N,C,fx,f1x,g1x,g2x,g3x)
    q3rdnew_free(fx) = q3rd_free(fx) - (u_free*dt/dx)*(q3rdhalf_free(f1x)-q3rdhalf_free(fx))
    q3rd_free=q3rdnew_free

    !FD/FV universal limiter 3rd order upwind
    call getqhalf(q3rdulhalf_free,q3rdul_free,3,N,C,fx,f1x,g1x,g2x,g3x)
    call universallimiter(q3rdulhalf_free,q3rdul_free,u_free,dt,dx,N)
    q3rdulnew_free(fx) = q3rdul_free(fx) - (u_free*dt/dx)*(q3rdulhalf_free(f1x)-q3rdulhalf_free(fx))
    q3rdul_free=q3rdulnew_free

    !FD/FV relaxed limiter 3rd order upwind
    call getqhalf(q3rdrlhalf_free,q3rdrl_free,3,N,C,fx,f1x,g1x,g2x,g3x)
    call relaxedlimiter(q3rdrlhalf_free,q3rdrl_free,u_free,dt,dx,N)
    q3rdrlnew_free(fx) = q3rdrl_free(fx) - (u_free*dt/dx)*(q3rdrlhalf_free(f1x)-q3rdrlhalf_free(fx))
    q3rdrl_free=q3rdrlnew_free

    !SLICE Scheme
    call slicescheme(qslicenew_free,qslice_free,C,x,dx,N)
    qslice_free=qslicenew_free

    !SLICE Scheme with BS limiter
    call sliceschemebs(qslicebsnew_free,qslicebs_free,C,x,dx,N)
    qslicebs_free=qslicebsnew_free

    !CALL THE (NON)LINEAR ADVECTION ROUTINES FOR THE REPLAY RUN
    !----------------------------------------------------------
    !PPM without Limiter
    call unlimitedppm1d(qppm_unlhalf_repl,qppm_unl_repl,u_repl,dt,dx,N)
    qppm_unlnew_repl(fx) = qppm_unl_repl(fx) - (u_repl*dt/dx)*(qppm_unlhalf_repl(f1x)-qppm_unlhalf_repl(fx))
    qppm_unl_repl=qppm_unlnew_repl

    !PPM with limiter
    call ppmlin1d(qppmhalf_lin_repl,qppm_lin_repl,u_repl,dt,dx,N)
    qppmnew_lin_repl(fx) = qppm_lin_repl(fx) - (u_repl*dt/dx)*(qppmhalf_lin_repl(f1x)-qppmhalf_lin_repl(fx))
    qppm_lin_repl=qppmnew_lin_repl

    !PPM with Colella-Woodward (JCP 1984) limiter
    call ppmoriginal1d(qppm_cwlhalf_repl,qppm_cwl_repl,u_repl,dt,dx,N)
    qppm_cwlnew_repl(fx) = qppm_cwl_repl(fx) - (u_repl*dt/dx)*(qppm_cwlhalf_repl(f1x)-qppm_cwlhalf_repl(fx))
    qppm_cwl_repl=qppm_cwlnew_repl

    !PPM with Colella-Sekora (JCP 2008) limiter
    call ppmcs1d(qppm_cslhalf_repl,qppm_csl_repl,u_repl,dt,dx,N)
    qppm_cslnew_repl(fx) = qppm_csl_repl(fx) - (u_repl*dt/dx)*(qppm_cslhalf_repl(f1x)-qppm_cslhalf_repl(fx))
    qppm_csl_repl=qppm_cslnew_repl

    !FD/FV unlimited 3rd order upwind
    call getqhalf(q3rdhalf_repl,q3rd_repl,3,N,C,fx,f1x,g1x,g2x,g3x)
    q3rdnew_repl(fx) = q3rd_repl(fx) - (u_repl*dt/dx)*(q3rdhalf_repl(f1x)-q3rdhalf_repl(fx))
    q3rd_repl=q3rdnew_repl

    !FD/FV universal limiter 3rd order upwind
    call getqhalf(q3rdulhalf_repl,q3rdul_repl,3,N,C,fx,f1x,g1x,g2x,g3x)
    call universallimiter(q3rdulhalf_repl,q3rdul_repl,u_repl,dt,dx,N)
    q3rdulnew_repl(fx) = q3rdul_repl(fx) - (u_repl*dt/dx)*(q3rdulhalf_repl(f1x)-q3rdulhalf_repl(fx))
    q3rdul_repl=q3rdulnew_repl

    !FD/FV relaxed limiter 3rd order upwind
    call getqhalf(q3rdrlhalf_repl,q3rdrl_repl,3,N,C,fx,f1x,g1x,g2x,g3x)
    call relaxedlimiter(q3rdrlhalf_repl,q3rdrl_repl,u_repl,dt,dx,N)
    q3rdrlnew_repl(fx) = q3rdrl_repl(fx) - (u_repl*dt/dx)*(q3rdrlhalf_repl(f1x)-q3rdrlhalf_repl(fx))
    q3rdrl_repl=q3rdrlnew_repl

    !SLICE Scheme
    call slicescheme(qslicenew_repl,qslice_repl,C,x,dx,N)
    qslice_repl=qslicenew_repl

    !SLICE Scheme with BS limiter
    call sliceschemebs(qslicebsnew_repl,qslicebs_repl,C,x,dx,N)
    qslicebs_repl=qslicebsnew_repl


    !CALL THE TANGENT LINEAR ADVECTION ROUTINES
    !------------------------------------------

    !PPM without limiter
    if (compute_j .eqv. .true.) then

       !Zero Jacobian
       JAC = 0.0
       JACDT = 0.0
       do k = 1,N

          qtraj = qppm_unl_traj
          utraj = u_traj

          qpert    = 0.0
          qpert(k) = 1.0 

          upert = 0.0

          call unlimitedppm1d_tlm(qhalftraj, qhalfpert, qtraj, qpert, utraj, upert, dt, dx, N)
          qnewpert(fx)  = qpert(fx)  - (utraj*dt/dx)*(qhalfpert (f1x)-qhalfpert (fx)) &
                                     - (upert *dt/dx)*(qhalftraj(f1x)-qhalftraj(fx))
          qpert=qnewpert

          JAC(:,k) = qpert

       enddo

       JACDT = (JAC - EYE)/DT

       if (tt == NSTOP) then
          do i=1,N
             do j = 1,N
                write(1,*), 'Jppmunl(',i,',',j,') = ', JACDT(i,j), ';'
             enddo
          enddo
       endif 

       !Move traj forward
       call unlimitedppm1d(qppm_unlhalf_traj, qppm_unl_traj, u_traj, dt, dx, N)
       qppm_unlnew_traj(fx) = qppm_unl_traj(fx) - (u_traj*dt/dx)*(qppm_unlhalf_traj(f1x)-qppm_unlhalf_traj(fx))
       qppm_unl_traj=qppm_unlnew_traj

       !Use Jacobian to advect pert
       qppm_unl_tlm_new = 0.0
       do i = 1,N
          do j = 1,N
             !qppm_lin_tlm_new(i) = qppm_lin_tlm_new(i) + JAC(i,j) * qppm_lin_tlm(j)
             qppm_unl_tlm_new(i) = qppm_unl_tlm_new(i) + EYE(i,j) * qppm_unl_tlm(j) + DT* JACDT(i,j) * qppm_unl_tlm(j)
          enddo
       enddo
       qppm_unl_tlm = qppm_unl_tlm_new

       !Compute the eigenvalues
       LWORK = -1
       work = 0.0
       info = 0
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       maxeval_ppm_unl(tt) = maxval(WR)

    else

       call unlimitedppm1d_tlm(qppm_unlhalf_traj, qppm_unlhalf_tlm, qppm_unl_traj, qppm_unl_tlm, u_traj, u_tlm, dt, dx, N)
       qppm_unlnew_tlm(fx)  = qppm_unl_tlm(fx)  - (u_traj*dt/dx)*(qppm_unlhalf_tlm (f1x)-qppm_unlhalf_tlm (fx)) &
                                                - (u_tlm *dt/dx)*(qppm_unlhalf_traj(f1x)-qppm_unlhalf_traj(fx))
       qppm_unlnew_traj(fx)  = qppm_unl_traj(fx)  - (u_traj*dt/dx)*(qppm_unlhalf_traj(f1x)-qppm_unlhalf_traj(fx))
       qppm_unl_tlm=qppm_unlnew_tlm
       qppm_unl_traj=qppm_unlnew_traj

    endif

    !PPM with limiter
    if (compute_j .eqv. .true.) then

       !Zero Jacobian
       JAC = 0.0
       JACDT = 0.0
       do k = 1,N

          qtraj = qppm_lin_traj
          utraj = u_traj

          qpert    = 0.0
          qpert(k) = 1.0 

          upert = 0.0

          call ppmlin1d_tlm(qhalftraj, qhalfpert, qtraj, qpert, utraj, upert, dt, dx, N)
          qnewpert(fx)  = qpert(fx)  - (utraj*dt/dx)*(qhalfpert(f1x)-qhalfpert(fx)) &
                                     - (upert*dt/dx)*(qhalftraj(f1x)-qhalftraj(fx))
          qpert=qnewpert

          JAC(:,k) = qpert

       enddo

       JACDT = (JAC - EYE)/DT

       if (tt == NSTOP) then
          do i=1,N
             do j = 1,N
                write(1,*), 'Jlin(',i,',',j,') = ', JACDT(i,j), ';'
             enddo
          enddo
       endif 

       !Move traj forward
       call ppmlin1d(qppmhalf_lin_traj,qppm_lin_traj,u_traj,dt,dx,N)
       qppmnew_lin_traj(fx) = qppm_lin_traj(fx) - (u_traj*dt/dx)*(qppmhalf_lin_traj(f1x)-qppmhalf_lin_traj(fx))
       qppm_lin_traj=qppmnew_lin_traj

       !Use Jacobian to advect pert
       qppm_lin_tlm_new = 0.0
       do i = 1,N
          do j = 1,N
             !qppm_lin_tlm_new(i) = qppm_lin_tlm_new(i) + JAC(i,j) * qppm_lin_tlm(j)
             qppm_lin_tlm_new(i) = qppm_lin_tlm_new(i) + EYE(i,j) * qppm_lin_tlm(j) + DT* JACDT(i,j) * qppm_lin_tlm(j)
          enddo
       enddo
       qppm_lin_tlm = qppm_lin_tlm_new

       !Compute the eigenvalues
       LWORK = -1
       work = 0.0
       info = 0
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       maxeval_ppm_lin(tt) = maxval(WR)

    else

       call ppmlin1d_tlm(qppmhalf_lin_traj, qppmhalf_lin_tlm, qppm_lin_traj, qppm_lin_tlm, u_traj, u_tlm, dt, dx, N)
       qppmnew_lin_tlm(fx)  = qppm_lin_tlm(fx)  - (u_traj*dt/dx)*(qppmhalf_lin_tlm (f1x)-qppmhalf_lin_tlm (fx)) &
                                                - (u_tlm *dt/dx)*(qppmhalf_lin_traj(f1x)-qppmhalf_lin_traj(fx))
       qppmnew_lin_traj(fx) = qppm_lin_traj(fx) - (u_repl*dt/dx)*(qppmhalf_lin_traj(f1x)-qppmhalf_lin_traj(fx))
       qppm_lin_tlm=qppmnew_lin_tlm
       qppm_lin_traj=qppmnew_lin_traj

    endif


    !PPM with Colella-Woodward (JCP 1984) limiter
    if (compute_j .eqv. .true.) then

       !Zero Jacobian
       JAC = 0.0
       JACDT = 0.0
       do k = 1,N

          qtraj = qppm_cwl_traj
          utraj = u_traj

          qpert    = 0.0
          qpert(k) = 1.0 

          upert = 0.0

          call ppmoriginal1d_tlm(qhalftraj, qhalfpert, qtraj, qpert, utraj, upert, dt, dx, N)
          qnewpert(fx)  = qpert(fx)  - (utraj*dt/dx)*(qhalfpert(f1x)-qhalfpert(fx)) &
                                     - (upert*dt/dx)*(qhalftraj(f1x)-qhalftraj(fx))
          qpert=qnewpert

          JAC(:,k) = qpert

       enddo

       JACDT = (JAC - EYE)/DT

       if (tt == NSTOP) then
          do i=1,N
             do j = 1,N
                write(1,*), 'Jcwl(',i,',',j,') = ', JACDT(i,j), ';'
             enddo
          enddo
       endif 

       !Move traj forward
       call ppmoriginal1d(qppm_cwlhalf_traj, qppm_cwl_traj, u_traj, dt, dx, N)
       qppm_cwlnew_traj(fx)  = qppm_cwl_traj(fx)  - (u_traj*dt/dx)*(qppm_cwlhalf_traj(f1x)-qppm_cwlhalf_traj(fx)) 
       qppm_cwl_traj=qppm_cwlnew_traj

       !Use Jacobian to advect pert
       qppm_cwl_tlm_new = 0.0
       do i = 1,N
          do j = 1,N
             qppm_cwl_tlm_new(i) = qppm_cwl_tlm_new(i) + EYE(i,j) * qppm_cwl_tlm(j) + DT* JACDT(i,j) * qppm_cwl_tlm(j)
          enddo
       enddo
       qppm_cwl_tlm = qppm_cwl_tlm_new

       !Compute the eigenvalues
       LWORK = -1
       work = 0.0
       info = 0
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       maxeval_ppm_cwl(tt) = maxval(WR)


    else

       call ppmoriginal1d_tlm(qppm_cwlhalf_traj, qppm_cwlhalf_tlm, qppm_cwl_traj, qppm_cwl_tlm, u_traj, u_tlm, dt, dx, N)
       qppm_cwlnew_tlm(fx)  = qppm_cwl_tlm(fx)  - (u_traj*dt/dx)*(qppm_cwlhalf_tlm (f1x)-qppm_cwlhalf_tlm (fx)) &
                                                - (u_tlm *dt/dx)*(qppm_cwlhalf_traj(f1x)-qppm_cwlhalf_traj(fx))
       qppm_cwlnew_traj(fx)  = qppm_cwl_traj(fx)  - (u_traj*dt/dx)*(qppm_cwlhalf_traj(f1x)-qppm_cwlhalf_traj(fx)) 
       qppm_cwl_tlm=qppm_cwlnew_tlm
       qppm_cwl_traj=qppm_cwlnew_traj

    endif

    !PPM with Colella-Sekora (JCP 2008) limiter
    if (compute_j .eqv. .true.) then

       !Zero Jacobian
       JAC = 0.0
       JACDT = 0.0
       do k = 1,N

          qtraj = qppm_csl_traj
          utraj = u_traj

          qpert    = 0.0
          qpert(k) = 1.0 

          upert = 0.0

          call ppmcs1d_tlm(qhalftraj, qhalfpert, qtraj, qpert, utraj, upert, dt, dx, N)
          qnewpert(fx)  = qpert(fx)  - (utraj*dt/dx)*(qhalfpert(f1x)-qhalfpert(fx)) &
                                     - (upert*dt/dx)*(qhalftraj(f1x)-qhalftraj(fx))
          qpert=qnewpert

          JAC(:,k) = qpert

       enddo

       JACDT = (JAC - EYE)/DT

       if (tt == NSTOP) then
          do i=1,N
             do j = 1,N
                write(1,*), 'Jcsl(',i,',',j,') = ', JACDT(i,j), ';'
             enddo
          enddo
       endif 

       !Move traj forward
       call ppmcs1d(qppm_cslhalf_traj, qppm_csl_traj, u_traj, dt, dx, N)
       qppm_cslnew_traj(fx)  = qppm_csl_traj(fx)  - (u_traj*dt/dx)*(qppm_cslhalf_traj(f1x)-qppm_cslhalf_traj(fx)) 
       qppm_csl_traj=qppm_cslnew_traj

       !Use Jacobian to advect pert
       qppm_csl_tlm_new = 0.0
       do i = 1,N
          do j = 1,N
             qppm_csl_tlm_new(i) = qppm_csl_tlm_new(i) + EYE(i,j) * qppm_csl_tlm(j) + DT* JACDT(i,j) * qppm_csl_tlm(j)
          enddo
       enddo
       qppm_csl_tlm = qppm_csl_tlm_new

       !Compute the eigenvalues
       LWORK = -1
       work = 0.0
       info = 0
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       maxeval_ppm_csl(tt) = maxval(WR)

    else

       call ppmcs1d_tlm(qppm_cslhalf_traj, qppm_cslhalf_tlm, qppm_csl_traj, qppm_csl_tlm, u_traj, u_tlm, dt, dx, N)
       qppm_cslnew_tlm(fx)  = qppm_csl_tlm(fx)  - (u_traj*dt/dx)*(qppm_cslhalf_tlm (f1x)-qppm_cslhalf_tlm (fx)) &
                                                - (u_tlm *dt/dx)*(qppm_cslhalf_traj(f1x)-qppm_cslhalf_traj(fx))
       qppm_cslnew_traj(fx)  = qppm_csl_traj(fx)  - (u_traj*dt/dx)*(qppm_cslhalf_traj(f1x)-qppm_cslhalf_traj(fx)) 
       qppm_csl_tlm=qppm_cslnew_tlm
       qppm_csl_traj=qppm_cslnew_traj

    endif

    !FD/FV unlimited 3rd order upwind
    if (compute_j .eqv. .true.) then

       !Zero Jacobian
       JAC = 0.0
       JACDT = 0.0
       do k = 1,N

          qtraj = q3rd_traj
          utraj = u_traj

          qpert    = 0.0
          qpert(k) = 1.0 

          upert = 0.0


          call getqhalf_tlm(qhalftraj,qhalfpert,qtraj,qpert,3,N,C,fx,f1x,g1x,g2x,g3x)
          qnewpert (fx) = qpert (fx) - (utraj*dt/dx)*(qhalfpert(f1x)-qhalfpert(fx)) &
                                               - (upert*dt/dx)*(qhalftraj(f1x)-qhalftraj(fx))
          qpert =qnewpert

          JAC(:,k) = qpert

       enddo

       JACDT = (JAC - EYE)/DT

       if (tt == NSTOP) then
          do i=1,N
             do j = 1,N
                write(1,*), 'J3rd(',i,',',j,') = ', JACDT(i,j), ';'
             enddo
          enddo
       endif 

       !Move traj forward
       call getqhalf(q3rdhalf_traj,q3rd_traj,3,N,C,fx,f1x,g1x,g2x,g3x)
       q3rdnew_traj(fx) = q3rd_traj(fx) - (u_traj*dt/dx)*(q3rdhalf_traj(f1x)-q3rdhalf_traj(fx))
       q3rd_traj=q3rdnew_traj

       !Use Jacobian to advect pert
       q3rd_tlm_new = 0.0
       do i = 1,N
          do j = 1,N
             q3rd_tlm_new(i) = q3rd_tlm_new(i) + EYE(i,j) * q3rd_tlm(j) + DT* JACDT(i,j) * q3rd_tlm(j)
          enddo
       enddo
       q3rd_tlm = q3rd_tlm_new

       !Compute the eigenvalues
       LWORK = -1
       work = 0.0
       info = 0
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       maxeval_3rd(tt) = maxval(WR)

    else

       call getqhalf_tlm(q3rdhalf_traj,q3rdhalf_tlm,q3rd_traj,q3rd_tlm,3,N,C,fx,f1x,g1x,g2x,g3x)
       q3rdnew_tlm(fx) = q3rd_tlm(fx) - (u_traj*dt/dx)*(q3rdhalf_tlm (f1x)-q3rdhalf_tlm (fx)) &
                                      - (u_tlm *dt/dx)*(q3rdhalf_traj(f1x)-q3rdhalf_traj(fx))
       q3rdnew_traj(fx) = q3rd_traj(fx) - (u_traj*dt/dx)*(q3rdhalf_traj(f1x)-q3rdhalf_traj(fx)) 
       q3rd_tlm=q3rdnew_tlm
       q3rd_traj=q3rdnew_traj

    endif


    !FD/FV universal limiter 3rd order upwind
    if (compute_j .eqv. .true.) then

       !Zero Jacobian
       JAC = 0.0
       JACDT = 0.0
       do k = 1,N

          qtraj = q3rdul_traj
          utraj = u_traj

          qpert    = 0.0
          qpert(k) = 1.0 

          upert = 0.0


          call getqhalf_tlm(qhalftraj,qhalfpert,qtraj,qpert,3,N,C,fx,f1x,g1x,g2x,g3x)
          call universallimiter_tlm(qhalftraj,qhalfpert,qtraj,qpert,utraj,upert,dt,dx,N)
          qnewpert (fx) = qpert (fx) - (utraj*dt/dx)*(qhalfpert(f1x)-qhalfpert(fx)) &
                                               - (upert*dt/dx)*(qhalftraj(f1x)-qhalftraj(fx))
          qpert =qnewpert

          JAC(:,k) = qpert

       enddo

       JACDT = (JAC - EYE)/DT

       if (tt == NSTOP) then
          do i=1,N
             do j = 1,N
                write(1,*), 'J3rdul(',i,',',j,') = ', JACDT(i,j), ';'
             enddo
          enddo
       endif 

       !Move traj forward
       call getqhalf(q3rdulhalf_traj,q3rdul_traj,3,N,C,fx,f1x,g1x,g2x,g3x)
       call universallimiter(q3rdulhalf_traj,q3rdul_traj,u_traj,dt,dx,N)
       q3rdulnew_traj(fx) = q3rdul_traj(fx) - (u_traj*dt/dx)*(q3rdulhalf_traj(f1x)-q3rdulhalf_traj(fx))
       q3rdul_traj=q3rdulnew_traj

       !Use Jacobian to advect pert
       q3rdul_tlm_new = 0.0
       do i = 1,N
          do j = 1,N
             q3rdul_tlm_new(i) = q3rdul_tlm_new(i) + EYE(i,j) * q3rdul_tlm(j) + DT* JACDT(i,j) * q3rdul_tlm(j)
          enddo
       enddo
       q3rdul_tlm = q3rdul_tlm_new

       !Compute the eigenvalues
       LWORK = -1
       work = 0.0
       info = 0
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       maxeval_3rdul(tt) = maxval(WR)

    else

       call getqhalf_tlm(q3rdulhalf_traj,q3rdulhalf_tlm,q3rdul_traj,q3rdul_tlm,3,N,C,fx,f1x,g1x,g2x,g3x)
       call universallimiter_tlm(q3rdulhalf_traj,q3rdulhalf_tlm,q3rdul_traj,q3rdul_tlm,u_traj,u_tlm,dt,dx,N)
       q3rdulnew_tlm (fx) = q3rdul_tlm (fx) - (u_traj*dt/dx)*(q3rdulhalf_tlm (f1x)-q3rdulhalf_tlm (fx)) &
                                            - (u_tlm *dt/dx)*(q3rdulhalf_traj(f1x)-q3rdulhalf_traj(fx))
       q3rdulnew_traj(fx) = q3rdul_traj(fx) - (u_traj*dt/dx)*(q3rdulhalf_traj(f1x)-q3rdulhalf_traj(fx))
       q3rdul_tlm =q3rdulnew_tlm
       q3rdul_traj=q3rdulnew_traj

    endif


    !FD/FV relaxed limiter 3rd order upwind
    if (compute_j .eqv. .true.) then

       !Zero Jacobian
       JAC = 0.0
       JACDT = 0.0
       do k = 1,N

          qtraj = q3rdrl_traj
          utraj = u_traj

          qpert    = 0.0
          qpert(k) = 1.0 

          upert = 0.0

          call getqhalf_tlm(qhalftraj,qhalfpert,qtraj,qpert,3,N,C,fx,f1x,g1x,g2x,g3x)
          call relaxedlimiter_tlm(qhalftraj,qhalfpert,qtraj,qpert,utraj,upert,dt,dx,N)
          qnewpert (fx) = qpert (fx) - (utraj*dt/dx)*(qhalfpert(f1x)-qhalfpert(fx)) &
                                     - (upert*dt/dx)*(qhalftraj(f1x)-qhalftraj(fx))
          qpert =qnewpert

          JAC(:,k) = qpert

       enddo

       JACDT = (JAC - EYE)/DT

       if (tt == NSTOP) then
          do i=1,N
             do j = 1,N
                write(1,*), 'J3rdrl(',i,',',j,') = ', JACDT(i,j), ';'
             enddo
          enddo
       endif 

       !Move traj forward
       call getqhalf(q3rdrlhalf_traj,q3rdrl_traj,3,N,C,fx,f1x,g1x,g2x,g3x)
       call relaxedlimiter(q3rdrlhalf_traj,q3rdrl_traj,u_traj,dt,dx,N)
       q3rdrlnew_traj(fx) = q3rdrl_traj(fx) - (u_traj*dt/dx)*(q3rdrlhalf_traj(f1x)-q3rdrlhalf_traj(fx))
       q3rdrl_traj=q3rdrlnew_traj

       !Use Jacobian to advect pert
       q3rdrl_tlm_new = 0.0
       do i = 1,N
          do j = 1,N
             q3rdrl_tlm_new(i) = q3rdrl_tlm_new(i) + EYE(i,j) * q3rdrl_tlm(j) + DT* JACDT(i,j) * q3rdrl_tlm(j)
          enddo
       enddo
       q3rdrl_tlm = q3rdrl_tlm_new

       !Compute the eigenvalues
       LWORK = -1
       work = 0.0
       info = 0
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       maxeval_3rdrl(tt) = maxval(WR)

    else

       call getqhalf_tlm(q3rdrlhalf_traj,q3rdrlhalf_tlm,q3rdrl_traj,q3rdrl_tlm,3,N,C,fx,f1x,g1x,g2x,g3x)
       call relaxedlimiter_tlm(q3rdrlhalf_traj,q3rdrlhalf_tlm,q3rdrl_traj,q3rdrl_tlm,u_traj,u_tlm,dt,dx,N)
       q3rdrlnew_tlm (fx) = q3rdrl_tlm (fx) - (u_traj*dt/dx)*(q3rdrlhalf_tlm (f1x)-q3rdrlhalf_tlm (fx)) &
                                            - (u_tlm *dt/dx)*(q3rdrlhalf_traj(f1x)-q3rdrlhalf_traj(fx))
       q3rdrlnew_traj(fx) = q3rdrl_traj(fx) - (u_traj*dt/dx)*(q3rdrlhalf_traj(f1x)-q3rdrlhalf_traj(fx))
       q3rdrl_tlm =q3rdrlnew_tlm
       q3rdrl_traj=q3rdrlnew_traj

    endif

    !SLICE scheme
    if (compute_j .eqv. .true.) then

       !Zero Jacobian
       JAC = 0.0
       JACDT = 0.0
       do k = 1,N

          qtraj = qslice_traj
          utraj = u_traj

          qpert    = 0.0
          qpert(k) = 1.0 

          upert = 0.0

          call slicescheme_tlm(qnewtraj,qnewpert,qtraj,qpert,C,x,dx,N)
          qpert=qnewpert

          JAC(:,k) = qpert

       enddo

       JACDT = (JAC - EYE)/DT

       if (tt == NSTOP) then
          do i=1,N
             do j = 1,N
                write(1,*), 'Jslice(',i,',',j,') = ', JACDT(i,j), ';'
             enddo
          enddo
       endif 
       

       !Move traj forward
       call slicescheme(qslicenew_traj,qslice_traj,C,x,dx,N)
       qslice_traj=qslicenew_traj

       !Use Jacobian to advect pert
       qslice_tlm_new = 0.0
       do i = 1,N
          do j = 1,N
             qslice_tlm_new(i) = qslice_tlm_new(i) + EYE(i,j) * qslice_tlm(j) + DT* JACDT(i,j) * qslice_tlm(j)
          enddo
       enddo
       qslice_tlm = qslice_tlm_new

       !Compute the eigenvalues
       LWORK = -1
       work = 0.0
       info = 0
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       maxeval_slice(tt) = maxval(WR)

    else

       call slicescheme_tlm(qslicenew_traj,qslicenew_tlm,qslice_traj,qslice_tlm,C,x,dx,N)
       qslice_traj=qslicenew_traj
       qslice_tlm=qslicenew_tlm

    endif


    !SLICE with BS Scheme
    if (compute_j .eqv. .true.) then

       !Zero Jacobian
       JAC = 0.0
       JACDT = 0.0
       do k = 1,N

          qtraj = qslicebs_traj
          utraj = u_traj

          qpert    = 0.0
          qpert(k) = 1.0 

          upert = 0.0

          call sliceschemebs_tlm(qnewtraj,qnewpert,qtraj,qpert,C,x,dx,N)
          qpert=qnewpert

          JAC(:,k) = qpert

       enddo

       JACDT = (JAC - EYE)/DT

       if (tt == NSTOP) then
          do i=1,N
             do j = 1,N
                write(1,*), 'Jslicebs(',i,',',j,') = ', JACDT(i,j), ';'
             enddo
          enddo

          do i=1,N
             write(1,*), 'qtraj(',i,') = ', qtraj(i), ';'
          enddo

       endif 

       !Move traj forward
       call sliceschemebs(qslicebsnew_traj,qslicebs_traj,C,x,dx,N)
       qslicebs_traj=qslicebsnew_traj

       !Use Jacobian to advect pert
       qslicebs_tlm_new = 0.0
       do i = 1,N
          do j = 1,N
             qslicebs_tlm_new(i) = qslicebs_tlm_new(i) + EYE(i,j) * qslicebs_tlm(j) + DT* JACDT(i,j) * qslicebs_tlm(j)
          enddo
       enddo
       qslicebs_tlm = qslicebs_tlm_new

       !Compute the eigenvalues
       LWORK = -1
       work = 0.0
       info = 0
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
       CALL DGEEV( 'Vectors', 'Vectors', N, JACDT, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       maxeval_slicebs(tt) = maxval(WR)

    else

       call sliceschemebs_tlm(qslicebsnew_traj,qslicebsnew_tlm,qslicebs_traj,qslicebs_tlm,C,x,dx,N)
       qslicebs_traj=qslicebsnew_traj
       qslicebs_tlm=qslicebsnew_tlm

    endif


 ENDDO


 !Compute difference between nonlinear runs at the endtime
 !--------------------------------------------------------
 nlm_pert_traj_qppm_unl = qppm_unl_repl - qppm_unl_free
 nlm_pert_traj_qppm     = qppm_lin_repl - qppm_lin_free
 nlm_pert_traj_qppm_cwl = qppm_cwl_repl - qppm_cwl_free
 nlm_pert_traj_qppm_csl = qppm_csl_repl - qppm_csl_free
 nlm_pert_traj_q3rd     = q3rd_repl     - q3rd_free
 nlm_pert_traj_q3rdul   = q3rdul_repl   - q3rdul_free
 nlm_pert_traj_q3rdrl   = q3rdrl_repl   - q3rdrl_free
 nlm_pert_traj_qslice   = qslice_repl   - qslice_free
 nlm_pert_traj_qslicebs = qslicebs_repl - qslicebs_free


 !Print out measures of the error between the nonlinea and linear perturbation trajectories
 !-----------------------------------------------------------------------------------------
 print*, ''
 print*, '==================================================='
 print*, 'STANDARD TLM GRADIENT TESTS'
 print*, '---------------------------'
 print*, 'TLM TEST PPM Unlimited       :', sum(nlm_pert_traj_qppm_unl/qppm_unl_tlm)/size(nlm_pert_traj_qppm_unl)
 print*, 'TLM TEST PPM Lin RMC etc     :', sum(nlm_pert_traj_qppm/qppm_lin_tlm)/size(nlm_pert_traj_qppm)
 print*, 'TLM TEST PPM with CW Lim     :', sum(nlm_pert_traj_qppm_cwl/qppm_cwl_tlm)/size(nlm_pert_traj_qppm_cwl)
 print*, 'TLM TEST PPM with CS Lim     :', sum(nlm_pert_traj_qppm_csl/qppm_csl_tlm)/size(nlm_pert_traj_qppm_csl)
 print*, 'TLM TEST 3rd Unlimited       :', sum(nlm_pert_traj_q3rd/q3rd_tlm)/size(nlm_pert_traj_q3rd)
 print*, 'TLM TEST 3rd Relaxed Lim     :', sum(nlm_pert_traj_q3rdrl/q3rdrl_tlm)/size(nlm_pert_traj_q3rdrl)
 print*, 'TLM TEST 3rd Universal Lim   :', sum(nlm_pert_traj_q3rdul/q3rdul_tlm)/size(nlm_pert_traj_q3rdul)
 print*, 'TLM TEST Semi-Lag (SLICE)    :', sum(nlm_pert_traj_qslice/qslice_tlm)/size(nlm_pert_traj_qslice)
 print*, 'TLM TEST Semi-Lag (SLICE BS) :', sum(nlm_pert_traj_qslicebs/qslicebs_tlm)/size(nlm_pert_traj_qslicebs)
 print*, '==================================================='
 print*, ''



 !Output for MATLAB
 !-----------------
 call dump1(q_free,'q_free',N)
 call dump1(qppm_unl_free,'qppm_unl_free',N)
 call dump1(qppm_lin_free,'qppm_lin_free',N)
 call dump1(qppm_cwl_free,'qppm_cwl_free',N)
 call dump1(qppm_csl_free,'qppm_csl_free',N)
 call dump1(q3rd_free,'q3rd_free',N)
 call dump1(q3rdul_free,'q3rdul_free',N)
 call dump1(q3rdrl_free,'q3rdrl_free',N)
 call dump1(qslice_free,'qslice_free',N)
 call dump1(qslicebs_free,'qslicebs_free',N)

 call dump1(q_repl,'q_repl',N)
 call dump1(qppm_unl_repl,'qppm_unl_repl',N)
 call dump1(qppm_lin_repl,'qppm_lin_repl',N)
 call dump1(qppm_cwl_repl,'qppm_cwl_repl',N)
 call dump1(qppm_csl_repl,'qppm_csl_repl',N)
 call dump1(q3rd_repl,'q3rd_repl',N)
 call dump1(q3rdul_repl,'q3rdul_repl',N)
 call dump1(q3rdrl_repl,'q3rdrl_repl',N)
 call dump1(qslice_repl,'qslice_repl',N)
 call dump1(qslicebs_repl,'qslicebs_repl',N)

 call dump1(q_tlm,'q_tlm',N)
 call dump1(qppm_lin_tlm,'qppm_unl_tlm',N)
 call dump1(qppm_lin_tlm,'qppm_lin_tlm',N)
 call dump1(qppm_cwl_tlm,'qppm_cwl_tlm',N)
 call dump1(qppm_csl_tlm,'qppm_csl_tlm',N)
 call dump1(q3rd_tlm,'q3rd_tlm',N)
 call dump1(q3rdul_tlm,'q3rdul_tlm',N)
 call dump1(q3rdrl_tlm,'q3rdrl_tlm',N)
 call dump1(qslice_tlm,'qslice_tlm',N)
 call dump1(qslicebs_tlm,'qslicebs_tlm',N)

 call dump2(maxeval_ppm_unl,'maxeval_ppm_unl',N1)
 call dump2(maxeval_ppm_lin,'maxeval_ppm_lin',N1)
 call dump2(maxeval_ppm_cwl,'maxeval_ppm_cwl',N1)
 call dump2(maxeval_ppm_csl,'maxeval_ppm_csl',N1)
 call dump2(maxeval_3rd,'maxeval_3rd',N1)
 call dump2(maxeval_3rdul,'maxeval_3rdul',N1)
 call dump2(maxeval_3rdrl,'maxeval_3rdrl',N1)
 call dump2(maxeval_slice,'maxeval_slice',N1)
 call dump2(maxeval_slicebs,'maxeval_slicebs',N1)

END PROGRAM advection
