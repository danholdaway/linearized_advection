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
DOUBLE PRECISION :: q_tlm_tmp(N)
DOUBLE PRECISION :: q_tlm_init(N)
DOUBLE PRECISION :: u_traj
DOUBLE PRECISION :: u_tlm
DOUBLE PRECISION :: qppm_lin_traj(N), qppmnew_lin_traj(N), qppmhalf_lin_traj(N)
DOUBLE PRECISION :: qppm_lin_tlm(N), qppmnew_lin_tlm(N), qppmhalf_lin_tlm(N)
DOUBLE PRECISION :: qppm_unl_traj(N), qppm_unlnew_traj(N), qppm_unlhalf_traj(N)
DOUBLE PRECISION :: qppm_unl_tlm(N), qppm_unlnew_tlm(N), qppm_unlhalf_tlm(N)
DOUBLE PRECISION :: qppm_cwl_traj(N), qppm_cwlnew_traj(N), qppm_cwlhalf_traj(N)
DOUBLE PRECISION :: qppm_cwl_tlm(N), qppm_cwlnew_tlm(N), qppm_cwlhalf_tlm(N)
DOUBLE PRECISION :: qppm_csl_traj(N), qppm_cslnew_traj(N), qppm_cslhalf_traj(N)
DOUBLE PRECISION :: qppm_csl_tlm(N), qppm_cslnew_tlm(N), qppm_cslhalf_tlm(N)
DOUBLE PRECISION :: q1st_traj(N), q1stnew_traj(N), q1sthalf_traj(N)
DOUBLE PRECISION :: q1st_tlm(N), q1stnew_tlm(N), q1sthalf_tlm(N)
DOUBLE PRECISION :: q2nd_traj(N), q2ndnew_traj(N), q2ndhalf_traj(N)
DOUBLE PRECISION :: q2nd_tlm(N), q2ndnew_tlm(N), q2ndhalf_tlm(N)
DOUBLE PRECISION :: q3rd_traj(N), q3rdnew_traj(N), q3rdhalf_traj(N)
DOUBLE PRECISION :: q3rd_tlm(N), q3rdnew_tlm(N), q3rdhalf_tlm(N)
DOUBLE PRECISION :: q3rdul_traj(N), q3rdulnew_traj(N), q3rdulhalf_traj(N)
DOUBLE PRECISION :: q3rdul_tlm(N), q3rdulnew_tlm(N), q3rdulhalf_tlm(N)
DOUBLE PRECISION :: q3rdrl_traj(N), q3rdrlnew_traj(N), q3rdrlhalf_traj(N)
DOUBLE PRECISION :: q3rdrl_tlm(N), q3rdrlnew_tlm(N), q3rdrlhalf_tlm(N)
DOUBLE PRECISION :: qsl_traj(N),      qslnew_traj(N)
DOUBLE PRECISION :: qsl_tlm(N),      qslnew_tlm(N)
DOUBLE PRECISION :: qslice_traj(N),   qslicenew_traj(N)
DOUBLE PRECISION :: qslice_tlm(N),   qslicenew_tlm(N)
DOUBLE PRECISION :: qslicebs_traj(N), qslicebsnew_traj(N)
DOUBLE PRECISION :: qslicebs_tlm(N), qslicebsnew_tlm(N)
DOUBLE PRECISION :: qspec_traj(N), qspecold_traj(N)
DOUBLE PRECISION :: qspec_tlm(N), qspecold_tlm(N) 

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

 q_free = q_free

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
 u_repl = u_free
 q_repl = q_free

 !call initialconditions(q_repl,x,dx,pi,fx,N,traceval)
 
 !Add a constant to the whole profile (testing only)
 !q_repl = q_repl+1e-4

 !Add a point perturbation
 q_repl(32) = q_repl(32) + 1e-4

 !Multiply by the whole profile to make a positive perturbation
 !q_repl = q_repl*(1+(1-0.9999))

 !Make a nonlinear perturbation by adding between two points
 !q_repl(25:39) = q_repl(25:39)+1e-4

 !Perturb one location
 !q_repl(10) = q_repl(10)+1e-3

! call initialconditions(q_tlm,x,dx,pi,fx,N,3)
! q_tlm = q_tlm * 1e-4 
! q_tlm = q_tlm + 1e-12
! q_repl = q_tlm + q_free


 !call initialconditions(q_tlm,x,dx,pi,fx,N,1)
 !q_tlm = q_tlm * 1e-4
 !q_tlm_tmp = q_tlm
 !q_tlm(1:N-5) = q_tlm_tmp(6:N)
 !q_tlm(N-4:N) = q_tlm_tmp(1:5)
 !q_repl = q_tlm + q_free

 q_tlm_init = q_repl-q_free

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

 !Run Advection
 !-------------
 DO tt=1, N1
   
    !Compute time step for spectral method
    if (tt .eq. 1) then
       twodt=dt
    else
       twodt=2*dt
    endif

    !CALL THE (NON)LINEAR ADVECTION ROUTINES FOR THE FREE RUN
    !--------------------------------------------------------
    !PPM with Lin limiter
    call ppmlin1d(qppmhalf_lin_free,qppm_lin_free,u_free,dt,dx,N)
    qppmnew_lin_free(fx) = qppm_lin_free(fx) - (u_free*dt/dx)*(qppmhalf_lin_free(f1x)-qppmhalf_lin_free(fx))
    qppm_lin_free=qppmnew_lin_free

    !PPM without Limiter
    call unlimitedppm1d(qppm_unlhalf_free,qppm_unl_free,u_free,dt,dx,N)
    qppm_unlnew_free(fx) = qppm_unl_free(fx) - (u_free*dt/dx)*(qppm_unlhalf_free(f1x)-qppm_unlhalf_free(fx))
    qppm_unl_free=qppm_unlnew_free

    !PPM with Colella-Woodward (JCP 1984) limiter
    call ppmoriginal1d(qppm_cwlhalf_free,qppm_cwl_free,u_free,dt,dx,N)
    qppm_cwlnew_free(fx) = qppm_cwl_free(fx) - (u_free*dt/dx)*(qppm_cwlhalf_free(f1x)-qppm_cwlhalf_free(fx))
    qppm_cwl_free=qppm_cwlnew_free

    !PPM with Colella-Sekora (JCP 2008) limiter
    call ppmcs1d(qppm_cslhalf_free,qppm_csl_free,u_free,dt,dx,N)
    qppm_cslnew_free(fx) = qppm_csl_free(fx) - (u_free*dt/dx)*(qppm_cslhalf_free(f1x)-qppm_cslhalf_free(fx))
    qppm_csl_free=qppm_cslnew_free
    
    !FD/FV unlimited 1st order upwind
    call getqhalf(q1sthalf_free,q1st_free,1,N,C,fx,f1x,g1x,g2x,g3x)
    q1stnew_free(fx) = q1st_free(fx) - (u_free*dt/dx)*(q1sthalf_free(f1x)-q1sthalf_free(fx))
    q1st_free=q1stnew_free

    !FD/FV unlimited 2nd order Lax-Wendroff
    call getqhalf(q2ndhalf_free,q2nd_free,2,N,C,fx,f1x,g1x,g2x,g3x)
    q2ndnew_free(fx) = q2nd_free(fx) - (u_free*dt/dx)*(q2ndhalf_free(f1x)-q2ndhalf_free(fx))
    q2nd_free=q2ndnew_free

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

    !Simple SL Scheme
    call simplesemilagrangianscheme(qslnew_free,qsl_free,C,N)
    qsl_free=qslnew_free

    !SLICE Scheme
    call slicescheme(qslicenew_free,qslice_free,C,x,dx,N)
    qslice_free=qslicenew_free

    !SLICE Scheme with BS limiter
    call sliceschemebs(qslicebsnew_free,qslicebs_free,C,x,dx,N)
    qslicebs_free=qslicebsnew_free

    !Spectral method
    call spectraladvection(qspec_free,qspecold_free,u_free,N,x,twodt)


    !CALL THE (NON)LINEAR ADVECTION ROUTINES FOR THE REPLAY RUN
    !----------------------------------------------------------
    !PPM with limiter
    call ppmlin1d(qppmhalf_lin_repl,qppm_lin_repl,u_repl,dt,dx,N)
    qppmnew_lin_repl(fx) = qppm_lin_repl(fx) - (u_repl*dt/dx)*(qppmhalf_lin_repl(f1x)-qppmhalf_lin_repl(fx))
    qppm_lin_repl=qppmnew_lin_repl

    !PPM without Limiter
    call unlimitedppm1d(qppm_unlhalf_repl,qppm_unl_repl,u_repl,dt,dx,N)
    qppm_unlnew_repl(fx) = qppm_unl_repl(fx) - (u_repl*dt/dx)*(qppm_unlhalf_repl(f1x)-qppm_unlhalf_repl(fx))
    qppm_unl_repl=qppm_unlnew_repl

    !PPM with Colella-Woodward (JCP 1984) limiter
    call ppmoriginal1d(qppm_cwlhalf_repl,qppm_cwl_repl,u_repl,dt,dx,N)
    qppm_cwlnew_repl(fx) = qppm_cwl_repl(fx) - (u_repl*dt/dx)*(qppm_cwlhalf_repl(f1x)-qppm_cwlhalf_repl(fx))
    qppm_cwl_repl=qppm_cwlnew_repl

    !PPM with Colella-Sekora (JCP 2008) limiter
    call ppmcs1d(qppm_cslhalf_repl,qppm_csl_repl,u_repl,dt,dx,N)
    qppm_cslnew_repl(fx) = qppm_csl_repl(fx) - (u_repl*dt/dx)*(qppm_cslhalf_repl(f1x)-qppm_cslhalf_repl(fx))
    qppm_csl_repl=qppm_cslnew_repl

    !FD/FV unlimited 1st order upwind
    call getqhalf(q1sthalf_repl,q1st_repl,1,N,C,fx,f1x,g1x,g2x,g3x)
    q1stnew_repl(fx) = q1st_repl(fx) - (u_repl*dt/dx)*(q1sthalf_repl(f1x)-q1sthalf_repl(fx))
    q1st_repl=q1stnew_repl

    !FD/FV unlimited 2nd order Lax-Wendroff
    call getqhalf(q2ndhalf_repl,q2nd_repl,2,N,C,fx,f1x,g1x,g2x,g3x)
    q2ndnew_repl(fx) = q2nd_repl(fx) - (u_repl*dt/dx)*(q2ndhalf_repl(f1x)-q2ndhalf_repl(fx))
    q2nd_repl=q2ndnew_repl

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

    !Simple SL Scheme
    call simplesemilagrangianscheme(qslnew_repl,qsl_repl,C,N)
    qsl_repl=qslnew_repl

    !SLICE Scheme
    call slicescheme(qslicenew_repl,qslice_repl,C,x,dx,N)
    qslice_repl=qslicenew_repl

    !SLICE Scheme with BS limiter
    call sliceschemebs(qslicebsnew_repl,qslicebs_repl,C,x,dx,N)
    qslicebs_repl=qslicebsnew_repl

    !Spectral method
    call spectraladvection(qspec_repl,qspecold_repl,u_repl,N,x,twodt)


    !CALL THE TANGENT LINEAR ADVECTION ROUTINES
    !------------------------------------------
    !PPM with limiter
    call ppmlin1d_tlm(qppmhalf_lin_traj, qppmhalf_lin_tlm, qppm_lin_traj, qppm_lin_tlm, u_traj, u_tlm, dt, dx, N)
    qppmnew_lin_tlm(fx)  = qppm_lin_tlm(fx)  - (u_traj*dt/dx)*(qppmhalf_lin_tlm (f1x)-qppmhalf_lin_tlm (fx)) &
                                             - (u_tlm *dt/dx)*(qppmhalf_lin_traj(f1x)-qppmhalf_lin_traj(fx))
    qppmnew_lin_traj(fx) = qppm_lin_traj(fx) - (u_repl*dt/dx)*(qppmhalf_lin_traj(f1x)-qppmhalf_lin_traj(fx))
    qppm_lin_tlm=qppmnew_lin_tlm
    qppm_lin_traj=qppmnew_lin_traj

    !PPM without Limiter
    call unlimitedppm1d_tlm(qppm_unlhalf_traj, qppm_unlhalf_tlm, qppm_unl_traj, qppm_unl_tlm, u_traj, u_tlm, dt, dx, N)
    qppm_unlnew_tlm(fx)  = qppm_unl_tlm(fx)  - (u_traj*dt/dx)*(qppm_unlhalf_tlm (f1x)-qppm_unlhalf_tlm (fx)) &
                                             - (u_tlm *dt/dx)*(qppm_unlhalf_traj(f1x)-qppm_unlhalf_traj(fx))
    qppm_unlnew_traj(fx)  = qppm_unl_traj(fx)  - (u_traj*dt/dx)*(qppm_unlhalf_traj(f1x)-qppm_unlhalf_traj(fx))
    qppm_unl_tlm=qppm_unlnew_tlm
    qppm_unl_traj=qppm_unlnew_traj

    !PPM with Colella-Woodward (JCP 1984) limiter
    call ppmoriginal1d_tlm(qppm_cwlhalf_traj, qppm_cwlhalf_tlm, qppm_cwl_traj, qppm_cwl_tlm, u_traj, u_tlm, dt, dx, N)
    qppm_cwlnew_tlm(fx)  = qppm_cwl_tlm(fx)  - (u_traj*dt/dx)*(qppm_cwlhalf_tlm (f1x)-qppm_cwlhalf_tlm (fx)) &
                                             - (u_tlm *dt/dx)*(qppm_cwlhalf_traj(f1x)-qppm_cwlhalf_traj(fx))
    qppm_cwlnew_traj(fx)  = qppm_cwl_traj(fx)  - (u_traj*dt/dx)*(qppm_cwlhalf_traj(f1x)-qppm_cwlhalf_traj(fx)) 
    qppm_cwl_tlm=qppm_cwlnew_tlm
    qppm_cwl_traj=qppm_cwlnew_traj

    !PPM with Colella-Sekora (JCP 2008) limiter
    call ppmcs1d_tlm(qppm_cslhalf_traj, qppm_cslhalf_tlm, qppm_csl_traj, qppm_csl_tlm, u_traj, u_tlm, dt, dx, N)
    qppm_cslnew_tlm(fx)  = qppm_csl_tlm(fx)  - (u_traj*dt/dx)*(qppm_cslhalf_tlm (f1x)-qppm_cslhalf_tlm (fx)) &
                                             - (u_tlm *dt/dx)*(qppm_cslhalf_traj(f1x)-qppm_cslhalf_traj(fx))
    qppm_cslnew_traj(fx)  = qppm_csl_traj(fx)  - (u_traj*dt/dx)*(qppm_cslhalf_traj(f1x)-qppm_cslhalf_traj(fx)) 
    qppm_csl_tlm=qppm_cslnew_tlm
    qppm_csl_traj=qppm_cslnew_traj


    !FD/FV unlimited 1st order upwind
    call getqhalf_tlm(q1sthalf_traj,q1sthalf_tlm,q1st_traj,q1st_tlm,1,N,C,fx,f1x,g1x,g2x,g3x)
    q1stnew_tlm(fx) = q1st_tlm(fx) - (u_traj*dt/dx)*(q1sthalf_tlm (f1x)-q1sthalf_tlm (fx)) &
                                   - (u_tlm *dt/dx)*(q1sthalf_traj(f1x)-q1sthalf_traj(fx))
    q1stnew_traj(fx) = q1st_traj(fx) - (u_traj*dt/dx)*(q1sthalf_traj(f1x)-q1sthalf_traj(fx)) 
    q1st_tlm=q1stnew_tlm
    q1st_traj=q1stnew_traj

    !FD/FV unlimited 2nd order Lax-Wendroff
    call getqhalf_tlm(q2ndhalf_traj,q2ndhalf_tlm,q2nd_traj,q2nd_tlm,2,N,C,fx,f1x,g1x,g2x,g3x)
    q2ndnew_tlm(fx) = q2nd_tlm(fx) - (u_traj*dt/dx)*(q2ndhalf_tlm (f1x)-q2ndhalf_tlm (fx)) &
                                   - (u_tlm *dt/dx)*(q2ndhalf_traj(f1x)-q2ndhalf_traj(fx))
    q2ndnew_traj(fx) = q2nd_traj(fx) - (u_traj*dt/dx)*(q2ndhalf_traj(f1x)-q2ndhalf_traj(fx)) 
    q2nd_tlm=q2ndnew_tlm
    q2nd_traj=q2ndnew_traj

    !FD/FV unlimited 3rd order upwind
    call getqhalf_tlm(q3rdhalf_traj,q3rdhalf_tlm,q3rd_traj,q3rd_tlm,3,N,C,fx,f1x,g1x,g2x,g3x)
    q3rdnew_tlm(fx) = q3rd_tlm(fx) - (u_traj*dt/dx)*(q3rdhalf_tlm (f1x)-q3rdhalf_tlm (fx)) &
                                   - (u_tlm *dt/dx)*(q3rdhalf_traj(f1x)-q3rdhalf_traj(fx))
    q3rdnew_traj(fx) = q3rd_traj(fx) - (u_traj*dt/dx)*(q3rdhalf_traj(f1x)-q3rdhalf_traj(fx)) 
    q3rd_tlm=q3rdnew_tlm
    q3rd_traj=q3rdnew_traj

    !FD/FV universal limiter 3rd order upwind
    call getqhalf_tlm(q3rdulhalf_traj,q3rdulhalf_tlm,q3rdul_traj,q3rdul_tlm,3,N,C,fx,f1x,g1x,g2x,g3x)
    call universallimiter_tlm(q3rdulhalf_traj,q3rdulhalf_tlm,q3rdul_traj,q3rdul_tlm,u_traj,u_tlm,dt,dx,N)
    q3rdulnew_tlm (fx) = q3rdul_tlm (fx) - (u_traj*dt/dx)*(q3rdulhalf_tlm (f1x)-q3rdulhalf_tlm (fx)) &
                                         - (u_tlm *dt/dx)*(q3rdulhalf_traj(f1x)-q3rdulhalf_traj(fx))
    q3rdulnew_traj(fx) = q3rdul_traj(fx) - (u_traj*dt/dx)*(q3rdulhalf_traj(f1x)-q3rdulhalf_traj(fx))
    q3rdul_tlm =q3rdulnew_tlm
    q3rdul_traj=q3rdulnew_traj

    !FD/FV relaxed limiter 3rd order upwind
    call getqhalf_tlm(q3rdrlhalf_traj,q3rdrlhalf_tlm,q3rdrl_traj,q3rdrl_tlm,3,N,C,fx,f1x,g1x,g2x,g3x)
    call relaxedlimiter_tlm(q3rdrlhalf_traj,q3rdrlhalf_tlm,q3rdrl_traj,q3rdrl_tlm,u_traj,u_tlm,dt,dx,N)
    q3rdrlnew_tlm (fx) = q3rdrl_tlm (fx) - (u_traj*dt/dx)*(q3rdrlhalf_tlm (f1x)-q3rdrlhalf_tlm (fx)) &
                                         - (u_tlm *dt/dx)*(q3rdrlhalf_traj(f1x)-q3rdrlhalf_traj(fx))
    q3rdrlnew_traj(fx) = q3rdrl_traj(fx) - (u_traj*dt/dx)*(q3rdrlhalf_traj(f1x)-q3rdrlhalf_traj(fx))
    q3rdrl_tlm =q3rdrlnew_tlm
    q3rdrl_traj=q3rdrlnew_traj

    !Simple SL Scheme
    call simplesemilagrangianscheme_tlm(qslnew_traj,qslnew_tlm,qsl_traj,qsl_tlm,C,N)
    qsl_traj=qslnew_traj
    qsl_tlm=qslnew_tlm

    !SLICE Scheme
    call slicescheme_tlm(qslicenew_traj,qslicenew_tlm,qslice_traj,qslice_tlm,C,x,dx,N)
    qslice_traj=qslicenew_traj
    qslice_tlm=qslicenew_tlm

    !SLICE with BS Scheme
    call sliceschemebs_tlm(qslicebsnew_traj,qslicebsnew_tlm,qslicebs_traj,qslicebs_tlm,C,x,dx,N)
    qslicebs_traj=qslicebsnew_traj
    qslicebs_tlm=qslicebsnew_tlm

    !Spectral method
    call spectraladvection_tlm(qspec_traj,qspec_tlm,qspecold_traj,qspecold_tlm,u_traj,u_tlm,N,x,twodt)

 ENDDO


 !Compute difference between nonlinear runs at the endtime
 !--------------------------------------------------------
 nlm_pert_traj_qppm     = qppm_lin_repl - qppm_lin_free
 nlm_pert_traj_qppm_unl = qppm_unl_repl - qppm_unl_free
 nlm_pert_traj_qppm_cwl = qppm_cwl_repl - qppm_cwl_free
 nlm_pert_traj_qppm_csl = qppm_csl_repl - qppm_csl_free
 nlm_pert_traj_q1st     = q1st_repl     - q1st_free
 nlm_pert_traj_q2nd     = q2nd_repl     - q2nd_free
 nlm_pert_traj_q3rd     = q3rd_repl     - q3rd_free
 nlm_pert_traj_q3rdul   = q3rdul_repl   - q3rdul_free
 nlm_pert_traj_q3rdrl   = q3rdrl_repl   - q3rdrl_free
 nlm_pert_traj_qsl      = qsl_repl      - qsl_free
 nlm_pert_traj_qslice   = qslice_repl   - qslice_free
 nlm_pert_traj_qslicebs = qslicebs_repl - qslicebs_free
 nlm_pert_traj_qspec    = qspec_repl    - qspec_free


 !Print out measures of the error between the nonlinear and linear perturbation trajectories
 !------------------------------------------------------------------------------------------
 print*, ''
 print*, '==================================================='
 print*, 'STANDARD TLM GRADIENT TESTS'
 print*, '---------------------------'
 print*, 'TLM TEST 1st Order FD        :', sum(nlm_pert_traj_q1st/q1st_tlm)/N
 print*, 'TLM TEST 2nd Order FD        :', sum(nlm_pert_traj_q2nd/q2nd_tlm)/N
 print*, 'TLM TEST 3rd Order FD        :', sum(nlm_pert_traj_q3rd/q3rd_tlm)/N
 print*, 'TLM TEST PPM Unlimited       :', sum(nlm_pert_traj_qppm_unl/qppm_unl_tlm)/N
 print*, 'TLM TEST PPM Lin RMC etc     :', sum(nlm_pert_traj_qppm/qppm_lin_tlm)/N
 print*, 'TLM TEST PPM with CW Lim     :', sum(nlm_pert_traj_qppm_cwl/qppm_cwl_tlm)/N
 print*, 'TLM TEST PPM with CS Lim     :', sum(nlm_pert_traj_qppm_csl/qppm_csl_tlm)/N
 print*, 'TLM TEST 3rd Relaxed Lim     :', sum(nlm_pert_traj_q3rdrl/q3rdrl_tlm)/N
 print*, 'TLM TEST 3rd Universal Lim   :', sum(nlm_pert_traj_q3rdul/q3rdul_tlm)/N
 print*, 'TLM TEST Semi-Lag (SIMP)     :', sum(nlm_pert_traj_qsl/qsl_tlm)/N
 print*, 'TLM TEST Semi-Lag (SLICE)    :', sum(nlm_pert_traj_qslice/qslice_tlm)/N
 print*, 'TLM TEST Semi-Lag (SLICE BS) :', sum(nlm_pert_traj_qslicebs/qslicebs_tlm)/N
 print*, 'TLM TEST Spectral method     :', sum(nlm_pert_traj_qspec/qspec_tlm)/N
 print*, '==================================================='
 print*, ''

 print*, '======================================================'
 print*, 'CORRELATION BETWEEN NONLINEAR AND LINEAR'
 print*, '----------------------------------------'
 print*, 'CORRELATION 1st Order FD        :', sum(nlm_pert_traj_q1st*q1st_tlm)/sqrt(sum(nlm_pert_traj_q1st**2)*sum(q1st_tlm**2))
 print*, 'CORRELATION 2nd Order FD        :', sum(nlm_pert_traj_q2nd*q2nd_tlm)/sqrt(sum(nlm_pert_traj_q2nd**2)*sum(q2nd_tlm**2))
 print*, 'CORRELATION 3rd Order FD        :', sum(nlm_pert_traj_q3rd*q3rd_tlm)/sqrt(sum(nlm_pert_traj_q3rd**2)*sum(q3rd_tlm**2))
 print*, 'CORRELATION PPM Unlimited       :', sum(nlm_pert_traj_qppm_unl*qppm_unl_tlm)&
                                              /sqrt(sum(nlm_pert_traj_qppm_unl**2)*sum(qppm_unl_tlm**2))
 print*, 'CORRELATION PPM Lin RMC etc     :', sum(nlm_pert_traj_qppm*qppm_lin_tlm)&
                                              /sqrt(sum(nlm_pert_traj_qppm**2)*sum(qppm_lin_tlm**2))
 print*, 'CORRELATION PPM with CW Lim     :', sum(nlm_pert_traj_qppm_cwl*qppm_cwl_tlm)&
                                              /sqrt(sum(nlm_pert_traj_qppm_cwl**2)*sum(qppm_cwl_tlm**2))
 print*, 'CORRELATION PPM with CS Lim     :', sum(nlm_pert_traj_qppm_csl*qppm_csl_tlm)&
                                              /sqrt(sum(nlm_pert_traj_qppm_csl**2)*sum(qppm_csl_tlm**2))
 print*, 'CORRELATION 3rd Relaxed Lim     :', sum(nlm_pert_traj_q3rdrl*q3rdrl_tlm)&
                                              /sqrt(sum(nlm_pert_traj_q3rdrl**2)*sum(q3rdrl_tlm**2))
 print*, 'CORRELATION 3rd Universal Lim   :', sum(nlm_pert_traj_q3rdul*q3rdul_tlm)&
                                              /sqrt(sum(nlm_pert_traj_q3rdul**2)*sum(q3rdul_tlm**2))
 print*, 'CORRELATION Semi-Lag (SIMP)     :', sum(nlm_pert_traj_qsl*qsl_tlm)&
                                              /sqrt(sum(nlm_pert_traj_qsl**2)*sum(qsl_tlm**2))
 print*, 'CORRELATION Semi-Lag (SLICE)    :', sum(nlm_pert_traj_qslice*qslice_tlm)&
                                              /sqrt(sum(nlm_pert_traj_qslice**2)*sum(qslice_tlm**2))
 print*, 'CORRELATION Semi-Lag (SLICE BS) :', sum(nlm_pert_traj_qslicebs*qslicebs_tlm)&
                                              /sqrt(sum(nlm_pert_traj_qslicebs**2)*sum(qslicebs_tlm**2))
 print*, 'CORRELATION Spectral method     :', sum(nlm_pert_traj_qspec*qspec_tlm)&
                                              /sqrt(sum(nlm_pert_traj_qspec**2)*sum(qspec_tlm**2))
 print*, '======================================================'

 print*, ''
 print*, '===================================================='
 print*, 'MEAN ERROR (L1-norm)'
 print*, '--------------------'
! print*, 'RMSE 1st Order FD        :', (1./N)*sum(abs(nlm_pert_traj_q1st-q1st_tlm))
! print*, 'RMSE 2nd Order FD        :', (1./N)*sum(abs(nlm_pert_traj_q2nd-q2nd_tlm))
! print*, 'RMSE 3rd Order FD        :', (1./N)*sum(abs(nlm_pert_traj_q3rd-q3rd_tlm))
! print*, 'RMSE PPM Unlimited       :', (1./N)*sum(abs(nlm_pert_traj_qppm_unl-qppm_unl_tlm))
! print*, 'RMSE PPM Lin RMC etc     :', (1./N)*sum(abs(nlm_pert_traj_qppm-qppm_lin_tlm))
! print*, 'RMSE PPM with CW Lim     :', (1./N)*sum(abs(nlm_pert_traj_qppm_cwl-qppm_cwl_tlm))
! print*, 'RMSE PPM with CS Lim     :', (1./N)*sum(abs(nlm_pert_traj_qppm_csl-qppm_csl_tlm))
! print*, 'RMSE 3rd Relaxed Lim     :', (1./N)*sum(abs(nlm_pert_traj_q3rdrl-q3rdrl_tlm))
! print*, 'RMSE 3rd Universal Lim   :', (1./N)*sum(abs(nlm_pert_traj_q3rdul-q3rdul_tlm))
! print*, 'RMSE Semi-Lag (SIMP)     :', (1./N)*sum(abs(nlm_pert_traj_qsl-qsl_tlm))
! print*, 'RMSE Semi-Lag (SLICE)    :', (1./N)*sum(abs(nlm_pert_traj_qslice-qslice_tlm))
! print*, 'RMSE Semi-Lag (SLICE BS) :', (1./N)*sum(abs(nlm_pert_traj_qslicebs-qslicebs_tlm))
! print*, 'RMSE Spectral method     :', (1./N)*sum(abs(nlm_pert_traj_qspec-qspec_tlm))
 print*, 'RMSE 1st Order FD        :', (1./N)*sum(abs(q_tlm_init-q1st_tlm))
 print*, 'RMSE 2nd Order FD        :', (1./N)*sum(abs(q_tlm_init-q2nd_tlm))
 print*, 'RMSE 3rd Order FD        :', (1./N)*sum(abs(q_tlm_init-q3rd_tlm))
 print*, 'RMSE PPM Unlimited       :', (1./N)*sum(abs(q_tlm_init-qppm_unl_tlm))
 print*, 'RMSE PPM Lin RMC etc     :', (1./N)*sum(abs(q_tlm_init-qppm_lin_tlm))
 print*, 'RMSE PPM with CW Lim     :', (1./N)*sum(abs(q_tlm_init-qppm_cwl_tlm))
 print*, 'RMSE PPM with CS Lim     :', (1./N)*sum(abs(q_tlm_init-qppm_csl_tlm))
 print*, 'RMSE 3rd Relaxed Lim     :', (1./N)*sum(abs(q_tlm_init-q3rdrl_tlm))
 print*, 'RMSE 3rd Universal Lim   :', (1./N)*sum(abs(q_tlm_init-q3rdul_tlm))
 print*, 'RMSE Semi-Lag (SIMP)     :', (1./N)*sum(abs(q_tlm_init-qsl_tlm))
 print*, 'RMSE Semi-Lag (SLICE)    :', (1./N)*sum(abs(q_tlm_init-qslice_tlm))
 print*, 'RMSE Semi-Lag (SLICE BS) :', (1./N)*sum(abs(q_tlm_init-qslicebs_tlm))
 print*, 'RMSE Spectral method     :', (1./N)*sum(abs(q_tlm_init-qspec_tlm))
 print*, '===================================================='
 print*, ''

 print*, '===================================================='
 print*, 'ROOT MEAN SQUARED ERROR (L2-norm)'
 print*, '---------------------------------'
! print*, 'RMSE 1st Order FD        :', sqrt((1./N)*sum((nlm_pert_traj_q1st-q1st_tlm)**2))
! print*, 'RMSE 2nd Order FD        :', sqrt((1./N)*sum((nlm_pert_traj_q2nd-q2nd_tlm)**2))
! print*, 'RMSE 3rd Order FD        :', sqrt((1./N)*sum((nlm_pert_traj_q3rd-q3rd_tlm)**2))
! print*, 'RMSE PPM Unlimited       :', sqrt((1./N)*sum((nlm_pert_traj_qppm_unl-qppm_unl_tlm)**2))
! print*, 'RMSE PPM Lin RMC etc     :', sqrt((1./N)*sum((nlm_pert_traj_qppm-qppm_lin_tlm)**2))
! print*, 'RMSE PPM with CW Lim     :', sqrt((1./N)*sum((nlm_pert_traj_qppm_cwl-qppm_cwl_tlm)**2))
! print*, 'RMSE PPM with CS Lim     :', sqrt((1./N)*sum((nlm_pert_traj_qppm_csl-qppm_csl_tlm)**2))
! print*, 'RMSE 3rd Relaxed Lim     :', sqrt((1./N)*sum((nlm_pert_traj_q3rdrl-q3rdrl_tlm)**2))
! print*, 'RMSE 3rd Universal Lim   :', sqrt((1./N)*sum((nlm_pert_traj_q3rdul-q3rdul_tlm)**2))
! print*, 'RMSE Semi-Lag (SIMP)     :', sqrt((1./N)*sum((nlm_pert_traj_qsl-qsl_tlm)**2))
! print*, 'RMSE Semi-Lag (SLICE)    :', sqrt((1./N)*sum((nlm_pert_traj_qslice-qslice_tlm)**2))
! print*, 'RMSE Semi-Lag (SLICE BS) :', sqrt((1./N)*sum((nlm_pert_traj_qslicebs-qslicebs_tlm)**2))
! print*, 'RMSE Spectral method     :', sqrt((1./N)*sum((nlm_pert_traj_qspec-qspec_tlm)**2))
 print*, 'RMSE 1st Order FD        :', sqrt((1./N)*sum((q_tlm_init-q1st_tlm)**2))
 print*, 'RMSE 2nd Order FD        :', sqrt((1./N)*sum((q_tlm_init-q2nd_tlm)**2))
 print*, 'RMSE 3rd Order FD        :', sqrt((1./N)*sum((q_tlm_init-q3rd_tlm)**2))
 print*, 'RMSE PPM Unlimited       :', sqrt((1./N)*sum((q_tlm_init-qppm_unl_tlm)**2))
 print*, 'RMSE PPM Lin RMC etc     :', sqrt((1./N)*sum((q_tlm_init-qppm_lin_tlm)**2))
 print*, 'RMSE PPM with CW Lim     :', sqrt((1./N)*sum((q_tlm_init-qppm_cwl_tlm)**2))
 print*, 'RMSE PPM with CS Lim     :', sqrt((1./N)*sum((q_tlm_init-qppm_csl_tlm)**2))
 print*, 'RMSE 3rd Relaxed Lim     :', sqrt((1./N)*sum((q_tlm_init-q3rdrl_tlm)**2))
 print*, 'RMSE 3rd Universal Lim   :', sqrt((1./N)*sum((q_tlm_init-q3rdul_tlm)**2))
 print*, 'RMSE Semi-Lag (SIMP)     :', sqrt((1./N)*sum((q_tlm_init-qsl_tlm)**2))
 print*, 'RMSE Semi-Lag (SLICE)    :', sqrt((1./N)*sum((q_tlm_init-qslice_tlm)**2))
 print*, 'RMSE Semi-Lag (SLICE BS) :', sqrt((1./N)*sum((q_tlm_init-qslicebs_tlm)**2))
 print*, 'RMSE Spectral method     :', sqrt((1./N)*sum((q_tlm_init-qspec_tlm)**2))
 print*, '===================================================='
 print*, ''


 print*, '===================================================='
 print*, 'ROOT MEAN SQUARED ERROR OF GRAD TEST (L2-norm)'
 print*, '---------------------------------'
 print*, 'RMSE 1st Order FD        :', sqrt((1./N)*sum((nlm_pert_traj_q1st-q1st_tlm)**2))
 print*, 'RMSE 2nd Order FD        :', sqrt((1./N)*sum((nlm_pert_traj_q2nd-q2nd_tlm)**2))
 print*, 'RMSE 3rd Order FD        :', sqrt((1./N)*sum((nlm_pert_traj_q3rd-q3rd_tlm)**2))
 print*, 'RMSE PPM Unlimited       :', sqrt((1./N)*sum((nlm_pert_traj_qppm_unl-qppm_unl_tlm)**2))
 print*, 'RMSE PPM Lin RMC etc     :', sqrt((1./N)*sum((nlm_pert_traj_qppm-qppm_lin_tlm)**2))/&
(maxval(qppm_lin_tlm)-minval(qppm_lin_tlm))
 print*, 'RMSE PPM with CW Lim     :', sqrt((1./N)*sum((nlm_pert_traj_qppm_cwl-qppm_cwl_tlm)**2))
 print*, 'RMSE PPM with CS Lim     :', sqrt((1./N)*sum((nlm_pert_traj_qppm_csl-qppm_csl_tlm)**2))
 print*, 'RMSE 3rd Relaxed Lim     :', sqrt((1./N)*sum((nlm_pert_traj_q3rdrl-q3rdrl_tlm)**2))
 print*, 'RMSE 3rd Universal Lim   :', sqrt((1./N)*sum((nlm_pert_traj_q3rdul-q3rdul_tlm)**2))
 print*, 'RMSE Semi-Lag (SIMP)     :', sqrt((1./N)*sum((nlm_pert_traj_qsl-qsl_tlm)**2))
 print*, 'RMSE Semi-Lag (SLICE)    :', sqrt((1./N)*sum((nlm_pert_traj_qslice-qslice_tlm)**2))
 print*, 'RMSE Semi-Lag (SLICE BS) :', sqrt((1./N)*sum((nlm_pert_traj_qslicebs-qslicebs_tlm)**2))
 print*, 'RMSE Spectral method     :', sqrt((1./N)*sum((nlm_pert_traj_qspec-qspec_tlm)**2))
 print*, '===================================================='
 print*, ''


 print*, '===================================================='
 print*, 'MAX ERROR (Linf-norm)'
 print*, '---------------------'
 print*, 'RMSE 1st Order FD        :', maxval(nlm_pert_traj_q1st-q1st_tlm)
 print*, 'RMSE 2nd Order FD        :', maxval(nlm_pert_traj_q2nd-q2nd_tlm)
 print*, 'RMSE 3rd Order FD        :', maxval(nlm_pert_traj_q3rd-q3rd_tlm)
 print*, 'RMSE PPM Unlimited       :', maxval(nlm_pert_traj_qppm_unl-qppm_unl_tlm)
 print*, 'RMSE PPM Lin RMC etc     :', maxval(nlm_pert_traj_qppm-qppm_lin_tlm)
 print*, 'RMSE PPM with CW Lim     :', maxval(nlm_pert_traj_qppm_cwl-qppm_cwl_tlm)
 print*, 'RMSE PPM with CS Lim     :', maxval(nlm_pert_traj_qppm_csl-qppm_csl_tlm)
 print*, 'RMSE 3rd Relaxed Lim     :', maxval(nlm_pert_traj_q3rdrl-q3rdrl_tlm)
 print*, 'RMSE 3rd Universal Lim   :', maxval(nlm_pert_traj_q3rdul-q3rdul_tlm)
 print*, 'RMSE Semi-Lag (SIMP)     :', maxval(nlm_pert_traj_qsl-qsl_tlm)
 print*, 'RMSE Semi-Lag (SLICE)    :', maxval(nlm_pert_traj_qslice-qslice_tlm)
 print*, 'RMSE Semi-Lag (SLICE BS) :', maxval(nlm_pert_traj_qslicebs-qslicebs_tlm)
 print*, 'RMSE Spectral method     :', maxval(nlm_pert_traj_qspec-qspec_tlm)
 print*, '===================================================='
 print*, ''



 !Output for MATLAB
 !-----------------
 call dump1(q_free,'q_free',N)
 call dump1(qppm_lin_free,'qppm_lin_free',N)
 call dump1(qppm_unl_free,'qppm_unl_free',N)
 call dump1(qppm_cwl_free,'qppm_cwl_free',N)
 call dump1(qppm_csl_free,'qppm_csl_free',N)
 call dump1(q1st_free,'q1st_free',N)
 call dump1(q2nd_free,'q2nd_free',N)
 call dump1(q3rd_free,'q3rd_free',N)
 call dump1(q3rdul_free,'q3rdul_free',N)
 call dump1(q3rdrl_free,'q3rdrl_free',N)
 call dump1(qsl_free,'qsl_free',N)
 call dump1(qslice_free,'qslice_free',N)
 call dump1(qslicebs_free,'qslicebs_free',N)
 call dump1(qspec_free,'qspec_free',N)

 call dump1(q_repl,'q_repl',N)
 call dump1(qppm_lin_repl,'qppm_lin_repl',N)
 call dump1(qppm_unl_repl,'qppm_unl_repl',N)
 call dump1(qppm_cwl_repl,'qppm_cwl_repl',N)
 call dump1(qppm_csl_repl,'qppm_csl_repl',N)
 call dump1(q1st_repl,'q1st_repl',N)
 call dump1(q2nd_repl,'q2nd_repl',N)
 call dump1(q3rd_repl,'q3rd_repl',N)
 call dump1(q3rdul_repl,'q3rdul_repl',N)
 call dump1(q3rdrl_repl,'q3rdrl_repl',N)
 call dump1(qsl_repl,'qsl_repl',N)
 call dump1(qslice_repl,'qslice_repl',N)
 call dump1(qslicebs_repl,'qslicebs_repl',N)
 call dump1(qspec_repl,'qspec_repl',N)

 call dump1(q_tlm,'q_tlm',N)
 call dump1(qppm_lin_tlm,'qppm_lin_tlm',N)
 call dump1(qppm_unl_tlm,'qppm_unl_tlm',N)
 call dump1(qppm_cwl_tlm,'qppm_cwl_tlm',N)
 call dump1(qppm_csl_tlm,'qppm_csl_tlm',N)
 call dump1(q1st_tlm,'q1st_tlm',N)
 call dump1(q2nd_tlm,'q2nd_tlm',N)
 call dump1(q3rd_tlm,'q3rd_tlm',N)
 call dump1(q3rdul_tlm,'q3rdul_tlm',N)
 call dump1(q3rdrl_tlm,'q3rdrl_tlm',N)
 call dump1(qsl_tlm,'qsl_tlm',N)
 call dump1(qslice_tlm,'qslice_tlm',N)
 call dump1(qslicebs_tlm,'qslicebs_tlm',N)
 call dump1(qspec_tlm,'qspec_tlm',N)

END PROGRAM advection



