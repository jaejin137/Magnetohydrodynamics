      program FASIP
      
      use datain

      use newtypes

      implicit none

      include "/usr/local/include/mpif.h"
      
      type (datain_type)::control
c      
c     
c     Static Variables
c     
      integer::
     $     il, jl, kl,          ! maximum cell number in 3D
     $     nl,                  ! equation number
     $     dim1, dim2,          ! maximum cell and face number in domain
     $     lhs_order,           ! 1st order MUSCL for LHS matrices
     $     rhs_order,           ! 1st order MUSCL for RHS
     $     idimen,              ! dimension number
     $     nstep,               ! No. of time steps
     $     ntotal,              ! No. of total steps
     $     vis_order,           ! order of viscous term
     $     limiter,
     $     integrate,           ! integrate method
     $     precondition,        ! precondition
     $     rhs_scheme,          ! rhside scheme
     $     lhs_scheme,
     $     iter_gs,             ! number of Gauss-Seidel Iterations
     $     turb,                ! turbulence model
     $     unidt,               ! using uniform dt in whole domain
     $     nn,                  ! iteration serial number
     $     loc(3),
     $     loc0(4),
     $     rank, np,            ! process number, total pro num
     $     intersteps, checksteps,
     $     moving,              ! 0, fixed, 1, forced moving, 2, self moving
     $     gcl,                 ! =1 goemetric conservation law
     $     strtp,               ! solid struture type
     $     nrbc_ex,             ! non-reflect boundary condition type
     $     recstep,             ! explicit unsteady record index (every tintvl)
     $     nfp,                 ! parameter for induced vibration of airfoil
     $     nmode,               ! total number of mode shapes
     $     main_dir,            ! main stream direction
     $     source,              ! source term
     $     wall_order,          ! wall accuracy
     $     rotide               ! (jaejin) for rotor problem

c
      logical::
     $     inviscid, suther,    ! viscousity and its computation method
     $     done = .False.,
     $     ready = .False.,
     $     pr_check = .False.,      ! save intermediate results
     $     ke = .False.,            !
     $     ex
c
      double precision::
     $     epsfactor,           ! epsilon used in linear reconstruction
     $     kfactor,             ! factor in linear reconstruction
     $     cfl,                 ! CFL number
     $     reynolds,            ! Reynolds number
     $     time,                ! marching time
     $     machinf,             ! Mach no. based on U_infinity and a_infinity
     $     gamma,               ! specific heat ratio
     $     prandtl,             ! Molecular Prandtl number
     $     prt,                 ! Turbulent Prandtl number
     $     tref,                ! Non-dimensional ref temperature in Sutherland
     $     resdulmax,
     $     resdulave,
     $     ptotal,              ! Inlet total pressure under bc 9
     $     ttotal,              ! Inlet total temperature under bc 9
     $     angl1, angl2,        ! Inlet angles
     $     poutlet,             ! Fixed subsonic outlet pressure
     $     resdl,               ! residual for solving structural eqn
     $     sigma(2),            ! coefficients for NRBC
     $     delta(5),            ! Entropy cut-off for rho, rhou, rhov, rhow,
     $                          ! rhoe
     $     tstart, tend,        ! job start time, end time
     $     systime,             ! systime fuction declarition
     $     k_prec,              ! constant for precondition
     $     rayleigh,
     $     froude,
     $     epsilon,
     $     rotalp               ! (jaejin) also for rotor problem
c
      double precision::
     $     iso_tw(10)           ! isothermal wall temperature values
c
      integer, parameter::
     $     snum = 20,           ! oscillating mesh number
     $     tag_turb = 1,
     $     tag_q = 2,
     $     tag_res = 3,
     $     tag_loc = 4,
     $     tag_max = 5,
     $     tag_qt = 6,
     $     tag_x = 7,
     $     tag_y = 8,
     $     tag_z = 9
c
      double precision::
     $     eps, dtemp
c     
c     variables for dual-time-stepping
c
      integer::
     $     dual_t,              ! 0. steady, local/uniform,
     $                          ! 1. unsteady, local,
     $                          ! 2. unsteady, uniform
     $     tsteps,              ! unsteady time step number
     $     dstep,               ! dual time stepping index
     $     dstart

      double precision::
     $     tintvl,              ! unsteady time step
     $     xsl(4,3),
     $     vbd(3,3), xyz(3,3),
     $     cs, mus, ub,         ! parameters for induced vibration of cylinder
     $     kc, wh, walf, uinf,  ! parameters for induced vibration of airfoil
     $     alf0,                ! parameters for induced vibration of airfoil
     $     xzl(2,3,6),          ! structural variables
     $     wpw(6)               ! (wj/wa)

c *** modified by baoyuan in 31 may,2005.
c
c *** the variables of boundary condition
c
      integer:: nb, bc_max, blen
      integer:: bcwake(3,2)
c
c *** used connect with namelist bcdef 
c
c     integer::
c    $     nb,block,iblock,start(3),end(3),istart(3),iend(3),bctype,
c    $     order(3),iorder(3),bc_max,blen,ic(3)

c
c *** the variables of multiblock and connection BC
c
      integer:: ix0,iq0,ivol,ivist,ibcxie,ibceta,ibczta,it0,
     $   iudi,iudj,iudk,iqt,ixs,
     $   ilow,iupp,jlow,jupp,klow,kupp
c
      integer, dimension(:), allocatable::
     $   bnum, ill, jll, kll
c
      double precision, dimension(:), allocatable::
     $   x0, y0, z0, q0, vist0, vol0,
     $   bcxie0, bcxie1, bceta0, bceta1, bczta0, bczta1,
     $   dt0, dist0
      double precision, dimension(:), allocatable::
     $   qi0, qii0, qiii0
c
      double precision, dimension(:), allocatable::
     $   udi0,vdi0,wdi0,udj0,vdj0,wdj0,udk0,vdk0,wdk0,
     $   xo0,yo0,zo0,qt0,xs0,ys0

      double precision, dimension(:,:,:), allocatable::
     $   hx, hy, hz           ! mode shape

      integer, dimension(:,:), allocatable::
     $   btb,bnb,bdir,rorder
      integer, dimension(:,:,:), allocatable::
     $   bstart,bend,iwake
c
      integer, parameter::
     $     fdatain = 15
c
      double precision, dimension(:), allocatable::
     $   ptemp, pttemp, tttemp, xw, yw, zw
      integer, dimension(:), allocatable::
     $   main_d
      integer:: i, j, k, ipos
      character(len=1):: choice
      double precision::
     $   velinit, theta, twpar, tkpar, ronum, varepsilon
      integer:: strm_dir
c
c *** SA 1eq turbulent model
c
      double precision:: 
     $     tko, cb1, cb2, cap_k, cw1, cw2, cw3,
     $     cv1, ct1, ct2, ct3, ct4, cdes
      integer:: iblnu, ipt, jpt, kpt, ic1, ic2, ic3, ides
c
c *** forced moving airfoill
c
      double precision:: lref, aref, tcl, tcd, tcm
      double precision:: xctr, yctr, xtor, ytor
      double precision:: xref, yref, xref0, yref0 
      double precision:: beta
      integer:: jfx
c
c *** following is statistical values for DES
c
      double precision, dimension(:), allocatable::
     $   rave, uave,vave,wave,eave
      integer :: nstatis  !the statistical times

      double precision, dimension(:), allocatable::
     $   diver
c *** the variables of boundary condition
      character(len=3)::
     $     bcdir

      integer::
     $   block,iblock,start(3),end(3),istart(3),iend(3),bctype,
     $   order(3)

c=================================jaejin(begin)=========================
c       MHD-related variables

        type(journaldata)::
     $  journal_info
        
        character::
     $  journal_caller*20,
     $  ionmodel,
     $  magtype*10,
     $  magref

        integer::
     $  realgas,        ! Whether real gas effect is used
     $  mflag,          ! mflag for TGAS1
     $  callperd,       ! calling period for TGAS1 i.t.o. steps
     $  mhdsim          ! Whether or not MHD simulation is carried out

        double precision::
     $  rho_inf,        ! Fluid density at infinity (in &mhd)
     $  p_inf,          ! Static pressure at infinity (in &mhd)
     $  T_inf,          ! Temperature at infinity (in &mhd)
     $  mu_inf,         ! Viscosity at infinity (in &mhd)
     $  c_sound,        ! Speed of sound (in &mhd)
     $  L_ref,          ! Characteristic length (in &mhd)
     $  T_onset,        ! threshold temperature for onset of nonzero
                        ! electrical conductivity
     $  sigma_e_ref,    ! Reference value of electrical conductivity
     $  B_ref,          ! Reference value for magnetic field
     $  I_in,           ! Reference value for input current (in &mhd)
     $  offset(3)      ! array for offset (see below)
        
        double precision, dimension(:,:,:), allocatable::
     $  sigma_e_star, B_x_star, B_y_star, B_z_star, B_norm

        double precision, dimension(:,:,:,:), allocatable::
     $  source_mhd
c=================================jaejin(end)===========================

c
c *** namelist definition
c
      namelist /consts/ suther, gamma, prandtl, prt, tref,
     $     delta, epsfactor
      namelist /flows/ blen, idimen, nl, machinf, reynolds, ptotal,
     $     ttotal, angl1, angl2, poutlet, tintvl, inviscid, turb


      namelist /comp/ nstep, eps, cfl, limiter, lhs_scheme, rhs_scheme,
     $     lhs_order, rhs_order, integrate, iter_gs, unidt, dual_t,
     $     intersteps, checksteps, tsteps, moving, gcl, precondition,
     $     k_prec, strtp, nrbc_ex, sigma, kfactor, vis_order,
     $     choice, velinit, strm_dir, theta, twpar, tkpar, source,
     $     rotide, varepsilon, wall_order
c
      namelist /rotor/ rotalp, ronum
      
c=================================jaejin(begin)=========================
c       Namelist for MHD Simulation (by Jaejin)
      namelist /mhd/ rho_inf, p_inf, T_inf, mu_inf, c_sound, L_ref,
c     $     realgas, mflag, callperd, mhdsim, sigma_e_ref, magref,
     $     mhdsim, ionmodel, T_onset, sigma_e_ref, magtype, magref,
     $     B_ref, I_in, offset
c=================================jaejin(end)===========================


      namelist /bcdef/ bcdir, block, bctype, start, end, iblock, 
     $     istart, iend, order
c
      namelist /iso_t/ iso_tw
c
      namelist /coef1eq/ iblnu, ipt, jpt, kpt, ic1, ic2, ic3, ides,
     $  tko, cb1, cb2, cap_k, cw2, cw3, cv1, ct1, ct2, ct3, ct4, cdes
c
c *** forced moving airfoil namelist definition
c
      namelist /mvgrid2/ lref, aref, xctr, yctr, xref, yref, jfx, beta
      namelist /strct_cyl/ cs, mus, ub
      namelist /strct_af/ kc, wh, walf, mus, uinf, nfp, alf0
      namelist /source_g/rayleigh, froude, epsilon
c    
c     Begin Program
c
c *** begin spawn the parallel processes
c
      call mpi_initial(np, rank)
c
c *** output headfile
c
      if (rank.eq.0) then
        write(6,*)'**************************************************'
        write(6,*)'*                                                *'
        write(6,*)'*         FASIP, Version 1.0                     *'
        write(6,*)'*  Fluid-Acoustics-Structure Interaction Package *'
        write(6,*)'*           04/22/2008                           *'
        write(6,*)'*                                                *'
        write(6,*)'*  Use of this software must  be authorized by   *'
        write(6,*)'*  Prof. Gecheng Zha  at University of Miami.    *'
        write(6,*)'*       www.eng.miami.edu/acfdlab                *'
        write(6,*)'**************************************************'
        write(6,*)'- ($Id: main.f 245 2008-09-10 15:46:52Z xychen $) -'
      end if
c
c *** read namelist
c
      open(fdatain, file='datain')
      read(fdatain, nml=consts)
      if (rank.eq.0) then
        write(*,*)'namelist consts have been read'
      end if
      read(fdatain, nml=flows)
      if (rank.eq.0) then
        write(*,*)'namelist flows have been read'
      end if
      read(fdatain, nml=comp)
      if (rank.eq.0) then
        write(*,*)'namelist comp have been read'
      end if
c
      if (nl.eq.6) then
        read(fdatain, nml=coef1eq)
        cw1 = cb1/cap_k**2+(1.d0+cb2)/tko
        if (rank.eq.0) then
          write(*,*)'namelist coef1eq have been read'
        end if
      end if
c
      if (dual_t.eq.1) then
        read(fdatain, nml=mvgrid2)
        if (rank.eq.0) then
          write(*,*)'namelist mvgrid2 have been read'
        end if
        xref0 = xref
        yref0 = yref
      end if
c
      if (moving.eq.2) then
        select case (strtp)
        case (1)
          read(fdatain, nml=strct_cyl)
          if (rank.eq.0) then
            write(*,*)'namelist strct_cyl have been read'
          end if
        case (3,5)
          read(fdatain, nml=strct_af)
          if (rank.eq.0) then
            write(*,*)'namelist strct_af have been read'
          end if
        end select
      end if
c
      if  (source.ge.1) then
          read(fdatain, nml=source_g)
          if (rank.eq.0) then
            write(*,*)'namelist source_g have been read'
          end if
          if (abs(reynolds-dsqrt(rayleigh*froude/prandtl)).ge.1e-3) then
             write(*,*)'Adjust Re and Ra to satisfy Re=Sqrt(Ra*Fr/Pr)'
             stop
          end if
      end if
c=================================jaejin(begin)=========================
c       Read namelist for rotor simulation (but not used)
      read(fdatain, nml=rotor)
      if (rank.eq.0) then
        write(*,*) 'namelist rotor have been read'
      end if
c=================================jaejin(end)===========================
      
c=================================jaejin(begin)=========================
c       Read namelist for MHD simulation
      read(fdatain, nml=mhd)
      if (rank.eq.0) then
        write(*,*) 'namelist mhd have been read'
      end if
c=================================jaejin(end)===========================

      bcdir='xie'
      iso_tw=0.
      do while (bcdir.ne.'end')
         read(fdatain,nml=bcdef)
         if (bctype.ge.101.and.bctype.le.110) then
            read(fdatain,nml=iso_t)
            control%t(:)=iso_tw(:)
            exit
         end if 
      end do
c
      close(fdatain)

c=================================jaejin(begin)=========================
c     Journaling (by Jaejin)
      if (choice.eq.'n') then
        journal_info%tag = ''
        journal_info%char1 = '1st'
        journal_info%char2 = ''
        journal_info%intg1 = 0
        journal_info%intg2 = 0
        journal_info%doub1 = 0.
        journal_info%doub2 = 0.
        journal_info%doub3 = 0.
        journal_info%doub4 = 0.
        journal_info%doub5 = 0.
        journal_info%doub6 = 0.
      else
        journal_info%tag = ''
        journal_info%char1 = '2nd'
      end if
      if (rank.eq.0) then
          journal_caller = 'main_1'
        call journal(journal_caller,journal_info)
      end if
c=================================jaejin(end)===========================
c
c *** check input data
c
c..   check turbulence
      if(turb.gt.1) then
        write(*,*) 'turb avaiables are 0,1'
        stop
      end if

c..   all moving case will enable dual time stepping
      if(moving.ge.1.and.dual_t.eq.0.and.moving.ne.5) then
        write(*,*) 'dual_t updated to 1'
        dual_t = 1
      end if

c..   make sure computation necessary
      if(nstep.lt.0) then
        write(*,*) "nstep is less than 1"
        stop
      end if

c..   oscillating cascade current only available in 2D
      if(strtp.eq.3.and.moving.eq.1) then
        write(*,*) 'normal frequency: ', kc
        if(idimen.ne.2) then
          idimen = 2
          write(*,*) 'idimen set to 2 when strtp = 3'
        end if
        if(kl.ne.1) then
          kl = 1
          write(*,*)
     $         'osci cascade is currently only available in 2D'
        end if
      end if

c..   turn off dual time stepping (set step to 1) for dual_t 0, 2
      if(dual_t.eq.0.or.dual_t.eq.2) tsteps = 1

c..   enforce universal time step for dual_t=2
      if(dual_t.eq.2) unidt = 1

c..   turn off turbulence for inviscid case
      if (inviscid) turb = 0 
c
c *** check complete
c
c *** convert angle to radiant unit
c
      angl1 = angl1*3.14159d0/180.d0
      angl2 = angl2*3.14159d0/180.d0
c
c *** read multiblock information
c
      open(fdatain, file='mblock.dat')
      read(fdatain,*)nb,bc_max
      read(fdatain,*)il,jl,kl
      allocate (ill(nb),jll(nb),kll(nb),bnum(nb))
      do i=1,nb
        read(fdatain,*)ill(i),jll(i),kll(i)
      end do
      bnum = 0
      if (bc_max.gt.0) then
        allocate(btb(nb,bc_max),bnb(nb,bc_max),bdir(nb,bc_max),
     $        bstart(3,nb,bc_max),bend(3,nb,bc_max),rorder(nb,bc_max))
        read(fdatain,*)(bnum(i),i=1,nb)
        do i=1,nb
          read(fdatain,*)(btb(i,j),bnb(i,j),bdir(i,j),rorder(i,j),
     $                           j=1,bnum(i))
          read(fdatain,*)((bstart(k,i,j),bend(k,i,j),
     $                                   k=1,3),j=1,bnum(i))
        end do
      end if
      allocate(main_d(nb),ptemp(nb),pttemp(nb),tttemp(nb))
      allocate(iwake(3,2,nb))
c
c *** read the poutlet condition
c
      read(fdatain,*)(main_d(i),ptemp(i),pttemp(i),tttemp(i),i=1,nb)
      main_dir = main_d(rank+1)
      poutlet = ptemp(rank+1)
      ptotal = pttemp(rank+1)
      ttotal = tttemp(rank+1)
      write(*,*) 'poutlet = ',poutlet
c
c *** read the wake BC
c
      read(fdatain,*)(((iwake(j,k,i),j=1,3),k=1,2),i=1,nb)
c
c *** read the wall points
c
      read(fdatain,*) ipos
      allocate(xw(ipos), yw(ipos), zw(ipos))
      read(fdatain,*) (xw(i),yw(i),zw(i), i=1,ipos)
c
      close(fdatain)
c
      do j=1,3
        do k=1,2
          bcwake(j,k) = iwake(j,k,rank+1)
        end do
      end do
      il=ill(rank+1)
      jl=jll(rank+1)
      kl=kll(rank+1)
c     
c *** compute boundary limit of variable arrays
c
      ilow=1-blen
      iupp=il+blen
      jlow=1-blen
      jupp=jl+blen
      klow=1-blen
      kupp=kl+blen
c
c *** define dim1,dim2
c
      dim1=max(il,jl,kl)
      dim2=max(iupp,jupp,kupp)
c
c *** compute index and allocate space for x0,y0,z0
c     
      ix0=(iupp-ilow)*(jupp-jlow)*
     $                          (kupp-klow)
      allocate (x0(ix0),y0(ix0),z0(ix0))
      if (dual_t.eq.1) allocate (xo0(ix0),yo0(ix0),zo0(ix0))
      if (moving.gt.0) allocate (hx(il+1,kl+1,6),
     $                 hy(il+1,kl+1,6),hz(il+1,kl+1,6))
c
c *** compute index and allocate space for q0
c     
      iq0=(iupp-ilow+1)*(jupp-jlow+1)*(kupp-klow+1)*nl
      allocate (q0(iq0))
c
      if (dual_t.eq.1) allocate (qi0(iq0),qii0(iq0))
      if (precondition.ge.1) allocate (qiii0(iq0),diver(iq0))
c
c *** for turbulence (DES)
c
      allocate (rave(iq0),uave(iq0),vave(iq0),wave(iq0),eave(iq0))
      rave=0.
      uave=0.
      vave=0.
      wave=0.
      eave=0.
      nstatis=0
c
c *** compute index and allocate space for vol0
c     
      ivol=(iupp-ilow-1)*(jupp-jlow-1)*
     $               (kupp-klow-1)
      allocate (vol0(ivol))
c
c *** compute index and allocate space for vist0
c     
      ivist=ivol
      allocate (vist0(ivist))
c
c *** compute index and allocate space for BC 
c     
      ibcxie=jl*kl
      allocate (bcxie0(ibcxie),bcxie1(ibcxie))

      ibceta=il*kl
      allocate (bceta0(ibceta),bceta1(ibceta))

      ibczta=il*jl
      allocate (bczta0(ibczta),bczta1(ibczta))
c
c *** compute index and allocate space for dt0
c     
      it0=il*jl*kl
      allocate (dt0(it0))
c
c *** allocate space for dist0
c     
      allocate (dist0(it0))
c
c *** compute the total number of the cells
c     
      it0 = 0
      do i=1,nb
        it0=it0+ill(i)*jll(i)*kll(i)
      end do
c
c *** compute index and allocate space for moving grid
c     
cc    if (moving.ge.1) then
        iudi=4*(jupp-jlow+1)*(kupp-klow+1)
        allocate (udi0(iudi),vdi0(iudi),wdi0(iudi))

        iudj=4*(iupp-ilow+1)*(kupp-klow+1)
        allocate (udj0(iudj),vdj0(iudj),wdj0(iudj))

        iudk=4*(iupp-ilow+1)*(jupp-jlow+1)
        allocate (udk0(iudk),vdk0(iudk),wdk0(iudk))

        iqt=(iupp-ilow+1)*(jupp-jlow+1)*(kupp-klow+1)*3
        allocate (qt0(iqt))
cc    end if
c
      if (strtp.eq.3.and.moving.eq.1) then
        ixs=(il+1)*(jl+1)*snum
        allocate (xs0(ixs),ys0(ixs))
      end if

c---------------------parameters for transfer

      control%blen=blen
      control%dual_t=dual_t
      control%gcl=gcl
      control%idimen=idimen
      control%integrate=integrate
      control%iter_gs=iter_gs
      control%lhs_order=lhs_order
      control%lhs_scheme=lhs_scheme
      control%limiter=limiter
      control%main_dir=main_dir
      control%moving=moving
      control%nl=nl
      control%nrbc_ex=nrbc_ex
      control%precondition=precondition
      control%rhs_order=rhs_order
      control%rhs_scheme=rhs_scheme
      control%source=source
      control%strtp=strtp
      control%tsteps=tsteps
      control%turb=turb
      control%unidt=unidt
      control%vis_order=vis_order
      control%wall_order=wall_order

      control%inviscid=inviscid
      control%suther=suther

      control%angl1=angl1
      control%angl2=angl2
      control%cfl=cfl
      control%epsilon=epsilon
      control%epsfactor=epsfactor
      control%froude=froude
      control%gamma=gamma
      control%kfactor=kfactor
      control%k_prec=k_prec
      control%machinf=machinf
      control%poutlet=poutlet
      control%ptotal=ptotal
      control%prandtl=prandtl
      control%prt=prt
      control%rayleigh=rayleigh
      control%reynolds=reynolds
      control%ronum=ronum
      control%tref=tref
      control%tintvl=tintvl
      control%ttotal=ttotal
      control%varepsilon=varepsilon
      control%velinit=velinit
c--------------------end parameters for transfer

c=================================jaejin(begin)=========================
c     Allocate space for MHD variables
c     if (mhdsim.eq.1) then
         allocate(sigma_e_star(ilow:iupp,jlow:jupp,klow:kupp),
     $     B_x_star(ilow:iupp,jlow:jupp,klow:kupp),
     $     B_y_star(ilow:iupp,jlow:jupp,klow:kupp),
     $     B_z_star(ilow:iupp,jlow:jupp,klow:kupp),
     $     B_norm(ilow:iupp,jlow:jupp,klow:kupp))
         allocate(source_mhd(il,jl,kl,2:5))
c     end if
c=================================jaejin(end)===========================

c
c *** initialization
c
      q0 = 1.d0
      if (rank.eq.0) then
        write(*,*)'initialization begin'
      end if
      if (precondition.ge.2) then
         write(*,*)'generate the medium data'
         call generate_medium(choice,rank+1,nl,
     $        ilow,iupp,jlow,jupp,klow,kupp,
     $        control)
      end if

      call initflow(rank+1, il, jl, kl, nl, ilow, iupp,
     $       jlow, jupp, klow, kupp,  
     $       snum, ntotal, time,
     $       bcxie0, bcxie1, bceta0,
     $       bceta1, bczta0, bczta1,
     $       xsl, xs0, ys0, vist0,
     $       x0, y0, z0, 
     $       xo0, yo0, zo0, 
     $       q0, qi0, qii0,
     $       qiii0,control)
c
c *** get average quantities for DES
c
      call initial_tbl(rank+1,il,jl,kl,nl,ilow,iupp,
     $       jlow,jupp,klow,kupp,nstatis,time,
     $       rave,uave,vave,wave,eave)
c
c *** compute the volume and/or udi,...,wdk
c
      if(moving.eq.0.or.moving.eq.5) then
        call volume(x0,y0,z0,
     $        vol0,ilow,iupp,jlow,jupp,
     $        klow, kupp, il, jl, kl,
     $        bcxie0,bcxie1,bceta0,
     $        bceta1,bczta0,bczta1)
        if ( moving.eq.5 ) then
          call wall_rtv(il,jl,kl,ilow,iupp,jlow,jupp,
     $     klow,kupp,x0,y0,z0,udi0,vdi0,wdi0,udj0,vdj0,wdj0,
     $     udk0,vdk0,wdk0,ronum)
          call rtvel(il, jl, kl, x0, y0, z0,
     $       ilow, iupp, jlow, jupp, klow, kupp,
     $       qt0, ronum)
        end if
      else
        call wall_vel(il,jl,kl,ilow,iupp,
     $        jlow,jupp,klow,kupp,
     $        x0,y0,z0, 
     $        xo0,yo0,zo0,tintvl, 
     $        udi0,vdi0,wdi0,
     $        udj0,vdj0,wdj0,
     $        udk0,vdk0,wdk0)
      end if
c
      if ( dabs(ronum).gt.1e-9 ) then
        call wall_rtv(il,jl,kl,ilow,iupp,jlow,jupp,
     >    klow,kupp,x0,y0,z0,udi0,vdi0,wdi0,udj0,vdj0,wdj0,
     >    udk0,vdk0,wdk0,ronum)
        call rtvel(il, jl, kl, x0, y0, z0,
     >       ilow, iupp, jlow, jupp, klow, kupp,
     >       qt0, ronum)
      end if
c
c *** update q0
c
      call exchange_q(rank+1,nb,blen,nl,
     $       bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $       ilow,iupp,jlow,jupp,klow,kupp,q0)

      if (precondition.ge.1) then
      call exchange_q(rank+1,nb,blen,nl,
     $       bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $       ilow,iupp,jlow,jupp,klow,kupp,qiii0)
      end if
c
c *** compute the nearest wall distance of the points for 1eq model
c
      dist0(:) = 0.d0
      if (nl.eq.6) then
        call wall_dist(il, jl, kl, ides, ipos,
     $     ilow, iupp, jlow, jupp, klow, kupp,
     $     bcxie0, bcxie1, bceta0, bceta1, bczta0, bczta1,
     $     x0, y0, z0, xw, yw, zw, cdes, dist0)
      end if
c
c *** show computation case
c
      if (rank.eq.0) then
        write(6,*) 'Start iteration ...'
        select case (integrate)
        case(1)
           write(6,*) 'Using AF method'
        case(2)
           write(6,*) 'Using RK method'
        case(3)
           write(6,*) 'Using EU method'
        case(4)
           write(6,*) 'Using GS method'
        end select

        select case (turb)
        case (0)
          write(6,*) 'Laminar flow'
        case (1)
          if (nl.eq.5) then
            write(6,*) 'Turbulent flow (BL model)'
          else if (nl.eq.6) then
            write(6,*) 'Turbulent flow (SA 1eq model)'
          end if
        case (2,3)
          write(6,*) 'Turbulent flow (K-Omega model)'
        end select
      end if
c
c *** set time to zero in steady state case
c
      if (dual_t.eq.0) time = 0.d0
c
c *** initial turb, time
c
      if (dual_t.eq.0) then
        dstart = ntotal+1
        nstep = ntotal + nstep
      else
        dstart = nint(time/tintvl)+1
      end if
c
c *** structure related info for strp=1,2
c
c     --- TIME INTEGRATION ---
c     
      if(dual_t.eq.2) then
         recstep = dstart
         if(rank.eq.0) write(*,*) 'start from step: ', recstep,
     $        ' at time: ', (recstep-1)*tintvl
      end if
      if (rank.eq.0) write(*,*)'time integration begin'

      tstart = systime()
c
      open(23, file = 'cdlt.his')
      open(24, file = 'mflow.plt')
      open(25, file = 'rtrslt.plt')
c     open(26, file = 'point_3.his')

      do dstep = dstart, nstep

c=================================jaejin(begin)=========================
c       Prepare the convergence history file
      if (choice.eq.'n') then
          open(27,file='convhist',status='replace')
        else
          open(27,file='convhist',status='old',access='append')
      end if
c=================================jaejin(end)===========================


        do nn = 1, tsteps
c *** computation according to each block
          if (moving.eq.2.or.(moving.eq.1.and.nn.eq.1)) then
c
c     --- update moving grid info ---
c
            call moving_grid(il,jl,kl,nl,
     $           ilow,iupp,jlow,jupp,klow,kupp,
     $           dim2,dstep,nn,nb, rank, snum,
     $           xsl, xzl, time,
     $           bcxie0,bcxie1, bceta0,bceta1, bczta0,bczta1,
     $           xs0,ys0, x0,y0,z0, xo0,yo0,zo0, qt0,q0,vol0,
     $           vbd, xyz,
     $           udi0,vdi0,wdi0, udj0,vdj0,wdj0, udk0,vdk0,wdk0,
     $           kc, cs, mus, ub, wh, walf, uinf, nfp, alf0,
     $           nmode, hx, hy, hz, wpw, bcwake, cv1, lref, aref,
     $           xctr, yctr, xref, yref, xtor, ytor, jfx, beta,
     $           xref0, yref0, control)
c
c *** update the velocities of cell and coordinates of boundary
c
            call exchange_q(rank+1,nb,blen,3,
     $          bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $          ilow,iupp,jlow,jupp,klow,kupp,qt0)

            call exchange_xyz(rank+1,nb,blen,
     $          bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $          ilow,iupp,jlow,jupp,klow,kupp,x0)

            call exchange_xyz(rank+1,nb,blen,
     $          bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $          ilow,iupp,jlow,jupp,klow,kupp,y0)

            if (idimen.gt.2) then
              call exchange_xyz(rank+1,nb,blen,
     $          bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $          ilow,iupp,jlow,jupp,klow,kupp,z0)
            end if

            call volume(x0,y0,z0,
     $           vol0,ilow,iupp,jlow,jupp,
     $           klow, kupp, il, jl, kl,
     $           bcxie0,bcxie1,bceta0,
     $           bceta1,bczta0,bczta1)

            call wall_vel(il,jl,kl,ilow,iupp,
     $           jlow,jupp,klow,kupp,
     $           x0,y0,z0, 
     $           xo0,yo0,zo0,tintvl, 
     $           udi0,vdi0,wdi0,
     $           udj0,vdj0,wdj0,
     $           udk0,vdk0,wdk0)
c
c *** compute the nearest wall distance of the points
c
c           if (nl.eq.6) then
c             call wall_dist(il, jl, kl, ides, ipos,
c    $           ilow, iupp, jlow, jupp, klow, kupp,
c    $           bcxie0, bcxie1, bceta0, bceta1, bczta0, bczta1,
c    $           x0, y0, z0, xw, yw, zw, cdes, dist0)
c           end if
          end if
c     
c     --- compute turbulence viscosity ---     
c
          if (turb.ne.0.and.(.not.inviscid)) then
c
c *** define the values of visturb 
c
            if ( dabs(ronum).gt.1e-9.and.nl.eq.5 ) then
            call turb_rot(x0,y0,z0,
     $           ilow,iupp,jlow,jupp,klow,kupp,
     $           il,jl,kl,nl,q0,
     $           bcxie0,bcxie1,
     $           bceta0,bceta1,
     $           bczta0,bczta1,
     $           vol0,vist0,
     $           udi0,vdi0,wdi0,udj0,vdj0,wdj0,udk0,vdk0,wdk0,
     $           qt0,control)
            else
            call turb_new(bcwake,x0,y0,z0, 
     $           ilow,iupp,jlow,jupp,klow,kupp,
     $           il, jl, kl, nl, q0,
     $           bcxie0,bcxie1,
     $           bceta0,bceta1,
     $           bczta0,bczta1,
     $           vol0, vist0, nb, rank+1,
     $           udi0,vdi0,wdi0,
     $           udj0,vdj0,wdj0,
     $           udk0,vdk0,wdk0,cv1,qt0,control)
            end if
c
c *** update the turbulence
c
            call exchange_vist(rank+1,nb,blen,
     $           bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $           ilow,iupp,jlow,jupp,klow,kupp,vist0)
          end if

          if (precondition.ge.2) then
             call diverg(x0,y0,z0,vol0,ilow,iupp,jlow,jupp,
     $            klow,kupp,il,jl,kl,nl,qiii0,diver)
          end if
c
c *** determines the allowable timestep
c
c
c --- compute time step interval ---
c
          call timestep(dt0,il,jl,kl,nl,
     $         x0,y0,z0, 
     $         q0,dim2,vol0,
     $         ilow,iupp,jlow,jupp,klow,kupp,
     $         bcxie0,bcxie1,
     $         bceta0,bceta1,
     $         bczta0,bczta1,
     $         qiii0,diver,
     $         qt0,control)
c
c
          if(nb.gt.1.and.unidt.gt.0) then
            dtemp=dt0(1)
            call exchange_dt(dtemp, rank, nb)
            dt0 = dtemp
          end if
c
c *** integrate one time step
c
          if(unidt.ne.0.and.dual_t.ne.1) time=time+dt0(1)
          ntotal=ntotal+1

c
          call integrate_ALL(nn,
     $         x0,y0,z0, 
     $         vol0,q0,
     $         bcxie0,bcxie1,
     $         bceta0,bceta1,
     $         bczta0,bczta1,
     $         delta, ntotal, time,
     $         il, jl, kl, nl, dim1, dim2,
     $         ilow,iupp,jlow,jupp,klow,kupp,
     $         resdulmax, loc, rank, nb, 
     $         vist0, dt0, resdulave,
     $         qi0, qii0,qiii0,
     $          dstep, ke,
     $         udi0,vdi0,wdi0,
     $         udj0,vdj0,wdj0,
     $         udk0,vdk0,wdk0,
     $         qt0,resdl,sigma,nrbc_ex,
     $         iblnu, ipt, jpt, kpt, ic1, ic2, ic3, tko, cb1, cb2,
     $         cap_k, cw1, cw2, cw3, cv1, ct1, ct2, ct3, ct4,
!     $         dist0,  diver, control) !original
     $         dist0, diver, control, !jaejin
c=================================jaejin(begin)=========================
c         Arguments for MHD simulation
     $         intersteps, rho_inf, T_inf, c_sound, L_ref, mhdsim,
     $         ionmodel, T_onset, sigma_e_ref, magtype, magref,
     $         B_ref, I_in, offset, sigma_e_star, B_x_star, B_y_star,
     $         B_z_star, B_norm, source_mhd)
c=================================jaejin(end)===========================

c      print *, 'Passed this point!'

c *** compute resduals and location
c
          resdulave = dfloat(il*jl*kl)*resdulave
          loc0(1:3) = loc(1:3)
          loc0(4) = rank+1
          call exchange_res(rank,np,it0,resdulave,resdulmax,loc0)
c
c *** update q0
c
          call exchange_q(rank+1,nb,blen,nl,
     $         bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $         ilow,iupp,jlow,jupp,klow,kupp,q0)


          if (precondition.ge.1) then
          call exchange_q(rank+1,nb,blen,nl,
     $         bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $         ilow,iupp,jlow,jupp,klow,kupp,qiii0)
          end if

c..   monitor info inside dual time stepping

          if(rank.eq.0.and.dual_t.eq.1) then
               if(nn.eq.1) write(*,'(/,"DTS, ", "Step: ",
     $              i4, ", Time: ", e12.4, ", Start: ", i10/)')
     $              dstep,  time+tintvl, ntotal
               if(mod(nn,checksteps).eq.0)
     $              print '(i6, d21.14, 4i4, d21.14)',
     $              nn, resdulmax, loc0, resdulave
          end if

          if(resdulmax.lt.eps) exit
            
        end do

c..   monitor info inside regular time marching

        if(dual_t.eq.0.and.mod(dstep,checksteps).eq.0) then
          if(rank.eq.0) then
            write(*,'(i6, d21.14, 4i5, d21.14)')
     $            dstep, resdulmax, loc0, resdulave
c=================================jaejin(begin)=========================
c       Record the convergence history to convhist
            write(27,'(i6, d21.14, 4i5, d21.14)')
     $            dstep, resdulmax, loc0, resdulave
c=================================jaejin(end)===========================

          end if
        end if
         
        if(dual_t.eq.2.and.mod(dstep,checksteps).eq.0) then
          if(rank.eq.0) then
            write(*,'(i10, d21.14, 4i5, d21.14, e12.4)')
     $            dstep, resdulmax, loc0, resdulave, time
          end if
        end if
         
        if(dual_t.eq.1.and.rank.eq.0) then
          if(resdulmax.lt.eps) then
            write(*,'(a,e12.5,/)') 'maximum residual reaches eps: ',
     $              resdulmax
          else
            write(*,'(a,i4,/)') 'maximum steps reaches tsteps: ',
     $              tsteps
          end if
        end if
c
c *** time marching
c
        if(dual_t.eq.1) then
          time = time + tintvl
          qii0 = qi0
          qi0 = q0
          if (moving.gt.0) then
            xref = xref0
            yref = yref0
          end if
        end if
c
c *** record time, xyz(1), xyz(2), cd, cl, ct each time step for
c *** case 1, 2 & 3
c
        if ( dual_t.eq.1) then
          if (turb.ne.0) then
            if ( dabs(ronum).gt.1e-9.and.nl.eq.5 ) then
            call turb_rot(x0,y0,z0,
     $           ilow,iupp,jlow,jupp,klow,kupp,
     $           il,jl,kl,nl,q0,
     $           bcxie0,bcxie1,
     $           bceta0,bceta1,
     $           bczta0,bczta1,
     $           vol0,vist0,
     $           udi0,vdi0,wdi0,udj0,vdj0,wdj0,udk0,vdk0,wdk0,
     $           qt0,control)
            else
            call turb_new(bcwake,x0,y0,z0, 
     $           ilow,iupp,jlow,jupp,klow,kupp,
     $           il, jl, kl, nl, q0,
     $           bcxie0,bcxie1,
     $           bceta0,bceta1,
     $           bczta0,bczta1,
     $           vol0, vist0, nb, rank+1,
     $           udi0,vdi0,wdi0,
     $           udj0,vdj0,wdj0,
     $           udk0,vdk0,wdk0,cv1,qt0,control)
            end if
c
c *** update the turbulence
c
            call exchange_vist(rank+1,nb,blen,
     $           bc_max,bnum,btb,bnb,bdir,bstart,bend,rorder,idimen,
     $           ilow,iupp,jlow,jupp,klow,kupp,vist0)
          end if
c
c *** compute the cl,cd and cm
c
          call cldm(rank, nb, il, jl, kl, nl,
     $    ilow, iupp, jlow, jupp, klow, kupp,
     $    bcxie0, bcxie1, bceta0, bceta1, bczta0, bczta1,
     $    x0, y0, z0, q0, vist0,
     $    udj0, vdj0, wdj0, lref, aref, xref, yref, tcl, tcd, tcm,
     $    control)
c
          if (rank.eq.0) then
            select case (strtp)
            case(0)
              write(23,102) time, tcl, tcd, tcm, xref, yref 
            case(1)
              write(23,103) time, xyz(1,1), xyz(2,1), tcl, tcd
            case(2)
              write(23,101) time, xyz(1,1), xyz(2,1),
     $                      tcl, tcd, tcm, xref, yref
            case(3)
              write(23,101) time, -2.d0*xyz(1,1), xyz(2,1),
     $                      tcl, tcd, tcm, xref, yref
            case(5)
              write(23,101) time, xyz(1,1), xyz(2,1),
     $                      tcl, tcd, tcm ,xref, yref
            end select
          end if
        end if                                                            
c
c *** record pressure history
c
c       if(dual_t.eq.1) then
c           call record_p(dstep,rank+1,time,il,jl,kl,nl,
c    $           ilow,iupp,jlow,jupp,klow,kupp,
c    $           x0,y0,q0,gamma,dual_t)
c       end if
c
        if(dual_t.eq.2.and.time.gt.tintvl*recstep) then
          if(rank.eq.0) write(*,*) 'write unsteady data record ',
     $           recstep, ' at time ', time

c
c *** record pressure history
c
          call record_p(recstep,rank+1,time,il,jl,kl,nl,
     $         ilow,iupp,jlow,jupp,klow,kupp,
     $         x0,y0,q0,gamma,dual_t)
c
c *** output results
c
          if(mod(recstep,intersteps).eq.0)
     $       call output(rank+1,dual_t,recstep,ntotal,il,jl,kl,
     $        nl,ilow,iupp,jlow,jupp,klow,kupp,
     $        x0,y0,z0, 
     $        xo0,yo0,zo0,xsl, 
     $        bcxie0,bcxie1,
     $        bceta0,bceta1,
     $        bczta0,bczta1,q0,
     $        qii0,time,tintvl,turb,vist0,blen,precondition,qiii0)
c=================================jaejin(begin)=========================
c         Output results for  MHD simulation
            if(mhdsim.eq.1)
     $        call output_mhd(rank+1,dual_t,recstep,ntotal,il,jl,kl,
     $          nl,ilow,iupp,jlow,jupp,klow,kupp,
     $          x0,y0,z0, 
     $          xo0,yo0,zo0,xsl, 
     $          bcxie0,bcxie1,
     $          bceta0,bceta1,
     $          bczta0,bczta1,q0,
     $          qii0,time,tintvl,turb,vist0,blen,precondition,qiii0,
     $          sigma_e_star, B_x_star, B_y_star, B_z_star, B_norm,
     $          source_mhd)
c=================================jaejin(end)===========================
            recstep = recstep+1
c
        end if

c..   record temporary or unsteady results

        if(dual_t.lt.2.and.mod(dstep,intersteps).eq.0) then
          call output(rank+1,dual_t,dstep,ntotal,il,jl,kl,
     $        nl,ilow,iupp,jlow,jupp,klow,kupp,
     $        x0,y0,z0, 
     $        xo0,yo0,zo0,xsl, 
     $        bcxie0,bcxie1,
     $        bceta0,bceta1,
     $        bczta0,bczta1,q0,
     $    qii0,time,tintvl,turb,vist0,blen,precondition,qiii0)

c=================================jaejin(begin)=========================
c       Record temporary or unsteady results for MHD simulation
          if (mhdsim.eq.1)
     $      call output_mhd(rank+1,dual_t,dstep,ntotal,il,jl,kl,
     $          nl,ilow,iupp,jlow,jupp,klow,kupp,
     $          x0,y0,z0, 
     $          xo0,yo0,zo0,xsl, 
     $          bcxie0,bcxie1,
     $          bceta0,bceta1,
     $          bczta0,bczta1,q0,
     $          qii0,time,tintvl,turb,vist0,blen,precondition,qiii0,
     $          sigma_e_star, B_x_star, B_y_star, B_z_star, B_norm,
     $          source_mhd)
c=================================jaejin(end)===========================

          if ( dabs(ronum).gt.1e-9 ) then
          call mass_flow(il, jl, kl, nl, x0, y0, z0, time, dstep,
     $                   ilow, iupp, jlow, jupp, klow, kupp, q0)
          call rotrslt(il, jl, kl, nl, x0, y0, z0, time, dstep,
     $                 ilow, iupp, jlow, jupp, klow, kupp, q0,
     $                 gamma, machinf, udj0, vdj0, wdj0, vol0)
          end if
        end if
        if(rank.eq.0.and.dual_t.eq.1) then
          if(resdulmax.lt.eps) print *, 'eps reached: ', resdulmax,
     $          ' at step: ', nn
          if(nn.ge.nstep) print *, 'step reached: ', nn,
     $          ' at residual: ', resdulmax
          tend = systime()
          write(*,'(a,e14.8,a)') 'time used in iteration: ',
     $           tend-tstart, ' s'
          tstart = tend
        end if

        inquire(file='stop', exist=ex)
        if(ex) then
          if(rank.eq.0) write(*,*)
     $           '"stop" signal found, terminate at step ', dstep
          exit
        end if
        
        if(resdulmax.lt.eps.and.dual_t.eq.0) exit
c
c *** computing and saving average quantities
c
             nstatis=nstatis+1
             call average(rank+1,nstatis,time,il,jl,kl,
     $        nl,ilow,iupp,jlow,jupp,klow,kupp,q0,
     $        x0,y0,z0,
     $        rave,uave,vave,wave,eave)

            if (mod(dstep,intersteps).eq.0) then
              call output_tbl(rank+1,il,jl,kl,
     $            nl,ilow,iupp,jlow,jupp,klow,kupp,nstatis,time,
     $            rave,uave,vave,wave,eave)
            end if

      end do

      if(dual_t.eq.0.or.dual_t.eq.1) then
         tend = systime()
         write(*,'(a,e14.8,a)') 'job computation time: ', tend-tstart,
     $        ' s'
      end if

c     
c     --- deallocate space ---
c     
      deallocate (ill,jll,kll,
     $     x0,y0,z0,q0,vol0,vist0,bnum,
     $     bcxie0,bcxie1,bceta0,bceta1,bczta0,bczta1,dt0,dist0)
      if (bc_max.gt.0) then
        deallocate (btb,bnb,bdir,bstart,bend,rorder)
      end if
      deallocate (main_d,ptemp,pttemp,tttemp,iwake,xw,yw,zw)

      if(dual_t.eq.1) deallocate(xo0, yo0, zo0, qi0, qii0)
      if (precondition.ge.1) deallocate(qiii0,diver)
      
      if(moving.ge.0)
     $     deallocate(udi0,vdi0,wdi0,udj0,vdj0,wdj0,
     $     udk0,vdk0,wdk0,qt0)

      if(strtp.eq.3.and.moving.eq.1) deallocate(xs0, ys0)

c=================================jaejin(begin)=========================
c       Deallocate MHD variables
c     if (mhdsim.eq.1)
      deallocate(sigma_e_star, B_x_star, B_y_star, B_z_star, B_norm,
     $  source_mhd)
c=================================jaejin(end)===========================


c     ... finalizing mpi ...
c     
      call mpi_final()

c=================================jaejin(begin)=========================
c       Journaling
      if (rank.eq.0) then
        if (choice.eq.'n') then
          journal_caller = 'main_2'
          journal_info%doub6 = tend - tstart
          call journal(journal_caller,journal_info)
        else
          journal_caller = 'main_add'
          journal_info%doub6 = tend - tstart
          call journal(journal_caller,journal_info)
        end if
      end if
c=================================jaejin(end)===========================

c
  101 format(8f12.6)
  102 format(6f12.6)
  103 format(5f12.6)
c     
      end program FASIP
C Build Time: Mon Jun  8 04:47:39 EDT 2009
