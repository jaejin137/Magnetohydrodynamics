        subroutine journal(journal_caller,journal_info)

        ! This subroutine journals useful information about each run so
        ! that the user can easily understand the setup for the
        ! simulation.

        use newtypes
        
        implicit none

        integer::
     $    blen,
     $    idimen,
     $    nl,
     $    dual_t,
     $    integrate,
     $    precondition,
     $    lhs_order,
     $    lhs_scheme,
     $    limiter,
     $    moving,
     $    rhs_order,
     $    rhs_scheme,
     $    turb,
     $    mhdsim,
     $    nstep,
     $    iter_gs,
     $    unidt,
     $    intersteps,
     $    checksteps,
     $    tsteps,
     $    gcl,
     $    source,
     $    strtp,
     $    nrbc_ex,
     $    vis_order,
     $    strm_dir,
     $    tkpar,
     $    twpar,
     $    wall_order,
     $    block,
     $    bctype,
     $    start(3),
     $    end(3),
     $    iblock,
     $    istart(3),
     $    iend(3),
     $    order(3),
     $    iso_tw_num,
     $    rotide

        logical::
     $    suther,
     $    inviscid

        double precision::
     $    delta(5),
     $    epsfactor,
     $    gamma,
     $    prandtl,
     $    prt,
     $    tref,
     $    angl1, angl2,
     $    machinf,
     $    poutlet, ptotal,
     $    ttotal,
     $    reynolds,        
     $    cfl,
     $    comptime,
     $    eps,
     $    rho_inf,
     $    p_inf,
     $    T_inf,
     $    mu_inf,
     $    c_sound,
     $    L_ref,
     $    T_onset,
     $    sigma_e_inf, sigma_e_ref, sigma_e_max,
     $    B_ref,
     $    I_in,
     $    offset(3),
     $    B_max,
     $    Re_m_inf,
     $    velinit,
     $    tintvl,
     $    kfactor,
     $    k_prec,
     $    sigma(2),
     $    kc,
     $    theta,
c     $    vbd(3,3),xyz(3,3),
     $    magref_val,
     $    ronum,
     $    rotalp,
     $    varepsilon,
     $    iso_tw(10),
     $    isotherm_temp

        character::
     $    choice,
     $    journal_caller*20,
     $    visc_model*16,
     $    viscous*3,
     $    turb_model*30,
     $    mhdswitch*3,
     $    lhs_ordr*22,
     $    lhs_schm*20,
     $    rhs_ordr*23,
     $    rhs_schm*21,
     $    integrate_method*20,
     $    flux_limiter*20,
     $    dual_time*3,
     $    ionmodel,
     $    ion_model*30,
     $    magtype*10,
     $    mag_type*30,
     $    magref,
     $    magref_mssg1*40,
     $    magref_mssg2*34,
     $    streamdir*20,
     $    bcdir*3,
     $    bc_type*70,
     $    interpret_bc_type*70

        integer, parameter::
     $    fdatain = 80,
     $    foutput = 81

        type(journaldata)::
     $    journal_info

        namelist /consts/ delta, epsfactor, gamma, prandtl, prt,
     $    suther, tref

        namelist /flows/ blen, idimen, nl, angl1, angl2, inviscid,
     $    machinf, poutlet, ptotal, reynolds, tintvl, ttotal, turb

        namelist /comp/ nstep, eps, cfl, limiter, lhs_scheme,
     $    rhs_scheme, lhs_order, rhs_order, integrate, iter_gs, unidt,
     $    dual_t, intersteps, checksteps, tsteps, moving, gcl,
     $    precondition, k_prec, strtp, nrbc_ex, sigma, kfactor, kc,
     $    vis_order, choice, velinit, strm_dir, theta, twpar, tkpar,
     $    source, rotide, varepsilon, wall_order

        namelist /rotor/ rotalp, ronum

        namelist /mhd/ rho_inf, p_inf, T_inf, mu_inf, c_sound, L_ref,
     $    mhdsim, ionmodel, T_onset, sigma_e_ref, magtype, magref,
     $    B_ref, I_in, offset

        namelist /bcdef/ bcdir, block, bctype, start, end,
     $    iblock, istart, iend, order
        
        namelist /iso_t/ iso_tw
        
        
        ! Subroutine Begins
        
        open(fdatain,file='datain')
        read(fdatain,nml=consts)
        read(fdatain,nml=flows)
        read(fdatain,nml=comp)
        read(fdatain,nml=rotor)
        read(fdatain,nml=mhd)
        read(fdatain,nml=iso_t)
        close(fdatain)

        ! Translate parameters to plain words

        sutherland: select case (suther)
        case(.true.)
          visc_model = "Sutherland's Law"
        case(.false.)
          visc_model = "Linear Model"
        case default
          print *, "*** Parameter error occurred: suther ***"
          stop
        end select sutherland
        
        viscosity: select case (inviscid)
        case(.true.)
           viscous = 'No'
        case(.false.)
           viscous = 'Yes'
        case default
           print *, "*** Parameter error occurred: inviscid ***"
           stop
        end select viscosity

        turbulence: select case (turb)
        case(0)
          turb_model = 'No -> Laminar Flow'
        case(1)
          turbulencemodel: select case(nl)
            case(5)
              turb_model = 'Baldwin-Lomax Model'
            case(6)
              turb_model = 'Spalart-Allmaras Model'
            case(7)
              turb_model = 'k-omega Model'
            case default
              print *, "*** Parameter error occurred: nl ***"
              stop
            end select turbulencemodel
        case default
          print *, "***Parameter error occurred: turb ***"
          stop
        end select turbulence

        ionization: select case(ionmodel)
        case('F')
          ion_model = 'Frozen Ionization Model'
        case('E')
          ion_model = 'Equilibrium Ionization Model'
        case default
          print *, "*** Parameter error occured: ionmodel ***"
          stop
        end select ionization

        mhdsimulation: select case (mhdsim)
        case(0)
          mhdswitch = 'Off'
        case(1)
          mhdswitch = 'On'
        case default
          print *, "*** Parameter error occurred: mhdsim ***"
          stop
        end select mhdsimulation

        magfldtype: select case (magtype)
        case('unifx')
          mag_type = 'Uniform in x-direction'
        case('unify')
          mag_type = 'Uniform in y-direction'
        case('unifz')
          mag_type = 'Uniform in z-direction'
        case('radial')
          mag_type = 'Radial'
        case('dipolex')
          mag_type = 'Dipole in x-direction'
        case('dipoley')
          mag_type = 'Dipole in y-direction'
        case('dipolez')
          mag_type = 'Dipole in z-direction'
        case('test')
          mag_type = 'Test Field'
        case default
          print *, "*** Parameter error occurred: magtype ***"
          stop
        end select magfldtype

        magreference: select case(magref)
        case('B')
          magref_mssg1 = 'Explicit Value for Magnetic Induction'
          magref_mssg2 = 'Reference Magnetic Induction (T) = '
          magref_val = B_ref
        case('I')
          magref_mssg1 = 'Input Current'
          magref_mssg2 = 'Reference Input Current (A) = '
          magref_val = I_in
        case default
          print *, "*** Parameter error occurred: magref ***"
          stop
        end select magreference

        lhsorder: select case(lhs_order)
        case(0)
          lhs_ordr = '1st Order MUSCL'
        case(1)
          lhs_ordr = '2nd or 3rd Order MUSCL'
        case default
          print *, "*** Parameter error occurred: lhs_order ***"
          stop
        end select lhsorder
        
        lhsscheme: select case(lhs_scheme)
        case(1)
          lhs_schm = 'Roe Scheme'
        case(2)
          lhs_schm = 'Zha Scheme'
        case(3)
          lhs_schm = 'None'
        case(4)
          lhs_schm = 'van Leer Scheme'
        case default
          print *, "*** Parameter error occurred: lhs_scheme ***"
          stop
        end select lhsscheme

        rhsorder: select case(rhs_order)
        case(0)
          rhs_ordr = '1st Order MUSCL'
        case(1)
          rhs_ordr = '2nd\/3rd Order MUSCL'
        case(4)
          rhs_ordr = '3rd Order WENO'
        case(5)
          rhs_ordr = '5th Order Fixed Stencil'
        case(6)
          rhs_ordr = '5th Order WENO'
        case(7)
          rhs_ordr = '7th Order Fixed Stencil'
        case(8)
          rhs_ordr = '7th Order WENO'
        case default
          print *, "*** Parameter error occurred: rhsorder ***"
          stop
        end select rhsorder
        
        rhsscheme: select case(rhs_scheme)
        case(1)
          rhs_schm = 'Roe Scheme'
        case(2)
          rhs_schm = 'Zha2 Scheme'
        case(3)
          rhs_schm = 'Zha ECUSPLD Scheme'
        case(4)
          rhs_schm = 'van Leer Scheme'
        case(5)
          rhs_schm = 'Edwards Scheme'
        case(6)
          rhs_schm = 'Zha6 Scheme'
        case(7)
          rhs_schm = 'AUSM+ Scheme'
        case(8)
          rhs_schm = 'Zha Scheme'
        case(9)
           rhs_schm = 'AUSMV Scheme'
        case(10)
          rhs_schm = 'AUSMD Scheme'
        case(11)
          rhs_schm = 'van Leer-Hanel'
        case default
          print *, "*** Parameter error occurred: rhs_scheme ***"
          stop
        end select rhsscheme

        integration: select case(integrate)
        case(1)
          integrate_method = 'AF Method'
        case(2)
          integrate_method = 'Runge-Kutta Method'
        case(3)
          integrate_method = 'Euler Method'
        case(4)
          integrate_method = 'Gauss-Seidel Method'
        case default
          print *, "*** Parameter error occurred: integrate ***"
          stop
        end select integration

        fluxlimiter: select case(limiter)
        case(0)
          flux_limiter = 'No Limiter'
        case(1)
          flux_limiter = 'MINMOD'
        case(2)
          flux_limiter = 'Super Bee'
        case default
          print *, "*** Parameter error occurred: limiter ***"
          stop
        end select fluxlimiter

        dualtime: select case(dual_t)
        case(0)
          dual_time = 'No'
        case(1)
          dual_time = 'Yes'
        case default
          print *, "*** Parameter error occurred: dual_t"
          stop
        end select dualtime

        streamdirect: select case(strm_dir)
        case(1)
          streamdir = 'Uniform Flow'
        case(4)
          streamdir = 'Blasius Profile'
        case default
          print *, "*** Parameter error occurred: strm_dir"
          stop
        end select streamdirect

8110    format(
     $    1x,"=============== MHD PARAMETERS ===============",/,
     $    1x,"Model: Frozen Ionization Model",/,
     $    1x,"Electrical Conductivity at Infinity (S/m) = ",es10.4,/,
     $    1x,"Electrical Conductivity of ROI(Max) (S/m) = ",es10.4,/,
     $    1x,"Type of Magnetic Induction: ",a30,/,
     $    1x,"Reference for Magnetic Induction: ",a27,/,
     $    1x,a34,es14.4,/,
     $    1x,"Maximum Value of Magnetic Induction (T) = ",es10.4,/,
     $    1x,"Magnetic Reynolds Number = ",es10.4,/)

8111    format(
     $    1x,"=============== MHD PARAMETERS ===============",/,
     $    1x,"Model: Equilibrium Ionization Model",/,
     $    1x,"Threshold Temperature for MHD Effects (K) = ",es10.4,/,
     $    1x,"Electrical Conductivity at Infinity (S/m) = ",es10.4,/,
     $    1x,"Electrical Conductivity of ROI(Max) (S/m) = ",es10.4,/,
     $    1x,"Type of Magnetic Induction: ",a30,/,
     $    1x,"Reference for Magnetic Induction: ",a40,/,
     $    1x,a34,es14.4,/,
     $    1x,"Maximum Value of Magnetic Induction (T) = ",es10.4,/,
     $    1x,"Magnetic Reynolds Number = ",es10.4,/)

8120    format(
     $    1x,"=============== FLUID PARAMETERS ===============",/,
     $    1x,"Density at Infinity (kg/m^3) = ",es10.4,/,
     $    1x,"Pressure at Infinity (Pa) = ",es10.4,/,
     $    1x,"Temperature at Infinity (K) = ",es10.4,/,
     $    1x,"Mach Number at Infinity = ",es10.4,/,
     $    1x,"Speed of Sound (m/s) = ",es10.4,/,
     $    1x,"Characteristic Length (m) = ",es10.4,/,
     $    1x,"Viscous: ",a3,/,
     $    1x,"Viscosity at Infinity (Pa*s) = ",es10.4,/,
     $    1x,"Viscosity Model: ",a16,/,
     $    1x,"Turbulence Model: ",a22,/,
     $    1x,"Molecular Prandtl Number = ",es10.4,/,
     $    1x,"Turbulent Prandtl Number = ",es10.4,/,
     $    1x,"Reynolds Number = ",es10.4,/,
     $    1x,"Total Pressure at Inlet (nondim.) = ",es10.4,/,
     $    1x,"Static Pressure at Outlet (nondim.) = ",es10.4,/,
     $    1x,"Total Temperature at Inlet (nondim.) = ",es10.4,/,
     $    1x,"Initial Velocity Profile:",a20,/)
        
8130    format(
     $    1x,"bcdir = ",a3,2x,"block = ",i2,2x,"start = ",i3,i3,i3,2x,
     $      "end = ",i3,i3,i3,/,
     $    1x,"---> ","bctype: ",a70)

8131    format(
     $    1x,"bcdir = ",a3,2x,"block = ",i2,2x,"start = ",i3,i3,i3,2x,
     $      "end=",i3,i3,i3,/,
     $    1x,"---> ","bctype: ",a20,2x,"wall temp=",es10.4)


8140    format(/,
     $    1x,"=============== COMPUTATION PARAMETERS ===============",/,
     $    1x,"Dimension = ",i1,/,
     $    1x,"LHS Order and Schemes : ",a22," and ",a15,/,
     $    1x,"RHS Order and Schemes : ",a22," and ",a15,/,
     $    1x,"Integration Method : ",a19,/,
     $    1x,"Flux Limiter : ",a10,/,
     $    1x,"Dual Time Steps : ",a3,/,
     $    1x,"Residual Limit(EPS) =",es10.4,/
     $    1x,"CFL for 1st Stage = ",es10.4,/
     $    1x,"Computation time = ",es16.9
     $  )

8150    format(
     $    1x,"CFL for Additional Stage =",es10.4
     $  )

8160    format(
     $    1x,"Computation Time for Additional Stage =",es16.9
     $  )

        
        select case (journal_caller)
        
        case ('main_1')
          if (journal_info%char1.eq.'1st') then
              open(foutput,file='journal.txt',status='replace')
              write(foutput,*) "This simulation was based on the
     $ following setup:"
              write(foutput,*) ""
              close(foutput)
            elseif (journal_info%char1.eq.'2nd') then
              open(foutput,file='journal.txt',status='old',
     $          access='append')
              write(foutput,8150) cfl
              close(foutput)
          end if
          
        case('source_mhd_lowrey')
          sigma_e_inf = journal_info%doub1
          sigma_e_max = journal_info%doub2
          B_max = journal_info%doub3
          Re_m_inf = journal_info%doub4
          open(foutput,file='journal.txt',status='old',
     $      access='append')
          select case(ionmodel)
            case('F')
              write(foutput,8110) sigma_e_inf, sigma_e_max, mag_type,
     $          magref_mssg1, magref_mssg2, magref_val, B_max, Re_m_inf
            case('E')
              write(foutput,8111) T_onset, sigma_e_inf, sigma_e_max,
     $          mag_type, magref_mssg1, magref_mssg2, magref_val,
     $          B_max, Re_m_inf
            end select
          
        case ('main_2')
          open(foutput,file='journal.txt',status='old',
     $      access='append')
          comptime = journal_info%doub6
          if (journal_info%char1.eq.'1st') then
            write(foutput,8120) rho_inf, p_inf, T_inf, machinf, c_sound,
     $        L_ref, viscous, mu_inf, visc_model, turb_model, prandtl,
     $        prt, reynolds, ptotal, poutlet, ttotal, streamdir
            write(foutput,*)
     $        "=============== BOUNDARY CONDITIONS ==============="
            open(fdatain,file='datain')
            do while(bcdir.ne.'end')
              read(fdatain,nml=bcdef)
              if(bcdir.eq.'end') then
                exit
              end if
              if (bctype.ge.101.and.bctype.le.110) then
                bc_type = 'Isothermal Wall B.C.'
                isotherm_temp = iso_tw(bctype-100)
                write(foutput,8131) bcdir, block, start(1), start(2),
     $            start(3), end(1), end(2), end(3), bc_type,
     $            isotherm_temp
              else
                bc_type = interpret_bc_type(bctype)
                write(foutput,8130) bcdir, block, start(1), start(2),
     $            start(3), end(1), end(2), end(3), bc_type
              end if
            end do
            close(fdatain)
            write(foutput,8140) idimen, lhs_ordr, lhs_schm, rhs_ordr,
     $        rhs_schm, integrate_method, flux_limiter, dual_time, eps,
     $        cfl, comptime
          elseif (journal_info%char1.eq.'2nd') then
            write(foutput,8160) comptime  
          end if
          close(foutput)
        
        end select
        return
            
        end subroutine journal

        function interpret_bc_type(bctype)
        implicit none
        integer:: bctype
        character:: interpret_bc_type*70
        select case(bctype)
        case(1)
          interpret_bc_type = 'Zero Gradient.'
        case(2)
          interpret_bc_type = 'Supersonic Inflow.'
        case(3)
          interpret_bc_type = 'No Slip Adiabatic Wall.'
        case(4)
          interpret_bc_type = 'Zero Gradient with w=0.'
        case(5)
          interpret_bc_type = 'Subsonic Outflow with Fixed Static
     $ Pressure (poutlet in datain).'
        case(6)
          interpret_bc_type = 'Subsonic Inflow with Fixed rho, u, v,
     $ and w at inlet.'
        case(7)
          interpret_bc_type = 'Inner Boundary for MPI.'
        case(8)
          interpret_bc_type = 'Symmetric B.C..'
        case(9)
          interpret_bc_type = 'Subsonic Inlet B.C. with Fixed Total
     $ Pressure and Temperature.'
        case(10)
          interpret_bc_type = 'Periodic B.C..'
        case(11)
          interpret_bc_type = 'Subsonic Outflow with Fixed Static
     $ Pressure (Computed in the Code).'
        case(19)
          interpret_bc_type = 'Isothermal Wall with Zero Gradient
     $ Presure.'
        case(20)
          interpret_bc_type = 'Periodic B.C. for Flow Variables Only.'
        case default
            print *, "*** Parameter error occurred: bctype ***"
            stop
        end select
        return
        end function interpret_bc_type
