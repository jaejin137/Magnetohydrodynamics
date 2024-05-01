      module datain

      type datain_type

      integer::
     $ blen,
     $ dual_t,
     $ gcl,
     $ idimen,
     $ integrate,
     $ iter_gs,
     $ lhs_order,
     $ lhs_scheme,
     $ limiter,
     $ main_dir,
     $ moving,
     $ nl,
     $ nrbc_ex,
     $ precondition,
     $ rhs_order,
     $ rhs_scheme,
     $ source,
     $ strtp,
     $ tsteps,
     $ turb,
     $ unidt,
     $ vis_order,
     $ wall_order 

      double precision:: 
     $ angl1,
     $ angl2,
     $ cfl,
     $ epsilon,
     $ epsfactor,
     $ froude,
     $ gamma,
     $ kfactor,
     $ k_prec,
     $ machinf,
     $ poutlet,
     $ ptotal,
     $ prandtl,
     $ prt,
     $ rayleigh,
     $ reynolds,
     $ ronum,
     $ tref,
     $ tintvl,
     $ ttotal,
     $ varepsilon,
     $ velinit


      double precision, dimension(10):: t  ! iso_tw

      logical:: 
     $ inviscid,
     $ suther




      end type datain_type
      
      end module datain
