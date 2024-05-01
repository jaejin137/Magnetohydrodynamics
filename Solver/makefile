objects = main.o invert.o rhside.o roe_matrix.o boundary.o \
	timestep.o negative.o initflow.o output.o volume.o \
	reconstruct.o updaterk.o invsid_flux.o lhs_matrix.o integrate_ALL.o\
	ursn.o mpi_exchange.o turb_new.o metric_all.o matrix_bnd.o\
        turb_line.o accessory.o zha_matrix.o van_leer_matrix.o \
        van_leer_flux.o vanleer_hanel_flux.o distance.o\
        compute_source.o kw_wall.o wall_func.o wall_vel.o\
        structure_para.o compute_qt.o spvel.o source_move.o cdltj.o\
	strmodel.o metric.o pbtrip.o block.o\
	save_mesh.o integ_nri_GS.o integ_nrj_GS.o\
        sr_nri.o isr_nri.o urs_nri.o urs_nrj.o sr_nrj.o isr_nrj.o record_p.o\
        moving_grid.o yplus.o\
        exchange_dt.o exchange_res.o exchange_xyz.o\
        exchange_q.o exchange_vist.o cldm.o lu_sgs.o ursn_1.o ursn_2.o\
        ivsnr.o weno_5.o \
        lu_ursn.o lu_ursn_gs.o \
        vorticity.o strmodelaf.o source_1eq.o wall_dist.o\
        initial_tbl.o average.o  output_tbl.o mvgrid_cyl.o mvgrid_af.o \
        gauss_seidel.o strafnlr.o roe_matrix_p.o \
        dwdq.o pre_matrix.o state.o diverg.o generate_medium.o source_gravity.o \
        rtvel.o wall_rtv.o vortm.o mass_flow.o turbulencej.o \
        turb_rot.o mcontrol.o \
	newtypes.o source_mhd_lowrey.o eleccond.o magfld.o journal.o \
	output_mhd.o


#
main.out: $(objects)
#
	sh utilities/addtimestamp.sh main.f
	mpif77 -O2 -o main.out $(objects)

.f.o:
#
	mpif77 -O2 -c $<


clean:
	rm *~ *.o
