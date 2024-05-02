	double precision::
     $       RC,ca1,ca2,cb1,cb2,cd1,
     $       cd2,ce1,ce2,cf1,cf2,tc,pc,aa,bb,mur,s,
     $       UW,TW,RHOWW,TW0

	integer istate
	double precision:: t_p(401,401),t_t(401),t_v(401),p_c(401)
     $                   ,t_cp(401),t_mu(401),t_h(401,401),v_c(401,2)
	common 
     $       /co/RC,ca1,ca2,cb1,cb2,cd1,cd2,ce1,ce2,
     $           cf1,cf2,tc,pc,aa,bb,mur,s,
     $           UW,TW,RHOWW,TW0,istate
	common /table/t_p,t_t,t_v,p_c,t_cp,t_mu,t_h,v_c

	double precision::
     $       sol,sog,svl,svg
	common /medi/sol,sog,svl,svg

