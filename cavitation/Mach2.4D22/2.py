#!/usr/bin/python

# Dependencies and Logistics ===================================================
# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path
import math

Ma=2.4
p0a=101325
rho_g=1.204
rho_l=1000
gam_g = 1.4
gam_l=6.12
pi_l=3.43E8
rho = 1


c_a = math.sqrt( gam_g * ( p0a  ) / rho_g )
psOp0a = ( Ma ** 2 -1 ) * 2 * gam_g / ( gam_g + 1 ) + 1
rhosOrho0a = ( 1 + ( gam_g + 1 ) / ( gam_g - 1) * psOp0a ) / ( ( gam_g + 1 ) / ( gam_g - 1) + psOp0a ) 
Ms = math.sqrt( ( gam_g + 1. ) / ( 2. * gam_g ) * ( psOp0a - 1. ) * ( p0a / ( p0a ) ) + 1. )
vel = c_a/gam_g * (psOp0a - 1.) * p0a / ( p0a  ) / Ms

ps=p0a*psOp0a
rho_post_g=rhosOrho0a*rho_g

c_l = math.sqrt( 1.4*ps/rho )

D = 0.022

MPD=100 #mesh per diameter
Ny = 12*MPD/2
Nx = 14*MPD
dx = D/MPD #8.3e-6

time_end = 0.0002#50us
cfl = 0.2

dt = cfl * dx/c_l #5.3E-9
Nt = int(time_end/dt)#10000

print 'Nt', Nt

# Navigating to script directory
if len(dirname(argv[0])) != 0: chdir(dirname(argv[0]))

# Adding master_scripts directory to module search path
mfc_dir = '../../src'; path[:0] = [mfc_dir + '/master_scripts']

# Command to execute the MFC components
from m_python_proxy import f_execute_mfc_component
# ==============================================================================

# Case Analysis Configuration ==================================================

# Selecting MFC component
comp_name = argv[1].strip()

# Serial or parallel computational engine:
engine = 'parallel'
if (comp_name=='pre_process'): engine='serial'
if (comp_name=='post_process'): engine='serial'
# Configuring case dictionary
ppn=64

case_dict =                                                                    \
    {                                                                          \
                    # Logistics ================================================
                    'case_dir'                     : '\'.\'',                  \
                    'run_time_info'                : 'T',                      \
                    'nodes'                        : 'node06',                        \
                    'ppn'                          : ppn,                        \
                    'queue'                        : 'node02',                 \
                    'walltime'                     : '24:00:00',               \
                    'mail_list'                    : '',
                    # ==========================================================
                                                                               \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 :  -7*D,                   \
                    'x_domain%end'                 :  7*D,                   \
                    'y_domain%beg'                 :  0.,                  \
                    'y_domain%end'                 :  6*D,  
                         
                    'm'                            : int(Nx),                      \
                    'n'                            : int(Ny),                        \
                    'p'                            : 0,                        \
                    'dt'                           : dt,                   \
                    't_step_start'                 : 0,                        \

                    't_step_stop'                  : 200,                     \
                    't_step_save'                  :1,#(Nt/1000),                     \
		    # ==========================================================
                                                                               \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 3,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'F',                      \
                    'num_fluids'                   : 2,                        \
		            'adv_alphan'                   : 'T',                      \
		            'mpp_lim'                      : 'T',                      \
		            'mixture_err'                  : 'T',                      \
		            'time_stepper'                 : 3,                        \
'weno_shock' : 'T',\
'mixture_err':'T',                     
'weno_vars'                    : 2,                        \
                    'weno_order'                   : 5,                        \
                    'weno_eps'                     : 1.E-16,                   \
                    'char_decomp'                  : 'F',                      \
                    'mapped_weno'                  : 'T',                      \
                    'null_weights'                 : 'F',                      \
                    'mp_weno'                      : 'F',                                          \
                    'weno_Re_flux'                 : 'F',                      \
		            'riemann_solver'               : 2,                        \
                    'wave_speeds'                  : 1,                        \
                    'avg_state'                    : 2,                        \
                    'commute_err'                  : 'F',                      \
                    'split_err'                    : 'F',                      \
                    'reg_eps'                      : 1.E+00,       
                    'bc_x%beg'                     : -6,#11,                       \
                    'bc_x%end'                     : -6,#12                       \
                    'bc_y%beg'                     : -2,                      \
                    'bc_y%end'                     : -6,                      \
                    # ==========================================================
                                                                               \
                    

# Formatted Database Files Structure Parameters ============
                    'lsq_deriv'      : 'T',  \
                    'We_src'                       : 'T',
'rho_wrt':'T',
'gamma_wrt':'T',
'heat_ratio_wrt':'T',
'pi_inf_wrt':'T',
'pres_inf_wrt':'T',                    
'c_wrt' :'T',
'E_wrt':'T',
'num_probes':3,
'probe_wrt':'T',
'probe(1)%x': D/2/4,
'probe(1)%y': 0,
'probe(1)%z': 0,
'probe(2)%x': D/2/4*2,
'probe(2)%y': 0,
'probe(2)%z': 0,
'probe(3)%x': D/2/4*3,
'probe(3)%y': 0,
'probe(3)%z': 0,
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
                    'schlieren_wrt'                 : 'T',                    \
                    'com_wrt(1)':'T',\
                    'cb_wrt(1)':'T',\
		    'threshold_mf(1)' : 1,\
		    'threshold_mf(2)' : 0.50,\
                    'threshold_mf(3)' : 0.75,\
                    'threshold_mf(4)' : 0.99,\
                    'schlieren_alpha(1)'                 : 400,                 \
                    'schlieren_alpha(2)'                 : 40,                    \
                    'fd_order'                 : 4,                    \
   
		    'parallel_io'                  :'T',                       \
		    # ==========================================================
                                                                                
		    # Patch 1: Background  ============================
                    'patch_icpp(1)%geometry'       : 3,                        \
                    'patch_icpp(1)%x_centroid'     : 0,                  \
                    'patch_icpp(1)%y_centroid'     : 3*D,          \
                    'patch_icpp(1)%length_x'       : 14*D,                   \
                    'patch_icpp(1)%length_y'       : 6*D,         \
                    'patch_icpp(1)%vel(1)'         : 0.,                   \
                    'patch_icpp(1)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(1)%pres'           : 101325.,                   \
                    'patch_icpp(1)%alpha_rho(1)'   : 0.,                \
                    'patch_icpp(1)%alpha_rho(2)'   : 1.204,                   \
                    'patch_icpp(1)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(1)%alpha(2)'       : 1.E+00,                   \
                    # ==========================================================

		    # Patch 2: Shocked state ============================
                    'patch_icpp(2)%geometry'       : 3,                        \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(2)%x_centroid'     : -11*D,                  \
                    'patch_icpp(2)%y_centroid'     : 3*D,          \
                    'patch_icpp(2)%length_x'       : 21*D,                   \
                    'patch_icpp(2)%length_y'       : 6*D,         \
                    'patch_icpp(2)%vel(1)'         : vel,                   \
                    'patch_icpp(2)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(2)%pres'           : ps,                   \
                    'patch_icpp(2)%alpha_rho(1)'   : 0.E+00,                \
                    'patch_icpp(2)%alpha_rho(2)'   : rho_post_g,                   \
                    'patch_icpp(2)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(2)%alpha(2)'       : 1.E+00,                   \
                    # ==========================================================
			
                    # Patch 3: Bubble  ======================================
                    'patch_icpp(3)%geometry'       : 2,                        \
                    'patch_icpp(3)%x_centroid'     : 0,                 \
                    'patch_icpp(3)%y_centroid'     : 0,                  \
                    'patch_icpp(3)%radius'         : D/2,                  \
                    'patch_icpp(3)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(3)%vel(1)'         : 0.,                   \
                    'patch_icpp(3)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(3)%pres'           : 101325.,                   \
                    'patch_icpp(3)%alpha_rho(1)'   : rho_l,                   \
                    'patch_icpp(3)%alpha_rho(2)'   : 0,#0.061,                  \
                    'patch_icpp(3)%alpha(1)'       : 1.,# 0.95                  \
                    'patch_icpp(3)%alpha(2)'       : 0.,#0.05,                   \
    
                    # ==========================================================
 
		    # Fluids Physical Parameters ==============================         \
                    'fluid_pp(1)%gamma'            : 1.E+00/(gam_l-1.E+00),  \
                    'fluid_pp(1)%pi_inf'           : pi_l*gam_l/(gam_l-1.E+00), \
                    'fluid_pp(2)%gamma'            : 1.E+00/(gam_g-1.E+00),  \
                    'fluid_pp(2)%pi_inf'           : 0.E+00,                   \
	            # ==========================================================

    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
