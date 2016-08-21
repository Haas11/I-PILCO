/* Include files */

#include <stddef.h>
#include "blas.h"
#include "IPILCO_relativeRPY_S_sfun.h"
#include "c1_IPILCO_relativeRPY_S.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "IPILCO_relativeRPY_S_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c1_debug_family_names[39] = { "R_e", "R_d", "Md", "p_e",
  "p_d", "delta_p", "dp_e", "prevp_d", "ddelta_p", "ddp_d", "a_t", "Re_d",
  "EUL_de", "T", "dT_e", "Ta", "omega_e", "omega_d", "deltaOmega_de", "dEUL_de",
  "domega_d", "a_o", "nargin", "nargout", "Kp", "Kd", "He", "xe", "dxe_n", "F",
  "Hd", "prevHd", "prevdp_d", "prevT_e", "dt", "a", "Fa", "dp_d", "T_e" };

static const char * c1_b_debug_family_names[3] = { "c", "nargin", "nargout" };

static const char * c1_c_debug_family_names[3] = { "r", "nargin", "nargout" };

static const char * c1_d_debug_family_names[6] = { "n", "t", "nargin", "nargout",
  "T", "R" };

static const char * c1_e_debug_family_names[4] = { "h", "d", "nargin", "nargout"
};

static const char * c1_f_debug_family_names[4] = { "nargin", "nargout", "x",
  "t1" };

static const char * c1_g_debug_family_names[11] = { "opt", "s", "sr", "cr", "sp",
  "cp", "type", "nargin", "nargout", "m", "rpy" };

static const char * c1_h_debug_family_names[21] = { "X", "Y", "Z", "alfa",
  "beta", "gamma", "sa", "sb", "sc", "ca", "cb", "cc", "sig", "l", "nargin",
  "nargout", "q", "i", "j", "k", "R" };

/* Function Declarations */
static void initialize_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void initialize_params_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void enable_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void disable_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void set_sim_state_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance, const mxArray *c1_st);
static void finalize_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void sf_gateway_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void mdl_start_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c1_chartstep_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void initSimStructsc1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c1_tr2rt(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                     real_T c1_T[16], real_T c1_R[9]);
static void c1_transl(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c1_x[16], real_T c1_t1[3]);
static void c1_my_tr2rpy(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_m[9], real_T c1_rpy[3]);
static void c1_cardatol(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_q[3], real_T c1_i, real_T c1_j, real_T c1_k, real_T c1_R[9]);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static void c1_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_b_T_e, const char_T *c1_identifier, real_T
  c1_y[9]);
static void c1_b_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[9]);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_c_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_b_dp_d, const char_T *c1_identifier, real_T
  c1_y[3]);
static void c1_d_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[3]);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_e_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_b_Fa, const char_T *c1_identifier, real_T
  c1_y[6]);
static void c1_f_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[6]);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_g_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static real_T c1_g_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_h_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_h_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[3]);
static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_i_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[36]);
static void c1_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_j_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[16]);
static void c1_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_i_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_j_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_k_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_k_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[3]);
static void c1_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_l_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_m_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static c1_snY5FzvLGhd6UMWwSEwxMX c1_l_emlrt_marshallIn
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance, const mxArray *c1_u,
   const emlrtMsgIdentifier *c1_parentId);
static boolean_T c1_m_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_info_helper(const mxArray **c1_info);
static const mxArray *c1_emlrt_marshallOut(const char * c1_u);
static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u);
static void c1_b_info_helper(const mxArray **c1_info);
static void c1_c_info_helper(const mxArray **c1_info);
static void c1_d_info_helper(const mxArray **c1_info);
static void c1_eml_switch_helper(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static void c1_rdivide(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_x[3], real_T c1_y, real_T c1_z[3]);
static void c1_inv(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                   real_T c1_x[9], real_T c1_y[9]);
static void c1_inv3x3(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c1_x[9], real_T c1_y[9]);
static real_T c1_norm(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c1_x[9]);
static void c1_eml_warning(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static void c1_eml_eps(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c1_b_eml_warning(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, char_T c1_varargin_2[14]);
static void c1_eml_scalar_eg(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static void c1_threshold(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c1_b_eml_scalar_eg(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static real_T c1_atan2(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_y, real_T c1_x);
static void c1_c_eml_scalar_eg(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static void c1_b_rdivide(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_x[9], real_T c1_y, real_T c1_z[9]);
static void c1_eye(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                   real_T c1_I[9]);
static void c1_mldivide(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_A[9], real_T c1_B[3], real_T c1_Y[3]);
static real_T c1_c_rdivide(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_x, real_T c1_y);
static void c1_b_inv(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                     real_T c1_x[36], real_T c1_y[36]);
static void c1_realmin(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c1_eml_matlab_zgetrf(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_A[36], real_T c1_b_A[36], int32_T c1_ipiv[6],
  int32_T *c1_info);
static int32_T c1_eml_ixamax(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, int32_T c1_n, real_T c1_x[36], int32_T c1_ix0);
static void c1_check_forloop_overflow_error
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance, boolean_T c1_overflow);
static void c1_b_eml_switch_helper(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static void c1_eml_xgeru(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  int32_T c1_m, int32_T c1_n, real_T c1_alpha1, int32_T c1_ix0, int32_T c1_iy0,
  real_T c1_A[36], int32_T c1_ia0, real_T c1_b_A[36]);
static void c1_eml_ipiv2perm(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, int32_T c1_ipiv[6], int32_T c1_perm[6]);
static void c1_eml_xtrsm(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_A[36], real_T c1_B[36], real_T c1_b_B[36]);
static void c1_c_eml_switch_helper(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static real_T c1_b_norm(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_x[36]);
static void c1_d_eml_scalar_eg(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static void c1_n_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_sprintf, const char_T *c1_identifier, char_T
  c1_y[14]);
static void c1_o_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  char_T c1_y[14]);
static const mxArray *c1_n_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_p_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_q_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_b_is_active_c1_IPILCO_relativeRPY_S, const
  char_T *c1_identifier);
static uint8_T c1_r_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_eml_matlab_zgetrf(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_A[36], int32_T c1_ipiv[6], int32_T *c1_info);
static void c1_b_eml_xgeru(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, int32_T c1_m, int32_T c1_n, real_T c1_alpha1, int32_T c1_ix0,
  int32_T c1_iy0, real_T c1_A[36], int32_T c1_ia0);
static void c1_b_eml_xtrsm(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_A[36], real_T c1_B[36]);
static void init_dsm_address_info(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static void init_simulink_io_address(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c1_is_active_c1_IPILCO_relativeRPY_S = 0U;
}

static void initialize_params_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c1_update_debugger_state_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  int32_T c1_i0;
  real_T c1_u[6];
  const mxArray *c1_b_y = NULL;
  int32_T c1_i1;
  real_T c1_b_u[9];
  const mxArray *c1_c_y = NULL;
  int32_T c1_i2;
  real_T c1_c_u[6];
  const mxArray *c1_d_y = NULL;
  int32_T c1_i3;
  real_T c1_d_u[3];
  const mxArray *c1_e_y = NULL;
  uint8_T c1_hoistedGlobal;
  uint8_T c1_e_u;
  const mxArray *c1_f_y = NULL;
  c1_st = NULL;
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellmatrix(5, 1), false);
  for (c1_i0 = 0; c1_i0 < 6; c1_i0++) {
    c1_u[c1_i0] = (*chartInstance->c1_Fa)[c1_i0];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 6, 1), false);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  for (c1_i1 = 0; c1_i1 < 9; c1_i1++) {
    c1_b_u[c1_i1] = (*chartInstance->c1_T_e)[c1_i1];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 3, 3),
                false);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  for (c1_i2 = 0; c1_i2 < 6; c1_i2++) {
    c1_c_u[c1_i2] = (*chartInstance->c1_a)[c1_i2];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_c_u, 0, 0U, 1U, 0U, 2, 6, 1),
                false);
  sf_mex_setcell(c1_y, 2, c1_d_y);
  for (c1_i3 = 0; c1_i3 < 3; c1_i3++) {
    c1_d_u[c1_i3] = (*chartInstance->c1_dp_d)[c1_i3];
  }

  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", c1_d_u, 0, 0U, 1U, 0U, 2, 3, 1),
                false);
  sf_mex_setcell(c1_y, 3, c1_e_y);
  c1_hoistedGlobal = chartInstance->c1_is_active_c1_IPILCO_relativeRPY_S;
  c1_e_u = c1_hoistedGlobal;
  c1_f_y = NULL;
  sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_e_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 4, c1_f_y);
  sf_mex_assign(&c1_st, c1_y, false);
  return c1_st;
}

static void set_sim_state_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance, const mxArray *c1_st)
{
  const mxArray *c1_u;
  real_T c1_dv0[6];
  int32_T c1_i4;
  real_T c1_dv1[9];
  int32_T c1_i5;
  real_T c1_dv2[6];
  int32_T c1_i6;
  real_T c1_dv3[3];
  int32_T c1_i7;
  chartInstance->c1_doneDoubleBufferReInit = true;
  c1_u = sf_mex_dup(c1_st);
  c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 0)), "Fa",
                        c1_dv0);
  for (c1_i4 = 0; c1_i4 < 6; c1_i4++) {
    (*chartInstance->c1_Fa)[c1_i4] = c1_dv0[c1_i4];
  }

  c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 1)), "T_e",
                      c1_dv1);
  for (c1_i5 = 0; c1_i5 < 9; c1_i5++) {
    (*chartInstance->c1_T_e)[c1_i5] = c1_dv1[c1_i5];
  }

  c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 2)), "a",
                        c1_dv2);
  for (c1_i6 = 0; c1_i6 < 6; c1_i6++) {
    (*chartInstance->c1_a)[c1_i6] = c1_dv2[c1_i6];
  }

  c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 3)),
                        "dp_d", c1_dv3);
  for (c1_i7 = 0; c1_i7 < 3; c1_i7++) {
    (*chartInstance->c1_dp_d)[c1_i7] = c1_dv3[c1_i7];
  }

  chartInstance->c1_is_active_c1_IPILCO_relativeRPY_S = c1_q_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 4)),
     "is_active_c1_IPILCO_relativeRPY_S");
  sf_mex_destroy(&c1_u);
  c1_update_debugger_state_c1_IPILCO_relativeRPY_S(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  int32_T c1_i8;
  int32_T c1_i9;
  int32_T c1_i10;
  int32_T c1_i11;
  int32_T c1_i12;
  int32_T c1_i13;
  int32_T c1_i14;
  int32_T c1_i15;
  int32_T c1_i16;
  int32_T c1_i17;
  int32_T c1_i18;
  int32_T c1_i19;
  int32_T c1_i20;
  int32_T c1_i21;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_IPILCO_relativeRPY_S(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_IPILCO_relativeRPY_SMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c1_i8 = 0; c1_i8 < 6; c1_i8++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_a)[c1_i8], 0U);
  }

  for (c1_i9 = 0; c1_i9 < 6; c1_i9++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_Fa)[c1_i9], 1U);
  }

  for (c1_i10 = 0; c1_i10 < 3; c1_i10++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_dp_d)[c1_i10], 2U);
  }

  for (c1_i11 = 0; c1_i11 < 36; c1_i11++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_Kp)[c1_i11], 3U);
  }

  for (c1_i12 = 0; c1_i12 < 36; c1_i12++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_Kd)[c1_i12], 4U);
  }

  for (c1_i13 = 0; c1_i13 < 16; c1_i13++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_He)[c1_i13], 5U);
  }

  for (c1_i14 = 0; c1_i14 < 6; c1_i14++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_xe)[c1_i14], 6U);
  }

  for (c1_i15 = 0; c1_i15 < 6; c1_i15++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_dxe_n)[c1_i15], 7U);
  }

  for (c1_i16 = 0; c1_i16 < 6; c1_i16++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_F)[c1_i16], 8U);
  }

  for (c1_i17 = 0; c1_i17 < 16; c1_i17++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_Hd)[c1_i17], 9U);
  }

  for (c1_i18 = 0; c1_i18 < 16; c1_i18++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_prevHd)[c1_i18], 10U);
  }

  for (c1_i19 = 0; c1_i19 < 3; c1_i19++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_prevdp_d)[c1_i19], 11U);
  }

  for (c1_i20 = 0; c1_i20 < 9; c1_i20++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_prevT_e)[c1_i20], 12U);
  }

  for (c1_i21 = 0; c1_i21 < 9; c1_i21++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_T_e)[c1_i21], 13U);
  }

  _SFD_DATA_RANGE_CHECK(*chartInstance->c1_dt, 14U);
}

static void mdl_start_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_chartstep_c1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  real_T c1_hoistedGlobal;
  int32_T c1_i22;
  real_T c1_b_Kp[36];
  int32_T c1_i23;
  real_T c1_b_Kd[36];
  int32_T c1_i24;
  real_T c1_b_He[16];
  int32_T c1_i25;
  real_T c1_b_xe[6];
  int32_T c1_i26;
  real_T c1_b_dxe_n[6];
  int32_T c1_i27;
  real_T c1_b_F[6];
  int32_T c1_i28;
  real_T c1_b_Hd[16];
  int32_T c1_i29;
  real_T c1_b_prevHd[16];
  int32_T c1_i30;
  real_T c1_b_prevdp_d[3];
  int32_T c1_i31;
  real_T c1_b_prevT_e[9];
  real_T c1_b_dt;
  uint32_T c1_debug_family_var_map[39];
  real_T c1_R_e[9];
  real_T c1_R_d[9];
  real_T c1_Md[36];
  real_T c1_p_e[3];
  real_T c1_p_d[3];
  real_T c1_delta_p[3];
  real_T c1_dp_e[3];
  real_T c1_prevp_d[3];
  real_T c1_ddelta_p[3];
  real_T c1_ddp_d[3];
  real_T c1_a_t[3];
  real_T c1_Re_d[9];
  real_T c1_EUL_de[3];
  real_T c1_T[9];
  real_T c1_dT_e[9];
  real_T c1_Ta[36];
  real_T c1_omega_e[3];
  real_T c1_omega_d[3];
  real_T c1_deltaOmega_de[3];
  real_T c1_dEUL_de[3];
  real_T c1_domega_d[3];
  real_T c1_a_o[3];
  real_T c1_nargin = 11.0;
  real_T c1_nargout = 4.0;
  real_T c1_b_a[6];
  real_T c1_b_Fa[6];
  real_T c1_b_dp_d[3];
  real_T c1_b_T_e[9];
  int32_T c1_i32;
  real_T c1_c_He[16];
  real_T c1_dv4[9];
  int32_T c1_i33;
  int32_T c1_i34;
  real_T c1_c_Hd[16];
  real_T c1_dv5[9];
  int32_T c1_i35;
  int32_T c1_i36;
  static real_T c1_dv6[36] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  int32_T c1_i37;
  int32_T c1_i38;
  real_T c1_d_Hd[16];
  real_T c1_dv7[3];
  int32_T c1_i39;
  int32_T c1_i40;
  int32_T c1_i41;
  int32_T c1_i42;
  real_T c1_c_prevHd[16];
  real_T c1_dv8[3];
  int32_T c1_i43;
  int32_T c1_i44;
  real_T c1_b_p_d[3];
  real_T c1_b[3];
  int32_T c1_i45;
  int32_T c1_i46;
  int32_T c1_i47;
  real_T c1_c_dp_d[3];
  real_T c1_dv9[3];
  int32_T c1_i48;
  int32_T c1_i49;
  int32_T c1_i50;
  int32_T c1_i51;
  int32_T c1_i52;
  real_T c1_c_a[9];
  int32_T c1_i53;
  int32_T c1_i54;
  real_T c1_y[3];
  int32_T c1_i55;
  int32_T c1_i56;
  int32_T c1_i57;
  int32_T c1_i58;
  int32_T c1_i59;
  int32_T c1_i60;
  int32_T c1_i61;
  int32_T c1_i62;
  real_T c1_b_y[3];
  int32_T c1_i63;
  int32_T c1_i64;
  int32_T c1_i65;
  static real_T c1_dv10[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c1_dv11[9];
  int32_T c1_i66;
  int32_T c1_i67;
  int32_T c1_i68;
  int32_T c1_i69;
  int32_T c1_i70;
  int32_T c1_i71;
  int32_T c1_i72;
  int32_T c1_i73;
  int32_T c1_i74;
  int32_T c1_i75;
  real_T c1_b_b[9];
  int32_T c1_i76;
  int32_T c1_i77;
  int32_T c1_i78;
  real_T c1_C[9];
  int32_T c1_i79;
  int32_T c1_i80;
  int32_T c1_i81;
  int32_T c1_i82;
  int32_T c1_i83;
  int32_T c1_i84;
  int32_T c1_i85;
  int32_T c1_i86;
  int32_T c1_i87;
  real_T c1_b_Re_d[9];
  real_T c1_dv12[3];
  int32_T c1_i88;
  int32_T c1_i89;
  real_T c1_b_EUL_de[3];
  real_T c1_dv13[9];
  int32_T c1_i90;
  int32_T c1_i91;
  int32_T c1_i92;
  int32_T c1_i93;
  int32_T c1_i94;
  int32_T c1_i95;
  int32_T c1_i96;
  int32_T c1_i97;
  int32_T c1_i98;
  int32_T c1_i99;
  int32_T c1_i100;
  int32_T c1_i101;
  int32_T c1_i102;
  int32_T c1_i103;
  int32_T c1_i104;
  real_T c1_c_T_e[9];
  real_T c1_dv14[9];
  int32_T c1_i105;
  int32_T c1_i106;
  int32_T c1_i107;
  int32_T c1_i108;
  int32_T c1_i109;
  int32_T c1_i110;
  int32_T c1_i111;
  int32_T c1_i112;
  int32_T c1_i113;
  int32_T c1_i114;
  int32_T c1_i115;
  int32_T c1_i116;
  int32_T c1_i117;
  int32_T c1_i118;
  int32_T c1_i119;
  int32_T c1_i120;
  int32_T c1_i121;
  int32_T c1_i122;
  int32_T c1_i123;
  real_T c1_d_T_e[9];
  int32_T c1_i124;
  real_T c1_b_deltaOmega_de[3];
  real_T c1_dv15[3];
  int32_T c1_i125;
  int32_T c1_i126;
  int32_T c1_i127;
  int32_T c1_i128;
  int32_T c1_i129;
  int32_T c1_i130;
  int32_T c1_i131;
  int32_T c1_i132;
  real_T c1_dv16[9];
  int32_T c1_i133;
  int32_T c1_i134;
  int32_T c1_i135;
  int32_T c1_i136;
  real_T c1_c_y[9];
  int32_T c1_i137;
  int32_T c1_i138;
  int32_T c1_i139;
  int32_T c1_i140;
  int32_T c1_i141;
  int32_T c1_i142;
  int32_T c1_i143;
  int32_T c1_i144;
  int32_T c1_i145;
  int32_T c1_i146;
  int32_T c1_i147;
  int32_T c1_i148;
  int32_T c1_i149;
  int32_T c1_i150;
  int32_T c1_i151;
  int32_T c1_i152;
  real_T c1_d_y[3];
  int32_T c1_i153;
  int32_T c1_i154;
  int32_T c1_i155;
  real_T c1_b_T[9];
  int32_T c1_i156;
  int32_T c1_i157;
  int32_T c1_i158;
  int32_T c1_i159;
  int32_T c1_i160;
  int32_T c1_i161;
  real_T c1_e_y[3];
  int32_T c1_i162;
  int32_T c1_i163;
  int32_T c1_i164;
  int32_T c1_i165;
  int32_T c1_i166;
  int32_T c1_i167;
  int32_T c1_i168;
  int32_T c1_i169;
  real_T c1_b_Ta[36];
  real_T c1_dv17[36];
  int32_T c1_i170;
  int32_T c1_i171;
  int32_T c1_i172;
  int32_T c1_i173;
  real_T c1_d_a[36];
  int32_T c1_i174;
  real_T c1_c_b[6];
  int32_T c1_i175;
  real_T c1_f_y[6];
  int32_T c1_i176;
  int32_T c1_i177;
  int32_T c1_i178;
  int32_T c1_i179;
  int32_T c1_i180;
  int32_T c1_i181;
  int32_T c1_i182;
  int32_T c1_i183;
  int32_T c1_i184;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  c1_hoistedGlobal = *chartInstance->c1_dt;
  for (c1_i22 = 0; c1_i22 < 36; c1_i22++) {
    c1_b_Kp[c1_i22] = (*chartInstance->c1_Kp)[c1_i22];
  }

  for (c1_i23 = 0; c1_i23 < 36; c1_i23++) {
    c1_b_Kd[c1_i23] = (*chartInstance->c1_Kd)[c1_i23];
  }

  for (c1_i24 = 0; c1_i24 < 16; c1_i24++) {
    c1_b_He[c1_i24] = (*chartInstance->c1_He)[c1_i24];
  }

  for (c1_i25 = 0; c1_i25 < 6; c1_i25++) {
    c1_b_xe[c1_i25] = (*chartInstance->c1_xe)[c1_i25];
  }

  for (c1_i26 = 0; c1_i26 < 6; c1_i26++) {
    c1_b_dxe_n[c1_i26] = (*chartInstance->c1_dxe_n)[c1_i26];
  }

  for (c1_i27 = 0; c1_i27 < 6; c1_i27++) {
    c1_b_F[c1_i27] = (*chartInstance->c1_F)[c1_i27];
  }

  for (c1_i28 = 0; c1_i28 < 16; c1_i28++) {
    c1_b_Hd[c1_i28] = (*chartInstance->c1_Hd)[c1_i28];
  }

  for (c1_i29 = 0; c1_i29 < 16; c1_i29++) {
    c1_b_prevHd[c1_i29] = (*chartInstance->c1_prevHd)[c1_i29];
  }

  for (c1_i30 = 0; c1_i30 < 3; c1_i30++) {
    c1_b_prevdp_d[c1_i30] = (*chartInstance->c1_prevdp_d)[c1_i30];
  }

  for (c1_i31 = 0; c1_i31 < 9; c1_i31++) {
    c1_b_prevT_e[c1_i31] = (*chartInstance->c1_prevT_e)[c1_i31];
  }

  c1_b_dt = c1_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 39U, 39U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_R_e, 0U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_R_d, 1U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_Md, 2U, c1_g_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_p_e, 3U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_p_d, 4U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_delta_p, 5U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_dp_e, 6U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_prevp_d, 7U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_ddelta_p, 8U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_ddp_d, 9U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_a_t, 10U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_Re_d, 11U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_EUL_de, 12U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_T, 13U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_dT_e, 14U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_Ta, 15U, c1_g_sf_marshallOut,
    c1_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_omega_e, 16U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_omega_d, 17U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_deltaOmega_de, 18U,
    c1_h_sf_marshallOut, c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_dEUL_de, 19U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_domega_d, 20U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_a_o, 21U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 22U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 23U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_Kp, 24U, c1_g_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_Kd, 25U, c1_g_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_He, 26U, c1_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_xe, 27U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_dxe_n, 28U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_F, 29U, c1_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_Hd, 30U, c1_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_prevHd, 31U, c1_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_prevdp_d, 32U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_prevT_e, 33U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_b_dt, 34U, c1_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_a, 35U, c1_c_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_Fa, 36U, c1_c_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_dp_d, 37U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_T_e, 38U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 28);
  for (c1_i32 = 0; c1_i32 < 16; c1_i32++) {
    c1_c_He[c1_i32] = c1_b_He[c1_i32];
  }

  c1_tr2rt(chartInstance, c1_c_He, c1_dv4);
  for (c1_i33 = 0; c1_i33 < 9; c1_i33++) {
    c1_R_e[c1_i33] = c1_dv4[c1_i33];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 29);
  for (c1_i34 = 0; c1_i34 < 16; c1_i34++) {
    c1_c_Hd[c1_i34] = c1_b_Hd[c1_i34];
  }

  c1_tr2rt(chartInstance, c1_c_Hd, c1_dv5);
  for (c1_i35 = 0; c1_i35 < 9; c1_i35++) {
    c1_R_d[c1_i35] = c1_dv5[c1_i35];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 31);
  for (c1_i36 = 0; c1_i36 < 36; c1_i36++) {
    c1_Md[c1_i36] = c1_dv6[c1_i36];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 34);
  for (c1_i37 = 0; c1_i37 < 3; c1_i37++) {
    c1_p_e[c1_i37] = c1_b_xe[c1_i37];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 35);
  for (c1_i38 = 0; c1_i38 < 16; c1_i38++) {
    c1_d_Hd[c1_i38] = c1_b_Hd[c1_i38];
  }

  c1_transl(chartInstance, c1_d_Hd, c1_dv7);
  for (c1_i39 = 0; c1_i39 < 3; c1_i39++) {
    c1_p_d[c1_i39] = c1_dv7[c1_i39];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 36);
  for (c1_i40 = 0; c1_i40 < 3; c1_i40++) {
    c1_delta_p[c1_i40] = c1_p_d[c1_i40] - c1_p_e[c1_i40];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 38);
  for (c1_i41 = 0; c1_i41 < 3; c1_i41++) {
    c1_dp_e[c1_i41] = c1_b_dxe_n[c1_i41];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 39);
  for (c1_i42 = 0; c1_i42 < 16; c1_i42++) {
    c1_c_prevHd[c1_i42] = c1_b_prevHd[c1_i42];
  }

  c1_transl(chartInstance, c1_c_prevHd, c1_dv8);
  for (c1_i43 = 0; c1_i43 < 3; c1_i43++) {
    c1_prevp_d[c1_i43] = c1_dv8[c1_i43];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 40);
  for (c1_i44 = 0; c1_i44 < 3; c1_i44++) {
    c1_b_p_d[c1_i44] = c1_p_d[c1_i44] - c1_prevp_d[c1_i44];
  }

  c1_rdivide(chartInstance, c1_b_p_d, c1_b_dt, c1_b);
  for (c1_i45 = 0; c1_i45 < 3; c1_i45++) {
    c1_b_dp_d[c1_i45] = c1_b[c1_i45];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 41);
  for (c1_i46 = 0; c1_i46 < 3; c1_i46++) {
    c1_ddelta_p[c1_i46] = c1_b_dp_d[c1_i46] - c1_dp_e[c1_i46];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 43);
  for (c1_i47 = 0; c1_i47 < 3; c1_i47++) {
    c1_c_dp_d[c1_i47] = c1_b_dp_d[c1_i47] - c1_b_prevdp_d[c1_i47];
  }

  c1_rdivide(chartInstance, c1_c_dp_d, c1_b_dt, c1_dv9);
  for (c1_i48 = 0; c1_i48 < 3; c1_i48++) {
    c1_ddp_d[c1_i48] = c1_dv9[c1_i48];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 45);
  c1_i49 = 0;
  c1_i50 = 0;
  for (c1_i51 = 0; c1_i51 < 3; c1_i51++) {
    for (c1_i52 = 0; c1_i52 < 3; c1_i52++) {
      c1_c_a[c1_i52 + c1_i49] = c1_b_Kd[c1_i52 + c1_i50];
    }

    c1_i49 += 3;
    c1_i50 += 6;
  }

  for (c1_i53 = 0; c1_i53 < 3; c1_i53++) {
    c1_b[c1_i53] = c1_ddelta_p[c1_i53];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i54 = 0; c1_i54 < 3; c1_i54++) {
    c1_y[c1_i54] = 0.0;
    c1_i55 = 0;
    for (c1_i56 = 0; c1_i56 < 3; c1_i56++) {
      c1_y[c1_i54] += c1_c_a[c1_i55 + c1_i54] * c1_b[c1_i56];
      c1_i55 += 3;
    }
  }

  c1_i57 = 0;
  c1_i58 = 0;
  for (c1_i59 = 0; c1_i59 < 3; c1_i59++) {
    for (c1_i60 = 0; c1_i60 < 3; c1_i60++) {
      c1_c_a[c1_i60 + c1_i57] = c1_b_Kp[c1_i60 + c1_i58];
    }

    c1_i57 += 3;
    c1_i58 += 6;
  }

  for (c1_i61 = 0; c1_i61 < 3; c1_i61++) {
    c1_b[c1_i61] = c1_delta_p[c1_i61];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i62 = 0; c1_i62 < 3; c1_i62++) {
    c1_b_y[c1_i62] = 0.0;
    c1_i63 = 0;
    for (c1_i64 = 0; c1_i64 < 3; c1_i64++) {
      c1_b_y[c1_i62] += c1_c_a[c1_i63 + c1_i62] * c1_b[c1_i64];
      c1_i63 += 3;
    }
  }

  for (c1_i65 = 0; c1_i65 < 9; c1_i65++) {
    c1_dv11[c1_i65] = c1_dv10[c1_i65];
  }

  c1_inv(chartInstance, c1_dv11, c1_c_a);
  for (c1_i66 = 0; c1_i66 < 3; c1_i66++) {
    c1_b[c1_i66] = (c1_y[c1_i66] + c1_b_y[c1_i66]) + c1_b_F[c1_i66];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i67 = 0; c1_i67 < 3; c1_i67++) {
    c1_y[c1_i67] = 0.0;
    c1_i68 = 0;
    for (c1_i69 = 0; c1_i69 < 3; c1_i69++) {
      c1_y[c1_i67] += c1_c_a[c1_i68 + c1_i67] * c1_b[c1_i69];
      c1_i68 += 3;
    }
  }

  for (c1_i70 = 0; c1_i70 < 3; c1_i70++) {
    c1_a_t[c1_i70] = c1_ddp_d[c1_i70] + c1_y[c1_i70];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 52);
  c1_i71 = 0;
  for (c1_i72 = 0; c1_i72 < 3; c1_i72++) {
    c1_i73 = 0;
    for (c1_i74 = 0; c1_i74 < 3; c1_i74++) {
      c1_c_a[c1_i74 + c1_i71] = c1_R_e[c1_i73 + c1_i72];
      c1_i73 += 3;
    }

    c1_i71 += 3;
  }

  for (c1_i75 = 0; c1_i75 < 9; c1_i75++) {
    c1_b_b[c1_i75] = c1_R_d[c1_i75];
  }

  c1_b_eml_scalar_eg(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  for (c1_i76 = 0; c1_i76 < 9; c1_i76++) {
    c1_Re_d[c1_i76] = 0.0;
  }

  for (c1_i77 = 0; c1_i77 < 9; c1_i77++) {
    c1_Re_d[c1_i77] = 0.0;
  }

  for (c1_i78 = 0; c1_i78 < 9; c1_i78++) {
    c1_C[c1_i78] = c1_Re_d[c1_i78];
  }

  for (c1_i79 = 0; c1_i79 < 9; c1_i79++) {
    c1_Re_d[c1_i79] = c1_C[c1_i79];
  }

  c1_threshold(chartInstance);
  for (c1_i80 = 0; c1_i80 < 9; c1_i80++) {
    c1_C[c1_i80] = c1_Re_d[c1_i80];
  }

  for (c1_i81 = 0; c1_i81 < 9; c1_i81++) {
    c1_Re_d[c1_i81] = c1_C[c1_i81];
  }

  for (c1_i82 = 0; c1_i82 < 3; c1_i82++) {
    c1_i83 = 0;
    for (c1_i84 = 0; c1_i84 < 3; c1_i84++) {
      c1_Re_d[c1_i83 + c1_i82] = 0.0;
      c1_i85 = 0;
      for (c1_i86 = 0; c1_i86 < 3; c1_i86++) {
        c1_Re_d[c1_i83 + c1_i82] += c1_c_a[c1_i85 + c1_i82] * c1_b_b[c1_i86 +
          c1_i83];
        c1_i85 += 3;
      }

      c1_i83 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 53);
  for (c1_i87 = 0; c1_i87 < 9; c1_i87++) {
    c1_b_Re_d[c1_i87] = c1_Re_d[c1_i87];
  }

  c1_my_tr2rpy(chartInstance, c1_b_Re_d, c1_dv12);
  for (c1_i88 = 0; c1_i88 < 3; c1_i88++) {
    c1_EUL_de[c1_i88] = c1_dv12[c1_i88];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 55);
  for (c1_i89 = 0; c1_i89 < 3; c1_i89++) {
    c1_b_EUL_de[c1_i89] = c1_EUL_de[c1_i89];
  }

  c1_cardatol(chartInstance, c1_b_EUL_de, 1.0, 2.0, 3.0, c1_dv13);
  for (c1_i90 = 0; c1_i90 < 9; c1_i90++) {
    c1_T[c1_i90] = c1_dv13[c1_i90];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 56);
  for (c1_i91 = 0; c1_i91 < 9; c1_i91++) {
    c1_c_a[c1_i91] = c1_R_e[c1_i91];
  }

  for (c1_i92 = 0; c1_i92 < 9; c1_i92++) {
    c1_b_b[c1_i92] = c1_T[c1_i92];
  }

  c1_b_eml_scalar_eg(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  for (c1_i93 = 0; c1_i93 < 9; c1_i93++) {
    c1_b_T_e[c1_i93] = 0.0;
  }

  for (c1_i94 = 0; c1_i94 < 9; c1_i94++) {
    c1_b_T_e[c1_i94] = 0.0;
  }

  for (c1_i95 = 0; c1_i95 < 9; c1_i95++) {
    c1_C[c1_i95] = c1_b_T_e[c1_i95];
  }

  for (c1_i96 = 0; c1_i96 < 9; c1_i96++) {
    c1_b_T_e[c1_i96] = c1_C[c1_i96];
  }

  c1_threshold(chartInstance);
  for (c1_i97 = 0; c1_i97 < 9; c1_i97++) {
    c1_C[c1_i97] = c1_b_T_e[c1_i97];
  }

  for (c1_i98 = 0; c1_i98 < 9; c1_i98++) {
    c1_b_T_e[c1_i98] = c1_C[c1_i98];
  }

  for (c1_i99 = 0; c1_i99 < 3; c1_i99++) {
    c1_i100 = 0;
    for (c1_i101 = 0; c1_i101 < 3; c1_i101++) {
      c1_b_T_e[c1_i100 + c1_i99] = 0.0;
      c1_i102 = 0;
      for (c1_i103 = 0; c1_i103 < 3; c1_i103++) {
        c1_b_T_e[c1_i100 + c1_i99] += c1_c_a[c1_i102 + c1_i99] * c1_b_b[c1_i103
          + c1_i100];
        c1_i102 += 3;
      }

      c1_i100 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 57);
  for (c1_i104 = 0; c1_i104 < 9; c1_i104++) {
    c1_c_T_e[c1_i104] = c1_b_T_e[c1_i104] - c1_b_prevT_e[c1_i104];
  }

  c1_b_rdivide(chartInstance, c1_c_T_e, c1_b_dt, c1_dv14);
  for (c1_i105 = 0; c1_i105 < 9; c1_i105++) {
    c1_dT_e[c1_i105] = c1_dv14[c1_i105];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 58);
  c1_eye(chartInstance, c1_C);
  c1_i106 = 0;
  c1_i107 = 0;
  for (c1_i108 = 0; c1_i108 < 3; c1_i108++) {
    for (c1_i109 = 0; c1_i109 < 3; c1_i109++) {
      c1_Ta[c1_i109 + c1_i106] = c1_C[c1_i109 + c1_i107];
    }

    c1_i106 += 6;
    c1_i107 += 3;
  }

  c1_i110 = 0;
  for (c1_i111 = 0; c1_i111 < 3; c1_i111++) {
    for (c1_i112 = 0; c1_i112 < 3; c1_i112++) {
      c1_Ta[(c1_i112 + c1_i110) + 18] = 0.0;
    }

    c1_i110 += 6;
  }

  c1_i113 = 0;
  for (c1_i114 = 0; c1_i114 < 3; c1_i114++) {
    for (c1_i115 = 0; c1_i115 < 3; c1_i115++) {
      c1_Ta[(c1_i115 + c1_i113) + 3] = 0.0;
    }

    c1_i113 += 6;
  }

  c1_i116 = 0;
  c1_i117 = 0;
  for (c1_i118 = 0; c1_i118 < 3; c1_i118++) {
    for (c1_i119 = 0; c1_i119 < 3; c1_i119++) {
      c1_Ta[(c1_i119 + c1_i116) + 21] = c1_T[c1_i119 + c1_i117];
    }

    c1_i116 += 6;
    c1_i117 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 60);
  for (c1_i120 = 0; c1_i120 < 3; c1_i120++) {
    c1_omega_e[c1_i120] = c1_b_dxe_n[c1_i120 + 3];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 61);
  for (c1_i121 = 0; c1_i121 < 3; c1_i121++) {
    c1_omega_d[c1_i121] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 62);
  for (c1_i122 = 0; c1_i122 < 3; c1_i122++) {
    c1_deltaOmega_de[c1_i122] = c1_omega_d[c1_i122] - c1_omega_e[c1_i122];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 63);
  for (c1_i123 = 0; c1_i123 < 9; c1_i123++) {
    c1_d_T_e[c1_i123] = c1_b_T_e[c1_i123];
  }

  for (c1_i124 = 0; c1_i124 < 3; c1_i124++) {
    c1_b_deltaOmega_de[c1_i124] = c1_deltaOmega_de[c1_i124];
  }

  c1_mldivide(chartInstance, c1_d_T_e, c1_b_deltaOmega_de, c1_dv15);
  for (c1_i125 = 0; c1_i125 < 3; c1_i125++) {
    c1_dEUL_de[c1_i125] = c1_dv15[c1_i125];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 65);
  for (c1_i126 = 0; c1_i126 < 3; c1_i126++) {
    c1_domega_d[c1_i126] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 66);
  for (c1_i127 = 0; c1_i127 < 9; c1_i127++) {
    c1_c_a[c1_i127] = c1_dT_e[c1_i127];
  }

  for (c1_i128 = 0; c1_i128 < 3; c1_i128++) {
    c1_b[c1_i128] = c1_dEUL_de[c1_i128];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i129 = 0; c1_i129 < 3; c1_i129++) {
    c1_y[c1_i129] = 0.0;
    c1_i130 = 0;
    for (c1_i131 = 0; c1_i131 < 3; c1_i131++) {
      c1_y[c1_i129] += c1_c_a[c1_i130 + c1_i129] * c1_b[c1_i131];
      c1_i130 += 3;
    }
  }

  for (c1_i132 = 0; c1_i132 < 9; c1_i132++) {
    c1_dv16[c1_i132] = c1_dv10[c1_i132];
  }

  c1_inv(chartInstance, c1_dv16, c1_c_a);
  for (c1_i133 = 0; c1_i133 < 9; c1_i133++) {
    c1_b_b[c1_i133] = c1_b_T_e[c1_i133];
  }

  c1_b_eml_scalar_eg(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i134 = 0; c1_i134 < 3; c1_i134++) {
    c1_i135 = 0;
    for (c1_i136 = 0; c1_i136 < 3; c1_i136++) {
      c1_c_y[c1_i135 + c1_i134] = 0.0;
      c1_i137 = 0;
      for (c1_i138 = 0; c1_i138 < 3; c1_i138++) {
        c1_c_y[c1_i135 + c1_i134] += c1_c_a[c1_i137 + c1_i134] * c1_b_b[c1_i138
          + c1_i135];
        c1_i137 += 3;
      }

      c1_i135 += 3;
    }
  }

  c1_i139 = 0;
  c1_i140 = 0;
  for (c1_i141 = 0; c1_i141 < 3; c1_i141++) {
    for (c1_i142 = 0; c1_i142 < 3; c1_i142++) {
      c1_c_a[c1_i142 + c1_i139] = c1_b_Kd[(c1_i142 + c1_i140) + 21];
    }

    c1_i139 += 3;
    c1_i140 += 6;
  }

  for (c1_i143 = 0; c1_i143 < 3; c1_i143++) {
    c1_b[c1_i143] = c1_dEUL_de[c1_i143];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i144 = 0; c1_i144 < 3; c1_i144++) {
    c1_b_y[c1_i144] = 0.0;
    c1_i145 = 0;
    for (c1_i146 = 0; c1_i146 < 3; c1_i146++) {
      c1_b_y[c1_i144] += c1_c_a[c1_i145 + c1_i144] * c1_b[c1_i146];
      c1_i145 += 3;
    }
  }

  c1_i147 = 0;
  c1_i148 = 0;
  for (c1_i149 = 0; c1_i149 < 3; c1_i149++) {
    for (c1_i150 = 0; c1_i150 < 3; c1_i150++) {
      c1_c_a[c1_i150 + c1_i147] = c1_b_Kp[(c1_i150 + c1_i148) + 21];
    }

    c1_i147 += 3;
    c1_i148 += 6;
  }

  for (c1_i151 = 0; c1_i151 < 3; c1_i151++) {
    c1_b[c1_i151] = c1_EUL_de[c1_i151];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i152 = 0; c1_i152 < 3; c1_i152++) {
    c1_d_y[c1_i152] = 0.0;
    c1_i153 = 0;
    for (c1_i154 = 0; c1_i154 < 3; c1_i154++) {
      c1_d_y[c1_i152] += c1_c_a[c1_i153 + c1_i152] * c1_b[c1_i154];
      c1_i153 += 3;
    }
  }

  for (c1_i155 = 0; c1_i155 < 9; c1_i155++) {
    c1_b_T[c1_i155] = c1_T[c1_i155];
  }

  c1_inv(chartInstance, c1_b_T, c1_C);
  c1_i156 = 0;
  for (c1_i157 = 0; c1_i157 < 3; c1_i157++) {
    c1_i158 = 0;
    for (c1_i159 = 0; c1_i159 < 3; c1_i159++) {
      c1_c_a[c1_i159 + c1_i156] = c1_C[c1_i158 + c1_i157];
      c1_i158 += 3;
    }

    c1_i156 += 3;
  }

  for (c1_i160 = 0; c1_i160 < 3; c1_i160++) {
    c1_b[c1_i160] = c1_b_F[c1_i160 + 3];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i161 = 0; c1_i161 < 3; c1_i161++) {
    c1_e_y[c1_i161] = 0.0;
    c1_i162 = 0;
    for (c1_i163 = 0; c1_i163 < 3; c1_i163++) {
      c1_e_y[c1_i161] += c1_c_a[c1_i162 + c1_i161] * c1_b[c1_i163];
      c1_i162 += 3;
    }
  }

  for (c1_i164 = 0; c1_i164 < 3; c1_i164++) {
    c1_b_y[c1_i164] = (c1_b_y[c1_i164] + c1_d_y[c1_i164]) + c1_e_y[c1_i164];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i165 = 0; c1_i165 < 3; c1_i165++) {
    c1_d_y[c1_i165] = 0.0;
    c1_i166 = 0;
    for (c1_i167 = 0; c1_i167 < 3; c1_i167++) {
      c1_d_y[c1_i165] += c1_c_y[c1_i166 + c1_i165] * c1_b_y[c1_i167];
      c1_i166 += 3;
    }
  }

  for (c1_i168 = 0; c1_i168 < 3; c1_i168++) {
    c1_a_o[c1_i168] = (c1_domega_d[c1_i168] + c1_y[c1_i168]) + c1_d_y[c1_i168];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 68);
  for (c1_i169 = 0; c1_i169 < 36; c1_i169++) {
    c1_b_Ta[c1_i169] = c1_Ta[c1_i169];
  }

  c1_b_inv(chartInstance, c1_b_Ta, c1_dv17);
  c1_i170 = 0;
  for (c1_i171 = 0; c1_i171 < 6; c1_i171++) {
    c1_i172 = 0;
    for (c1_i173 = 0; c1_i173 < 6; c1_i173++) {
      c1_d_a[c1_i173 + c1_i170] = c1_dv17[c1_i172 + c1_i171];
      c1_i172 += 6;
    }

    c1_i170 += 6;
  }

  for (c1_i174 = 0; c1_i174 < 6; c1_i174++) {
    c1_c_b[c1_i174] = c1_b_F[c1_i174];
  }

  c1_d_eml_scalar_eg(chartInstance);
  c1_d_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i175 = 0; c1_i175 < 6; c1_i175++) {
    c1_f_y[c1_i175] = 0.0;
    c1_i176 = 0;
    for (c1_i177 = 0; c1_i177 < 6; c1_i177++) {
      c1_f_y[c1_i175] += c1_d_a[c1_i176 + c1_i175] * c1_c_b[c1_i177];
      c1_i176 += 6;
    }
  }

  for (c1_i178 = 0; c1_i178 < 6; c1_i178++) {
    c1_b_Fa[c1_i178] = c1_f_y[c1_i178];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 69);
  for (c1_i179 = 0; c1_i179 < 3; c1_i179++) {
    c1_b_a[c1_i179] = c1_a_t[c1_i179];
  }

  for (c1_i180 = 0; c1_i180 < 3; c1_i180++) {
    c1_b_a[c1_i180 + 3] = c1_a_o[c1_i180];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -69);
  _SFD_SYMBOL_SCOPE_POP();
  for (c1_i181 = 0; c1_i181 < 6; c1_i181++) {
    (*chartInstance->c1_a)[c1_i181] = c1_b_a[c1_i181];
  }

  for (c1_i182 = 0; c1_i182 < 6; c1_i182++) {
    (*chartInstance->c1_Fa)[c1_i182] = c1_b_Fa[c1_i182];
  }

  for (c1_i183 = 0; c1_i183 < 3; c1_i183++) {
    (*chartInstance->c1_dp_d)[c1_i183] = c1_b_dp_d[c1_i183];
  }

  for (c1_i184 = 0; c1_i184 < 9; c1_i184++) {
    (*chartInstance->c1_T_e)[c1_i184] = c1_b_T_e[c1_i184];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
}

static void initSimStructsc1_IPILCO_relativeRPY_S
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_tr2rt(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                     real_T c1_T[16], real_T c1_R[9])
{
  uint32_T c1_debug_family_var_map[6];
  real_T c1_n;
  real_T c1_t[3];
  real_T c1_nargin = 1.0;
  real_T c1_nargout = 1.0;
  uint32_T c1_b_debug_family_var_map[3];
  real_T c1_c;
  real_T c1_b_nargin = 1.0;
  real_T c1_b_nargout = 1.0;
  real_T c1_r;
  real_T c1_c_nargin = 1.0;
  real_T c1_c_nargout = 1.0;
  real_T c1_b_c;
  real_T c1_d_nargin = 1.0;
  real_T c1_d_nargout = 1.0;
  int32_T c1_i185;
  int32_T c1_i186;
  int32_T c1_i187;
  int32_T c1_i188;
  int32_T c1_i189;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c1_d_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_n, 0U, c1_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_t, 1U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 2U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 3U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_T, 4U, c1_e_sf_marshallOut,
    c1_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_R, 5U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 41);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c1_b_debug_family_names,
    c1_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_c, 0U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargin, 1U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargout, 2U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 31);
  c1_c = 4.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c1_c_debug_family_names,
    c1_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_r, 0U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_c_nargin, 1U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_c_nargout, 2U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  CV_SCRIPT_FCN(2, 0);
  _SFD_SCRIPT_CALL(2U, chartInstance->c1_sfEvent, 33);
  c1_r = 4.0;
  _SFD_SCRIPT_CALL(2U, chartInstance->c1_sfEvent, -33);
  _SFD_SYMBOL_SCOPE_POP();
  CV_SCRIPT_IF(0, 0, CV_RELATIONAL_EVAL(14U, 0U, 0, 4.0, 4.0, -1, 1U, 0));
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 45);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c1_b_debug_family_names,
    c1_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_c, 0U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_d_nargin, 1U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_d_nargout, 2U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, 31);
  c1_b_c = 4.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c1_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
  c1_n = 4.0;
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 47);
  CV_SCRIPT_IF(0, 1, CV_RELATIONAL_EVAL(14U, 0U, 1, 1.0, 1.0, -1, 4U, 0));
  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 55);
  c1_i185 = 0;
  c1_i186 = 0;
  for (c1_i187 = 0; c1_i187 < 3; c1_i187++) {
    for (c1_i188 = 0; c1_i188 < 3; c1_i188++) {
      c1_R[c1_i188 + c1_i185] = c1_T[c1_i188 + c1_i186];
    }

    c1_i185 += 3;
    c1_i186 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 56);
  for (c1_i189 = 0; c1_i189 < 3; c1_i189++) {
    c1_t[c1_i189] = c1_T[c1_i189 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, -56);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c1_transl(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c1_x[16], real_T c1_t1[3])
{
  uint32_T c1_debug_family_var_map[4];
  real_T c1_nargin = 1.0;
  real_T c1_nargout = 1.0;
  boolean_T c1_h;
  real_T c1_d[2];
  real_T c1_b_nargin = 1.0;
  real_T c1_b_nargout = 1.0;
  int32_T c1_i190;
  int32_T c1_i191;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c1_f_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 0U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 1U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_x, 2U, c1_e_sf_marshallOut,
    c1_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_t1, 3U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  CV_SCRIPT_FCN(3, 0);
  _SFD_SCRIPT_CALL(3U, chartInstance->c1_sfEvent, 53);
  CV_SCRIPT_IF(3, 0, CV_RELATIONAL_EVAL(14U, 3U, 0, 1.0, 1.0, -1, 0U, 1));
  _SFD_SCRIPT_CALL(3U, chartInstance->c1_sfEvent, 54);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c1_e_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_h, 0U, c1_j_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_d, 1U, c1_i_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargin, 2U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargout, 3U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  CV_SCRIPT_FCN(4, 0);
  _SFD_SCRIPT_CALL(4U, chartInstance->c1_sfEvent, 37);
  for (c1_i190 = 0; c1_i190 < 2; c1_i190++) {
    c1_d[c1_i190] = 4.0;
  }

  _SFD_SCRIPT_CALL(4U, chartInstance->c1_sfEvent, 38);
  CV_SCRIPT_IF(4, 0, CV_RELATIONAL_EVAL(14U, 4U, 0, 2.0, 2.0, -1, 5U, 1));
  _SFD_SCRIPT_CALL(4U, chartInstance->c1_sfEvent, 39);
  c1_h = true;
  _SFD_SCRIPT_CALL(4U, chartInstance->c1_sfEvent, 41);
  CV_SCRIPT_COND(4, 0, c1_h);
  CV_SCRIPT_COND(4, 1, CV_RELATIONAL_EVAL(14U, 4U, 1, 1.0, 1.0, -1, 4U, 0));
  CV_SCRIPT_MCDC(4, 0, false);
  CV_SCRIPT_IF(4, 1, false);
  _SFD_SCRIPT_CALL(4U, chartInstance->c1_sfEvent, -46);
  _SFD_SYMBOL_SCOPE_POP();
  CV_SCRIPT_IF(3, 1, true);
  _SFD_SCRIPT_CALL(3U, chartInstance->c1_sfEvent, 55);
  CV_SCRIPT_IF(3, 2, CV_RELATIONAL_EVAL(14U, 3U, 1, 2.0, 3.0, -1, 0U, 0));
  _SFD_SCRIPT_CALL(3U, chartInstance->c1_sfEvent, 66);
  CV_SCRIPT_COND(3, 0, CV_RELATIONAL_EVAL(14U, 3U, 4, 1.0, 1.0, -1, 0U, 1));
  CV_SCRIPT_MCDC(3, 0, true);
  CV_SCRIPT_IF(3, 5, true);
  _SFD_SCRIPT_CALL(3U, chartInstance->c1_sfEvent, 67);
  for (c1_i191 = 0; c1_i191 < 3; c1_i191++) {
    c1_t1[c1_i191] = c1_x[c1_i191 + 12];
  }

  _SFD_SCRIPT_CALL(3U, chartInstance->c1_sfEvent, -89);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c1_my_tr2rpy(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_m[9], real_T c1_rpy[3])
{
  uint32_T c1_debug_family_var_map[11];
  c1_snY5FzvLGhd6UMWwSEwxMX c1_opt;
  real_T c1_s[2];
  real_T c1_sr;
  real_T c1_cr;
  real_T c1_sp;
  real_T c1_cp;
  char_T c1_type[3];
  real_T c1_nargin = 2.0;
  real_T c1_nargout = 1.0;
  int32_T c1_i192;
  static char_T c1_cv0[3] = { 'x', 'y', 'z' };

  int32_T c1_i193;
  int32_T c1_i194;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_d0;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_d1;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_d2;
  real_T c1_k_x;
  real_T c1_l_x;
  real_T c1_d3;
  real_T c1_m_x;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_p_x;
  int32_T c1_i195;
  real_T c1_b_a[3];
  int32_T c1_i196;
  int32_T c1_i197;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 11U, 11U, c1_g_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_opt, 0U, c1_m_sf_marshallOut,
    c1_i_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_s, 1U, c1_i_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_sr, 2U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_cr, 3U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_sp, 4U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_cp, 5U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_type, 6U, c1_l_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 7U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 8U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_m, 9U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_rpy, 10U, c1_k_sf_marshallOut,
    c1_h_sf_marshallIn);
  for (c1_i192 = 0; c1_i192 < 3; c1_i192++) {
    c1_type[c1_i192] = c1_cv0[c1_i192];
  }

  CV_SCRIPT_FCN(5, 0);
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 49);
  c1_opt.deg = false;
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 50);
  c1_opt.zyx = false;
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 51);
  CV_SCRIPT_IF(5, 0, false);
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 56);
  for (c1_i193 = 0; c1_i193 < 2; c1_i193++) {
    c1_s[c1_i193] = 3.0;
  }

  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 57);
  CV_SCRIPT_IF(5, 1, CV_RELATIONAL_EVAL(14U, 5U, 0, 2.0, 2.0, -1, 4U, 0));
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 64);
  for (c1_i194 = 0; c1_i194 < 3; c1_i194++) {
    c1_rpy[c1_i194] = 0.0;
  }

  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 66);
  if (CV_SCRIPT_IF(5, 2, CV_SCRIPT_MCDC(5, 0, !CV_SCRIPT_COND(5, 0, c1_opt.zyx))))
  {
    _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 68);
    c1_eml_eps(chartInstance);
    c1_x = c1_m[8];
    c1_b_x = c1_x;
    c1_d0 = muDoubleScalarAbs(c1_b_x);
    guard2 = false;
    if (CV_SCRIPT_COND(5, 1, CV_RELATIONAL_EVAL(14U, 5U, 1, c1_d0,
          2.2204460492503131E-16, -1, 2U, c1_d0 < 2.2204460492503131E-16))) {
      c1_eml_eps(chartInstance);
      c1_c_x = c1_m[7];
      c1_d_x = c1_c_x;
      c1_d1 = muDoubleScalarAbs(c1_d_x);
      if (CV_SCRIPT_COND(5, 2, CV_RELATIONAL_EVAL(14U, 5U, 2, c1_d1,
            2.2204460492503131E-16, -1, 2U, c1_d1 < 2.2204460492503131E-16))) {
        CV_SCRIPT_MCDC(5, 1, true);
        CV_SCRIPT_IF(5, 3, true);
        _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 70);
        c1_rpy[0] = 0.0;
        _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 71);
        c1_rpy[1] = c1_atan2(chartInstance, c1_m[6], c1_m[8]);
        _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 72);
        c1_rpy[2] = c1_atan2(chartInstance, c1_m[1], c1_m[4]);
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2 == true) {
      CV_SCRIPT_MCDC(5, 1, false);
      CV_SCRIPT_IF(5, 3, false);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 74);
      c1_rpy[0] = c1_atan2(chartInstance, -c1_m[7], c1_m[8]);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 76);
      c1_e_x = c1_rpy[0];
      c1_sr = c1_e_x;
      c1_f_x = c1_sr;
      c1_sr = c1_f_x;
      c1_sr = muDoubleScalarSin(c1_sr);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 77);
      c1_g_x = c1_rpy[0];
      c1_cr = c1_g_x;
      c1_h_x = c1_cr;
      c1_cr = c1_h_x;
      c1_cr = muDoubleScalarCos(c1_cr);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 78);
      c1_rpy[1] = c1_atan2(chartInstance, c1_m[6], c1_cr * c1_m[8] - c1_sr *
                           c1_m[7]);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 79);
      c1_rpy[2] = c1_atan2(chartInstance, -c1_m[3], c1_m[0]);
    }
  } else {
    _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 83);
    c1_eml_eps(chartInstance);
    c1_i_x = c1_m[0];
    c1_j_x = c1_i_x;
    c1_d2 = muDoubleScalarAbs(c1_j_x);
    guard1 = false;
    if (CV_SCRIPT_COND(5, 3, CV_RELATIONAL_EVAL(14U, 5U, 3, c1_d2,
          2.2204460492503131E-16, -1, 2U, c1_d2 < 2.2204460492503131E-16))) {
      c1_eml_eps(chartInstance);
      c1_k_x = c1_m[1];
      c1_l_x = c1_k_x;
      c1_d3 = muDoubleScalarAbs(c1_l_x);
      if (CV_SCRIPT_COND(5, 4, CV_RELATIONAL_EVAL(14U, 5U, 4, c1_d3,
            2.2204460492503131E-16, -1, 2U, c1_d3 < 2.2204460492503131E-16))) {
        CV_SCRIPT_MCDC(5, 2, true);
        CV_SCRIPT_IF(5, 4, true);
        _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 85);
        c1_rpy[0] = 0.0;
        _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 86);
        c1_rpy[1] = c1_atan2(chartInstance, -c1_m[2], c1_m[0]);
        _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 87);
        c1_rpy[2] = c1_atan2(chartInstance, -c1_m[7], c1_m[4]);
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 == true) {
      CV_SCRIPT_MCDC(5, 2, false);
      CV_SCRIPT_IF(5, 4, false);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 89);
      c1_rpy[0] = c1_atan2(chartInstance, c1_m[1], c1_m[0]);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 90);
      c1_m_x = c1_rpy[0];
      c1_sp = c1_m_x;
      c1_n_x = c1_sp;
      c1_sp = c1_n_x;
      c1_sp = muDoubleScalarSin(c1_sp);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 91);
      c1_o_x = c1_rpy[0];
      c1_cp = c1_o_x;
      c1_p_x = c1_cp;
      c1_cp = c1_p_x;
      c1_cp = muDoubleScalarCos(c1_cp);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 92);
      c1_rpy[1] = c1_atan2(chartInstance, -c1_m[2], c1_cp * c1_m[0] + c1_sp *
                           c1_m[1]);
      _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 93);
      c1_rpy[2] = c1_atan2(chartInstance, c1_sp * c1_m[6] - c1_cp * c1_m[7],
                           c1_cp * c1_m[4] - c1_sp * c1_m[3]);
    }
  }

  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 96);
  if (CV_SCRIPT_IF(5, 5, c1_opt.deg)) {
    _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 97);
    for (c1_i195 = 0; c1_i195 < 3; c1_i195++) {
      c1_b_a[c1_i195] = c1_rpy[c1_i195];
    }

    for (c1_i196 = 0; c1_i196 < 3; c1_i196++) {
      c1_b_a[c1_i196] *= 180.0;
    }

    for (c1_i197 = 0; c1_i197 < 3; c1_i197++) {
      c1_rpy[c1_i197] = c1_b_a[c1_i197] / 3.1415926535897931;
    }
  }

  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, -97);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c1_cardatol(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_q[3], real_T c1_i, real_T c1_j, real_T c1_k, real_T c1_R[9])
{
  uint32_T c1_debug_family_var_map[21];
  real_T c1_X;
  real_T c1_Y;
  real_T c1_Z;
  real_T c1_alfa;
  real_T c1_beta;
  real_T c1_gamma;
  real_T c1_sa;
  real_T c1_sb;
  real_T c1_sc;
  real_T c1_ca;
  real_T c1_cb;
  real_T c1_cc;
  real_T c1_sig;
  real_T c1_l;
  real_T c1_nargin = 4.0;
  real_T c1_nargout = 1.0;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_l_x;
  int32_T c1_i198;
  static char_T c1_varargin_1[49] = { '	', '*', '*', '*', '*', ' ', '	', 'S',
    'p', 'a', 'c', 'e', 'l', 'i', 'b', ' ', 'E', 'r', 'r', 'o', 'r', ':', ' ',
    'I', 'l', 'l', 'e', 'g', 'a', 'l', ' ', 'r', 'o', 't', 'a', 't', 'i', 'o',
    'n', ' ', 'a', 'x', 'i', 's', '	', '*', '*', '*', '*' };

  char_T c1_u[49];
  const mxArray *c1_y = NULL;
  real_T c1_m_x;
  real_T c1_xk;
  real_T c1_d4;
  int32_T c1_i199;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T guard5 = false;
  boolean_T guard6 = false;
  boolean_T guard7 = false;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 21U, 21U, c1_h_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_X, 0U, c1_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_Y, 1U, c1_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_Z, 2U, c1_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_alfa, 3U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_beta, 4U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_gamma, 5U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_sa, 6U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_sb, 7U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_sc, 8U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_ca, 9U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_cb, 10U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_cc, 11U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_sig, 12U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_l, 13U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 14U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 15U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_q, 16U, c1_h_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_i, 17U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_j, 18U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_k, 19U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_R, 20U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_SCRIPT_FCN(6, 0);
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 23);
  c1_X = 1.0;
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 23);
  c1_Y = 2.0;
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 23);
  c1_Z = 3.0;
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 25);
  c1_alfa = c1_q[0];
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 26);
  c1_beta = c1_q[1];
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 27);
  c1_gamma = c1_q[2];
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 29);
  c1_x = c1_alfa;
  c1_sa = c1_x;
  c1_b_x = c1_sa;
  c1_sa = c1_b_x;
  c1_sa = muDoubleScalarSin(c1_sa);
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 29);
  c1_c_x = c1_beta;
  c1_sb = c1_c_x;
  c1_d_x = c1_sb;
  c1_sb = c1_d_x;
  c1_sb = muDoubleScalarSin(c1_sb);
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 29);
  c1_e_x = c1_gamma;
  c1_sc = c1_e_x;
  c1_f_x = c1_sc;
  c1_sc = c1_f_x;
  c1_sc = muDoubleScalarSin(c1_sc);
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 30);
  c1_g_x = c1_alfa;
  c1_ca = c1_g_x;
  c1_h_x = c1_ca;
  c1_ca = c1_h_x;
  c1_ca = muDoubleScalarCos(c1_ca);
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 30);
  c1_i_x = c1_beta;
  c1_cb = c1_i_x;
  c1_j_x = c1_cb;
  c1_cb = c1_j_x;
  c1_cb = muDoubleScalarCos(c1_cb);
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 30);
  c1_k_x = c1_gamma;
  c1_cc = c1_k_x;
  c1_l_x = c1_cc;
  c1_cc = c1_l_x;
  c1_cc = muDoubleScalarCos(c1_cc);
  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 32);
  guard1 = false;
  guard2 = false;
  guard3 = false;
  guard4 = false;
  guard5 = false;
  guard6 = false;
  guard7 = false;
  if (CV_SCRIPT_COND(6, 0, CV_RELATIONAL_EVAL(14U, 6U, 0, c1_i, c1_X, -1, 2U,
        c1_i < c1_X))) {
    guard7 = true;
  } else if (CV_SCRIPT_COND(6, 1, CV_RELATIONAL_EVAL(14U, 6U, 1, c1_i, c1_Z, -1,
    4U, c1_i > c1_Z))) {
    guard7 = true;
  } else if (CV_SCRIPT_COND(6, 2, CV_RELATIONAL_EVAL(14U, 6U, 2, c1_j, c1_X, -1,
    2U, c1_j < c1_X))) {
    guard6 = true;
  } else if (CV_SCRIPT_COND(6, 3, CV_RELATIONAL_EVAL(14U, 6U, 3, c1_j, c1_Z, -1,
    4U, c1_j > c1_Z))) {
    guard5 = true;
  } else if (CV_SCRIPT_COND(6, 4, CV_RELATIONAL_EVAL(14U, 6U, 4, c1_k, c1_X, -1,
    2U, c1_k < c1_X))) {
    guard4 = true;
  } else if (CV_SCRIPT_COND(6, 5, CV_RELATIONAL_EVAL(14U, 6U, 5, c1_k, c1_Z, -1,
    4U, c1_k > c1_Z))) {
    guard3 = true;
  } else if (CV_SCRIPT_COND(6, 6, CV_RELATIONAL_EVAL(14U, 6U, 6, c1_i, c1_j, -1,
    0U, c1_i == c1_j))) {
    guard2 = true;
  } else if (CV_SCRIPT_COND(6, 7, CV_RELATIONAL_EVAL(14U, 6U, 7, c1_j, c1_k, -1,
    0U, c1_j == c1_k))) {
    guard1 = true;
  } else {
    CV_SCRIPT_MCDC(6, 0, false);
    CV_SCRIPT_IF(6, 0, false);
  }

  if (guard7 == true) {
    guard6 = true;
  }

  if (guard6 == true) {
    guard5 = true;
  }

  if (guard5 == true) {
    guard4 = true;
  }

  if (guard4 == true) {
    guard3 = true;
  }

  if (guard3 == true) {
    guard2 = true;
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    CV_SCRIPT_MCDC(6, 0, true);
    CV_SCRIPT_IF(6, 0, true);
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 33);
    for (c1_i198 = 0; c1_i198 < 49; c1_i198++) {
      c1_u[c1_i198] = c1_varargin_1[c1_i198];
    }

    c1_y = NULL;
    sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 49),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14, c1_y);
  }

  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 36);
  c1_m_x = (c1_j - c1_i) + 3.0;
  c1_c_eml_scalar_eg(chartInstance);
  c1_xk = c1_m_x;
  c1_d4 = muDoubleScalarRem(c1_xk, 3.0);
  if (CV_SCRIPT_IF(6, 1, CV_RELATIONAL_EVAL(14U, 6U, 8, c1_d4, 1.0, -1, 0U,
        c1_d4 == 1.0))) {
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 36);
    c1_sig = 1.0;
  } else {
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 37);
    c1_sig = -1.0;
  }

  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 40);
  for (c1_i199 = 0; c1_i199 < 9; c1_i199++) {
    c1_R[c1_i199] = 0.0;
  }

  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 42);
  if (CV_SCRIPT_IF(6, 2, CV_RELATIONAL_EVAL(14U, 6U, 9, c1_i, c1_k, -1, 1U, c1_i
        != c1_k))) {
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 43);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("i", c1_i),
      1, 3, 1, 0) - 1] = 1.0;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 44);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("i", c1_i),
      1, 3, 1, 0) + 2] = 0.0;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 45);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("i", c1_i),
      1, 3, 1, 0) + 5] = c1_sig * c1_sb;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 47);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("j", c1_j),
      1, 3, 1, 0) - 1] = 0.0;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 48);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("j", c1_j),
      1, 3, 1, 0) + 2] = c1_ca;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 49);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("j", c1_j),
      1, 3, 1, 0) + 5] = -c1_sig * c1_sa * c1_cb;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 51);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("k", c1_k),
      1, 3, 1, 0) - 1] = 0.0;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 52);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("k", c1_k),
      1, 3, 1, 0) + 2] = c1_sig * c1_sa;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 53);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("k", c1_k),
      1, 3, 1, 0) + 5] = c1_ca * c1_cb;
  } else {
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 57);
    c1_l = (6.0 - c1_i) - c1_j;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 59);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("i", c1_i),
      1, 3, 1, 0) - 1] = 1.0;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 60);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("i", c1_i),
      1, 3, 1, 0) + 2] = 0.0;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 61);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("i", c1_i),
      1, 3, 1, 0) + 5] = c1_cb;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 63);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("j", c1_j),
      1, 3, 1, 0) - 1] = 0.0;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 64);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("j", c1_j),
      1, 3, 1, 0) + 2] = c1_ca;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 65);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("j", c1_j),
      1, 3, 1, 0) + 5] = c1_sa * c1_sb;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 67);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("l", c1_l),
      1, 3, 1, 0) - 1] = 0.0;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 68);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("l", c1_l),
      1, 3, 1, 0) + 2] = c1_sig * c1_sa;
    _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, 69);
    c1_R[_SFD_EML_ARRAY_BOUNDS_CHECK("R", (int32_T)_SFD_INTEGER_CHECK("l", c1_l),
      1, 3, 1, 0) + 5] = -c1_sig * c1_ca * c1_sb;
  }

  _SFD_SCRIPT_CALL(6U, chartInstance->c1_sfEvent, -69);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber)
{
  (void)c1_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 0U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\robot\\tr2rt.m"));
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 1U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\common\\numcols.m"));
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 2U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\common\\numrows.m"));
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 3U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\robot\\transl.m"));
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 4U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\common\\ishomog.m"));
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 5U,
    sf_debug_get_script_id("D:\\victor\\MATLAB\\I-PILCO\\util\\my_tr2rpy.m"));
  _SFD_SCRIPT_TRANSLATION(c1_chartNumber, c1_instanceNumber, 6U,
    sf_debug_get_script_id("D:\\victor\\MATLAB\\I-PILCO\\util\\cardatol.m"));
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i200;
  int32_T c1_i201;
  int32_T c1_i202;
  real_T c1_b_inData[9];
  int32_T c1_i203;
  int32_T c1_i204;
  int32_T c1_i205;
  real_T c1_u[9];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i200 = 0;
  for (c1_i201 = 0; c1_i201 < 3; c1_i201++) {
    for (c1_i202 = 0; c1_i202 < 3; c1_i202++) {
      c1_b_inData[c1_i202 + c1_i200] = (*(real_T (*)[9])c1_inData)[c1_i202 +
        c1_i200];
    }

    c1_i200 += 3;
  }

  c1_i203 = 0;
  for (c1_i204 = 0; c1_i204 < 3; c1_i204++) {
    for (c1_i205 = 0; c1_i205 < 3; c1_i205++) {
      c1_u[c1_i205 + c1_i203] = c1_b_inData[c1_i205 + c1_i203];
    }

    c1_i203 += 3;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_b_T_e, const char_T *c1_identifier, real_T
  c1_y[9])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_T_e), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_T_e);
}

static void c1_b_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[9])
{
  real_T c1_dv18[9];
  int32_T c1_i206;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv18, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c1_i206 = 0; c1_i206 < 9; c1_i206++) {
    c1_y[c1_i206] = c1_dv18[c1_i206];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_T_e;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[9];
  int32_T c1_i207;
  int32_T c1_i208;
  int32_T c1_i209;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_b_T_e = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_T_e), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_T_e);
  c1_i207 = 0;
  for (c1_i208 = 0; c1_i208 < 3; c1_i208++) {
    for (c1_i209 = 0; c1_i209 < 3; c1_i209++) {
      (*(real_T (*)[9])c1_outData)[c1_i209 + c1_i207] = c1_y[c1_i209 + c1_i207];
    }

    c1_i207 += 3;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i210;
  real_T c1_b_inData[3];
  int32_T c1_i211;
  real_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i210 = 0; c1_i210 < 3; c1_i210++) {
    c1_b_inData[c1_i210] = (*(real_T (*)[3])c1_inData)[c1_i210];
  }

  for (c1_i211 = 0; c1_i211 < 3; c1_i211++) {
    c1_u[c1_i211] = c1_b_inData[c1_i211];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 3, 1), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_c_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_b_dp_d, const char_T *c1_identifier, real_T
  c1_y[3])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_dp_d), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_dp_d);
}

static void c1_d_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[3])
{
  real_T c1_dv19[3];
  int32_T c1_i212;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv19, 1, 0, 0U, 1, 0U, 2, 3, 1);
  for (c1_i212 = 0; c1_i212 < 3; c1_i212++) {
    c1_y[c1_i212] = c1_dv19[c1_i212];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_dp_d;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[3];
  int32_T c1_i213;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_b_dp_d = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_dp_d), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_dp_d);
  for (c1_i213 = 0; c1_i213 < 3; c1_i213++) {
    (*(real_T (*)[3])c1_outData)[c1_i213] = c1_y[c1_i213];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i214;
  real_T c1_b_inData[6];
  int32_T c1_i215;
  real_T c1_u[6];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i214 = 0; c1_i214 < 6; c1_i214++) {
    c1_b_inData[c1_i214] = (*(real_T (*)[6])c1_inData)[c1_i214];
  }

  for (c1_i215 = 0; c1_i215 < 6; c1_i215++) {
    c1_u[c1_i215] = c1_b_inData[c1_i215];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 6, 1), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_e_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_b_Fa, const char_T *c1_identifier, real_T
  c1_y[6])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Fa), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Fa);
}

static void c1_f_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[6])
{
  real_T c1_dv20[6];
  int32_T c1_i216;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv20, 1, 0, 0U, 1, 0U, 2, 6, 1);
  for (c1_i216 = 0; c1_i216 < 6; c1_i216++) {
    c1_y[c1_i216] = c1_dv20[c1_i216];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_Fa;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[6];
  int32_T c1_i217;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_b_Fa = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_Fa), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_Fa);
  for (c1_i217 = 0; c1_i217 < 6; c1_i217++) {
    (*(real_T (*)[6])c1_outData)[c1_i217] = c1_y[c1_i217];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  real_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(real_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i218;
  int32_T c1_i219;
  int32_T c1_i220;
  real_T c1_b_inData[16];
  int32_T c1_i221;
  int32_T c1_i222;
  int32_T c1_i223;
  real_T c1_u[16];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i218 = 0;
  for (c1_i219 = 0; c1_i219 < 4; c1_i219++) {
    for (c1_i220 = 0; c1_i220 < 4; c1_i220++) {
      c1_b_inData[c1_i220 + c1_i218] = (*(real_T (*)[16])c1_inData)[c1_i220 +
        c1_i218];
    }

    c1_i218 += 4;
  }

  c1_i221 = 0;
  for (c1_i222 = 0; c1_i222 < 4; c1_i222++) {
    for (c1_i223 = 0; c1_i223 < 4; c1_i223++) {
      c1_u[c1_i223 + c1_i221] = c1_b_inData[c1_i223 + c1_i221];
    }

    c1_i221 += 4;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i224;
  real_T c1_b_inData[6];
  int32_T c1_i225;
  real_T c1_u[6];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i224 = 0; c1_i224 < 6; c1_i224++) {
    c1_b_inData[c1_i224] = (*(real_T (*)[6])c1_inData)[c1_i224];
  }

  for (c1_i225 = 0; c1_i225 < 6; c1_i225++) {
    c1_u[c1_i225] = c1_b_inData[c1_i225];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_g_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i226;
  int32_T c1_i227;
  int32_T c1_i228;
  real_T c1_b_inData[36];
  int32_T c1_i229;
  int32_T c1_i230;
  int32_T c1_i231;
  real_T c1_u[36];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i226 = 0;
  for (c1_i227 = 0; c1_i227 < 6; c1_i227++) {
    for (c1_i228 = 0; c1_i228 < 6; c1_i228++) {
      c1_b_inData[c1_i228 + c1_i226] = (*(real_T (*)[36])c1_inData)[c1_i228 +
        c1_i226];
    }

    c1_i226 += 6;
  }

  c1_i229 = 0;
  for (c1_i230 = 0; c1_i230 < 6; c1_i230++) {
    for (c1_i231 = 0; c1_i231 < 6; c1_i231++) {
      c1_u[c1_i231 + c1_i229] = c1_b_inData[c1_i231 + c1_i229];
    }

    c1_i229 += 6;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 6, 6), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static real_T c1_g_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d5;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d5, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d5;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_nargout;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_nargout = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_nargout), &c1_thisId);
  sf_mex_destroy(&c1_nargout);
  *(real_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_h_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i232;
  real_T c1_b_inData[3];
  int32_T c1_i233;
  real_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i232 = 0; c1_i232 < 3; c1_i232++) {
    c1_b_inData[c1_i232] = (*(real_T (*)[3])c1_inData)[c1_i232];
  }

  for (c1_i233 = 0; c1_i233 < 3; c1_i233++) {
    c1_u[c1_i233] = c1_b_inData[c1_i233];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_h_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[3])
{
  real_T c1_dv21[3];
  int32_T c1_i234;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv21, 1, 0, 0U, 1, 0U, 1, 3);
  for (c1_i234 = 0; c1_i234 < 3; c1_i234++) {
    c1_y[c1_i234] = c1_dv21[c1_i234];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_a_o;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[3];
  int32_T c1_i235;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_a_o = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_a_o), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_a_o);
  for (c1_i235 = 0; c1_i235 < 3; c1_i235++) {
    (*(real_T (*)[3])c1_outData)[c1_i235] = c1_y[c1_i235];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static void c1_i_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[36])
{
  real_T c1_dv22[36];
  int32_T c1_i236;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv22, 1, 0, 0U, 1, 0U, 2, 6, 6);
  for (c1_i236 = 0; c1_i236 < 36; c1_i236++) {
    c1_y[c1_i236] = c1_dv22[c1_i236];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Ta;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[36];
  int32_T c1_i237;
  int32_T c1_i238;
  int32_T c1_i239;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_Ta = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Ta), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_Ta);
  c1_i237 = 0;
  for (c1_i238 = 0; c1_i238 < 6; c1_i238++) {
    for (c1_i239 = 0; c1_i239 < 6; c1_i239++) {
      (*(real_T (*)[36])c1_outData)[c1_i239 + c1_i237] = c1_y[c1_i239 + c1_i237];
    }

    c1_i237 += 6;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static void c1_j_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[16])
{
  real_T c1_dv23[16];
  int32_T c1_i240;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv23, 1, 0, 0U, 1, 0U, 2, 4, 4);
  for (c1_i240 = 0; c1_i240 < 16; c1_i240++) {
    c1_y[c1_i240] = c1_dv23[c1_i240];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_T;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[16];
  int32_T c1_i241;
  int32_T c1_i242;
  int32_T c1_i243;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_T = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_T), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_T);
  c1_i241 = 0;
  for (c1_i242 = 0; c1_i242 < 4; c1_i242++) {
    for (c1_i243 = 0; c1_i243 < 4; c1_i243++) {
      (*(real_T (*)[16])c1_outData)[c1_i243 + c1_i241] = c1_y[c1_i243 + c1_i241];
    }

    c1_i241 += 4;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_i_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i244;
  real_T c1_b_inData[2];
  int32_T c1_i245;
  real_T c1_u[2];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i244 = 0; c1_i244 < 2; c1_i244++) {
    c1_b_inData[c1_i244] = (*(real_T (*)[2])c1_inData)[c1_i244];
  }

  for (c1_i245 = 0; c1_i245 < 2; c1_i245++) {
    c1_u[c1_i245] = c1_b_inData[c1_i245];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 1, 2), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_j_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  boolean_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(boolean_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 11, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_k_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i246;
  real_T c1_b_inData[3];
  int32_T c1_i247;
  real_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i246 = 0; c1_i246 < 3; c1_i246++) {
    c1_b_inData[c1_i246] = (*(real_T (*)[3])c1_inData)[c1_i246];
  }

  for (c1_i247 = 0; c1_i247 < 3; c1_i247++) {
    c1_u[c1_i247] = c1_b_inData[c1_i247];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 1, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_k_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  real_T c1_y[3])
{
  real_T c1_dv24[3];
  int32_T c1_i248;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv24, 1, 0, 0U, 1, 0U, 2, 1, 3);
  for (c1_i248 = 0; c1_i248 < 3; c1_i248++) {
    c1_y[c1_i248] = c1_dv24[c1_i248];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_rpy;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[3];
  int32_T c1_i249;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_rpy = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_rpy), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_rpy);
  for (c1_i249 = 0; c1_i249 < 3; c1_i249++) {
    (*(real_T (*)[3])c1_outData)[c1_i249] = c1_y[c1_i249];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_l_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i250;
  char_T c1_b_inData[3];
  int32_T c1_i251;
  char_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i250 = 0; c1_i250 < 3; c1_i250++) {
    c1_b_inData[c1_i250] = (*(char_T (*)[3])c1_inData)[c1_i250];
  }

  for (c1_i251 = 0; c1_i251 < 3; c1_i251++) {
    c1_u[c1_i251] = c1_b_inData[c1_i251];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_m_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  c1_snY5FzvLGhd6UMWwSEwxMX c1_u;
  const mxArray *c1_y = NULL;
  boolean_T c1_b_u;
  const mxArray *c1_b_y = NULL;
  boolean_T c1_c_u;
  const mxArray *c1_c_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(c1_snY5FzvLGhd6UMWwSEwxMX *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  c1_b_u = c1_u.deg;
  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", &c1_b_u, 11, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_b_y, "deg", "deg", 0);
  c1_c_u = c1_u.zyx;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_c_u, 11, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c1_y, c1_c_y, "zyx", "zyx", 0);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static c1_snY5FzvLGhd6UMWwSEwxMX c1_l_emlrt_marshallIn
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance, const mxArray *c1_u,
   const emlrtMsgIdentifier *c1_parentId)
{
  c1_snY5FzvLGhd6UMWwSEwxMX c1_y;
  emlrtMsgIdentifier c1_thisId;
  static const char * c1_fieldNames[2] = { "deg", "zyx" };

  c1_thisId.fParent = c1_parentId;
  sf_mex_check_struct(c1_parentId, c1_u, 2, c1_fieldNames, 0U, NULL);
  c1_thisId.fIdentifier = "deg";
  c1_y.deg = c1_m_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "deg", "deg", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "zyx";
  c1_y.zyx = c1_m_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "zyx", "zyx", 0)), &c1_thisId);
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static boolean_T c1_m_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  boolean_T c1_y;
  boolean_T c1_b0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_b0, 1, 11, 0U, 0, 0U, 0);
  c1_y = c1_b0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_opt;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  c1_snY5FzvLGhd6UMWwSEwxMX c1_y;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_opt = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_opt), &c1_thisId);
  sf_mex_destroy(&c1_opt);
  *(c1_snY5FzvLGhd6UMWwSEwxMX *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_IPILCO_relativeRPY_S_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  sf_mex_assign(&c1_nameCaptureInfo, sf_mex_createstruct("structure", 2, 224, 1),
                false);
  c1_info_helper(&c1_nameCaptureInfo);
  c1_b_info_helper(&c1_nameCaptureInfo);
  c1_c_info_helper(&c1_nameCaptureInfo);
  c1_d_info_helper(&c1_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c1_nameCaptureInfo);
  return c1_nameCaptureInfo;
}

static void c1_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs0 = NULL;
  const mxArray *c1_lhs0 = NULL;
  const mxArray *c1_rhs1 = NULL;
  const mxArray *c1_lhs1 = NULL;
  const mxArray *c1_rhs2 = NULL;
  const mxArray *c1_lhs2 = NULL;
  const mxArray *c1_rhs3 = NULL;
  const mxArray *c1_lhs3 = NULL;
  const mxArray *c1_rhs4 = NULL;
  const mxArray *c1_lhs4 = NULL;
  const mxArray *c1_rhs5 = NULL;
  const mxArray *c1_lhs5 = NULL;
  const mxArray *c1_rhs6 = NULL;
  const mxArray *c1_lhs6 = NULL;
  const mxArray *c1_rhs7 = NULL;
  const mxArray *c1_lhs7 = NULL;
  const mxArray *c1_rhs8 = NULL;
  const mxArray *c1_lhs8 = NULL;
  const mxArray *c1_rhs9 = NULL;
  const mxArray *c1_lhs9 = NULL;
  const mxArray *c1_rhs10 = NULL;
  const mxArray *c1_lhs10 = NULL;
  const mxArray *c1_rhs11 = NULL;
  const mxArray *c1_lhs11 = NULL;
  const mxArray *c1_rhs12 = NULL;
  const mxArray *c1_lhs12 = NULL;
  const mxArray *c1_rhs13 = NULL;
  const mxArray *c1_lhs13 = NULL;
  const mxArray *c1_rhs14 = NULL;
  const mxArray *c1_lhs14 = NULL;
  const mxArray *c1_rhs15 = NULL;
  const mxArray *c1_lhs15 = NULL;
  const mxArray *c1_rhs16 = NULL;
  const mxArray *c1_lhs16 = NULL;
  const mxArray *c1_rhs17 = NULL;
  const mxArray *c1_lhs17 = NULL;
  const mxArray *c1_rhs18 = NULL;
  const mxArray *c1_lhs18 = NULL;
  const mxArray *c1_rhs19 = NULL;
  const mxArray *c1_lhs19 = NULL;
  const mxArray *c1_rhs20 = NULL;
  const mxArray *c1_lhs20 = NULL;
  const mxArray *c1_rhs21 = NULL;
  const mxArray *c1_lhs21 = NULL;
  const mxArray *c1_rhs22 = NULL;
  const mxArray *c1_lhs22 = NULL;
  const mxArray *c1_rhs23 = NULL;
  const mxArray *c1_lhs23 = NULL;
  const mxArray *c1_rhs24 = NULL;
  const mxArray *c1_lhs24 = NULL;
  const mxArray *c1_rhs25 = NULL;
  const mxArray *c1_lhs25 = NULL;
  const mxArray *c1_rhs26 = NULL;
  const mxArray *c1_lhs26 = NULL;
  const mxArray *c1_rhs27 = NULL;
  const mxArray *c1_lhs27 = NULL;
  const mxArray *c1_rhs28 = NULL;
  const mxArray *c1_lhs28 = NULL;
  const mxArray *c1_rhs29 = NULL;
  const mxArray *c1_lhs29 = NULL;
  const mxArray *c1_rhs30 = NULL;
  const mxArray *c1_lhs30 = NULL;
  const mxArray *c1_rhs31 = NULL;
  const mxArray *c1_lhs31 = NULL;
  const mxArray *c1_rhs32 = NULL;
  const mxArray *c1_lhs32 = NULL;
  const mxArray *c1_rhs33 = NULL;
  const mxArray *c1_lhs33 = NULL;
  const mxArray *c1_rhs34 = NULL;
  const mxArray *c1_lhs34 = NULL;
  const mxArray *c1_rhs35 = NULL;
  const mxArray *c1_lhs35 = NULL;
  const mxArray *c1_rhs36 = NULL;
  const mxArray *c1_lhs36 = NULL;
  const mxArray *c1_rhs37 = NULL;
  const mxArray *c1_lhs37 = NULL;
  const mxArray *c1_rhs38 = NULL;
  const mxArray *c1_lhs38 = NULL;
  const mxArray *c1_rhs39 = NULL;
  const mxArray *c1_lhs39 = NULL;
  const mxArray *c1_rhs40 = NULL;
  const mxArray *c1_lhs40 = NULL;
  const mxArray *c1_rhs41 = NULL;
  const mxArray *c1_lhs41 = NULL;
  const mxArray *c1_rhs42 = NULL;
  const mxArray *c1_lhs42 = NULL;
  const mxArray *c1_rhs43 = NULL;
  const mxArray *c1_lhs43 = NULL;
  const mxArray *c1_rhs44 = NULL;
  const mxArray *c1_lhs44 = NULL;
  const mxArray *c1_rhs45 = NULL;
  const mxArray *c1_lhs45 = NULL;
  const mxArray *c1_rhs46 = NULL;
  const mxArray *c1_lhs46 = NULL;
  const mxArray *c1_rhs47 = NULL;
  const mxArray *c1_lhs47 = NULL;
  const mxArray *c1_rhs48 = NULL;
  const mxArray *c1_lhs48 = NULL;
  const mxArray *c1_rhs49 = NULL;
  const mxArray *c1_lhs49 = NULL;
  const mxArray *c1_rhs50 = NULL;
  const mxArray *c1_lhs50 = NULL;
  const mxArray *c1_rhs51 = NULL;
  const mxArray *c1_lhs51 = NULL;
  const mxArray *c1_rhs52 = NULL;
  const mxArray *c1_lhs52 = NULL;
  const mxArray *c1_rhs53 = NULL;
  const mxArray *c1_lhs53 = NULL;
  const mxArray *c1_rhs54 = NULL;
  const mxArray *c1_lhs54 = NULL;
  const mxArray *c1_rhs55 = NULL;
  const mxArray *c1_lhs55 = NULL;
  const mxArray *c1_rhs56 = NULL;
  const mxArray *c1_lhs56 = NULL;
  const mxArray *c1_rhs57 = NULL;
  const mxArray *c1_lhs57 = NULL;
  const mxArray *c1_rhs58 = NULL;
  const mxArray *c1_lhs58 = NULL;
  const mxArray *c1_rhs59 = NULL;
  const mxArray *c1_lhs59 = NULL;
  const mxArray *c1_rhs60 = NULL;
  const mxArray *c1_lhs60 = NULL;
  const mxArray *c1_rhs61 = NULL;
  const mxArray *c1_lhs61 = NULL;
  const mxArray *c1_rhs62 = NULL;
  const mxArray *c1_lhs62 = NULL;
  const mxArray *c1_rhs63 = NULL;
  const mxArray *c1_lhs63 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("tr2rt"), "name", "name", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/tr2rt.m"),
                  "resolved", "resolved", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464793910U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c1_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/tr2rt.m"),
                  "context", "context", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("numcols"), "name", "name", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/numcols.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464793908U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c1_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/tr2rt.m"),
                  "context", "context", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("numrows"), "name", "name", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/numrows.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464793908U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c1_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("diag"), "name", "name", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713854U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c1_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("ismatrix"), "name", "name", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1331304858U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c1_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c1_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c1_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c1_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/diag.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c1_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isfi"), "name", "name", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/fixedpoint/isfi.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1346510358U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c1_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/fixedpoint/isfi.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnumerictype"), "name",
                  "name", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/fixedpoint/isnumerictype.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1398875598U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c1_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c1_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c1_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c1_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c1_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("transl"), "name", "name", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/transl.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464793910U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c1_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/transl.m"),
                  "context", "context", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("ishomog"), "name", "name", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/ishomog.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464793908U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c1_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/ishomog.m"),
                  "context", "context", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("all"), "name", "name", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582414U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c1_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "context", "context",
                  18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c1_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "context", "context",
                  19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c1_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "context", "context",
                  20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.allOrAny"),
                  "name", "name", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583158U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c1_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c1_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isequal"), "name", "name", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818758U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c1_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_isequal_core"), "name",
                  "name", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818786U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c1_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.constNonSingletonDim"), "name", "name", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/constNonSingletonDim.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c1_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713880U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c1_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c1_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818796U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c1_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1386423952U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c1_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c1_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("inv"), "name", "name", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m"), "resolved",
                  "resolved", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1305318000U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c1_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c1_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713852U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c1_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c1_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818712U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c1_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1386423952U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c1_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3"), "context",
                  "context", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c1_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c1_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("norm"), "name", "name", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713868U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c1_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "context",
                  "context", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c1_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713852U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c1_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713858U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c1_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c1_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm"),
                  "context", "context", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_guarded_nan"), "name",
                  "name", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818776U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c1_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m"),
                  "context", "context", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818782U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c1_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_warning"), "name", "name",
                  45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818802U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c1_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isnan"), "name", "name", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713858U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c1_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c1_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_float_class"), "name",
                  "name", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818782U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c1_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_eps"), "name", "name", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c1_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c1_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond"),
                  "context", "context", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_flt2str"), "name", "name",
                  51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "resolved",
                  "resolved", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1360282350U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c1_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m"), "context",
                  "context", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "name", "name", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m"), "resolved",
                  "resolved", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1319729968U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c1_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c1_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c1_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c1_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c1_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs56), "lhs", "lhs",
                  56);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980690U), "fileTimeLo",
                  "fileTimeLo", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 57);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 57);
  sf_mex_assign(&c1_rhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs57, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs57), "rhs", "rhs",
                  57);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs57), "lhs", "lhs",
                  57);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 58);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 58);
  sf_mex_assign(&c1_rhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs58, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs58), "rhs", "rhs",
                  58);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs58), "lhs", "lhs",
                  58);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 59);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 59);
  sf_mex_assign(&c1_rhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs59, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs59), "rhs", "rhs",
                  59);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs59), "lhs", "lhs",
                  59);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 60);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 60);
  sf_mex_assign(&c1_rhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs60, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs60), "rhs", "rhs",
                  60);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs60), "lhs", "lhs",
                  60);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 61);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 61);
  sf_mex_assign(&c1_rhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs61, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs61), "rhs", "rhs",
                  61);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs61), "lhs", "lhs",
                  61);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 62);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 62);
  sf_mex_assign(&c1_rhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs62, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs62), "rhs", "rhs",
                  62);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs62), "lhs", "lhs",
                  62);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 63);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 63);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 63);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 63);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 63);
  sf_mex_assign(&c1_rhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs63, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs63), "rhs", "rhs",
                  63);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs63), "lhs", "lhs",
                  63);
  sf_mex_destroy(&c1_rhs0);
  sf_mex_destroy(&c1_lhs0);
  sf_mex_destroy(&c1_rhs1);
  sf_mex_destroy(&c1_lhs1);
  sf_mex_destroy(&c1_rhs2);
  sf_mex_destroy(&c1_lhs2);
  sf_mex_destroy(&c1_rhs3);
  sf_mex_destroy(&c1_lhs3);
  sf_mex_destroy(&c1_rhs4);
  sf_mex_destroy(&c1_lhs4);
  sf_mex_destroy(&c1_rhs5);
  sf_mex_destroy(&c1_lhs5);
  sf_mex_destroy(&c1_rhs6);
  sf_mex_destroy(&c1_lhs6);
  sf_mex_destroy(&c1_rhs7);
  sf_mex_destroy(&c1_lhs7);
  sf_mex_destroy(&c1_rhs8);
  sf_mex_destroy(&c1_lhs8);
  sf_mex_destroy(&c1_rhs9);
  sf_mex_destroy(&c1_lhs9);
  sf_mex_destroy(&c1_rhs10);
  sf_mex_destroy(&c1_lhs10);
  sf_mex_destroy(&c1_rhs11);
  sf_mex_destroy(&c1_lhs11);
  sf_mex_destroy(&c1_rhs12);
  sf_mex_destroy(&c1_lhs12);
  sf_mex_destroy(&c1_rhs13);
  sf_mex_destroy(&c1_lhs13);
  sf_mex_destroy(&c1_rhs14);
  sf_mex_destroy(&c1_lhs14);
  sf_mex_destroy(&c1_rhs15);
  sf_mex_destroy(&c1_lhs15);
  sf_mex_destroy(&c1_rhs16);
  sf_mex_destroy(&c1_lhs16);
  sf_mex_destroy(&c1_rhs17);
  sf_mex_destroy(&c1_lhs17);
  sf_mex_destroy(&c1_rhs18);
  sf_mex_destroy(&c1_lhs18);
  sf_mex_destroy(&c1_rhs19);
  sf_mex_destroy(&c1_lhs19);
  sf_mex_destroy(&c1_rhs20);
  sf_mex_destroy(&c1_lhs20);
  sf_mex_destroy(&c1_rhs21);
  sf_mex_destroy(&c1_lhs21);
  sf_mex_destroy(&c1_rhs22);
  sf_mex_destroy(&c1_lhs22);
  sf_mex_destroy(&c1_rhs23);
  sf_mex_destroy(&c1_lhs23);
  sf_mex_destroy(&c1_rhs24);
  sf_mex_destroy(&c1_lhs24);
  sf_mex_destroy(&c1_rhs25);
  sf_mex_destroy(&c1_lhs25);
  sf_mex_destroy(&c1_rhs26);
  sf_mex_destroy(&c1_lhs26);
  sf_mex_destroy(&c1_rhs27);
  sf_mex_destroy(&c1_lhs27);
  sf_mex_destroy(&c1_rhs28);
  sf_mex_destroy(&c1_lhs28);
  sf_mex_destroy(&c1_rhs29);
  sf_mex_destroy(&c1_lhs29);
  sf_mex_destroy(&c1_rhs30);
  sf_mex_destroy(&c1_lhs30);
  sf_mex_destroy(&c1_rhs31);
  sf_mex_destroy(&c1_lhs31);
  sf_mex_destroy(&c1_rhs32);
  sf_mex_destroy(&c1_lhs32);
  sf_mex_destroy(&c1_rhs33);
  sf_mex_destroy(&c1_lhs33);
  sf_mex_destroy(&c1_rhs34);
  sf_mex_destroy(&c1_lhs34);
  sf_mex_destroy(&c1_rhs35);
  sf_mex_destroy(&c1_lhs35);
  sf_mex_destroy(&c1_rhs36);
  sf_mex_destroy(&c1_lhs36);
  sf_mex_destroy(&c1_rhs37);
  sf_mex_destroy(&c1_lhs37);
  sf_mex_destroy(&c1_rhs38);
  sf_mex_destroy(&c1_lhs38);
  sf_mex_destroy(&c1_rhs39);
  sf_mex_destroy(&c1_lhs39);
  sf_mex_destroy(&c1_rhs40);
  sf_mex_destroy(&c1_lhs40);
  sf_mex_destroy(&c1_rhs41);
  sf_mex_destroy(&c1_lhs41);
  sf_mex_destroy(&c1_rhs42);
  sf_mex_destroy(&c1_lhs42);
  sf_mex_destroy(&c1_rhs43);
  sf_mex_destroy(&c1_lhs43);
  sf_mex_destroy(&c1_rhs44);
  sf_mex_destroy(&c1_lhs44);
  sf_mex_destroy(&c1_rhs45);
  sf_mex_destroy(&c1_lhs45);
  sf_mex_destroy(&c1_rhs46);
  sf_mex_destroy(&c1_lhs46);
  sf_mex_destroy(&c1_rhs47);
  sf_mex_destroy(&c1_lhs47);
  sf_mex_destroy(&c1_rhs48);
  sf_mex_destroy(&c1_lhs48);
  sf_mex_destroy(&c1_rhs49);
  sf_mex_destroy(&c1_lhs49);
  sf_mex_destroy(&c1_rhs50);
  sf_mex_destroy(&c1_lhs50);
  sf_mex_destroy(&c1_rhs51);
  sf_mex_destroy(&c1_lhs51);
  sf_mex_destroy(&c1_rhs52);
  sf_mex_destroy(&c1_lhs52);
  sf_mex_destroy(&c1_rhs53);
  sf_mex_destroy(&c1_lhs53);
  sf_mex_destroy(&c1_rhs54);
  sf_mex_destroy(&c1_lhs54);
  sf_mex_destroy(&c1_rhs55);
  sf_mex_destroy(&c1_lhs55);
  sf_mex_destroy(&c1_rhs56);
  sf_mex_destroy(&c1_lhs56);
  sf_mex_destroy(&c1_rhs57);
  sf_mex_destroy(&c1_lhs57);
  sf_mex_destroy(&c1_rhs58);
  sf_mex_destroy(&c1_lhs58);
  sf_mex_destroy(&c1_rhs59);
  sf_mex_destroy(&c1_lhs59);
  sf_mex_destroy(&c1_rhs60);
  sf_mex_destroy(&c1_lhs60);
  sf_mex_destroy(&c1_rhs61);
  sf_mex_destroy(&c1_lhs61);
  sf_mex_destroy(&c1_rhs62);
  sf_mex_destroy(&c1_lhs62);
  sf_mex_destroy(&c1_rhs63);
  sf_mex_destroy(&c1_lhs63);
}

static const mxArray *c1_emlrt_marshallOut(const char * c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c1_u)), false);
  return c1_y;
}

static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u)
{
  const mxArray *c1_y = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 7, 0U, 0U, 0U, 0), false);
  return c1_y;
}

static void c1_b_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs64 = NULL;
  const mxArray *c1_lhs64 = NULL;
  const mxArray *c1_rhs65 = NULL;
  const mxArray *c1_lhs65 = NULL;
  const mxArray *c1_rhs66 = NULL;
  const mxArray *c1_lhs66 = NULL;
  const mxArray *c1_rhs67 = NULL;
  const mxArray *c1_lhs67 = NULL;
  const mxArray *c1_rhs68 = NULL;
  const mxArray *c1_lhs68 = NULL;
  const mxArray *c1_rhs69 = NULL;
  const mxArray *c1_lhs69 = NULL;
  const mxArray *c1_rhs70 = NULL;
  const mxArray *c1_lhs70 = NULL;
  const mxArray *c1_rhs71 = NULL;
  const mxArray *c1_lhs71 = NULL;
  const mxArray *c1_rhs72 = NULL;
  const mxArray *c1_lhs72 = NULL;
  const mxArray *c1_rhs73 = NULL;
  const mxArray *c1_lhs73 = NULL;
  const mxArray *c1_rhs74 = NULL;
  const mxArray *c1_lhs74 = NULL;
  const mxArray *c1_rhs75 = NULL;
  const mxArray *c1_lhs75 = NULL;
  const mxArray *c1_rhs76 = NULL;
  const mxArray *c1_lhs76 = NULL;
  const mxArray *c1_rhs77 = NULL;
  const mxArray *c1_lhs77 = NULL;
  const mxArray *c1_rhs78 = NULL;
  const mxArray *c1_lhs78 = NULL;
  const mxArray *c1_rhs79 = NULL;
  const mxArray *c1_lhs79 = NULL;
  const mxArray *c1_rhs80 = NULL;
  const mxArray *c1_lhs80 = NULL;
  const mxArray *c1_rhs81 = NULL;
  const mxArray *c1_lhs81 = NULL;
  const mxArray *c1_rhs82 = NULL;
  const mxArray *c1_lhs82 = NULL;
  const mxArray *c1_rhs83 = NULL;
  const mxArray *c1_lhs83 = NULL;
  const mxArray *c1_rhs84 = NULL;
  const mxArray *c1_lhs84 = NULL;
  const mxArray *c1_rhs85 = NULL;
  const mxArray *c1_lhs85 = NULL;
  const mxArray *c1_rhs86 = NULL;
  const mxArray *c1_lhs86 = NULL;
  const mxArray *c1_rhs87 = NULL;
  const mxArray *c1_lhs87 = NULL;
  const mxArray *c1_rhs88 = NULL;
  const mxArray *c1_lhs88 = NULL;
  const mxArray *c1_rhs89 = NULL;
  const mxArray *c1_lhs89 = NULL;
  const mxArray *c1_rhs90 = NULL;
  const mxArray *c1_lhs90 = NULL;
  const mxArray *c1_rhs91 = NULL;
  const mxArray *c1_lhs91 = NULL;
  const mxArray *c1_rhs92 = NULL;
  const mxArray *c1_lhs92 = NULL;
  const mxArray *c1_rhs93 = NULL;
  const mxArray *c1_lhs93 = NULL;
  const mxArray *c1_rhs94 = NULL;
  const mxArray *c1_lhs94 = NULL;
  const mxArray *c1_rhs95 = NULL;
  const mxArray *c1_lhs95 = NULL;
  const mxArray *c1_rhs96 = NULL;
  const mxArray *c1_lhs96 = NULL;
  const mxArray *c1_rhs97 = NULL;
  const mxArray *c1_lhs97 = NULL;
  const mxArray *c1_rhs98 = NULL;
  const mxArray *c1_lhs98 = NULL;
  const mxArray *c1_rhs99 = NULL;
  const mxArray *c1_lhs99 = NULL;
  const mxArray *c1_rhs100 = NULL;
  const mxArray *c1_lhs100 = NULL;
  const mxArray *c1_rhs101 = NULL;
  const mxArray *c1_lhs101 = NULL;
  const mxArray *c1_rhs102 = NULL;
  const mxArray *c1_lhs102 = NULL;
  const mxArray *c1_rhs103 = NULL;
  const mxArray *c1_lhs103 = NULL;
  const mxArray *c1_rhs104 = NULL;
  const mxArray *c1_lhs104 = NULL;
  const mxArray *c1_rhs105 = NULL;
  const mxArray *c1_lhs105 = NULL;
  const mxArray *c1_rhs106 = NULL;
  const mxArray *c1_lhs106 = NULL;
  const mxArray *c1_rhs107 = NULL;
  const mxArray *c1_lhs107 = NULL;
  const mxArray *c1_rhs108 = NULL;
  const mxArray *c1_lhs108 = NULL;
  const mxArray *c1_rhs109 = NULL;
  const mxArray *c1_lhs109 = NULL;
  const mxArray *c1_rhs110 = NULL;
  const mxArray *c1_lhs110 = NULL;
  const mxArray *c1_rhs111 = NULL;
  const mxArray *c1_lhs111 = NULL;
  const mxArray *c1_rhs112 = NULL;
  const mxArray *c1_lhs112 = NULL;
  const mxArray *c1_rhs113 = NULL;
  const mxArray *c1_lhs113 = NULL;
  const mxArray *c1_rhs114 = NULL;
  const mxArray *c1_lhs114 = NULL;
  const mxArray *c1_rhs115 = NULL;
  const mxArray *c1_lhs115 = NULL;
  const mxArray *c1_rhs116 = NULL;
  const mxArray *c1_lhs116 = NULL;
  const mxArray *c1_rhs117 = NULL;
  const mxArray *c1_lhs117 = NULL;
  const mxArray *c1_rhs118 = NULL;
  const mxArray *c1_lhs118 = NULL;
  const mxArray *c1_rhs119 = NULL;
  const mxArray *c1_lhs119 = NULL;
  const mxArray *c1_rhs120 = NULL;
  const mxArray *c1_lhs120 = NULL;
  const mxArray *c1_rhs121 = NULL;
  const mxArray *c1_lhs121 = NULL;
  const mxArray *c1_rhs122 = NULL;
  const mxArray *c1_lhs122 = NULL;
  const mxArray *c1_rhs123 = NULL;
  const mxArray *c1_lhs123 = NULL;
  const mxArray *c1_rhs124 = NULL;
  const mxArray *c1_lhs124 = NULL;
  const mxArray *c1_rhs125 = NULL;
  const mxArray *c1_lhs125 = NULL;
  const mxArray *c1_rhs126 = NULL;
  const mxArray *c1_lhs126 = NULL;
  const mxArray *c1_rhs127 = NULL;
  const mxArray *c1_lhs127 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 64);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 64);
  sf_mex_assign(&c1_rhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs64, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs64), "rhs", "rhs",
                  64);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs64), "lhs", "lhs",
                  64);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("my_tr2rpy"), "name", "name",
                  65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "resolved", "resolved", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464800136U), "fileTimeLo",
                  "fileTimeLo", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 65);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 65);
  sf_mex_assign(&c1_rhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs65, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs65), "rhs", "rhs",
                  65);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs65), "lhs", "lhs",
                  65);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("length"), "name", "name", 66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1303146206U), "fileTimeLo",
                  "fileTimeLo", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 66);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 66);
  sf_mex_assign(&c1_rhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs66, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs66), "rhs", "rhs",
                  66);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs66), "lhs", "lhs",
                  66);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713852U), "fileTimeLo",
                  "fileTimeLo", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 67);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 67);
  sf_mex_assign(&c1_rhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs67, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs67), "rhs", "rhs",
                  67);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs67), "lhs", "lhs",
                  67);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 68);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 68);
  sf_mex_assign(&c1_rhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs68, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs68), "rhs", "rhs",
                  68);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs68), "lhs", "lhs",
                  68);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("atan2"), "name", "name", 69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m"), "resolved",
                  "resolved", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
                  "fileTimeLo", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 69);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 69);
  sf_mex_assign(&c1_rhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs69, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs69), "rhs", "rhs",
                  69);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs69), "lhs", "lhs",
                  69);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m"), "context",
                  "context", 70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 70);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 70);
  sf_mex_assign(&c1_rhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs70, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs70), "rhs", "rhs",
                  70);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs70), "lhs", "lhs",
                  70);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m"), "context",
                  "context", 71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 71);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 71);
  sf_mex_assign(&c1_rhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs71, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs71), "rhs", "rhs",
                  71);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs71), "lhs", "lhs",
                  71);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 72);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 72);
  sf_mex_assign(&c1_rhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs72, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs72), "rhs", "rhs",
                  72);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs72), "lhs", "lhs",
                  72);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m"), "context",
                  "context", 73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_atan2"), "name",
                  "name", 73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m"),
                  "resolved", "resolved", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818720U), "fileTimeLo",
                  "fileTimeLo", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 73);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 73);
  sf_mex_assign(&c1_rhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs73, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs73), "rhs", "rhs",
                  73);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs73), "lhs", "lhs",
                  73);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sin"), "name", "name", 74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328504U), "fileTimeLo",
                  "fileTimeLo", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 74);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 74);
  sf_mex_assign(&c1_rhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs74, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs74), "rhs", "rhs",
                  74);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs74), "lhs", "lhs",
                  74);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818736U), "fileTimeLo",
                  "fileTimeLo", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 75);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 75);
  sf_mex_assign(&c1_rhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs75, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs75), "rhs", "rhs",
                  75);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs75), "lhs", "lhs",
                  75);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cos"), "name", "name", 76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
                  "fileTimeLo", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 76);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 76);
  sf_mex_assign(&c1_rhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs76, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs76), "rhs", "rhs",
                  76);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs76), "lhs", "lhs",
                  76);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818722U), "fileTimeLo",
                  "fileTimeLo", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 77);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 77);
  sf_mex_assign(&c1_rhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs77, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs77), "rhs", "rhs",
                  77);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs77), "lhs", "lhs",
                  77);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 78);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 78);
  sf_mex_assign(&c1_rhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs78, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs78), "rhs", "rhs",
                  78);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs78), "lhs", "lhs",
                  78);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name", 79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807648U), "fileTimeLo",
                  "fileTimeLo", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370009886U), "mFileTimeLo",
                  "mFileTimeLo", 79);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 79);
  sf_mex_assign(&c1_rhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs79, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs79), "rhs", "rhs",
                  79);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs79), "lhs", "lhs",
                  79);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 80);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 80);
  sf_mex_assign(&c1_rhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs80, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs80), "rhs", "rhs",
                  80);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs80), "lhs", "lhs",
                  80);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713880U), "fileTimeLo",
                  "fileTimeLo", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 81);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 81);
  sf_mex_assign(&c1_rhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs81, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs81), "rhs", "rhs",
                  81);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs81), "lhs", "lhs",
                  81);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cardatol"), "name", "name", 82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/cardatol.m"), "resolved", "resolved", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464800136U), "fileTimeLo",
                  "fileTimeLo", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 82);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 82);
  sf_mex_assign(&c1_rhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs82, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs82), "rhs", "rhs",
                  82);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs82), "lhs", "lhs",
                  82);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/cardatol.m"), "context", "context", 83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sin"), "name", "name", 83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328504U), "fileTimeLo",
                  "fileTimeLo", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 83);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 83);
  sf_mex_assign(&c1_rhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs83, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs83), "rhs", "rhs",
                  83);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs83), "lhs", "lhs",
                  83);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/cardatol.m"), "context", "context", 84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cos"), "name", "name", 84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
                  "fileTimeLo", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 84);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 84);
  sf_mex_assign(&c1_rhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs84, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs84), "rhs", "rhs",
                  84);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs84), "lhs", "lhs",
                  84);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/cardatol.m"), "context", "context", 85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("error"), "name", "name", 85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/lang/error.m"), "resolved",
                  "resolved", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1319729966U), "fileTimeLo",
                  "fileTimeLo", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 85);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 85);
  sf_mex_assign(&c1_rhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs85, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs85), "rhs", "rhs",
                  85);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs85), "lhs", "lhs",
                  85);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/cardatol.m"), "context", "context", 86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rem"), "name", "name", 86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/rem.m"), "resolved",
                  "resolved", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713854U), "fileTimeLo",
                  "fileTimeLo", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 86);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 86);
  sf_mex_assign(&c1_rhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs86, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs86), "rhs", "rhs",
                  86);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs86), "lhs", "lhs",
                  86);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/rem.m"), "context",
                  "context", 87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 87);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 87);
  sf_mex_assign(&c1_rhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs87, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs87), "rhs", "rhs",
                  87);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs87), "lhs", "lhs",
                  87);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/rem.m"), "context",
                  "context", 88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 88);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 88);
  sf_mex_assign(&c1_rhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs88, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs88), "rhs", "rhs",
                  88);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs88), "lhs", "lhs",
                  88);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/rem.m"), "context",
                  "context", 89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 89);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 89);
  sf_mex_assign(&c1_rhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs89, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs89), "rhs", "rhs",
                  89);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs89), "lhs", "lhs",
                  89);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eye"), "name", "name", 90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "resolved",
                  "resolved", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1406813148U), "fileTimeLo",
                  "fileTimeLo", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 90);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 90);
  sf_mex_assign(&c1_rhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs90, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs90), "rhs", "rhs",
                  90);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs90), "lhs", "lhs",
                  90);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_assert_valid_size_arg"),
                  "name", "name", 91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "resolved", "resolved", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1368183030U), "fileTimeLo",
                  "fileTimeLo", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 91);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 91);
  sf_mex_assign(&c1_rhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs91, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs91), "rhs", "rhs",
                  91);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs91), "lhs", "lhs",
                  91);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 92);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 92);
  sf_mex_assign(&c1_rhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs92, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs92), "rhs", "rhs",
                  92);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs92), "lhs", "lhs",
                  92);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral"),
                  "context", "context", 93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isinf"), "name", "name", 93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "resolved",
                  "resolved", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713856U), "fileTimeLo",
                  "fileTimeLo", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 93);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 93);
  sf_mex_assign(&c1_rhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs93, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs93), "rhs", "rhs",
                  93);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs93), "lhs", "lhs",
                  93);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m"), "context",
                  "context", 94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 94);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 94);
  sf_mex_assign(&c1_rhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs94, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs94), "rhs", "rhs",
                  94);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs94), "lhs", "lhs",
                  94);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_is_integer_class"), "name",
                  "name", 95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_integer_class.m"),
                  "resolved", "resolved", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818782U), "fileTimeLo",
                  "fileTimeLo", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 95);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 95);
  sf_mex_assign(&c1_rhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs95, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs95), "rhs", "rhs",
                  95);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs95), "lhs", "lhs",
                  95);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 96);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 96);
  sf_mex_assign(&c1_rhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs96, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs96), "rhs", "rhs",
                  96);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs96), "lhs", "lhs",
                  96);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 97);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 97);
  sf_mex_assign(&c1_rhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs97, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs97), "rhs", "rhs",
                  97);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs97), "lhs", "lhs",
                  97);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isinbounds"),
                  "context", "context", 98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326728322U), "fileTimeLo",
                  "fileTimeLo", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 98);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 98);
  sf_mex_assign(&c1_rhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs98, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs98), "rhs", "rhs",
                  98);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs98), "lhs", "lhs",
                  98);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!apply_float_relop"),
                  "context", "context", 99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 99);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 99);
  sf_mex_assign(&c1_rhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs99, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs99), "rhs", "rhs",
                  99);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs99), "lhs", "lhs",
                  99);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!float_class_contains_indexIntClass"),
                  "context", "context", 100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 100);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 100);
  sf_mex_assign(&c1_rhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs100, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs100), "rhs", "rhs",
                  100);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs100), "lhs", "lhs",
                  100);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m!is_signed_indexIntClass"),
                  "context", "context", 101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 101);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 101);
  sf_mex_assign(&c1_rhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs101, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs101), "rhs", "rhs",
                  101);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs101), "lhs", "lhs",
                  101);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 102);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 102);
  sf_mex_assign(&c1_rhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs102, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs102), "rhs", "rhs",
                  102);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs102), "lhs", "lhs",
                  102);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"),
                  "context", "context", 103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 103);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 103);
  sf_mex_assign(&c1_rhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs103, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs103), "rhs", "rhs",
                  103);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs103), "lhs", "lhs",
                  103);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m"), "context",
                  "context", 104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 104);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 104);
  sf_mex_assign(&c1_rhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs104, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs104), "rhs", "rhs",
                  104);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs104), "lhs", "lhs",
                  104);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mldivide"), "name", "name",
                  105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "resolved",
                  "resolved", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807648U), "fileTimeLo",
                  "fileTimeLo", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1319729966U), "mFileTimeLo",
                  "mFileTimeLo", 105);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 105);
  sf_mex_assign(&c1_rhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs105, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs105), "rhs", "rhs",
                  105);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs105), "lhs", "lhs",
                  105);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p"), "context",
                  "context", 106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_lusolve"), "name", "name",
                  106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m"), "resolved",
                  "resolved", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370009886U), "fileTimeLo",
                  "fileTimeLo", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 106);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 106);
  sf_mex_assign(&c1_rhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs106, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs106), "rhs", "rhs",
                  106);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs106), "lhs", "lhs",
                  106);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xcabs1"), "name", "name",
                  107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "resolved", "resolved", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 107);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 107);
  sf_mex_assign(&c1_rhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs107, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs107), "rhs", "rhs",
                  107);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs107), "lhs", "lhs",
                  107);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"),
                  "context", "context", 108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xcabs1"),
                  "name", "name", 108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "resolved", "resolved", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 108);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 108);
  sf_mex_assign(&c1_rhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs108, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs108), "rhs", "rhs",
                  108);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs108), "lhs", "lhs",
                  108);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "context", "context", 109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713852U), "fileTimeLo",
                  "fileTimeLo", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 109);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 109);
  sf_mex_assign(&c1_rhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs109, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs109), "rhs", "rhs",
                  109);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs109), "lhs", "lhs",
                  109);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713880U), "fileTimeLo",
                  "fileTimeLo", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 110);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 110);
  sf_mex_assign(&c1_rhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs110, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs110), "rhs", "rhs",
                  110);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs110), "lhs", "lhs",
                  110);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular"),
                  "context", "context", 111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_warning"), "name", "name",
                  111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m"), "resolved",
                  "resolved", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818802U), "fileTimeLo",
                  "fileTimeLo", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 111);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 111);
  sf_mex_assign(&c1_rhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs111, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs111), "rhs", "rhs",
                  111);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs111), "lhs", "lhs",
                  111);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolve3x3"),
                  "context", "context", 112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 112);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 112);
  sf_mex_assign(&c1_rhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs112, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs112), "rhs", "rhs",
                  112);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs112), "lhs", "lhs",
                  112);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 113);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 113);
  sf_mex_assign(&c1_rhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs113, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs113), "rhs", "rhs",
                  113);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs113), "lhs", "lhs",
                  113);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgetrf"), "name", "name",
                  114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "resolved", "resolved", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818806U), "fileTimeLo",
                  "fileTimeLo", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 114);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 114);
  sf_mex_assign(&c1_rhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs114, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs114), "rhs", "rhs",
                  114);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs114), "lhs", "lhs",
                  114);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"),
                  "context", "context", 115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_lapack_xgetrf"), "name",
                  "name", 115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "resolved", "resolved", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818810U), "fileTimeLo",
                  "fileTimeLo", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 115);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 115);
  sf_mex_assign(&c1_rhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs115, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs115), "rhs", "rhs",
                  115);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs115), "lhs", "lhs",
                  115);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"),
                  "context", "context", 116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_matlab_zgetrf"), "name",
                  "name", 116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "resolved", "resolved", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1302688994U), "fileTimeLo",
                  "fileTimeLo", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 116);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 116);
  sf_mex_assign(&c1_rhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs116, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs116), "rhs", "rhs",
                  116);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs116), "lhs", "lhs",
                  116);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("realmin"), "name", "name", 117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307651242U), "fileTimeLo",
                  "fileTimeLo", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 117);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 117);
  sf_mex_assign(&c1_rhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs117, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs117), "rhs", "rhs",
                  117);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs117), "lhs", "lhs",
                  117);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_realmin"), "name", "name",
                  118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1307651244U), "fileTimeLo",
                  "fileTimeLo", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 118);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 118);
  sf_mex_assign(&c1_rhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs118, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs118), "rhs", "rhs",
                  118);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs118), "lhs", "lhs",
                  118);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 119);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 119);
  sf_mex_assign(&c1_rhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs119, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs119), "rhs", "rhs",
                  119);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs119), "lhs", "lhs",
                  119);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 120);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 120);
  sf_mex_assign(&c1_rhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs120, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs120), "rhs", "rhs",
                  120);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs120), "lhs", "lhs",
                  120);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("min"), "name", "name", 121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1311255318U), "fileTimeLo",
                  "fileTimeLo", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 121);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 121);
  sf_mex_assign(&c1_rhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs121, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs121), "rhs", "rhs",
                  121);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs121), "lhs", "lhs",
                  121);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "context",
                  "context", 122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_min_or_max"), "name",
                  "name", 122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"),
                  "resolved", "resolved", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1378295984U), "fileTimeLo",
                  "fileTimeLo", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 122);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 122);
  sf_mex_assign(&c1_rhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs122, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs122), "rhs", "rhs",
                  122);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs122), "lhs", "lhs",
                  122);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 123);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 123);
  sf_mex_assign(&c1_rhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs123, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs123), "rhs", "rhs",
                  123);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs123), "lhs", "lhs",
                  123);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 124);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 124);
  sf_mex_assign(&c1_rhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs124, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs124), "rhs", "rhs",
                  124);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs124), "lhs", "lhs",
                  124);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 125);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 125);
  sf_mex_assign(&c1_rhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs125, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs125), "rhs", "rhs",
                  125);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs125), "lhs", "lhs",
                  125);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 126);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 126);
  sf_mex_assign(&c1_rhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs126, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs126), "rhs", "rhs",
                  126);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs126), "lhs", "lhs",
                  126);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 127);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 127);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 127);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 127);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 127);
  sf_mex_assign(&c1_rhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs127, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs127), "rhs", "rhs",
                  127);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs127), "lhs", "lhs",
                  127);
  sf_mex_destroy(&c1_rhs64);
  sf_mex_destroy(&c1_lhs64);
  sf_mex_destroy(&c1_rhs65);
  sf_mex_destroy(&c1_lhs65);
  sf_mex_destroy(&c1_rhs66);
  sf_mex_destroy(&c1_lhs66);
  sf_mex_destroy(&c1_rhs67);
  sf_mex_destroy(&c1_lhs67);
  sf_mex_destroy(&c1_rhs68);
  sf_mex_destroy(&c1_lhs68);
  sf_mex_destroy(&c1_rhs69);
  sf_mex_destroy(&c1_lhs69);
  sf_mex_destroy(&c1_rhs70);
  sf_mex_destroy(&c1_lhs70);
  sf_mex_destroy(&c1_rhs71);
  sf_mex_destroy(&c1_lhs71);
  sf_mex_destroy(&c1_rhs72);
  sf_mex_destroy(&c1_lhs72);
  sf_mex_destroy(&c1_rhs73);
  sf_mex_destroy(&c1_lhs73);
  sf_mex_destroy(&c1_rhs74);
  sf_mex_destroy(&c1_lhs74);
  sf_mex_destroy(&c1_rhs75);
  sf_mex_destroy(&c1_lhs75);
  sf_mex_destroy(&c1_rhs76);
  sf_mex_destroy(&c1_lhs76);
  sf_mex_destroy(&c1_rhs77);
  sf_mex_destroy(&c1_lhs77);
  sf_mex_destroy(&c1_rhs78);
  sf_mex_destroy(&c1_lhs78);
  sf_mex_destroy(&c1_rhs79);
  sf_mex_destroy(&c1_lhs79);
  sf_mex_destroy(&c1_rhs80);
  sf_mex_destroy(&c1_lhs80);
  sf_mex_destroy(&c1_rhs81);
  sf_mex_destroy(&c1_lhs81);
  sf_mex_destroy(&c1_rhs82);
  sf_mex_destroy(&c1_lhs82);
  sf_mex_destroy(&c1_rhs83);
  sf_mex_destroy(&c1_lhs83);
  sf_mex_destroy(&c1_rhs84);
  sf_mex_destroy(&c1_lhs84);
  sf_mex_destroy(&c1_rhs85);
  sf_mex_destroy(&c1_lhs85);
  sf_mex_destroy(&c1_rhs86);
  sf_mex_destroy(&c1_lhs86);
  sf_mex_destroy(&c1_rhs87);
  sf_mex_destroy(&c1_lhs87);
  sf_mex_destroy(&c1_rhs88);
  sf_mex_destroy(&c1_lhs88);
  sf_mex_destroy(&c1_rhs89);
  sf_mex_destroy(&c1_lhs89);
  sf_mex_destroy(&c1_rhs90);
  sf_mex_destroy(&c1_lhs90);
  sf_mex_destroy(&c1_rhs91);
  sf_mex_destroy(&c1_lhs91);
  sf_mex_destroy(&c1_rhs92);
  sf_mex_destroy(&c1_lhs92);
  sf_mex_destroy(&c1_rhs93);
  sf_mex_destroy(&c1_lhs93);
  sf_mex_destroy(&c1_rhs94);
  sf_mex_destroy(&c1_lhs94);
  sf_mex_destroy(&c1_rhs95);
  sf_mex_destroy(&c1_lhs95);
  sf_mex_destroy(&c1_rhs96);
  sf_mex_destroy(&c1_lhs96);
  sf_mex_destroy(&c1_rhs97);
  sf_mex_destroy(&c1_lhs97);
  sf_mex_destroy(&c1_rhs98);
  sf_mex_destroy(&c1_lhs98);
  sf_mex_destroy(&c1_rhs99);
  sf_mex_destroy(&c1_lhs99);
  sf_mex_destroy(&c1_rhs100);
  sf_mex_destroy(&c1_lhs100);
  sf_mex_destroy(&c1_rhs101);
  sf_mex_destroy(&c1_lhs101);
  sf_mex_destroy(&c1_rhs102);
  sf_mex_destroy(&c1_lhs102);
  sf_mex_destroy(&c1_rhs103);
  sf_mex_destroy(&c1_lhs103);
  sf_mex_destroy(&c1_rhs104);
  sf_mex_destroy(&c1_lhs104);
  sf_mex_destroy(&c1_rhs105);
  sf_mex_destroy(&c1_lhs105);
  sf_mex_destroy(&c1_rhs106);
  sf_mex_destroy(&c1_lhs106);
  sf_mex_destroy(&c1_rhs107);
  sf_mex_destroy(&c1_lhs107);
  sf_mex_destroy(&c1_rhs108);
  sf_mex_destroy(&c1_lhs108);
  sf_mex_destroy(&c1_rhs109);
  sf_mex_destroy(&c1_lhs109);
  sf_mex_destroy(&c1_rhs110);
  sf_mex_destroy(&c1_lhs110);
  sf_mex_destroy(&c1_rhs111);
  sf_mex_destroy(&c1_lhs111);
  sf_mex_destroy(&c1_rhs112);
  sf_mex_destroy(&c1_lhs112);
  sf_mex_destroy(&c1_rhs113);
  sf_mex_destroy(&c1_lhs113);
  sf_mex_destroy(&c1_rhs114);
  sf_mex_destroy(&c1_lhs114);
  sf_mex_destroy(&c1_rhs115);
  sf_mex_destroy(&c1_lhs115);
  sf_mex_destroy(&c1_rhs116);
  sf_mex_destroy(&c1_lhs116);
  sf_mex_destroy(&c1_rhs117);
  sf_mex_destroy(&c1_lhs117);
  sf_mex_destroy(&c1_rhs118);
  sf_mex_destroy(&c1_lhs118);
  sf_mex_destroy(&c1_rhs119);
  sf_mex_destroy(&c1_lhs119);
  sf_mex_destroy(&c1_rhs120);
  sf_mex_destroy(&c1_lhs120);
  sf_mex_destroy(&c1_rhs121);
  sf_mex_destroy(&c1_lhs121);
  sf_mex_destroy(&c1_rhs122);
  sf_mex_destroy(&c1_lhs122);
  sf_mex_destroy(&c1_rhs123);
  sf_mex_destroy(&c1_lhs123);
  sf_mex_destroy(&c1_rhs124);
  sf_mex_destroy(&c1_lhs124);
  sf_mex_destroy(&c1_rhs125);
  sf_mex_destroy(&c1_lhs125);
  sf_mex_destroy(&c1_rhs126);
  sf_mex_destroy(&c1_lhs126);
  sf_mex_destroy(&c1_rhs127);
  sf_mex_destroy(&c1_lhs127);
}

static void c1_c_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs128 = NULL;
  const mxArray *c1_lhs128 = NULL;
  const mxArray *c1_rhs129 = NULL;
  const mxArray *c1_lhs129 = NULL;
  const mxArray *c1_rhs130 = NULL;
  const mxArray *c1_lhs130 = NULL;
  const mxArray *c1_rhs131 = NULL;
  const mxArray *c1_lhs131 = NULL;
  const mxArray *c1_rhs132 = NULL;
  const mxArray *c1_lhs132 = NULL;
  const mxArray *c1_rhs133 = NULL;
  const mxArray *c1_lhs133 = NULL;
  const mxArray *c1_rhs134 = NULL;
  const mxArray *c1_lhs134 = NULL;
  const mxArray *c1_rhs135 = NULL;
  const mxArray *c1_lhs135 = NULL;
  const mxArray *c1_rhs136 = NULL;
  const mxArray *c1_lhs136 = NULL;
  const mxArray *c1_rhs137 = NULL;
  const mxArray *c1_lhs137 = NULL;
  const mxArray *c1_rhs138 = NULL;
  const mxArray *c1_lhs138 = NULL;
  const mxArray *c1_rhs139 = NULL;
  const mxArray *c1_lhs139 = NULL;
  const mxArray *c1_rhs140 = NULL;
  const mxArray *c1_lhs140 = NULL;
  const mxArray *c1_rhs141 = NULL;
  const mxArray *c1_lhs141 = NULL;
  const mxArray *c1_rhs142 = NULL;
  const mxArray *c1_lhs142 = NULL;
  const mxArray *c1_rhs143 = NULL;
  const mxArray *c1_lhs143 = NULL;
  const mxArray *c1_rhs144 = NULL;
  const mxArray *c1_lhs144 = NULL;
  const mxArray *c1_rhs145 = NULL;
  const mxArray *c1_lhs145 = NULL;
  const mxArray *c1_rhs146 = NULL;
  const mxArray *c1_lhs146 = NULL;
  const mxArray *c1_rhs147 = NULL;
  const mxArray *c1_lhs147 = NULL;
  const mxArray *c1_rhs148 = NULL;
  const mxArray *c1_lhs148 = NULL;
  const mxArray *c1_rhs149 = NULL;
  const mxArray *c1_lhs149 = NULL;
  const mxArray *c1_rhs150 = NULL;
  const mxArray *c1_lhs150 = NULL;
  const mxArray *c1_rhs151 = NULL;
  const mxArray *c1_lhs151 = NULL;
  const mxArray *c1_rhs152 = NULL;
  const mxArray *c1_lhs152 = NULL;
  const mxArray *c1_rhs153 = NULL;
  const mxArray *c1_lhs153 = NULL;
  const mxArray *c1_rhs154 = NULL;
  const mxArray *c1_lhs154 = NULL;
  const mxArray *c1_rhs155 = NULL;
  const mxArray *c1_lhs155 = NULL;
  const mxArray *c1_rhs156 = NULL;
  const mxArray *c1_lhs156 = NULL;
  const mxArray *c1_rhs157 = NULL;
  const mxArray *c1_lhs157 = NULL;
  const mxArray *c1_rhs158 = NULL;
  const mxArray *c1_lhs158 = NULL;
  const mxArray *c1_rhs159 = NULL;
  const mxArray *c1_lhs159 = NULL;
  const mxArray *c1_rhs160 = NULL;
  const mxArray *c1_lhs160 = NULL;
  const mxArray *c1_rhs161 = NULL;
  const mxArray *c1_lhs161 = NULL;
  const mxArray *c1_rhs162 = NULL;
  const mxArray *c1_lhs162 = NULL;
  const mxArray *c1_rhs163 = NULL;
  const mxArray *c1_lhs163 = NULL;
  const mxArray *c1_rhs164 = NULL;
  const mxArray *c1_lhs164 = NULL;
  const mxArray *c1_rhs165 = NULL;
  const mxArray *c1_lhs165 = NULL;
  const mxArray *c1_rhs166 = NULL;
  const mxArray *c1_lhs166 = NULL;
  const mxArray *c1_rhs167 = NULL;
  const mxArray *c1_lhs167 = NULL;
  const mxArray *c1_rhs168 = NULL;
  const mxArray *c1_lhs168 = NULL;
  const mxArray *c1_rhs169 = NULL;
  const mxArray *c1_lhs169 = NULL;
  const mxArray *c1_rhs170 = NULL;
  const mxArray *c1_lhs170 = NULL;
  const mxArray *c1_rhs171 = NULL;
  const mxArray *c1_lhs171 = NULL;
  const mxArray *c1_rhs172 = NULL;
  const mxArray *c1_lhs172 = NULL;
  const mxArray *c1_rhs173 = NULL;
  const mxArray *c1_lhs173 = NULL;
  const mxArray *c1_rhs174 = NULL;
  const mxArray *c1_lhs174 = NULL;
  const mxArray *c1_rhs175 = NULL;
  const mxArray *c1_lhs175 = NULL;
  const mxArray *c1_rhs176 = NULL;
  const mxArray *c1_lhs176 = NULL;
  const mxArray *c1_rhs177 = NULL;
  const mxArray *c1_lhs177 = NULL;
  const mxArray *c1_rhs178 = NULL;
  const mxArray *c1_lhs178 = NULL;
  const mxArray *c1_rhs179 = NULL;
  const mxArray *c1_lhs179 = NULL;
  const mxArray *c1_rhs180 = NULL;
  const mxArray *c1_lhs180 = NULL;
  const mxArray *c1_rhs181 = NULL;
  const mxArray *c1_lhs181 = NULL;
  const mxArray *c1_rhs182 = NULL;
  const mxArray *c1_lhs182 = NULL;
  const mxArray *c1_rhs183 = NULL;
  const mxArray *c1_lhs183 = NULL;
  const mxArray *c1_rhs184 = NULL;
  const mxArray *c1_lhs184 = NULL;
  const mxArray *c1_rhs185 = NULL;
  const mxArray *c1_lhs185 = NULL;
  const mxArray *c1_rhs186 = NULL;
  const mxArray *c1_lhs186 = NULL;
  const mxArray *c1_rhs187 = NULL;
  const mxArray *c1_lhs187 = NULL;
  const mxArray *c1_rhs188 = NULL;
  const mxArray *c1_lhs188 = NULL;
  const mxArray *c1_rhs189 = NULL;
  const mxArray *c1_lhs189 = NULL;
  const mxArray *c1_rhs190 = NULL;
  const mxArray *c1_lhs190 = NULL;
  const mxArray *c1_rhs191 = NULL;
  const mxArray *c1_lhs191 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 128);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 128);
  sf_mex_assign(&c1_rhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs128, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs128), "rhs", "rhs",
                  128);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs128), "lhs", "lhs",
                  128);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 129);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 129);
  sf_mex_assign(&c1_rhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs129, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs129), "rhs", "rhs",
                  129);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs129), "lhs", "lhs",
                  129);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("colon"), "name", "name", 130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1378295988U), "fileTimeLo",
                  "fileTimeLo", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 130);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 130);
  sf_mex_assign(&c1_rhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs130, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs130), "rhs", "rhs",
                  130);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs130), "lhs", "lhs",
                  130);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("colon"), "name", "name", 131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1378295988U), "fileTimeLo",
                  "fileTimeLo", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 131);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 131);
  sf_mex_assign(&c1_rhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs131, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs131), "rhs", "rhs",
                  131);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs131), "lhs", "lhs",
                  131);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 132);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 132);
  sf_mex_assign(&c1_rhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs132, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs132), "rhs", "rhs",
                  132);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs132), "lhs", "lhs",
                  132);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 133);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 133);
  sf_mex_assign(&c1_rhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs133, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs133), "rhs", "rhs",
                  133);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs133), "lhs", "lhs",
                  133);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "context",
                  "context", 134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("floor"), "name", "name", 134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713854U), "fileTimeLo",
                  "fileTimeLo", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 134);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 134);
  sf_mex_assign(&c1_rhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs134, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs134), "rhs", "rhs",
                  134);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs134), "lhs", "lhs",
                  134);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 135);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 135);
  sf_mex_assign(&c1_rhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs135, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs135), "rhs", "rhs",
                  135);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs135), "lhs", "lhs",
                  135);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818726U), "fileTimeLo",
                  "fileTimeLo", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 136);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 136);
  sf_mex_assign(&c1_rhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs136, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs136), "rhs", "rhs",
                  136);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs136), "lhs", "lhs",
                  136);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 137);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 137);
  sf_mex_assign(&c1_rhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs137, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs137), "rhs", "rhs",
                  137);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs137), "lhs", "lhs",
                  137);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange"),
                  "context", "context", 138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 138);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 138);
  sf_mex_assign(&c1_rhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs138, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs138), "rhs", "rhs",
                  138);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs138), "lhs", "lhs",
                  138);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmin"), "name", "name", 139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 139);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 139);
  sf_mex_assign(&c1_rhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs139, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs139), "rhs", "rhs",
                  139);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs139), "lhs", "lhs",
                  139);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 140);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 140);
  sf_mex_assign(&c1_rhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs140, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs140), "rhs", "rhs",
                  140);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs140), "lhs", "lhs",
                  140);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher"),
                  "context", "context", 141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_isa_uint"), "name", "name",
                  141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 141);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 141);
  sf_mex_assign(&c1_rhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs141, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs141), "rhs", "rhs",
                  141);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs141), "lhs", "lhs",
                  141);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "context",
                  "context", 142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.isaUint"),
                  "name", "name", 142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/isaUint.p"),
                  "resolved", "resolved", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 142);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 142);
  sf_mex_assign(&c1_rhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs142, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs142), "rhs", "rhs",
                  142);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs142), "lhs", "lhs",
                  142);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_unsigned_class"), "name",
                  "name", 143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "resolved", "resolved", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 143);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 143);
  sf_mex_assign(&c1_rhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs143, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs143), "rhs", "rhs",
                  143);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs143), "lhs", "lhs",
                  143);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m"),
                  "context", "context", 144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.unsignedClass"),
                  "name", "name", 144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "resolved", "resolved", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 144);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 144);
  sf_mex_assign(&c1_rhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs144, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs144), "rhs", "rhs",
                  144);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs144), "lhs", "lhs",
                  144);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "context", "context", 145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 145);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 145);
  sf_mex_assign(&c1_rhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs145, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs145), "rhs", "rhs",
                  145);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs145), "lhs", "lhs",
                  145);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/unsignedClass.p"),
                  "context", "context", 146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 146);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 146);
  sf_mex_assign(&c1_rhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs146, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs146), "rhs", "rhs",
                  146);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs146), "lhs", "lhs",
                  146);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 147);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 147);
  sf_mex_assign(&c1_rhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs147, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs147), "rhs", "rhs",
                  147);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs147), "lhs", "lhs",
                  147);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 148);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 148);
  sf_mex_assign(&c1_rhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs148, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs148), "rhs", "rhs",
                  148);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs148), "lhs", "lhs",
                  148);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_isa_uint"), "name", "name",
                  149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m"), "resolved",
                  "resolved", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 149);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 149);
  sf_mex_assign(&c1_rhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs149, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs149), "rhs", "rhs",
                  149);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs149), "lhs", "lhs",
                  149);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd"),
                  "context", "context", 150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 150);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 150);
  sf_mex_assign(&c1_rhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs150, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs150), "rhs", "rhs",
                  150);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs150), "lhs", "lhs",
                  150);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_signed_integer_colon"),
                  "context", "context", 151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 151);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 151);
  sf_mex_assign(&c1_rhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs151, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs151), "rhs", "rhs",
                  151);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs151), "lhs", "lhs",
                  151);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 152);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 152);
  sf_mex_assign(&c1_rhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs152, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs152), "rhs", "rhs",
                  152);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs152), "lhs", "lhs",
                  152);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 153);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 153);
  sf_mex_assign(&c1_rhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs153, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs153), "rhs", "rhs",
                  153);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs153), "lhs", "lhs",
                  153);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 154);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 154);
  sf_mex_assign(&c1_rhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs154, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs154), "rhs", "rhs",
                  154);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs154), "lhs", "lhs",
                  154);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 155);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 155);
  sf_mex_assign(&c1_rhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs155, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs155), "rhs", "rhs",
                  155);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs155), "lhs", "lhs",
                  155);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 156);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 156);
  sf_mex_assign(&c1_rhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs156, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs156), "rhs", "rhs",
                  156);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs156), "lhs", "lhs",
                  156);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 157);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 157);
  sf_mex_assign(&c1_rhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs157, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs157), "rhs", "rhs",
                  157);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs157), "lhs", "lhs",
                  157);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 158);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 158);
  sf_mex_assign(&c1_rhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs158, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs158), "rhs", "rhs",
                  158);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs158), "lhs", "lhs",
                  158);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 159);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 159);
  sf_mex_assign(&c1_rhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs159, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs159), "rhs", "rhs",
                  159);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs159), "lhs", "lhs",
                  159);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 160);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 160);
  sf_mex_assign(&c1_rhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs160, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs160), "rhs", "rhs",
                  160);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs160), "lhs", "lhs",
                  160);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 161);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 161);
  sf_mex_assign(&c1_rhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs161, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs161), "rhs", "rhs",
                  161);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs161), "lhs", "lhs",
                  161);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 162);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 162);
  sf_mex_assign(&c1_rhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs162, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs162), "rhs", "rhs",
                  162);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs162), "lhs", "lhs",
                  162);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_ixamax"), "name", "name",
                  163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "resolved", "resolved", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 163);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 163);
  sf_mex_assign(&c1_rhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs163, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs163), "rhs", "rhs",
                  163);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs163), "lhs", "lhs",
                  163);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 164);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 164);
  sf_mex_assign(&c1_rhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs164, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs164), "rhs", "rhs",
                  164);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs164), "lhs", "lhs",
                  164);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"),
                  "context", "context", 165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.ixamax"),
                  "name", "name", 165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "resolved", "resolved", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 165);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 165);
  sf_mex_assign(&c1_rhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs165, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs165), "rhs", "rhs",
                  165);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs165), "lhs", "lhs",
                  165);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "context", "context", 166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 166);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 166);
  sf_mex_assign(&c1_rhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs166, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs166), "rhs", "rhs",
                  166);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs166), "lhs", "lhs",
                  166);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p!below_threshold"),
                  "context", "context", 167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 167);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 167);
  sf_mex_assign(&c1_rhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs167, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs167), "rhs", "rhs",
                  167);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs167), "lhs", "lhs",
                  167);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p!below_threshold"),
                  "context", "context", 168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("length"), "name", "name", 168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1303146206U), "fileTimeLo",
                  "fileTimeLo", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 168);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 168);
  sf_mex_assign(&c1_rhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs168, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs168), "rhs", "rhs",
                  168);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs168), "lhs", "lhs",
                  168);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength"),
                  "context", "context", 169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 169);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 169);
  sf_mex_assign(&c1_rhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs169, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs169), "rhs", "rhs",
                  169);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs169), "lhs", "lhs",
                  169);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/ixamax.p"),
                  "context", "context", 170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.ixamax"),
                  "name", "name", 170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "resolved", "resolved", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 170);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 170);
  sf_mex_assign(&c1_rhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs170, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs170), "rhs", "rhs",
                  170);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs170), "lhs", "lhs",
                  170);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xcabs1"),
                  "name", "name", 171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xcabs1.p"),
                  "resolved", "resolved", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 171);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 171);
  sf_mex_assign(&c1_rhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs171, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs171), "rhs", "rhs",
                  171);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs171), "lhs", "lhs",
                  171);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 172);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 172);
  sf_mex_assign(&c1_rhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs172, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs172), "rhs", "rhs",
                  172);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs172), "lhs", "lhs",
                  172);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/ixamax.p"),
                  "context", "context", 173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 173);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 173);
  sf_mex_assign(&c1_rhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs173, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs173), "rhs", "rhs",
                  173);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs173), "lhs", "lhs",
                  173);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xswap"), "name", "name",
                  174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"),
                  "resolved", "resolved", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980692U), "fileTimeLo",
                  "fileTimeLo", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 174);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 174);
  sf_mex_assign(&c1_rhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs174, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs174), "rhs", "rhs",
                  174);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs174), "lhs", "lhs",
                  174);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 175);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 175);
  sf_mex_assign(&c1_rhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs175, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs175), "rhs", "rhs",
                  175);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs175), "lhs", "lhs",
                  175);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"), "context",
                  "context", 176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xswap"),
                  "name", "name", 176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "resolved", "resolved", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 176);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 176);
  sf_mex_assign(&c1_rhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs176, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs176), "rhs", "rhs",
                  176);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs176), "lhs", "lhs",
                  176);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 177);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 177);
  sf_mex_assign(&c1_rhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs177, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs177), "rhs", "rhs",
                  177);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs177), "lhs", "lhs",
                  177);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p!below_threshold"),
                  "context", "context", 178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 178);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 178);
  sf_mex_assign(&c1_rhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs178, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs178), "rhs", "rhs",
                  178);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs178), "lhs", "lhs",
                  178);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xswap.p"),
                  "context", "context", 179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xswap"),
                  "name", "name", 179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "resolved", "resolved", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 179);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 179);
  sf_mex_assign(&c1_rhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs179, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs179), "rhs", "rhs",
                  179);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs179), "lhs", "lhs",
                  179);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713852U), "fileTimeLo",
                  "fileTimeLo", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 180);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 180);
  sf_mex_assign(&c1_rhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs180, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs180), "rhs", "rhs",
                  180);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs180), "lhs", "lhs",
                  180);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 181);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 181);
  sf_mex_assign(&c1_rhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs181, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs181), "rhs", "rhs",
                  181);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs181), "lhs", "lhs",
                  181);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818712U), "fileTimeLo",
                  "fileTimeLo", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 182);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 182);
  sf_mex_assign(&c1_rhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs182, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs182), "rhs", "rhs",
                  182);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs182), "lhs", "lhs",
                  182);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 183);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 183);
  sf_mex_assign(&c1_rhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs183, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs183), "rhs", "rhs",
                  183);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs183), "lhs", "lhs",
                  183);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xswap.p"),
                  "context", "context", 184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 184);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 184);
  sf_mex_assign(&c1_rhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs184, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs184), "rhs", "rhs",
                  184);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs184), "lhs", "lhs",
                  184);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1386423952U), "fileTimeLo",
                  "fileTimeLo", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 185);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 185);
  sf_mex_assign(&c1_rhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs185, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs185), "rhs", "rhs",
                  185);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs185), "lhs", "lhs",
                  185);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"),
                  "context", "context", 186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgeru"), "name", "name",
                  186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"),
                  "resolved", "resolved", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980690U), "fileTimeLo",
                  "fileTimeLo", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 186);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 186);
  sf_mex_assign(&c1_rhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs186, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs186), "rhs", "rhs",
                  186);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs186), "lhs", "lhs",
                  186);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 187);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 187);
  sf_mex_assign(&c1_rhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs187, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs187), "rhs", "rhs",
                  187);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs187), "lhs", "lhs",
                  187);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m"), "context",
                  "context", 188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xgeru"),
                  "name", "name", 188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgeru.p"),
                  "resolved", "resolved", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 188);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 188);
  sf_mex_assign(&c1_rhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs188, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs188), "rhs", "rhs",
                  188);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs188), "lhs", "lhs",
                  188);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgeru.p"),
                  "context", "context", 189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xger"),
                  "name", "name", 189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "resolved", "resolved", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 189);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 189);
  sf_mex_assign(&c1_rhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs189, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs189), "rhs", "rhs",
                  189);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs189), "lhs", "lhs",
                  189);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "context", "context", 190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 190);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 190);
  sf_mex_assign(&c1_rhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs190, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs190), "rhs", "rhs",
                  190);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs190), "lhs", "lhs",
                  190);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 191);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 191);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 191);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 191);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 191);
  sf_mex_assign(&c1_rhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs191, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs191), "rhs", "rhs",
                  191);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs191), "lhs", "lhs",
                  191);
  sf_mex_destroy(&c1_rhs128);
  sf_mex_destroy(&c1_lhs128);
  sf_mex_destroy(&c1_rhs129);
  sf_mex_destroy(&c1_lhs129);
  sf_mex_destroy(&c1_rhs130);
  sf_mex_destroy(&c1_lhs130);
  sf_mex_destroy(&c1_rhs131);
  sf_mex_destroy(&c1_lhs131);
  sf_mex_destroy(&c1_rhs132);
  sf_mex_destroy(&c1_lhs132);
  sf_mex_destroy(&c1_rhs133);
  sf_mex_destroy(&c1_lhs133);
  sf_mex_destroy(&c1_rhs134);
  sf_mex_destroy(&c1_lhs134);
  sf_mex_destroy(&c1_rhs135);
  sf_mex_destroy(&c1_lhs135);
  sf_mex_destroy(&c1_rhs136);
  sf_mex_destroy(&c1_lhs136);
  sf_mex_destroy(&c1_rhs137);
  sf_mex_destroy(&c1_lhs137);
  sf_mex_destroy(&c1_rhs138);
  sf_mex_destroy(&c1_lhs138);
  sf_mex_destroy(&c1_rhs139);
  sf_mex_destroy(&c1_lhs139);
  sf_mex_destroy(&c1_rhs140);
  sf_mex_destroy(&c1_lhs140);
  sf_mex_destroy(&c1_rhs141);
  sf_mex_destroy(&c1_lhs141);
  sf_mex_destroy(&c1_rhs142);
  sf_mex_destroy(&c1_lhs142);
  sf_mex_destroy(&c1_rhs143);
  sf_mex_destroy(&c1_lhs143);
  sf_mex_destroy(&c1_rhs144);
  sf_mex_destroy(&c1_lhs144);
  sf_mex_destroy(&c1_rhs145);
  sf_mex_destroy(&c1_lhs145);
  sf_mex_destroy(&c1_rhs146);
  sf_mex_destroy(&c1_lhs146);
  sf_mex_destroy(&c1_rhs147);
  sf_mex_destroy(&c1_lhs147);
  sf_mex_destroy(&c1_rhs148);
  sf_mex_destroy(&c1_lhs148);
  sf_mex_destroy(&c1_rhs149);
  sf_mex_destroy(&c1_lhs149);
  sf_mex_destroy(&c1_rhs150);
  sf_mex_destroy(&c1_lhs150);
  sf_mex_destroy(&c1_rhs151);
  sf_mex_destroy(&c1_lhs151);
  sf_mex_destroy(&c1_rhs152);
  sf_mex_destroy(&c1_lhs152);
  sf_mex_destroy(&c1_rhs153);
  sf_mex_destroy(&c1_lhs153);
  sf_mex_destroy(&c1_rhs154);
  sf_mex_destroy(&c1_lhs154);
  sf_mex_destroy(&c1_rhs155);
  sf_mex_destroy(&c1_lhs155);
  sf_mex_destroy(&c1_rhs156);
  sf_mex_destroy(&c1_lhs156);
  sf_mex_destroy(&c1_rhs157);
  sf_mex_destroy(&c1_lhs157);
  sf_mex_destroy(&c1_rhs158);
  sf_mex_destroy(&c1_lhs158);
  sf_mex_destroy(&c1_rhs159);
  sf_mex_destroy(&c1_lhs159);
  sf_mex_destroy(&c1_rhs160);
  sf_mex_destroy(&c1_lhs160);
  sf_mex_destroy(&c1_rhs161);
  sf_mex_destroy(&c1_lhs161);
  sf_mex_destroy(&c1_rhs162);
  sf_mex_destroy(&c1_lhs162);
  sf_mex_destroy(&c1_rhs163);
  sf_mex_destroy(&c1_lhs163);
  sf_mex_destroy(&c1_rhs164);
  sf_mex_destroy(&c1_lhs164);
  sf_mex_destroy(&c1_rhs165);
  sf_mex_destroy(&c1_lhs165);
  sf_mex_destroy(&c1_rhs166);
  sf_mex_destroy(&c1_lhs166);
  sf_mex_destroy(&c1_rhs167);
  sf_mex_destroy(&c1_lhs167);
  sf_mex_destroy(&c1_rhs168);
  sf_mex_destroy(&c1_lhs168);
  sf_mex_destroy(&c1_rhs169);
  sf_mex_destroy(&c1_lhs169);
  sf_mex_destroy(&c1_rhs170);
  sf_mex_destroy(&c1_lhs170);
  sf_mex_destroy(&c1_rhs171);
  sf_mex_destroy(&c1_lhs171);
  sf_mex_destroy(&c1_rhs172);
  sf_mex_destroy(&c1_lhs172);
  sf_mex_destroy(&c1_rhs173);
  sf_mex_destroy(&c1_lhs173);
  sf_mex_destroy(&c1_rhs174);
  sf_mex_destroy(&c1_lhs174);
  sf_mex_destroy(&c1_rhs175);
  sf_mex_destroy(&c1_lhs175);
  sf_mex_destroy(&c1_rhs176);
  sf_mex_destroy(&c1_lhs176);
  sf_mex_destroy(&c1_rhs177);
  sf_mex_destroy(&c1_lhs177);
  sf_mex_destroy(&c1_rhs178);
  sf_mex_destroy(&c1_lhs178);
  sf_mex_destroy(&c1_rhs179);
  sf_mex_destroy(&c1_lhs179);
  sf_mex_destroy(&c1_rhs180);
  sf_mex_destroy(&c1_lhs180);
  sf_mex_destroy(&c1_rhs181);
  sf_mex_destroy(&c1_lhs181);
  sf_mex_destroy(&c1_rhs182);
  sf_mex_destroy(&c1_lhs182);
  sf_mex_destroy(&c1_rhs183);
  sf_mex_destroy(&c1_lhs183);
  sf_mex_destroy(&c1_rhs184);
  sf_mex_destroy(&c1_lhs184);
  sf_mex_destroy(&c1_rhs185);
  sf_mex_destroy(&c1_lhs185);
  sf_mex_destroy(&c1_rhs186);
  sf_mex_destroy(&c1_lhs186);
  sf_mex_destroy(&c1_rhs187);
  sf_mex_destroy(&c1_lhs187);
  sf_mex_destroy(&c1_rhs188);
  sf_mex_destroy(&c1_lhs188);
  sf_mex_destroy(&c1_rhs189);
  sf_mex_destroy(&c1_lhs189);
  sf_mex_destroy(&c1_rhs190);
  sf_mex_destroy(&c1_lhs190);
  sf_mex_destroy(&c1_rhs191);
  sf_mex_destroy(&c1_lhs191);
}

static void c1_d_info_helper(const mxArray **c1_info)
{
  const mxArray *c1_rhs192 = NULL;
  const mxArray *c1_lhs192 = NULL;
  const mxArray *c1_rhs193 = NULL;
  const mxArray *c1_lhs193 = NULL;
  const mxArray *c1_rhs194 = NULL;
  const mxArray *c1_lhs194 = NULL;
  const mxArray *c1_rhs195 = NULL;
  const mxArray *c1_lhs195 = NULL;
  const mxArray *c1_rhs196 = NULL;
  const mxArray *c1_lhs196 = NULL;
  const mxArray *c1_rhs197 = NULL;
  const mxArray *c1_lhs197 = NULL;
  const mxArray *c1_rhs198 = NULL;
  const mxArray *c1_lhs198 = NULL;
  const mxArray *c1_rhs199 = NULL;
  const mxArray *c1_lhs199 = NULL;
  const mxArray *c1_rhs200 = NULL;
  const mxArray *c1_lhs200 = NULL;
  const mxArray *c1_rhs201 = NULL;
  const mxArray *c1_lhs201 = NULL;
  const mxArray *c1_rhs202 = NULL;
  const mxArray *c1_lhs202 = NULL;
  const mxArray *c1_rhs203 = NULL;
  const mxArray *c1_lhs203 = NULL;
  const mxArray *c1_rhs204 = NULL;
  const mxArray *c1_lhs204 = NULL;
  const mxArray *c1_rhs205 = NULL;
  const mxArray *c1_lhs205 = NULL;
  const mxArray *c1_rhs206 = NULL;
  const mxArray *c1_lhs206 = NULL;
  const mxArray *c1_rhs207 = NULL;
  const mxArray *c1_lhs207 = NULL;
  const mxArray *c1_rhs208 = NULL;
  const mxArray *c1_lhs208 = NULL;
  const mxArray *c1_rhs209 = NULL;
  const mxArray *c1_lhs209 = NULL;
  const mxArray *c1_rhs210 = NULL;
  const mxArray *c1_lhs210 = NULL;
  const mxArray *c1_rhs211 = NULL;
  const mxArray *c1_lhs211 = NULL;
  const mxArray *c1_rhs212 = NULL;
  const mxArray *c1_lhs212 = NULL;
  const mxArray *c1_rhs213 = NULL;
  const mxArray *c1_lhs213 = NULL;
  const mxArray *c1_rhs214 = NULL;
  const mxArray *c1_lhs214 = NULL;
  const mxArray *c1_rhs215 = NULL;
  const mxArray *c1_lhs215 = NULL;
  const mxArray *c1_rhs216 = NULL;
  const mxArray *c1_lhs216 = NULL;
  const mxArray *c1_rhs217 = NULL;
  const mxArray *c1_lhs217 = NULL;
  const mxArray *c1_rhs218 = NULL;
  const mxArray *c1_lhs218 = NULL;
  const mxArray *c1_rhs219 = NULL;
  const mxArray *c1_lhs219 = NULL;
  const mxArray *c1_rhs220 = NULL;
  const mxArray *c1_lhs220 = NULL;
  const mxArray *c1_rhs221 = NULL;
  const mxArray *c1_lhs221 = NULL;
  const mxArray *c1_rhs222 = NULL;
  const mxArray *c1_lhs222 = NULL;
  const mxArray *c1_rhs223 = NULL;
  const mxArray *c1_lhs223 = NULL;
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.int"),
                  "name", "name", 192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/int.p"),
                  "resolved", "resolved", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 192);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 192);
  sf_mex_assign(&c1_rhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs192, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs192), "rhs", "rhs",
                  192);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs192), "lhs", "lhs",
                  192);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("intmax"), "name", "name", 193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 193);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 193);
  sf_mex_assign(&c1_rhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs193, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs193), "rhs", "rhs",
                  193);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs193), "lhs", "lhs",
                  193);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p!below_threshold"),
                  "context", "context", 194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("min"), "name", "name", 194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m"), "resolved",
                  "resolved", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1311255318U), "fileTimeLo",
                  "fileTimeLo", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 194);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 194);
  sf_mex_assign(&c1_rhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs194, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs194), "rhs", "rhs",
                  194);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs194), "lhs", "lhs",
                  194);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 195);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 195);
  sf_mex_assign(&c1_rhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs195, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs195), "rhs", "rhs",
                  195);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs195), "lhs", "lhs",
                  195);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum"),
                  "context", "context", 196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 196);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 196);
  sf_mex_assign(&c1_rhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs196, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs196), "rhs", "rhs",
                  196);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs196), "lhs", "lhs",
                  196);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 197);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 197);
  sf_mex_assign(&c1_rhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs197, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs197), "rhs", "rhs",
                  197);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs197), "lhs", "lhs",
                  197);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum"),
                  "context", "context", 198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 198);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 198);
  sf_mex_assign(&c1_rhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs198, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs198), "rhs", "rhs",
                  198);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs198), "lhs", "lhs",
                  198);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xger.p"),
                  "context", "context", 199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xger"),
                  "name", "name", 199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xger.p"),
                  "resolved", "resolved", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 199);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 199);
  sf_mex_assign(&c1_rhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs199, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs199), "rhs", "rhs",
                  199);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs199), "lhs", "lhs",
                  199);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xger.p"),
                  "context", "context", 200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xgerx"),
                  "name", "name", 200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "resolved", "resolved", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 200);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 200);
  sf_mex_assign(&c1_rhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs200, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs200), "rhs", "rhs",
                  200);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs200), "lhs", "lhs",
                  200);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713852U), "fileTimeLo",
                  "fileTimeLo", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 201);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 201);
  sf_mex_assign(&c1_rhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs201, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs201), "rhs", "rhs",
                  201);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs201), "lhs", "lhs",
                  201);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 202);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 202);
  sf_mex_assign(&c1_rhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs202, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs202), "rhs", "rhs",
                  202);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs202), "lhs", "lhs",
                  202);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 203);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 203);
  sf_mex_assign(&c1_rhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs203, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs203), "rhs", "rhs",
                  203);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs203), "lhs", "lhs",
                  203);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 204);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 204);
  sf_mex_assign(&c1_rhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs204, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs204), "rhs", "rhs",
                  204);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs204), "lhs", "lhs",
                  204);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgerx.p"),
                  "context", "context", 205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 205);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 205);
  sf_mex_assign(&c1_rhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs205, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs205), "rhs", "rhs",
                  205);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs205), "lhs", "lhs",
                  205);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_ipiv2perm"), "name",
                  "name", 206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_ipiv2perm.m"), "resolved",
                  "resolved", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818782U), "fileTimeLo",
                  "fileTimeLo", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 206);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 206);
  sf_mex_assign(&c1_rhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs206, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs206), "rhs", "rhs",
                  206);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs206), "lhs", "lhs",
                  206);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_ipiv2perm.m"), "context",
                  "context", 207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("colon"), "name", "name", 207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m"), "resolved",
                  "resolved", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1378295988U), "fileTimeLo",
                  "fileTimeLo", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 207);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 207);
  sf_mex_assign(&c1_rhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs207, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs207), "rhs", "rhs",
                  207);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs207), "lhs", "lhs",
                  207);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_ipiv2perm.m"), "context",
                  "context", 208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 208);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 208);
  sf_mex_assign(&c1_rhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs208, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs208), "rhs", "rhs",
                  208);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs208), "lhs", "lhs",
                  208);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_ipiv2perm.m"), "context",
                  "context", 209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.indexIntRelop"),
                  "name", "name", 209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"),
                  "resolved", "resolved", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326728322U), "fileTimeLo",
                  "fileTimeLo", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 209);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 209);
  sf_mex_assign(&c1_rhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs209, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs209), "rhs", "rhs",
                  209);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs209), "lhs", "lhs",
                  209);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 210);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 210);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 210);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 210);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 210);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 210);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 210);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 210);
  sf_mex_assign(&c1_rhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs210, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs210), "rhs", "rhs",
                  210);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs210), "lhs", "lhs",
                  210);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 211);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 211);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 211);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 211);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 211);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 211);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 211);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 211);
  sf_mex_assign(&c1_rhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs211, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs211), "rhs", "rhs",
                  211);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs211), "lhs", "lhs",
                  211);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 212);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 212);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 212);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 212);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 212);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 212);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 212);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 212);
  sf_mex_assign(&c1_rhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs212, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs212), "rhs", "rhs",
                  212);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs212), "lhs", "lhs",
                  212);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!invNxN"), "context",
                  "context", 213);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xtrsm"), "name", "name",
                  213);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 213);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"),
                  "resolved", "resolved", 213);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980692U), "fileTimeLo",
                  "fileTimeLo", 213);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 213);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 213);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 213);
  sf_mex_assign(&c1_rhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs213, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs213), "rhs", "rhs",
                  213);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs213), "lhs", "lhs",
                  213);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 214);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 214);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 214);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 214);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 214);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 214);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 214);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 214);
  sf_mex_assign(&c1_rhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs214, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs214), "rhs", "rhs",
                  214);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs214), "lhs", "lhs",
                  214);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"), "context",
                  "context", 215);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xtrsm"),
                  "name", "name", 215);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 215);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "resolved", "resolved", 215);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 215);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 215);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 215);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 215);
  sf_mex_assign(&c1_rhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs215, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs215), "rhs", "rhs",
                  215);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs215), "lhs", "lhs",
                  215);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 216);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 216);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 216);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 216);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 216);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 216);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 216);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 216);
  sf_mex_assign(&c1_rhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs216, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs216), "rhs", "rhs",
                  216);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs216), "lhs", "lhs",
                  216);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p!below_threshold"),
                  "context", "context", 217);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 217);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 217);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 217);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 217);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 217);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 217);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 217);
  sf_mex_assign(&c1_rhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs217, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs217), "rhs", "rhs",
                  217);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs217), "lhs", "lhs",
                  217);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 218);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 218);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 218);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 218);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 218);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 218);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 218);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 218);
  sf_mex_assign(&c1_rhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs218, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs218), "rhs", "rhs",
                  218);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs218), "lhs", "lhs",
                  218);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xtrsm.p"),
                  "context", "context", 219);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xtrsm"),
                  "name", "name", 219);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 219);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "resolved", "resolved", 219);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 219);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 219);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 219);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 219);
  sf_mex_assign(&c1_rhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs219, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs219), "rhs", "rhs",
                  219);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs219), "lhs", "lhs",
                  219);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 220);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 220);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 220);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 220);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 220);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 220);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 220);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 220);
  sf_mex_assign(&c1_rhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs220, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs220), "rhs", "rhs",
                  220);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs220), "lhs", "lhs",
                  220);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 221);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 221);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 221);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 221);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 221);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 221);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 221);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 221);
  sf_mex_assign(&c1_rhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs221, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs221), "rhs", "rhs",
                  221);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs221), "lhs", "lhs",
                  221);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 222);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 222);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 222);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 222);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 222);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 222);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 222);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 222);
  sf_mex_assign(&c1_rhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs222, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs222), "rhs", "rhs",
                  222);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs222), "lhs", "lhs",
                  222);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xtrsm.p"),
                  "context", "context", 223);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 223);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 223);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 223);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713880U), "fileTimeLo",
                  "fileTimeLo", 223);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 223);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 223);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 223);
  sf_mex_assign(&c1_rhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs223, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs223), "rhs", "rhs",
                  223);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs223), "lhs", "lhs",
                  223);
  sf_mex_destroy(&c1_rhs192);
  sf_mex_destroy(&c1_lhs192);
  sf_mex_destroy(&c1_rhs193);
  sf_mex_destroy(&c1_lhs193);
  sf_mex_destroy(&c1_rhs194);
  sf_mex_destroy(&c1_lhs194);
  sf_mex_destroy(&c1_rhs195);
  sf_mex_destroy(&c1_lhs195);
  sf_mex_destroy(&c1_rhs196);
  sf_mex_destroy(&c1_lhs196);
  sf_mex_destroy(&c1_rhs197);
  sf_mex_destroy(&c1_lhs197);
  sf_mex_destroy(&c1_rhs198);
  sf_mex_destroy(&c1_lhs198);
  sf_mex_destroy(&c1_rhs199);
  sf_mex_destroy(&c1_lhs199);
  sf_mex_destroy(&c1_rhs200);
  sf_mex_destroy(&c1_lhs200);
  sf_mex_destroy(&c1_rhs201);
  sf_mex_destroy(&c1_lhs201);
  sf_mex_destroy(&c1_rhs202);
  sf_mex_destroy(&c1_lhs202);
  sf_mex_destroy(&c1_rhs203);
  sf_mex_destroy(&c1_lhs203);
  sf_mex_destroy(&c1_rhs204);
  sf_mex_destroy(&c1_lhs204);
  sf_mex_destroy(&c1_rhs205);
  sf_mex_destroy(&c1_lhs205);
  sf_mex_destroy(&c1_rhs206);
  sf_mex_destroy(&c1_lhs206);
  sf_mex_destroy(&c1_rhs207);
  sf_mex_destroy(&c1_lhs207);
  sf_mex_destroy(&c1_rhs208);
  sf_mex_destroy(&c1_lhs208);
  sf_mex_destroy(&c1_rhs209);
  sf_mex_destroy(&c1_lhs209);
  sf_mex_destroy(&c1_rhs210);
  sf_mex_destroy(&c1_lhs210);
  sf_mex_destroy(&c1_rhs211);
  sf_mex_destroy(&c1_lhs211);
  sf_mex_destroy(&c1_rhs212);
  sf_mex_destroy(&c1_lhs212);
  sf_mex_destroy(&c1_rhs213);
  sf_mex_destroy(&c1_lhs213);
  sf_mex_destroy(&c1_rhs214);
  sf_mex_destroy(&c1_lhs214);
  sf_mex_destroy(&c1_rhs215);
  sf_mex_destroy(&c1_lhs215);
  sf_mex_destroy(&c1_rhs216);
  sf_mex_destroy(&c1_lhs216);
  sf_mex_destroy(&c1_rhs217);
  sf_mex_destroy(&c1_lhs217);
  sf_mex_destroy(&c1_rhs218);
  sf_mex_destroy(&c1_lhs218);
  sf_mex_destroy(&c1_rhs219);
  sf_mex_destroy(&c1_lhs219);
  sf_mex_destroy(&c1_rhs220);
  sf_mex_destroy(&c1_lhs220);
  sf_mex_destroy(&c1_rhs221);
  sf_mex_destroy(&c1_lhs221);
  sf_mex_destroy(&c1_rhs222);
  sf_mex_destroy(&c1_lhs222);
  sf_mex_destroy(&c1_rhs223);
  sf_mex_destroy(&c1_lhs223);
}

static void c1_eml_switch_helper(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_rdivide(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_x[3], real_T c1_y, real_T c1_z[3])
{
  real_T c1_b_y;
  real_T c1_c_y;
  int32_T c1_i252;
  (void)chartInstance;
  c1_b_y = c1_y;
  c1_c_y = c1_b_y;
  for (c1_i252 = 0; c1_i252 < 3; c1_i252++) {
    c1_z[c1_i252] = c1_x[c1_i252] / c1_c_y;
  }
}

static void c1_inv(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                   real_T c1_x[9], real_T c1_y[9])
{
  int32_T c1_i253;
  real_T c1_b_x[9];
  int32_T c1_i254;
  real_T c1_c_x[9];
  real_T c1_n1x;
  int32_T c1_i255;
  real_T c1_b_y[9];
  real_T c1_n1xinv;
  real_T c1_rc;
  real_T c1_d_x;
  boolean_T c1_b;
  real_T c1_e_x;
  int32_T c1_i256;
  static char_T c1_cv1[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T c1_u[8];
  const mxArray *c1_c_y = NULL;
  real_T c1_b_u;
  const mxArray *c1_d_y = NULL;
  real_T c1_c_u;
  const mxArray *c1_e_y = NULL;
  real_T c1_d_u;
  const mxArray *c1_f_y = NULL;
  char_T c1_str[14];
  int32_T c1_i257;
  char_T c1_b_str[14];
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  for (c1_i253 = 0; c1_i253 < 9; c1_i253++) {
    c1_b_x[c1_i253] = c1_x[c1_i253];
  }

  c1_inv3x3(chartInstance, c1_b_x, c1_y);
  for (c1_i254 = 0; c1_i254 < 9; c1_i254++) {
    c1_c_x[c1_i254] = c1_x[c1_i254];
  }

  c1_n1x = c1_norm(chartInstance, c1_c_x);
  for (c1_i255 = 0; c1_i255 < 9; c1_i255++) {
    c1_b_y[c1_i255] = c1_y[c1_i255];
  }

  c1_n1xinv = c1_norm(chartInstance, c1_b_y);
  c1_rc = 1.0 / (c1_n1x * c1_n1xinv);
  guard1 = false;
  guard2 = false;
  if (c1_n1x == 0.0) {
    guard2 = true;
  } else if (c1_n1xinv == 0.0) {
    guard2 = true;
  } else if (c1_rc == 0.0) {
    guard1 = true;
  } else {
    c1_d_x = c1_rc;
    c1_b = muDoubleScalarIsNaN(c1_d_x);
    guard3 = false;
    if (c1_b) {
      guard3 = true;
    } else {
      c1_eml_eps(chartInstance);
      if (c1_rc < 2.2204460492503131E-16) {
        guard3 = true;
      }
    }

    if (guard3 == true) {
      c1_e_x = c1_rc;
      for (c1_i256 = 0; c1_i256 < 8; c1_i256++) {
        c1_u[c1_i256] = c1_cv1[c1_i256];
      }

      c1_c_y = NULL;
      sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 8),
                    false);
      c1_b_u = 14.0;
      c1_d_y = NULL;
      sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0),
                    false);
      c1_c_u = 6.0;
      c1_e_y = NULL;
      sf_mex_assign(&c1_e_y, sf_mex_create("y", &c1_c_u, 0, 0U, 0U, 0U, 0),
                    false);
      c1_d_u = c1_e_x;
      c1_f_y = NULL;
      sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_d_u, 0, 0U, 0U, 0U, 0),
                    false);
      c1_n_emlrt_marshallIn(chartInstance, sf_mex_call_debug
                            (sfGlobalDebugInstanceStruct, "sprintf", 1U, 2U, 14,
        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "sprintf", 1U, 3U, 14,
                          c1_c_y, 14, c1_d_y, 14, c1_e_y), 14, c1_f_y),
                            "sprintf", c1_str);
      for (c1_i257 = 0; c1_i257 < 14; c1_i257++) {
        c1_b_str[c1_i257] = c1_str[c1_i257];
      }

      c1_b_eml_warning(chartInstance, c1_b_str);
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c1_eml_warning(chartInstance);
  }
}

static void c1_inv3x3(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c1_x[9], real_T c1_y[9])
{
  int32_T c1_p1;
  int32_T c1_p2;
  int32_T c1_p3;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_absx11;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_absx21;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_absx31;
  real_T c1_t1;
  real_T c1_h_x;
  real_T c1_b_y;
  real_T c1_i_x;
  real_T c1_c_y;
  real_T c1_z;
  real_T c1_j_x;
  real_T c1_d_y;
  real_T c1_k_x;
  real_T c1_e_y;
  real_T c1_b_z;
  real_T c1_l_x;
  real_T c1_m_x;
  real_T c1_f_y;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_g_y;
  int32_T c1_itmp;
  real_T c1_p_x;
  real_T c1_h_y;
  real_T c1_q_x;
  real_T c1_i_y;
  real_T c1_c_z;
  real_T c1_r_x;
  real_T c1_j_y;
  real_T c1_s_x;
  real_T c1_k_y;
  real_T c1_t3;
  real_T c1_t_x;
  real_T c1_l_y;
  real_T c1_u_x;
  real_T c1_m_y;
  real_T c1_t2;
  int32_T c1_b_a;
  int32_T c1_c_a;
  int32_T c1_c;
  real_T c1_v_x;
  real_T c1_n_y;
  real_T c1_w_x;
  real_T c1_o_y;
  real_T c1_d_z;
  int32_T c1_d_a;
  int32_T c1_e_a;
  int32_T c1_b_c;
  int32_T c1_f_a;
  int32_T c1_g_a;
  int32_T c1_c_c;
  real_T c1_x_x;
  real_T c1_p_y;
  real_T c1_y_x;
  real_T c1_q_y;
  real_T c1_ab_x;
  real_T c1_r_y;
  real_T c1_bb_x;
  real_T c1_s_y;
  int32_T c1_h_a;
  int32_T c1_i_a;
  int32_T c1_d_c;
  real_T c1_cb_x;
  real_T c1_t_y;
  real_T c1_db_x;
  real_T c1_u_y;
  real_T c1_e_z;
  int32_T c1_j_a;
  int32_T c1_k_a;
  int32_T c1_e_c;
  int32_T c1_l_a;
  int32_T c1_m_a;
  int32_T c1_f_c;
  real_T c1_v_y;
  real_T c1_w_y;
  real_T c1_eb_x;
  real_T c1_x_y;
  real_T c1_fb_x;
  real_T c1_y_y;
  int32_T c1_n_a;
  int32_T c1_o_a;
  int32_T c1_g_c;
  real_T c1_gb_x;
  real_T c1_ab_y;
  real_T c1_hb_x;
  real_T c1_bb_y;
  real_T c1_f_z;
  int32_T c1_p_a;
  int32_T c1_q_a;
  int32_T c1_h_c;
  int32_T c1_r_a;
  int32_T c1_s_a;
  int32_T c1_i_c;
  boolean_T guard1 = false;
  (void)chartInstance;
  c1_p1 = 0;
  c1_p2 = 3;
  c1_p3 = 6;
  c1_b_x = c1_x[0];
  c1_c_x = c1_b_x;
  c1_absx11 = muDoubleScalarAbs(c1_c_x);
  c1_d_x = c1_x[1];
  c1_e_x = c1_d_x;
  c1_absx21 = muDoubleScalarAbs(c1_e_x);
  c1_f_x = c1_x[2];
  c1_g_x = c1_f_x;
  c1_absx31 = muDoubleScalarAbs(c1_g_x);
  guard1 = false;
  if (c1_absx21 > c1_absx11) {
    if (c1_absx21 > c1_absx31) {
      c1_p1 = 3;
      c1_p2 = 0;
      c1_t1 = c1_x[0];
      c1_x[0] = c1_x[1];
      c1_x[1] = c1_t1;
      c1_t1 = c1_x[3];
      c1_x[3] = c1_x[4];
      c1_x[4] = c1_t1;
      c1_t1 = c1_x[6];
      c1_x[6] = c1_x[7];
      c1_x[7] = c1_t1;
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1 == true) {
    if (c1_absx31 > c1_absx11) {
      c1_p1 = 6;
      c1_p3 = 0;
      c1_t1 = c1_x[0];
      c1_x[0] = c1_x[2];
      c1_x[2] = c1_t1;
      c1_t1 = c1_x[3];
      c1_x[3] = c1_x[5];
      c1_x[5] = c1_t1;
      c1_t1 = c1_x[6];
      c1_x[6] = c1_x[8];
      c1_x[8] = c1_t1;
    }
  }

  c1_h_x = c1_x[1];
  c1_b_y = c1_x[0];
  c1_i_x = c1_h_x;
  c1_c_y = c1_b_y;
  c1_z = c1_i_x / c1_c_y;
  c1_x[1] = c1_z;
  c1_j_x = c1_x[2];
  c1_d_y = c1_x[0];
  c1_k_x = c1_j_x;
  c1_e_y = c1_d_y;
  c1_b_z = c1_k_x / c1_e_y;
  c1_x[2] = c1_b_z;
  c1_x[4] -= c1_x[1] * c1_x[3];
  c1_x[5] -= c1_x[2] * c1_x[3];
  c1_x[7] -= c1_x[1] * c1_x[6];
  c1_x[8] -= c1_x[2] * c1_x[6];
  c1_l_x = c1_x[5];
  c1_m_x = c1_l_x;
  c1_f_y = muDoubleScalarAbs(c1_m_x);
  c1_n_x = c1_x[4];
  c1_o_x = c1_n_x;
  c1_g_y = muDoubleScalarAbs(c1_o_x);
  if (c1_f_y > c1_g_y) {
    c1_itmp = c1_p2;
    c1_p2 = c1_p3;
    c1_p3 = c1_itmp;
    c1_t1 = c1_x[1];
    c1_x[1] = c1_x[2];
    c1_x[2] = c1_t1;
    c1_t1 = c1_x[4];
    c1_x[4] = c1_x[5];
    c1_x[5] = c1_t1;
    c1_t1 = c1_x[7];
    c1_x[7] = c1_x[8];
    c1_x[8] = c1_t1;
  }

  c1_p_x = c1_x[5];
  c1_h_y = c1_x[4];
  c1_q_x = c1_p_x;
  c1_i_y = c1_h_y;
  c1_c_z = c1_q_x / c1_i_y;
  c1_x[5] = c1_c_z;
  c1_x[8] -= c1_x[5] * c1_x[7];
  c1_r_x = c1_x[5] * c1_x[1] - c1_x[2];
  c1_j_y = c1_x[8];
  c1_s_x = c1_r_x;
  c1_k_y = c1_j_y;
  c1_t3 = c1_s_x / c1_k_y;
  c1_t_x = -(c1_x[1] + c1_x[7] * c1_t3);
  c1_l_y = c1_x[4];
  c1_u_x = c1_t_x;
  c1_m_y = c1_l_y;
  c1_t2 = c1_u_x / c1_m_y;
  c1_b_a = c1_p1;
  c1_c_a = c1_b_a + 1;
  c1_c = c1_c_a;
  c1_v_x = (1.0 - c1_x[3] * c1_t2) - c1_x[6] * c1_t3;
  c1_n_y = c1_x[0];
  c1_w_x = c1_v_x;
  c1_o_y = c1_n_y;
  c1_d_z = c1_w_x / c1_o_y;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_c), 1, 9, 1, 0) - 1] = c1_d_z;
  c1_d_a = c1_p1;
  c1_e_a = c1_d_a + 2;
  c1_b_c = c1_e_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_b_c), 1, 9, 1, 0) - 1] = c1_t2;
  c1_f_a = c1_p1;
  c1_g_a = c1_f_a + 3;
  c1_c_c = c1_g_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_c_c), 1, 9, 1, 0) - 1] = c1_t3;
  c1_x_x = -c1_x[5];
  c1_p_y = c1_x[8];
  c1_y_x = c1_x_x;
  c1_q_y = c1_p_y;
  c1_t3 = c1_y_x / c1_q_y;
  c1_ab_x = 1.0 - c1_x[7] * c1_t3;
  c1_r_y = c1_x[4];
  c1_bb_x = c1_ab_x;
  c1_s_y = c1_r_y;
  c1_t2 = c1_bb_x / c1_s_y;
  c1_h_a = c1_p2;
  c1_i_a = c1_h_a + 1;
  c1_d_c = c1_i_a;
  c1_cb_x = -(c1_x[3] * c1_t2 + c1_x[6] * c1_t3);
  c1_t_y = c1_x[0];
  c1_db_x = c1_cb_x;
  c1_u_y = c1_t_y;
  c1_e_z = c1_db_x / c1_u_y;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_d_c), 1, 9, 1, 0) - 1] = c1_e_z;
  c1_j_a = c1_p2;
  c1_k_a = c1_j_a + 2;
  c1_e_c = c1_k_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_e_c), 1, 9, 1, 0) - 1] = c1_t2;
  c1_l_a = c1_p2;
  c1_m_a = c1_l_a + 3;
  c1_f_c = c1_m_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_f_c), 1, 9, 1, 0) - 1] = c1_t3;
  c1_v_y = c1_x[8];
  c1_w_y = c1_v_y;
  c1_t3 = 1.0 / c1_w_y;
  c1_eb_x = -c1_x[7] * c1_t3;
  c1_x_y = c1_x[4];
  c1_fb_x = c1_eb_x;
  c1_y_y = c1_x_y;
  c1_t2 = c1_fb_x / c1_y_y;
  c1_n_a = c1_p3;
  c1_o_a = c1_n_a + 1;
  c1_g_c = c1_o_a;
  c1_gb_x = -(c1_x[3] * c1_t2 + c1_x[6] * c1_t3);
  c1_ab_y = c1_x[0];
  c1_hb_x = c1_gb_x;
  c1_bb_y = c1_ab_y;
  c1_f_z = c1_hb_x / c1_bb_y;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_g_c), 1, 9, 1, 0) - 1] = c1_f_z;
  c1_p_a = c1_p3;
  c1_q_a = c1_p_a + 2;
  c1_h_c = c1_q_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_h_c), 1, 9, 1, 0) - 1] = c1_t2;
  c1_r_a = c1_p3;
  c1_s_a = c1_r_a + 3;
  c1_i_c = c1_s_a;
  c1_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_i_c), 1, 9, 1, 0) - 1] = c1_t3;
}

static real_T c1_norm(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c1_x[9])
{
  real_T c1_y;
  int32_T c1_j;
  real_T c1_b_j;
  real_T c1_s;
  int32_T c1_i;
  real_T c1_b_i;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_b_y;
  real_T c1_d_x;
  boolean_T c1_b;
  boolean_T exitg1;
  (void)chartInstance;
  c1_y = 0.0;
  c1_j = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c1_j < 3)) {
    c1_b_j = 1.0 + (real_T)c1_j;
    c1_s = 0.0;
    for (c1_i = 0; c1_i < 3; c1_i++) {
      c1_b_i = 1.0 + (real_T)c1_i;
      c1_b_x = c1_x[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_i), 1, 3, 1, 0) + 3 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 3, 2, 0) - 1)) - 1];
      c1_c_x = c1_b_x;
      c1_b_y = muDoubleScalarAbs(c1_c_x);
      c1_s += c1_b_y;
    }

    c1_d_x = c1_s;
    c1_b = muDoubleScalarIsNaN(c1_d_x);
    if (c1_b) {
      c1_y = rtNaN;
      exitg1 = true;
    } else {
      if (c1_s > c1_y) {
        c1_y = c1_s;
      }

      c1_j++;
    }
  }

  return c1_y;
}

static void c1_eml_warning(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  int32_T c1_i258;
  static char_T c1_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c1_u[27];
  const mxArray *c1_y = NULL;
  (void)chartInstance;
  for (c1_i258 = 0; c1_i258 < 27; c1_i258++) {
    c1_u[c1_i258] = c1_varargin_1[c1_i258];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 27), false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    1U, 14, c1_y));
}

static void c1_eml_eps(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_b_eml_warning(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, char_T c1_varargin_2[14])
{
  int32_T c1_i259;
  static char_T c1_varargin_1[33] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 'i', 'l', 'l', 'C', 'o', 'n', 'd', 'i', 't', 'i',
    'o', 'n', 'e', 'd', 'M', 'a', 't', 'r', 'i', 'x' };

  char_T c1_u[33];
  const mxArray *c1_y = NULL;
  int32_T c1_i260;
  char_T c1_b_u[14];
  const mxArray *c1_b_y = NULL;
  (void)chartInstance;
  for (c1_i259 = 0; c1_i259 < 33; c1_i259++) {
    c1_u[c1_i259] = c1_varargin_1[c1_i259];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 33), false);
  for (c1_i260 = 0; c1_i260 < 14; c1_i260++) {
    c1_b_u[c1_i260] = c1_varargin_2[c1_i260];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 14),
                false);
  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "warning", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c1_y, 14, c1_b_y));
}

static void c1_eml_scalar_eg(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_threshold(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_b_eml_scalar_eg(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static real_T c1_atan2(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_y, real_T c1_x)
{
  real_T c1_b_y;
  real_T c1_b_x;
  c1_c_eml_scalar_eg(chartInstance);
  c1_b_y = c1_y;
  c1_b_x = c1_x;
  return muDoubleScalarAtan2(c1_b_y, c1_b_x);
}

static void c1_c_eml_scalar_eg(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_b_rdivide(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_x[9], real_T c1_y, real_T c1_z[9])
{
  real_T c1_b_y;
  real_T c1_c_y;
  int32_T c1_i261;
  (void)chartInstance;
  c1_b_y = c1_y;
  c1_c_y = c1_b_y;
  for (c1_i261 = 0; c1_i261 < 9; c1_i261++) {
    c1_z[c1_i261] = c1_x[c1_i261] / c1_c_y;
  }
}

static void c1_eye(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                   real_T c1_I[9])
{
  int32_T c1_i262;
  int32_T c1_k;
  int32_T c1_b_k;
  for (c1_i262 = 0; c1_i262 < 9; c1_i262++) {
    c1_I[c1_i262] = 0.0;
  }

  c1_eml_switch_helper(chartInstance);
  for (c1_k = 1; c1_k < 4; c1_k++) {
    c1_b_k = c1_k;
    c1_I[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 3, 1, 0) + 3 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_k), 1, 3, 2, 0) - 1)) -
      1] = 1.0;
  }
}

static void c1_mldivide(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_A[9], real_T c1_B[3], real_T c1_Y[3])
{
  int32_T c1_i263;
  real_T c1_b_A[9];
  int32_T c1_r1;
  int32_T c1_r2;
  int32_T c1_r3;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_y;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_b_y;
  real_T c1_maxval;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_c_y;
  real_T c1_k_x;
  real_T c1_l_x;
  real_T c1_d_y;
  real_T c1_a21;
  real_T c1_m_x;
  real_T c1_n_x;
  real_T c1_o_x;
  real_T c1_p_x;
  real_T c1_e_y;
  real_T c1_q_x;
  real_T c1_r_x;
  real_T c1_f_y;
  real_T c1_d;
  real_T c1_s_x;
  real_T c1_t_x;
  real_T c1_u_x;
  real_T c1_v_x;
  real_T c1_g_y;
  real_T c1_w_x;
  real_T c1_x_x;
  real_T c1_h_y;
  real_T c1_b_d;
  real_T c1_y_x;
  real_T c1_ab_x;
  real_T c1_bb_x;
  real_T c1_cb_x;
  real_T c1_i_y;
  real_T c1_db_x;
  real_T c1_eb_x;
  real_T c1_j_y;
  real_T c1_c_d;
  int32_T c1_rtemp;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  for (c1_i263 = 0; c1_i263 < 9; c1_i263++) {
    c1_b_A[c1_i263] = c1_A[c1_i263];
  }

  c1_r1 = 1;
  c1_r2 = 2;
  c1_r3 = 3;
  c1_x = c1_b_A[0];
  c1_b_x = c1_x;
  c1_c_x = c1_b_x;
  c1_d_x = c1_c_x;
  c1_y = muDoubleScalarAbs(c1_d_x);
  c1_e_x = 0.0;
  c1_f_x = c1_e_x;
  c1_b_y = muDoubleScalarAbs(c1_f_x);
  c1_maxval = c1_y + c1_b_y;
  c1_g_x = c1_b_A[1];
  c1_h_x = c1_g_x;
  c1_i_x = c1_h_x;
  c1_j_x = c1_i_x;
  c1_c_y = muDoubleScalarAbs(c1_j_x);
  c1_k_x = 0.0;
  c1_l_x = c1_k_x;
  c1_d_y = muDoubleScalarAbs(c1_l_x);
  c1_a21 = c1_c_y + c1_d_y;
  if (c1_a21 > c1_maxval) {
    c1_maxval = c1_a21;
    c1_r1 = 2;
    c1_r2 = 1;
  }

  c1_m_x = c1_b_A[2];
  c1_n_x = c1_m_x;
  c1_o_x = c1_n_x;
  c1_p_x = c1_o_x;
  c1_e_y = muDoubleScalarAbs(c1_p_x);
  c1_q_x = 0.0;
  c1_r_x = c1_q_x;
  c1_f_y = muDoubleScalarAbs(c1_r_x);
  c1_d = c1_e_y + c1_f_y;
  if (c1_d > c1_maxval) {
    c1_r1 = 3;
    c1_r2 = 2;
    c1_r3 = 1;
  }

  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r2), 1, 3, 1, 0) - 1] = c1_c_rdivide(chartInstance,
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) - 1], c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) - 1]);
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) - 1] = c1_c_rdivide(chartInstance,
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r3), 1, 3, 1, 0) - 1], c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) - 1]);
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r2), 1, 3, 1, 0) + 2] = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 2] -
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) - 1] * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 2];
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 2] = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 2] -
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r3), 1, 3, 1, 0) - 1] * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 2];
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r2), 1, 3, 1, 0) + 5] = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 5] -
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) - 1] * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 5];
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 5] = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5] -
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r3), 1, 3, 1, 0) - 1] * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 5];
  c1_s_x = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r3), 1, 3, 1, 0) + 2];
  c1_t_x = c1_s_x;
  c1_u_x = c1_t_x;
  c1_v_x = c1_u_x;
  c1_g_y = muDoubleScalarAbs(c1_v_x);
  c1_w_x = 0.0;
  c1_x_x = c1_w_x;
  c1_h_y = muDoubleScalarAbs(c1_x_x);
  c1_b_d = c1_g_y + c1_h_y;
  c1_y_x = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) + 2];
  c1_ab_x = c1_y_x;
  c1_bb_x = c1_ab_x;
  c1_cb_x = c1_bb_x;
  c1_i_y = muDoubleScalarAbs(c1_cb_x);
  c1_db_x = 0.0;
  c1_eb_x = c1_db_x;
  c1_j_y = muDoubleScalarAbs(c1_eb_x);
  c1_c_d = c1_i_y + c1_j_y;
  if (c1_b_d > c1_c_d) {
    c1_rtemp = c1_r2;
    c1_r2 = c1_r3;
    c1_r3 = c1_rtemp;
  }

  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 2] = c1_c_rdivide(chartInstance,
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r3), 1, 3, 1, 0) + 2], c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 2]);
  c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)
    c1_r3), 1, 3, 1, 0) + 5] = c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c1_r3), 1, 3, 1, 0) + 5] -
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r3), 1, 3, 1, 0) + 2] * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("",
    (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 5];
  guard1 = false;
  guard2 = false;
  if (c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_r1), 1, 3, 1, 0) - 1] == 0.0) {
    guard2 = true;
  } else if (c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
               "", (real_T)c1_r2), 1, 3, 1, 0) + 2] == 0.0) {
    guard2 = true;
  } else {
    if (c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_r3), 1, 3, 1, 0) + 5] == 0.0) {
      guard1 = true;
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c1_eml_warning(chartInstance);
  }

  c1_eml_scalar_eg(chartInstance);
  c1_Y[0] = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) - 1];
  c1_Y[1] = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) - 1] - c1_Y[0] *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) - 1];
  c1_Y[2] = (c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
               (real_T)c1_r3), 1, 3, 1, 0) - 1] - c1_Y[0] *
             c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
    "", (real_T)c1_r3), 1, 3, 1, 0) - 1]) - c1_Y[1] *
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r3), 1, 3, 1, 0) + 2];
  c1_Y[2] = c1_c_rdivide(chartInstance, c1_Y[2],
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r3), 1, 3, 1, 0) + 5]);
  c1_Y[0] -= c1_Y[2] * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 5];
  c1_Y[1] -= c1_Y[2] * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c1_r2), 1, 3, 1, 0) + 5];
  c1_Y[1] = c1_c_rdivide(chartInstance, c1_Y[1],
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r2), 1, 3, 1, 0) + 2]);
  c1_Y[0] -= c1_Y[1] * c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
    _SFD_INTEGER_CHECK("", (real_T)c1_r1), 1, 3, 1, 0) + 2];
  c1_Y[0] = c1_c_rdivide(chartInstance, c1_Y[0],
    c1_b_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
    (real_T)c1_r1), 1, 3, 1, 0) - 1]);
}

static real_T c1_c_rdivide(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_x, real_T c1_y)
{
  real_T c1_b_x;
  real_T c1_b_y;
  real_T c1_c_x;
  real_T c1_c_y;
  (void)chartInstance;
  c1_b_x = c1_x;
  c1_b_y = c1_y;
  c1_c_x = c1_b_x;
  c1_c_y = c1_b_y;
  return c1_c_x / c1_c_y;
}

static void c1_b_inv(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                     real_T c1_x[36], real_T c1_y[36])
{
  int32_T c1_i264;
  real_T c1_b_x[36];
  int32_T c1_i265;
  int32_T c1_info;
  int32_T c1_ipiv[6];
  int32_T c1_i266;
  int32_T c1_b_ipiv[6];
  int32_T c1_p[6];
  int32_T c1_k;
  int32_T c1_b_k;
  int32_T c1_c;
  int32_T c1_c_k;
  int32_T c1_b_a;
  int32_T c1_c_a;
  boolean_T c1_overflow;
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_d_a;
  int32_T c1_e_a;
  int32_T c1_i267;
  int32_T c1_f_a;
  int32_T c1_g_a;
  boolean_T c1_b_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  int32_T c1_i268;
  real_T c1_c_x[36];
  int32_T c1_i269;
  real_T c1_d_x[36];
  real_T c1_n1x;
  int32_T c1_i270;
  real_T c1_b_y[36];
  real_T c1_n1xinv;
  real_T c1_rc;
  real_T c1_e_x;
  boolean_T c1_b;
  real_T c1_f_x;
  int32_T c1_i271;
  static char_T c1_cv2[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T c1_u[8];
  const mxArray *c1_c_y = NULL;
  real_T c1_b_u;
  const mxArray *c1_d_y = NULL;
  real_T c1_c_u;
  const mxArray *c1_e_y = NULL;
  real_T c1_d_u;
  const mxArray *c1_f_y = NULL;
  char_T c1_str[14];
  int32_T c1_i272;
  char_T c1_b_str[14];
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  for (c1_i264 = 0; c1_i264 < 36; c1_i264++) {
    c1_b_x[c1_i264] = c1_x[c1_i264];
  }

  for (c1_i265 = 0; c1_i265 < 36; c1_i265++) {
    c1_y[c1_i265] = 0.0;
  }

  c1_b_eml_matlab_zgetrf(chartInstance, c1_b_x, c1_ipiv, &c1_info);
  for (c1_i266 = 0; c1_i266 < 6; c1_i266++) {
    c1_b_ipiv[c1_i266] = c1_ipiv[c1_i266];
  }

  c1_eml_ipiv2perm(chartInstance, c1_b_ipiv, c1_p);
  c1_eml_switch_helper(chartInstance);
  for (c1_k = 1; c1_k < 7; c1_k++) {
    c1_b_k = c1_k;
    c1_c = c1_p[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c1_b_k), 1, 6, 1, 0) - 1];
    c1_y[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_k), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
            (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_c), 1, 6, 2, 0) - 1)) - 1]
      = 1.0;
    c1_c_k = c1_b_k;
    c1_b_a = c1_c_k;
    c1_c_a = c1_b_a;
    if (c1_c_a > 6) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_eml_switch_helper(chartInstance);
      c1_overflow = false;
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_j = c1_c_k; c1_j < 7; c1_j++) {
      c1_b_j = c1_j;
      if (c1_y[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
             (real_T)c1_b_j), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK("",
             (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_c), 1, 6, 2, 0) - 1)) -
          1] != 0.0) {
        c1_d_a = c1_b_j;
        c1_e_a = c1_d_a + 1;
        c1_i267 = c1_e_a;
        c1_f_a = c1_i267;
        c1_g_a = c1_f_a;
        if (c1_g_a > 6) {
          c1_b_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_eml_switch_helper(chartInstance);
          c1_b_overflow = false;
        }

        if (c1_b_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
        }

        for (c1_i = c1_i267; c1_i < 7; c1_i++) {
          c1_b_i = c1_i;
          c1_y[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_i), 1, 6, 1, 0) + 6 *
                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                   (real_T)c1_c), 1, 6, 2, 0) - 1)) - 1] = c1_y
            [(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                (real_T)c1_b_i), 1, 6, 1, 0) + 6 * (_SFD_EML_ARRAY_BOUNDS_CHECK(
                "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_c), 1, 6, 2, 0) -
               1)) - 1] - c1_y[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 6, 1, 0) + 6 *
                                (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_c), 1, 6, 2, 0) - 1)) - 1] *
            c1_b_x[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
                      "", (real_T)c1_b_i), 1, 6, 1, 0) + 6 *
                    (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", (real_T)c1_b_j), 1, 6, 2, 0) - 1)) - 1];
        }
      }
    }
  }

  for (c1_i268 = 0; c1_i268 < 36; c1_i268++) {
    c1_c_x[c1_i268] = c1_b_x[c1_i268];
  }

  c1_b_eml_xtrsm(chartInstance, c1_c_x, c1_y);
  for (c1_i269 = 0; c1_i269 < 36; c1_i269++) {
    c1_d_x[c1_i269] = c1_x[c1_i269];
  }

  c1_n1x = c1_b_norm(chartInstance, c1_d_x);
  for (c1_i270 = 0; c1_i270 < 36; c1_i270++) {
    c1_b_y[c1_i270] = c1_y[c1_i270];
  }

  c1_n1xinv = c1_b_norm(chartInstance, c1_b_y);
  c1_rc = 1.0 / (c1_n1x * c1_n1xinv);
  guard1 = false;
  guard2 = false;
  if (c1_n1x == 0.0) {
    guard2 = true;
  } else if (c1_n1xinv == 0.0) {
    guard2 = true;
  } else if (c1_rc == 0.0) {
    guard1 = true;
  } else {
    c1_e_x = c1_rc;
    c1_b = muDoubleScalarIsNaN(c1_e_x);
    guard3 = false;
    if (c1_b) {
      guard3 = true;
    } else {
      c1_eml_eps(chartInstance);
      if (c1_rc < 2.2204460492503131E-16) {
        guard3 = true;
      }
    }

    if (guard3 == true) {
      c1_f_x = c1_rc;
      for (c1_i271 = 0; c1_i271 < 8; c1_i271++) {
        c1_u[c1_i271] = c1_cv2[c1_i271];
      }

      c1_c_y = NULL;
      sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 8),
                    false);
      c1_b_u = 14.0;
      c1_d_y = NULL;
      sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0),
                    false);
      c1_c_u = 6.0;
      c1_e_y = NULL;
      sf_mex_assign(&c1_e_y, sf_mex_create("y", &c1_c_u, 0, 0U, 0U, 0U, 0),
                    false);
      c1_d_u = c1_f_x;
      c1_f_y = NULL;
      sf_mex_assign(&c1_f_y, sf_mex_create("y", &c1_d_u, 0, 0U, 0U, 0U, 0),
                    false);
      c1_n_emlrt_marshallIn(chartInstance, sf_mex_call_debug
                            (sfGlobalDebugInstanceStruct, "sprintf", 1U, 2U, 14,
        sf_mex_call_debug(sfGlobalDebugInstanceStruct, "sprintf", 1U, 3U, 14,
                          c1_c_y, 14, c1_d_y, 14, c1_e_y), 14, c1_f_y),
                            "sprintf", c1_str);
      for (c1_i272 = 0; c1_i272 < 14; c1_i272++) {
        c1_b_str[c1_i272] = c1_str[c1_i272];
      }

      c1_b_eml_warning(chartInstance, c1_b_str);
    }
  }

  if (guard2 == true) {
    guard1 = true;
  }

  if (guard1 == true) {
    c1_eml_warning(chartInstance);
  }
}

static void c1_realmin(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_eml_matlab_zgetrf(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_A[36], real_T c1_b_A[36], int32_T c1_ipiv[6],
  int32_T *c1_info)
{
  int32_T c1_i273;
  for (c1_i273 = 0; c1_i273 < 36; c1_i273++) {
    c1_b_A[c1_i273] = c1_A[c1_i273];
  }

  c1_b_eml_matlab_zgetrf(chartInstance, c1_b_A, c1_ipiv, c1_info);
}

static int32_T c1_eml_ixamax(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, int32_T c1_n, real_T c1_x[36], int32_T c1_ix0)
{
  int32_T c1_idxmax;
  int32_T c1_b_n;
  int32_T c1_b_ix0;
  int32_T c1_c_n;
  int32_T c1_c_ix0;
  int32_T c1_ix;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_y;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_b_y;
  real_T c1_smax;
  int32_T c1_d_n;
  int32_T c1_b;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  int32_T c1_b_a;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_i_x;
  real_T c1_c_y;
  real_T c1_j_x;
  real_T c1_k_x;
  real_T c1_d_y;
  real_T c1_s;
  c1_b_n = c1_n;
  c1_b_ix0 = c1_ix0;
  c1_c_n = c1_b_n;
  c1_c_ix0 = c1_b_ix0;
  if (c1_c_n < 1) {
    c1_idxmax = 0;
  } else {
    c1_idxmax = 1;
    if (c1_c_n > 1) {
      c1_ix = c1_c_ix0;
      c1_b_x = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c1_ix), 1, 36, 1, 0) - 1];
      c1_c_x = c1_b_x;
      c1_d_x = c1_c_x;
      c1_y = muDoubleScalarAbs(c1_d_x);
      c1_e_x = 0.0;
      c1_f_x = c1_e_x;
      c1_b_y = muDoubleScalarAbs(c1_f_x);
      c1_smax = c1_y + c1_b_y;
      c1_d_n = c1_c_n;
      c1_b = c1_d_n;
      c1_b_b = c1_b;
      if (2 > c1_b_b) {
        c1_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_eml_switch_helper(chartInstance);
        c1_overflow = (c1_b_b > 2147483646);
      }

      if (c1_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_overflow);
      }

      for (c1_k = 2; c1_k <= c1_d_n; c1_k++) {
        c1_b_k = c1_k;
        c1_b_a = c1_ix + 1;
        c1_ix = c1_b_a;
        c1_g_x = c1_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_ix), 1, 36, 1, 0) - 1];
        c1_h_x = c1_g_x;
        c1_i_x = c1_h_x;
        c1_c_y = muDoubleScalarAbs(c1_i_x);
        c1_j_x = 0.0;
        c1_k_x = c1_j_x;
        c1_d_y = muDoubleScalarAbs(c1_k_x);
        c1_s = c1_c_y + c1_d_y;
        if (c1_s > c1_smax) {
          c1_idxmax = c1_b_k;
          c1_smax = c1_s;
        }
      }
    }
  }

  return c1_idxmax;
}

static void c1_check_forloop_overflow_error
  (SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance, boolean_T c1_overflow)
{
  int32_T c1_i274;
  static char_T c1_cv3[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c1_u[34];
  const mxArray *c1_y = NULL;
  int32_T c1_i275;
  static char_T c1_cv4[5] = { 'i', 'n', 't', '3', '2' };

  char_T c1_b_u[5];
  const mxArray *c1_b_y = NULL;
  (void)chartInstance;
  if (!c1_overflow) {
  } else {
    for (c1_i274 = 0; c1_i274 < 34; c1_i274++) {
      c1_u[c1_i274] = c1_cv3[c1_i274];
    }

    c1_y = NULL;
    sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  false);
    for (c1_i275 = 0; c1_i275 < 5; c1_i275++) {
      c1_b_u[c1_i275] = c1_cv4[c1_i275];
    }

    c1_b_y = NULL;
    sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 5),
                  false);
    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message",
      1U, 2U, 14, c1_y, 14, c1_b_y));
  }
}

static void c1_b_eml_switch_helper(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_eml_xgeru(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  int32_T c1_m, int32_T c1_n, real_T c1_alpha1, int32_T c1_ix0, int32_T c1_iy0,
  real_T c1_A[36], int32_T c1_ia0, real_T c1_b_A[36])
{
  int32_T c1_i276;
  for (c1_i276 = 0; c1_i276 < 36; c1_i276++) {
    c1_b_A[c1_i276] = c1_A[c1_i276];
  }

  c1_b_eml_xgeru(chartInstance, c1_m, c1_n, c1_alpha1, c1_ix0, c1_iy0, c1_b_A,
                 c1_ia0);
}

static void c1_eml_ipiv2perm(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, int32_T c1_ipiv[6], int32_T c1_perm[6])
{
  int32_T c1_i277;
  int32_T c1_k;
  real_T c1_b_k;
  int32_T c1_ipk;
  int32_T c1_b_a;
  real_T c1_b;
  int32_T c1_c_a;
  real_T c1_b_b;
  int32_T c1_idx;
  real_T c1_flt;
  boolean_T c1_p;
  int32_T c1_pipk;
  (void)chartInstance;
  for (c1_i277 = 0; c1_i277 < 6; c1_i277++) {
    c1_perm[c1_i277] = 1 + c1_i277;
  }

  for (c1_k = 0; c1_k < 5; c1_k++) {
    c1_b_k = 1.0 + (real_T)c1_k;
    c1_ipk = c1_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", c1_b_k), 1, 6, 1, 0) - 1];
    c1_b_a = c1_ipk;
    c1_b = c1_b_k;
    c1_c_a = c1_b_a;
    c1_b_b = c1_b;
    c1_idx = c1_c_a;
    c1_flt = c1_b_b;
    c1_p = ((real_T)c1_idx > c1_flt);
    if (c1_p) {
      c1_pipk = c1_perm[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c1_ipk), 1, 6, 1, 0) - 1];
      c1_perm[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c1_ipk), 1, 6, 1, 0) - 1] = c1_perm[_SFD_EML_ARRAY_BOUNDS_CHECK(
        "", (int32_T)_SFD_INTEGER_CHECK("", c1_b_k), 1, 6, 1, 0) - 1];
      c1_perm[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c1_b_k), 1, 6, 1, 0) - 1] = c1_pipk;
    }
  }
}

static void c1_eml_xtrsm(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_A[36], real_T c1_B[36], real_T c1_b_B[36])
{
  int32_T c1_i278;
  int32_T c1_i279;
  real_T c1_b_A[36];
  for (c1_i278 = 0; c1_i278 < 36; c1_i278++) {
    c1_b_B[c1_i278] = c1_B[c1_i278];
  }

  for (c1_i279 = 0; c1_i279 < 36; c1_i279++) {
    c1_b_A[c1_i279] = c1_A[c1_i279];
  }

  c1_b_eml_xtrsm(chartInstance, c1_b_A, c1_b_B);
}

static void c1_c_eml_switch_helper(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static real_T c1_b_norm(SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c1_x[36])
{
  real_T c1_y;
  int32_T c1_j;
  real_T c1_b_j;
  real_T c1_s;
  int32_T c1_i;
  real_T c1_b_i;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_b_y;
  real_T c1_d_x;
  boolean_T c1_b;
  boolean_T exitg1;
  (void)chartInstance;
  c1_y = 0.0;
  c1_j = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c1_j < 6)) {
    c1_b_j = 1.0 + (real_T)c1_j;
    c1_s = 0.0;
    for (c1_i = 0; c1_i < 6; c1_i++) {
      c1_b_i = 1.0 + (real_T)c1_i;
      c1_b_x = c1_x[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
                      ("", c1_b_i), 1, 6, 1, 0) + 6 *
                     (_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", c1_b_j), 1, 6, 2, 0) - 1)) - 1];
      c1_c_x = c1_b_x;
      c1_b_y = muDoubleScalarAbs(c1_c_x);
      c1_s += c1_b_y;
    }

    c1_d_x = c1_s;
    c1_b = muDoubleScalarIsNaN(c1_d_x);
    if (c1_b) {
      c1_y = rtNaN;
      exitg1 = true;
    } else {
      if (c1_s > c1_y) {
        c1_y = c1_s;
      }

      c1_j++;
    }
  }

  return c1_y;
}

static void c1_d_eml_scalar_eg(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_n_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_sprintf, const char_T *c1_identifier, char_T
  c1_y[14])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_sprintf), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_sprintf);
}

static void c1_o_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId,
  char_T c1_y[14])
{
  char_T c1_cv5[14];
  int32_T c1_i280;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_cv5, 1, 10, 0U, 1, 0U, 2, 1,
                14);
  for (c1_i280 = 0; c1_i280 < 14; c1_i280++) {
    c1_y[c1_i280] = c1_cv5[c1_i280];
  }

  sf_mex_destroy(&c1_u);
}

static const mxArray *c1_n_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(int32_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static int32_T c1_p_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_y;
  int32_T c1_i281;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_i281, 1, 6, 0U, 0, 0U, 0);
  c1_y = c1_i281;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_y;
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_p_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_q_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_b_is_active_c1_IPILCO_relativeRPY_S, const
  char_T *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_r_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_IPILCO_relativeRPY_S), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_IPILCO_relativeRPY_S);
  return c1_y;
}

static uint8_T c1_r_emlrt_marshallIn(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_eml_matlab_zgetrf(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_A[36], int32_T c1_ipiv[6], int32_T *c1_info)
{
  int32_T c1_i282;
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_b_a;
  int32_T c1_c_a;
  int32_T c1_jm1;
  int32_T c1_b;
  int32_T c1_b_b;
  int32_T c1_mmj;
  int32_T c1_d_a;
  int32_T c1_e_a;
  int32_T c1_c;
  int32_T c1_c_b;
  int32_T c1_d_b;
  int32_T c1_jj;
  int32_T c1_f_a;
  int32_T c1_g_a;
  int32_T c1_jp1j;
  int32_T c1_h_a;
  int32_T c1_i_a;
  int32_T c1_b_c;
  int32_T c1_i283;
  real_T c1_b_A[36];
  int32_T c1_j_a;
  int32_T c1_k_a;
  int32_T c1_jpiv_offset;
  int32_T c1_l_a;
  int32_T c1_e_b;
  int32_T c1_m_a;
  int32_T c1_f_b;
  int32_T c1_jpiv;
  int32_T c1_n_a;
  int32_T c1_g_b;
  int32_T c1_o_a;
  int32_T c1_h_b;
  int32_T c1_c_c;
  int32_T c1_i_b;
  int32_T c1_j_b;
  int32_T c1_jrow;
  int32_T c1_p_a;
  int32_T c1_k_b;
  int32_T c1_q_a;
  int32_T c1_l_b;
  int32_T c1_jprow;
  int32_T c1_ix0;
  int32_T c1_iy0;
  int32_T c1_b_ix0;
  int32_T c1_b_iy0;
  int32_T c1_c_ix0;
  int32_T c1_c_iy0;
  int32_T c1_ix;
  int32_T c1_iy;
  int32_T c1_k;
  real_T c1_temp;
  int32_T c1_r_a;
  int32_T c1_s_a;
  int32_T c1_b_jp1j;
  int32_T c1_t_a;
  int32_T c1_u_a;
  int32_T c1_d_c;
  int32_T c1_v_a;
  int32_T c1_m_b;
  int32_T c1_w_a;
  int32_T c1_n_b;
  int32_T c1_i284;
  int32_T c1_x_a;
  int32_T c1_o_b;
  int32_T c1_y_a;
  int32_T c1_p_b;
  boolean_T c1_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  real_T c1_x;
  real_T c1_y;
  real_T c1_b_x;
  real_T c1_b_y;
  real_T c1_z;
  int32_T c1_q_b;
  int32_T c1_r_b;
  int32_T c1_e_c;
  int32_T c1_ab_a;
  int32_T c1_bb_a;
  int32_T c1_f_c;
  int32_T c1_cb_a;
  int32_T c1_db_a;
  int32_T c1_g_c;
  c1_realmin(chartInstance);
  c1_eml_eps(chartInstance);
  for (c1_i282 = 0; c1_i282 < 6; c1_i282++) {
    c1_ipiv[c1_i282] = 1 + c1_i282;
  }

  *c1_info = 0;
  c1_eml_switch_helper(chartInstance);
  for (c1_j = 1; c1_j < 6; c1_j++) {
    c1_b_j = c1_j;
    c1_b_a = c1_b_j;
    c1_c_a = c1_b_a - 1;
    c1_jm1 = c1_c_a;
    c1_b = c1_b_j;
    c1_b_b = c1_b;
    c1_mmj = 6 - c1_b_b;
    c1_d_a = c1_jm1;
    c1_e_a = c1_d_a;
    c1_c = c1_e_a * 7;
    c1_c_b = c1_c;
    c1_d_b = c1_c_b + 1;
    c1_jj = c1_d_b;
    c1_f_a = c1_jj;
    c1_g_a = c1_f_a + 1;
    c1_jp1j = c1_g_a;
    c1_h_a = c1_mmj;
    c1_i_a = c1_h_a;
    c1_b_c = c1_i_a;
    for (c1_i283 = 0; c1_i283 < 36; c1_i283++) {
      c1_b_A[c1_i283] = c1_A[c1_i283];
    }

    c1_j_a = c1_eml_ixamax(chartInstance, c1_b_c + 1, c1_b_A, c1_jj);
    c1_k_a = c1_j_a - 1;
    c1_jpiv_offset = c1_k_a;
    c1_l_a = c1_jj;
    c1_e_b = c1_jpiv_offset;
    c1_m_a = c1_l_a;
    c1_f_b = c1_e_b;
    c1_jpiv = c1_m_a + c1_f_b;
    if (c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_jpiv), 1, 36, 1, 0) - 1] != 0.0) {
      if (c1_jpiv_offset != 0) {
        c1_n_a = c1_b_j;
        c1_g_b = c1_jpiv_offset;
        c1_o_a = c1_n_a;
        c1_h_b = c1_g_b;
        c1_c_c = c1_o_a + c1_h_b;
        c1_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_j), 1, 6, 1, 0) - 1] = c1_c_c;
        c1_i_b = c1_jm1;
        c1_j_b = c1_i_b + 1;
        c1_jrow = c1_j_b;
        c1_p_a = c1_jrow;
        c1_k_b = c1_jpiv_offset;
        c1_q_a = c1_p_a;
        c1_l_b = c1_k_b;
        c1_jprow = c1_q_a + c1_l_b;
        c1_ix0 = c1_jrow;
        c1_iy0 = c1_jprow;
        c1_b_ix0 = c1_ix0;
        c1_b_iy0 = c1_iy0;
        c1_b_eml_switch_helper(chartInstance);
        c1_c_ix0 = c1_b_ix0;
        c1_c_iy0 = c1_b_iy0;
        c1_ix = c1_c_ix0;
        c1_iy = c1_c_iy0;
        c1_eml_switch_helper(chartInstance);
        for (c1_k = 1; c1_k < 7; c1_k++) {
          c1_temp = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ix), 1, 36, 1, 0) - 1];
          c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_ix), 1, 36, 1, 0) - 1] = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK
            ("", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_iy), 1, 36, 1, 0) -
            1];
          c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_iy), 1, 36, 1, 0) - 1] = c1_temp;
          c1_r_a = c1_ix + 6;
          c1_ix = c1_r_a;
          c1_s_a = c1_iy + 6;
          c1_iy = c1_s_a;
        }
      }

      c1_b_jp1j = c1_jp1j;
      c1_t_a = c1_mmj;
      c1_u_a = c1_t_a;
      c1_d_c = c1_u_a;
      c1_v_a = c1_jp1j;
      c1_m_b = c1_d_c - 1;
      c1_w_a = c1_v_a;
      c1_n_b = c1_m_b;
      c1_i284 = c1_w_a + c1_n_b;
      c1_x_a = c1_b_jp1j;
      c1_o_b = c1_i284;
      c1_y_a = c1_x_a;
      c1_p_b = c1_o_b;
      if (c1_y_a > c1_p_b) {
        c1_overflow = false;
      } else {
        c1_eml_switch_helper(chartInstance);
        c1_eml_switch_helper(chartInstance);
        c1_overflow = (c1_p_b > 2147483646);
      }

      if (c1_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_overflow);
      }

      for (c1_i = c1_b_jp1j; c1_i <= c1_i284; c1_i++) {
        c1_b_i = c1_i;
        c1_x = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)c1_b_i), 1, 36, 1, 0) - 1];
        c1_y = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)c1_jj), 1, 36, 1, 0) - 1];
        c1_b_x = c1_x;
        c1_b_y = c1_y;
        c1_z = c1_b_x / c1_b_y;
        c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_i), 1, 36, 1, 0) - 1] = c1_z;
      }
    } else {
      *c1_info = c1_b_j;
    }

    c1_q_b = c1_b_j;
    c1_r_b = c1_q_b;
    c1_e_c = 6 - c1_r_b;
    c1_ab_a = c1_jj;
    c1_bb_a = c1_ab_a;
    c1_f_c = c1_bb_a;
    c1_cb_a = c1_jj;
    c1_db_a = c1_cb_a;
    c1_g_c = c1_db_a;
    c1_b_eml_xgeru(chartInstance, c1_mmj, c1_e_c, -1.0, c1_jp1j, c1_f_c + 6,
                   c1_A, c1_g_c + 7);
  }

  if (*c1_info == 0) {
    if (!(c1_A[35] != 0.0)) {
      *c1_info = 6;
    }
  }
}

static void c1_b_eml_xgeru(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, int32_T c1_m, int32_T c1_n, real_T c1_alpha1, int32_T c1_ix0,
  int32_T c1_iy0, real_T c1_A[36], int32_T c1_ia0)
{
  int32_T c1_b_m;
  int32_T c1_b_n;
  real_T c1_b_alpha1;
  int32_T c1_b_ix0;
  int32_T c1_b_iy0;
  int32_T c1_b_ia0;
  int32_T c1_c_m;
  int32_T c1_c_n;
  real_T c1_c_alpha1;
  int32_T c1_c_ix0;
  int32_T c1_c_iy0;
  int32_T c1_c_ia0;
  int32_T c1_d_m;
  int32_T c1_d_n;
  real_T c1_d_alpha1;
  int32_T c1_d_ix0;
  int32_T c1_d_iy0;
  int32_T c1_d_ia0;
  int32_T c1_e_m;
  int32_T c1_e_n;
  real_T c1_e_alpha1;
  int32_T c1_e_ix0;
  int32_T c1_e_iy0;
  int32_T c1_e_ia0;
  int32_T c1_ixstart;
  int32_T c1_b_a;
  int32_T c1_jA;
  int32_T c1_jy;
  int32_T c1_f_n;
  int32_T c1_b;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_j;
  real_T c1_yjy;
  real_T c1_temp;
  int32_T c1_ix;
  int32_T c1_c_b;
  int32_T c1_i285;
  int32_T c1_c_a;
  int32_T c1_d_b;
  int32_T c1_i286;
  int32_T c1_d_a;
  int32_T c1_e_b;
  int32_T c1_e_a;
  int32_T c1_f_b;
  boolean_T c1_b_overflow;
  int32_T c1_ijA;
  int32_T c1_b_ijA;
  int32_T c1_f_a;
  int32_T c1_g_a;
  int32_T c1_h_a;
  c1_b_m = c1_m;
  c1_b_n = c1_n;
  c1_b_alpha1 = c1_alpha1;
  c1_b_ix0 = c1_ix0;
  c1_b_iy0 = c1_iy0;
  c1_b_ia0 = c1_ia0;
  c1_c_m = c1_b_m;
  c1_c_n = c1_b_n;
  c1_c_alpha1 = c1_b_alpha1;
  c1_c_ix0 = c1_b_ix0;
  c1_c_iy0 = c1_b_iy0;
  c1_c_ia0 = c1_b_ia0;
  c1_d_m = c1_c_m;
  c1_d_n = c1_c_n;
  c1_d_alpha1 = c1_c_alpha1;
  c1_d_ix0 = c1_c_ix0;
  c1_d_iy0 = c1_c_iy0;
  c1_d_ia0 = c1_c_ia0;
  c1_e_m = c1_d_m;
  c1_e_n = c1_d_n;
  c1_e_alpha1 = c1_d_alpha1;
  c1_e_ix0 = c1_d_ix0;
  c1_e_iy0 = c1_d_iy0;
  c1_e_ia0 = c1_d_ia0;
  if (c1_e_alpha1 == 0.0) {
  } else {
    c1_ixstart = c1_e_ix0;
    c1_b_a = c1_e_ia0 - 1;
    c1_jA = c1_b_a;
    c1_jy = c1_e_iy0;
    c1_f_n = c1_e_n;
    c1_b = c1_f_n;
    c1_b_b = c1_b;
    if (1 > c1_b_b) {
      c1_overflow = false;
    } else {
      c1_eml_switch_helper(chartInstance);
      c1_eml_switch_helper(chartInstance);
      c1_overflow = (c1_b_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_j = 1; c1_j <= c1_f_n; c1_j++) {
      c1_yjy = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c1_jy), 1, 36, 1, 0) - 1];
      if (c1_yjy != 0.0) {
        c1_temp = c1_yjy * c1_e_alpha1;
        c1_ix = c1_ixstart;
        c1_c_b = c1_jA + 1;
        c1_i285 = c1_c_b;
        c1_c_a = c1_e_m;
        c1_d_b = c1_jA;
        c1_i286 = c1_c_a + c1_d_b;
        c1_d_a = c1_i285;
        c1_e_b = c1_i286;
        c1_e_a = c1_d_a;
        c1_f_b = c1_e_b;
        if (c1_e_a > c1_f_b) {
          c1_b_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_eml_switch_helper(chartInstance);
          c1_b_overflow = (c1_f_b > 2147483646);
        }

        if (c1_b_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
        }

        for (c1_ijA = c1_i285; c1_ijA <= c1_i286; c1_ijA++) {
          c1_b_ijA = c1_ijA;
          c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ijA), 1, 36, 1, 0) - 1] =
            c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ijA), 1, 36, 1, 0) - 1] +
            c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_ix), 1, 36, 1, 0) - 1] * c1_temp;
          c1_f_a = c1_ix + 1;
          c1_ix = c1_f_a;
        }
      }

      c1_g_a = c1_jy + 6;
      c1_jy = c1_g_a;
      c1_h_a = c1_jA + 6;
      c1_jA = c1_h_a;
    }
  }
}

static void c1_b_eml_xtrsm(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_A[36], real_T c1_B[36])
{
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_jBcol;
  int32_T c1_k;
  int32_T c1_b_k;
  int32_T c1_kAcol;
  int32_T c1_i287;
  int32_T c1_b;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  c1_c_eml_switch_helper(chartInstance);
  c1_eml_switch_helper(chartInstance);
  for (c1_j = 1; c1_j < 7; c1_j++) {
    c1_b_j = c1_j - 1;
    c1_jBcol = 6 * c1_b_j;
    c1_eml_switch_helper(chartInstance);
    for (c1_k = 6; c1_k > 0; c1_k--) {
      c1_b_k = c1_k;
      c1_kAcol = 6 * (c1_b_k - 1);
      if (c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c1_b_k + c1_jBcol)), 1, 36, 1, 0) - 1] != 0.0) {
        c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)(c1_b_k + c1_jBcol)), 1, 36, 1, 0) - 1] = c1_c_rdivide
          (chartInstance, c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)(c1_b_k + c1_jBcol)), 1, 36, 1, 0) -
           1], c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
             "", (real_T)(c1_b_k + c1_kAcol)), 1, 36, 1, 0) - 1]);
        c1_i287 = c1_b_k - 1;
        c1_b = c1_i287;
        c1_b_b = c1_b;
        if (1 > c1_b_b) {
          c1_overflow = false;
        } else {
          c1_eml_switch_helper(chartInstance);
          c1_eml_switch_helper(chartInstance);
          c1_overflow = (c1_b_b > 2147483646);
        }

        if (c1_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_overflow);
        }

        for (c1_i = 1; c1_i <= c1_i287; c1_i++) {
          c1_b_i = c1_i;
          c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c1_b_i + c1_jBcol)), 1, 36, 1, 0) - 1] =
            c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c1_b_i + c1_jBcol)), 1, 36, 1, 0) - 1] -
            c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c1_b_k + c1_jBcol)), 1, 36, 1, 0) - 1] *
            c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)(c1_b_i + c1_kAcol)), 1, 36, 1, 0) - 1];
        }
      }
    }
  }
}

static void init_dsm_address_info(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc1_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  chartInstance->c1_a = (real_T (*)[6])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c1_Fa = (real_T (*)[6])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 2);
  chartInstance->c1_dp_d = (real_T (*)[3])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 3);
  chartInstance->c1_Kp = (real_T (*)[36])ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c1_Kd = (real_T (*)[36])ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c1_He = (real_T (*)[16])ssGetInputPortSignal_wrapper
    (chartInstance->S, 2);
  chartInstance->c1_xe = (real_T (*)[6])ssGetInputPortSignal_wrapper
    (chartInstance->S, 3);
  chartInstance->c1_dxe_n = (real_T (*)[6])ssGetInputPortSignal_wrapper
    (chartInstance->S, 4);
  chartInstance->c1_F = (real_T (*)[6])ssGetInputPortSignal_wrapper
    (chartInstance->S, 5);
  chartInstance->c1_Hd = (real_T (*)[16])ssGetInputPortSignal_wrapper
    (chartInstance->S, 6);
  chartInstance->c1_prevHd = (real_T (*)[16])ssGetInputPortSignal_wrapper
    (chartInstance->S, 7);
  chartInstance->c1_prevdp_d = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 8);
  chartInstance->c1_prevT_e = (real_T (*)[9])ssGetInputPortSignal_wrapper
    (chartInstance->S, 9);
  chartInstance->c1_T_e = (real_T (*)[9])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 4);
  chartInstance->c1_dt = (real_T *)ssGetInputPortSignal_wrapper(chartInstance->S,
    10);
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c1_IPILCO_relativeRPY_S_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4214587420U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3338436263U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3323366938U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1518274374U);
}

mxArray* sf_c1_IPILCO_relativeRPY_S_get_post_codegen_info(void);
mxArray *sf_c1_IPILCO_relativeRPY_S_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("hJ2PtGyDDYKE3q6nLrSOPF");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,11,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(6);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(6);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(4);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(4);
      mxSetField(mxData,6,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,6,"type",mxType);
    }

    mxSetField(mxData,6,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(4);
      mxSetField(mxData,7,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,7,"type",mxType);
    }

    mxSetField(mxData,7,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,8,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,8,"type",mxType);
    }

    mxSetField(mxData,8,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(3);
      mxSetField(mxData,9,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,9,"type",mxType);
    }

    mxSetField(mxData,9,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,10,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,10,"type",mxType);
    }

    mxSetField(mxData,10,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(3);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxArray* mxPostCodegenInfo =
      sf_c1_IPILCO_relativeRPY_S_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c1_IPILCO_relativeRPY_S_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c1_IPILCO_relativeRPY_S_jit_fallback_info(void)
{
  const char *infoFields[] = { "fallbackType", "fallbackReason",
    "incompatibleSymbol", };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 3, infoFields);
  mxArray *fallbackReason = mxCreateString("feature_off");
  mxArray *incompatibleSymbol = mxCreateString("");
  mxArray *fallbackType = mxCreateString("early");
  mxSetField(mxInfo, 0, infoFields[0], fallbackType);
  mxSetField(mxInfo, 0, infoFields[1], fallbackReason);
  mxSetField(mxInfo, 0, infoFields[2], incompatibleSymbol);
  return mxInfo;
}

mxArray *sf_c1_IPILCO_relativeRPY_S_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c1_IPILCO_relativeRPY_S_get_post_codegen_info(void)
{
  const char* fieldNames[] = { "exportedFunctionsUsedByThisChart",
    "exportedFunctionsChecksum" };

  mwSize dims[2] = { 1, 1 };

  mxArray* mxPostCodegenInfo = mxCreateStructArray(2, dims, sizeof(fieldNames)/
    sizeof(fieldNames[0]), fieldNames);

  {
    mxArray* mxExportedFunctionsChecksum = mxCreateString("");
    mwSize exp_dims[2] = { 0, 1 };

    mxArray* mxExportedFunctionsUsedByThisChart = mxCreateCellArray(2, exp_dims);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsUsedByThisChart",
               mxExportedFunctionsUsedByThisChart);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsChecksum",
               mxExportedFunctionsChecksum);
  }

  return mxPostCodegenInfo;
}

static const mxArray *sf_get_sim_state_info_c1_IPILCO_relativeRPY_S(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x5'type','srcId','name','auxInfo'{{M[1],M[22],T\"Fa\",},{M[1],M[11],T\"T_e\",},{M[1],M[16],T\"a\",},{M[1],M[21],T\"dp_d\",},{M[8],M[0],T\"is_active_c1_IPILCO_relativeRPY_S\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 5, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_IPILCO_relativeRPY_S_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _IPILCO_relativeRPY_SMachineNumber_,
           1,
           1,
           1,
           0,
           15,
           0,
           0,
           0,
           0,
           7,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize its own list of scripts */
        init_script_number_translation(_IPILCO_relativeRPY_SMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_IPILCO_relativeRPY_SMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _IPILCO_relativeRPY_SMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,2,0,1,"a");
          _SFD_SET_DATA_PROPS(1,2,0,1,"Fa");
          _SFD_SET_DATA_PROPS(2,2,0,1,"dp_d");
          _SFD_SET_DATA_PROPS(3,1,1,0,"Kp");
          _SFD_SET_DATA_PROPS(4,1,1,0,"Kd");
          _SFD_SET_DATA_PROPS(5,1,1,0,"He");
          _SFD_SET_DATA_PROPS(6,1,1,0,"xe");
          _SFD_SET_DATA_PROPS(7,1,1,0,"dxe_n");
          _SFD_SET_DATA_PROPS(8,1,1,0,"F");
          _SFD_SET_DATA_PROPS(9,1,1,0,"Hd");
          _SFD_SET_DATA_PROPS(10,1,1,0,"prevHd");
          _SFD_SET_DATA_PROPS(11,1,1,0,"prevdp_d");
          _SFD_SET_DATA_PROPS(12,1,1,0,"prevT_e");
          _SFD_SET_DATA_PROPS(13,2,0,1,"T_e");
          _SFD_SET_DATA_PROPS(14,1,1,0,"dt");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,3493);
        _SFD_CV_INIT_SCRIPT(0,1,2,0,0,0,1,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"tr2rt",1354,-1,1749);
        _SFD_CV_INIT_SCRIPT_IF(0,0,1384,1411,-1,1454);
        _SFD_CV_INIT_SCRIPT_IF(0,1,1481,1497,1684,1748);
        _SFD_CV_INIT_SCRIPT_FOR(0,0,1573,1591,1679);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(0,0,1387,1411,-1,1);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(0,1,1484,1497,-1,4);
        _SFD_CV_INIT_SCRIPT(1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(1,0,"numcols",951,-1,991);
        _SFD_CV_INIT_SCRIPT(2,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(2,0,"numrows",946,-1,988);
        _SFD_CV_INIT_SCRIPT(3,1,9,0,0,0,0,0,2,1);
        _SFD_CV_INIT_SCRIPT_FCN(3,0,"transl",1953,-1,3159);
        _SFD_CV_INIT_SCRIPT_IF(3,0,1995,2009,3045,3063);
        _SFD_CV_INIT_SCRIPT_IF(3,1,2018,2031,2716,3036);
        _SFD_CV_INIT_SCRIPT_IF(3,2,2044,2060,2394,2707);
        _SFD_CV_INIT_SCRIPT_IF(3,3,2128,2143,2207,2226);
        _SFD_CV_INIT_SCRIPT_IF(3,4,2207,2226,-1,2226);
        _SFD_CV_INIT_SCRIPT_IF(3,5,2449,2480,2532,2551);
        _SFD_CV_INIT_SCRIPT_IF(3,6,2532,2551,-1,2551);
        _SFD_CV_INIT_SCRIPT_IF(3,7,2716,2737,2872,3036);
        _SFD_CV_INIT_SCRIPT_IF(3,8,3045,3063,-1,3063);

        {
          static int condStart[] = { 2452, 2468 };

          static int condEnd[] = { 2464, 2480 };

          static int pfixExpr[] = { 0, 1, -2 };

          _SFD_CV_INIT_SCRIPT_MCDC(3,0,2452,2480,2,0,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(3,0,1998,2009,-1,0);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(3,1,2047,2060,-1,0);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(3,4,2452,2464,-1,0);
        _SFD_CV_INIT_SCRIPT(4,1,2,0,0,0,0,0,2,1);
        _SFD_CV_INIT_SCRIPT_FCN(4,0,"ishomog",1166,-1,1399);
        _SFD_CV_INIT_SCRIPT_IF(4,0,1220,1237,1367,1398);
        _SFD_CV_INIT_SCRIPT_IF(4,1,1282,1300,-1,1361);

        {
          static int condStart[] = { 1285, 1290 };

          static int condEnd[] = { 1286, 1300 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_SCRIPT_MCDC(4,0,1285,1300,2,0,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(4,0,1223,1237,-1,5);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(4,1,1290,1300,-1,4);
        _SFD_CV_INIT_SCRIPT(5,1,6,0,0,0,1,0,5,3);
        _SFD_CV_INIT_SCRIPT_FCN(5,0,"my_tr2rpy",1704,-1,3253);
        _SFD_CV_INIT_SCRIPT_IF(5,0,1786,1807,-1,1852);
        _SFD_CV_INIT_SCRIPT_IF(5,1,1899,1915,2031,3198);
        _SFD_CV_INIT_SCRIPT_IF(5,2,2031,2042,2624,3198);
        _SFD_CV_INIT_SCRIPT_IF(5,3,2071,2112,2310,2619);
        _SFD_CV_INIT_SCRIPT_IF(5,4,2680,2721,2928,3190);
        _SFD_CV_INIT_SCRIPT_IF(5,5,3203,3213,-1,3249);
        _SFD_CV_INIT_SCRIPT_FOR(5,0,1942,1955,1992);

        {
          static int condStart[] = { 2035 };

          static int condEnd[] = { 2042 };

          static int pfixExpr[] = { 0, -1 };

          _SFD_CV_INIT_SCRIPT_MCDC(5,0,2034,2042,1,0,&(condStart[0]),&(condEnd[0]),
            2,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 2074, 2095 };

          static int condEnd[] = { 2091, 2112 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_SCRIPT_MCDC(5,1,2074,2112,2,1,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 2683, 2704 };

          static int condEnd[] = { 2700, 2721 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_SCRIPT_MCDC(5,2,2683,2721,2,3,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(5,0,1902,1915,-1,4);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(5,1,2074,2091,-1,2);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(5,2,2095,2112,-1,2);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(5,3,2683,2700,-1,2);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(5,4,2704,2721,-1,2);
        _SFD_CV_INIT_SCRIPT(6,1,3,0,0,0,0,0,8,1);
        _SFD_CV_INIT_SCRIPT_FCN(6,0,"cardatol",1,-1,1455);
        _SFD_CV_INIT_SCRIPT_IF(6,0,888,940,-1,1007);
        _SFD_CV_INIT_SCRIPT_IF(6,1,1009,1029,1051,1089);
        _SFD_CV_INIT_SCRIPT_IF(6,2,1130,1139,1305,1454);

        {
          static int condStart[] = { 893, 899, 905, 911, 917, 923, 929, 936 };

          static int condEnd[] = { 896, 902, 908, 914, 920, 926, 933, 940 };

          static int pfixExpr[] = { 0, 1, -2, 2, -2, 3, -2, 4, -2, 5, -2, 6, -2,
            7, -2 };

          _SFD_CV_INIT_SCRIPT_MCDC(6,0,893,940,8,0,&(condStart[0]),&(condEnd[0]),
            15,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,0,893,896,-1,2);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,1,899,902,-1,4);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,2,905,908,-1,2);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,3,911,914,-1,4);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,4,917,920,-1,2);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,5,923,926,-1,4);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,6,929,933,-1,0);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,7,936,940,-1,0);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,8,1014,1029,-1,0);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,9,1134,1138,-1,1);

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)
            c1_c_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)
            c1_c_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)
            c1_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_g_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_g_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_f_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(12,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(13,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)
            c1_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(14,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_d_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_VALUE_PTR(0U, *chartInstance->c1_a);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c1_Fa);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c1_dp_d);
        _SFD_SET_DATA_VALUE_PTR(3U, *chartInstance->c1_Kp);
        _SFD_SET_DATA_VALUE_PTR(4U, *chartInstance->c1_Kd);
        _SFD_SET_DATA_VALUE_PTR(5U, *chartInstance->c1_He);
        _SFD_SET_DATA_VALUE_PTR(6U, *chartInstance->c1_xe);
        _SFD_SET_DATA_VALUE_PTR(7U, *chartInstance->c1_dxe_n);
        _SFD_SET_DATA_VALUE_PTR(8U, *chartInstance->c1_F);
        _SFD_SET_DATA_VALUE_PTR(9U, *chartInstance->c1_Hd);
        _SFD_SET_DATA_VALUE_PTR(10U, *chartInstance->c1_prevHd);
        _SFD_SET_DATA_VALUE_PTR(11U, *chartInstance->c1_prevdp_d);
        _SFD_SET_DATA_VALUE_PTR(12U, *chartInstance->c1_prevT_e);
        _SFD_SET_DATA_VALUE_PTR(13U, *chartInstance->c1_T_e);
        _SFD_SET_DATA_VALUE_PTR(14U, chartInstance->c1_dt);
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _IPILCO_relativeRPY_SMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "MQjIm31OZ1w7wMIbr85yjH";
}

static void sf_opaque_initialize_c1_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c1_IPILCO_relativeRPY_S
    ((SFc1_IPILCO_relativeRPY_SInstanceStruct*) chartInstanceVar);
  initialize_c1_IPILCO_relativeRPY_S((SFc1_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c1_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  enable_c1_IPILCO_relativeRPY_S((SFc1_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_disable_c1_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  disable_c1_IPILCO_relativeRPY_S((SFc1_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c1_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  sf_gateway_c1_IPILCO_relativeRPY_S((SFc1_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c1_IPILCO_relativeRPY_S(SimStruct*
  S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c1_IPILCO_relativeRPY_S
    ((SFc1_IPILCO_relativeRPY_SInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
}

static void sf_opaque_set_sim_state_c1_IPILCO_relativeRPY_S(SimStruct* S, const
  mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c1_IPILCO_relativeRPY_S((SFc1_IPILCO_relativeRPY_SInstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c1_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_IPILCO_relativeRPY_SInstanceStruct*) chartInstanceVar
      )->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_IPILCO_relativeRPY_S_optimization_info();
    }

    finalize_c1_IPILCO_relativeRPY_S((SFc1_IPILCO_relativeRPY_SInstanceStruct*)
      chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_IPILCO_relativeRPY_S((SFc1_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_IPILCO_relativeRPY_S(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c1_IPILCO_relativeRPY_S
      ((SFc1_IPILCO_relativeRPY_SInstanceStruct*)(chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c1_IPILCO_relativeRPY_S(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_IPILCO_relativeRPY_S_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,1,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,1,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,1);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 7, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 8, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 9, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 10, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,11);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,4);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=4; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 11; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2606862326U));
  ssSetChecksum1(S,(1361738556U));
  ssSetChecksum2(S,(830630678U));
  ssSetChecksum3(S,(3398029618U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c1_IPILCO_relativeRPY_S(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_IPILCO_relativeRPY_S(SimStruct *S)
{
  SFc1_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc1_IPILCO_relativeRPY_SInstanceStruct *)utMalloc(sizeof
    (SFc1_IPILCO_relativeRPY_SInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc1_IPILCO_relativeRPY_SInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.enableChart =
    sf_opaque_enable_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c1_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.callAtomicSubchartUserFcn = NULL;
  chartInstance->chartInfo.callAtomicSubchartAutoFcn = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->checksum = SF_RUNTIME_INFO_CHECKSUM;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  crtInfo->compiledInfo = NULL;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  init_simulink_io_address(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c1_IPILCO_relativeRPY_S_method_dispatcher(SimStruct *S, int_T method, void *
  data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_IPILCO_relativeRPY_S(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_IPILCO_relativeRPY_S(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_IPILCO_relativeRPY_S(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_IPILCO_relativeRPY_S_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
