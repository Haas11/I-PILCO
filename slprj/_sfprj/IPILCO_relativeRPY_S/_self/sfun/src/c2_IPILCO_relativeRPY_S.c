/* Include files */

#include <stddef.h>
#include "blas.h"
#include "IPILCO_relativeRPY_S_sfun.h"
#include "c2_IPILCO_relativeRPY_S.h"
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
static const char * c2_debug_family_names[24] = { "Fext", "contact", "cf_v",
  "insertStiffness", "r_hole", "ywall_l", "ywall_r", "R0_n", "p0_n", "H0_n",
  "W0_n", "Wn_n", "nargin", "nargout", "Kp_env", "Kd_env", "xe", "dxe", "xc",
  "x_hole", "unused", "y", "inHole", "INSERTED" };

static const char * c2_b_debug_family_names[3] = { "c", "nargin", "nargout" };

static const char * c2_c_debug_family_names[3] = { "r", "nargin", "nargout" };

static const char * c2_d_debug_family_names[6] = { "ct", "st", "nargin",
  "nargout", "t", "R" };

static const char * c2_e_debug_family_names[6] = { "ct", "st", "nargin",
  "nargout", "t", "R" };

static const char * c2_f_debug_family_names[6] = { "ct", "st", "nargin",
  "nargout", "t", "R" };

static const char * c2_g_debug_family_names[10] = { "opt", "pitch", "yaw",
  "roll", "d2r", "type", "nargin", "nargout", "rpy_in", "R" };

static const char * c2_h_debug_family_names[3] = { "c", "nargin", "nargout" };

static const char * c2_i_debug_family_names[3] = { "r", "nargin", "nargout" };

static const char * c2_j_debug_family_names[3] = { "r", "nargin", "nargout" };

static const char * c2_k_debug_family_names[3] = { "c", "nargin", "nargout" };

static const char * c2_l_debug_family_names[5] = { "nargin", "nargout", "R", "t",
  "T" };

static const char * c2_m_debug_family_names[4] = { "h", "d", "nargin", "nargout"
};

static const char * c2_n_debug_family_names[4] = { "nargin", "nargout", "x",
  "t1" };

static const char * c2_o_debug_family_names[6] = { "d", "n", "nargin", "nargout",
  "T", "R" };

static const char * c2_p_debug_family_names[10] = { "f", "m", "k", "ft", "mt",
  "nargin", "nargout", "T", "W", "Wt" };

/* Function Declarations */
static void initialize_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void initialize_params_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void enable_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void disable_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c2_update_debugger_state_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void set_sim_state_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance, const mxArray *c2_st);
static void finalize_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void sf_gateway_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void mdl_start_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c2_chartstep_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void initSimStructsc2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c2_my_rpy2r(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c2_rpy_in[3], real_T c2_R[9]);
static void c2_numcols(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c2_numrows(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c2_rotx(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                    real_T c2_t, real_T c2_R[9]);
static void c2_roty(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                    real_T c2_t, real_T c2_R[9]);
static void c2_rotz(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                    real_T c2_t, real_T c2_R[9]);
static void c2_wtrans(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c2_T[16], real_T c2_W[6], real_T c2_Wt[6]);
static void c2_transl(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c2_x[16], real_T c2_t1[3]);
static void c2_t2r(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                   real_T c2_T[16], real_T c2_R[9]);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber);
static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData);
static boolean_T c2_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_b_INSERTED, const char_T *c2_identifier);
static boolean_T c2_b_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static real_T c2_c_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_b_inHole, const char_T *c2_identifier);
static real_T c2_d_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_e_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_b_y, const char_T *c2_identifier, real_T
  c2_c_y[12]);
static void c2_f_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[12]);
static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static real_T c2_g_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_h_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[6]);
static void c2_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_g_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_i_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[16]);
static void c2_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_j_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[3]);
static void c2_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_h_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_k_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[9]);
static void c2_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_i_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_l_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[3]);
static void c2_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_j_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_k_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static c2_s3TYGqYEZfhjEunCONSnN4C c2_m_emlrt_marshallIn
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance, const mxArray *c2_u,
   const emlrtMsgIdentifier *c2_parentId);
static boolean_T c2_n_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_l_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_m_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_info_helper(const mxArray **c2_info);
static const mxArray *c2_emlrt_marshallOut(const char * c2_u);
static const mxArray *c2_b_emlrt_marshallOut(const uint32_T c2_u);
static void c2_eml_scalar_eg(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static void c2_threshold(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance);
static void c2_cross(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                     real_T c2_a[3], real_T c2_b[3], real_T c2_c[3]);
static void c2_b_eml_scalar_eg(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static const mxArray *c2_n_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static int32_T c2_o_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static uint8_T c2_p_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_IPILCO_relativeRPY_S, const
  char_T *c2_identifier);
static uint8_T c2_q_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void init_dsm_address_info(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);
static void init_simulink_io_address(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  chartInstance->c2_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c2_inHole_not_empty = false;
  chartInstance->c2_INSERTED_not_empty = false;
  chartInstance->c2_is_active_c2_IPILCO_relativeRPY_S = 0U;
}

static void initialize_params_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c2_update_debugger_state_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  const mxArray *c2_st;
  const mxArray *c2_b_y = NULL;
  int32_T c2_i0;
  real_T c2_u[12];
  const mxArray *c2_c_y = NULL;
  boolean_T c2_hoistedGlobal;
  boolean_T c2_b_u;
  const mxArray *c2_d_y = NULL;
  real_T c2_b_hoistedGlobal;
  real_T c2_c_u;
  const mxArray *c2_e_y = NULL;
  uint8_T c2_c_hoistedGlobal;
  uint8_T c2_d_u;
  const mxArray *c2_f_y = NULL;
  c2_st = NULL;
  c2_st = NULL;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_createcellmatrix(4, 1), false);
  for (c2_i0 = 0; c2_i0 < 12; c2_i0++) {
    c2_u[c2_i0] = (*chartInstance->c2_y)[c2_i0];
  }

  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 12), false);
  sf_mex_setcell(c2_b_y, 0, c2_c_y);
  c2_hoistedGlobal = chartInstance->c2_INSERTED;
  c2_b_u = c2_hoistedGlobal;
  c2_d_y = NULL;
  if (!chartInstance->c2_INSERTED_not_empty) {
    sf_mex_assign(&c2_d_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c2_d_y, sf_mex_create("y", &c2_b_u, 11, 0U, 0U, 0U, 0), false);
  }

  sf_mex_setcell(c2_b_y, 1, c2_d_y);
  c2_b_hoistedGlobal = chartInstance->c2_inHole;
  c2_c_u = c2_b_hoistedGlobal;
  c2_e_y = NULL;
  if (!chartInstance->c2_inHole_not_empty) {
    sf_mex_assign(&c2_e_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c2_e_y, sf_mex_create("y", &c2_c_u, 0, 0U, 0U, 0U, 0), false);
  }

  sf_mex_setcell(c2_b_y, 2, c2_e_y);
  c2_c_hoistedGlobal = chartInstance->c2_is_active_c2_IPILCO_relativeRPY_S;
  c2_d_u = c2_c_hoistedGlobal;
  c2_f_y = NULL;
  sf_mex_assign(&c2_f_y, sf_mex_create("y", &c2_d_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c2_b_y, 3, c2_f_y);
  sf_mex_assign(&c2_st, c2_b_y, false);
  return c2_st;
}

static void set_sim_state_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance, const mxArray *c2_st)
{
  const mxArray *c2_u;
  real_T c2_dv0[12];
  int32_T c2_i1;
  chartInstance->c2_doneDoubleBufferReInit = true;
  c2_u = sf_mex_dup(c2_st);
  c2_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 0)), "y",
                        c2_dv0);
  for (c2_i1 = 0; c2_i1 < 12; c2_i1++) {
    (*chartInstance->c2_y)[c2_i1] = c2_dv0[c2_i1];
  }

  chartInstance->c2_INSERTED = c2_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c2_u, 1)), "INSERTED");
  chartInstance->c2_inHole = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c2_u, 2)), "inHole");
  chartInstance->c2_is_active_c2_IPILCO_relativeRPY_S = c2_p_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 3)),
     "is_active_c2_IPILCO_relativeRPY_S");
  sf_mex_destroy(&c2_u);
  c2_update_debugger_state_c2_IPILCO_relativeRPY_S(chartInstance);
  sf_mex_destroy(&c2_st);
}

static void finalize_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  int32_T c2_i2;
  int32_T c2_i3;
  int32_T c2_i4;
  int32_T c2_i5;
  int32_T c2_i6;
  int32_T c2_i7;
  int32_T c2_i8;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  for (c2_i2 = 0; c2_i2 < 6; c2_i2++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_Kp_env)[c2_i2], 0U);
  }

  chartInstance->c2_sfEvent = CALL_EVENT;
  c2_chartstep_c2_IPILCO_relativeRPY_S(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_IPILCO_relativeRPY_SMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c2_i3 = 0; c2_i3 < 12; c2_i3++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_y)[c2_i3], 1U);
  }

  for (c2_i4 = 0; c2_i4 < 6; c2_i4++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_Kd_env)[c2_i4], 2U);
  }

  for (c2_i5 = 0; c2_i5 < 6; c2_i5++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_xe)[c2_i5], 3U);
  }

  for (c2_i6 = 0; c2_i6 < 6; c2_i6++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_dxe)[c2_i6], 4U);
  }

  for (c2_i7 = 0; c2_i7 < 6; c2_i7++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_xc)[c2_i7], 5U);
  }

  for (c2_i8 = 0; c2_i8 < 3; c2_i8++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c2_x_hole)[c2_i8], 6U);
  }

  _SFD_DATA_RANGE_CHECK(*chartInstance->c2_unused, 7U);
}

static void mdl_start_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c2_chartstep_c2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  real_T c2_hoistedGlobal;
  int32_T c2_i9;
  real_T c2_b_Kp_env[6];
  int32_T c2_i10;
  real_T c2_b_Kd_env[6];
  int32_T c2_i11;
  real_T c2_b_xe[6];
  int32_T c2_i12;
  real_T c2_b_dxe[6];
  int32_T c2_i13;
  real_T c2_b_xc[6];
  int32_T c2_i14;
  real_T c2_b_x_hole[3];
  real_T c2_b_unused;
  uint32_T c2_debug_family_var_map[24];
  real_T c2_Fext[6];
  real_T c2_contact;
  real_T c2_cf_v;
  real_T c2_insertStiffness;
  real_T c2_r_hole;
  real_T c2_ywall_l;
  real_T c2_ywall_r;
  real_T c2_R0_n[9];
  real_T c2_p0_n[3];
  real_T c2_H0_n[16];
  real_T c2_W0_n[6];
  real_T c2_Wn_n[6];
  real_T c2_nargin = 7.0;
  real_T c2_nargout = 1.0;
  real_T c2_b_y[12];
  int32_T c2_i15;
  int32_T c2_i16;
  real_T c2_c_xe[3];
  real_T c2_dv1[9];
  int32_T c2_i17;
  int32_T c2_i18;
  int32_T c2_i19;
  real_T c2_R[9];
  int32_T c2_i20;
  real_T c2_t[3];
  uint32_T c2_b_debug_family_var_map[5];
  real_T c2_b_nargin = 2.0;
  real_T c2_b_nargout = 1.0;
  uint32_T c2_c_debug_family_var_map[3];
  real_T c2_c;
  real_T c2_c_nargin = 1.0;
  real_T c2_c_nargout = 1.0;
  real_T c2_r;
  real_T c2_d_nargin = 1.0;
  real_T c2_d_nargout = 1.0;
  real_T c2_b_r;
  real_T c2_e_nargin = 1.0;
  real_T c2_e_nargout = 1.0;
  real_T c2_c_r;
  real_T c2_f_nargin = 1.0;
  real_T c2_f_nargout = 1.0;
  real_T c2_b_c;
  real_T c2_g_nargin = 1.0;
  real_T c2_g_nargout = 1.0;
  real_T c2_c_c;
  real_T c2_h_nargin = 1.0;
  real_T c2_h_nargout = 1.0;
  int32_T c2_i21;
  int32_T c2_i22;
  int32_T c2_i23;
  int32_T c2_i24;
  int32_T c2_i25;
  int32_T c2_i26;
  int32_T c2_i27;
  static real_T c2_dv2[4] = { 0.0, 0.0, 0.0, 1.0 };

  int32_T c2_i28;
  int32_T c2_i29;
  real_T c2_b_H0_n[16];
  int32_T c2_i30;
  real_T c2_b_W0_n[6];
  real_T c2_dv3[6];
  int32_T c2_i31;
  int32_T c2_i32;
  int32_T c2_i33;
  int32_T c2_i34;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T guard5 = false;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  c2_hoistedGlobal = *chartInstance->c2_unused;
  for (c2_i9 = 0; c2_i9 < 6; c2_i9++) {
    c2_b_Kp_env[c2_i9] = (*chartInstance->c2_Kp_env)[c2_i9];
  }

  for (c2_i10 = 0; c2_i10 < 6; c2_i10++) {
    c2_b_Kd_env[c2_i10] = (*chartInstance->c2_Kd_env)[c2_i10];
  }

  for (c2_i11 = 0; c2_i11 < 6; c2_i11++) {
    c2_b_xe[c2_i11] = (*chartInstance->c2_xe)[c2_i11];
  }

  for (c2_i12 = 0; c2_i12 < 6; c2_i12++) {
    c2_b_dxe[c2_i12] = (*chartInstance->c2_dxe)[c2_i12];
  }

  for (c2_i13 = 0; c2_i13 < 6; c2_i13++) {
    c2_b_xc[c2_i13] = (*chartInstance->c2_xc)[c2_i13];
  }

  for (c2_i14 = 0; c2_i14 < 3; c2_i14++) {
    c2_b_x_hole[c2_i14] = (*chartInstance->c2_x_hole)[c2_i14];
  }

  c2_b_unused = c2_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 24U, 24U, c2_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_Fext, 0U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_contact, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_cf_v, 2U, c2_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_insertStiffness, 3U, c2_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_r_hole, 4U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_ywall_l, 5U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_ywall_r, 6U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_R0_n, 7U, c2_h_sf_marshallOut,
    c2_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_p0_n, 8U, c2_e_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_H0_n, 9U, c2_g_sf_marshallOut,
    c2_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_W0_n, 10U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_Wn_n, 11U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 12U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 13U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_Kp_env, 14U, c2_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_Kd_env, 15U, c2_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_xe, 16U, c2_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_dxe, 17U, c2_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_xc, 18U, c2_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_b_x_hole, 19U, c2_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_b_unused, 20U, c2_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_b_y, 21U, c2_c_sf_marshallOut,
    c2_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&chartInstance->c2_inHole, 22U,
    c2_b_sf_marshallOut, c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&chartInstance->c2_INSERTED, 23U,
    c2_sf_marshallOut, c2_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 10);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 11);
  if (CV_EML_IF(0, 1, 0, !chartInstance->c2_inHole_not_empty)) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 12);
    chartInstance->c2_inHole = 0.0;
    chartInstance->c2_inHole_not_empty = true;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 14);
  if (CV_EML_IF(0, 1, 1, !chartInstance->c2_INSERTED_not_empty)) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 15);
    chartInstance->c2_INSERTED = false;
    chartInstance->c2_INSERTED_not_empty = true;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 17);
  for (c2_i15 = 0; c2_i15 < 6; c2_i15++) {
    c2_Fext[c2_i15] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 18);
  c2_contact = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 21);
  c2_cf_v = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 22);
  c2_insertStiffness = 500.0;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 23);
  c2_r_hole = 0.01;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 24);
  c2_ywall_l = c2_b_x_hole[1] + c2_r_hole;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 24);
  c2_ywall_r = c2_b_x_hole[1] - c2_r_hole;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 27);
  for (c2_i16 = 0; c2_i16 < 3; c2_i16++) {
    c2_c_xe[c2_i16] = c2_b_xe[c2_i16 + 3];
  }

  c2_my_rpy2r(chartInstance, c2_c_xe, c2_dv1);
  for (c2_i17 = 0; c2_i17 < 9; c2_i17++) {
    c2_R0_n[c2_i17] = c2_dv1[c2_i17];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 28);
  for (c2_i18 = 0; c2_i18 < 3; c2_i18++) {
    c2_p0_n[c2_i18] = c2_b_xe[c2_i18];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 29);
  for (c2_i19 = 0; c2_i19 < 9; c2_i19++) {
    c2_R[c2_i19] = c2_R0_n[c2_i19];
  }

  for (c2_i20 = 0; c2_i20 < 3; c2_i20++) {
    c2_t[c2_i20] = c2_p0_n[c2_i20];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 5U, 5U, c2_l_debug_family_names,
    c2_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_nargin, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_nargout, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_R, 2U, c2_h_sf_marshallOut,
    c2_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_t, 3U, c2_e_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_H0_n, 4U, c2_g_sf_marshallOut,
    c2_f_sf_marshallIn);
  CV_SCRIPT_FCN(6, 0);
  _SFD_SCRIPT_CALL(6U, chartInstance->c2_sfEvent, 38);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c2_h_debug_family_names,
    c2_c_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_c, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_c_nargin, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_c_nargout, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, 31);
  c2_c = 3.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c2_i_debug_family_names,
    c2_c_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_r, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_d_nargin, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_d_nargout, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_SCRIPT_FCN(2, 0);
  _SFD_SCRIPT_CALL(2U, chartInstance->c2_sfEvent, 33);
  c2_r = 3.0;
  _SFD_SCRIPT_CALL(2U, chartInstance->c2_sfEvent, -33);
  _SFD_SYMBOL_SCOPE_POP();
  CV_SCRIPT_IF(6, 0, CV_RELATIONAL_EVAL(14U, 6U, 0, 3.0, 3.0, -1, 1U, 0));
  _SFD_SCRIPT_CALL(6U, chartInstance->c2_sfEvent, 41);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c2_i_debug_family_names,
    c2_c_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_r, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_e_nargin, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_e_nargout, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_SCRIPT_FCN(2, 0);
  _SFD_SCRIPT_CALL(2U, chartInstance->c2_sfEvent, 33);
  c2_b_r = 3.0;
  _SFD_SCRIPT_CALL(2U, chartInstance->c2_sfEvent, -33);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c2_j_debug_family_names,
    c2_c_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_c_r, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_f_nargin, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_f_nargout, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_SCRIPT_FCN(2, 0);
  _SFD_SCRIPT_CALL(2U, chartInstance->c2_sfEvent, 33);
  c2_c_r = 3.0;
  _SFD_SCRIPT_CALL(2U, chartInstance->c2_sfEvent, -33);
  _SFD_SYMBOL_SCOPE_POP();
  CV_SCRIPT_IF(6, 1, CV_RELATIONAL_EVAL(14U, 6U, 1, 3.0, 3.0, -1, 1U, 0));
  _SFD_SCRIPT_CALL(6U, chartInstance->c2_sfEvent, 45);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c2_k_debug_family_names,
    c2_c_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_c, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_g_nargin, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_g_nargout, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, 31);
  c2_b_c = 1.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
  CV_SCRIPT_IF(6, 2, CV_RELATIONAL_EVAL(14U, 6U, 2, 1.0, 1.0, -1, 1U, 0));
  _SFD_SCRIPT_CALL(6U, chartInstance->c2_sfEvent, 49);
  CV_SCRIPT_IF(6, 3, CV_RELATIONAL_EVAL(14U, 6U, 3, 1.0, 1.0, -1, 4U, 0));
  _SFD_SCRIPT_CALL(6U, chartInstance->c2_sfEvent, 57);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c2_h_debug_family_names,
    c2_c_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_c_c, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_h_nargin, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_h_nargout, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, 31);
  c2_c_c = 3.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
  c2_i21 = 0;
  c2_i22 = 0;
  for (c2_i23 = 0; c2_i23 < 3; c2_i23++) {
    for (c2_i24 = 0; c2_i24 < 3; c2_i24++) {
      c2_H0_n[c2_i24 + c2_i21] = c2_R[c2_i24 + c2_i22];
    }

    c2_i21 += 4;
    c2_i22 += 3;
  }

  for (c2_i25 = 0; c2_i25 < 3; c2_i25++) {
    c2_H0_n[c2_i25 + 12] = c2_t[c2_i25];
  }

  c2_i26 = 0;
  for (c2_i27 = 0; c2_i27 < 4; c2_i27++) {
    c2_H0_n[c2_i26 + 3] = c2_dv2[c2_i27];
    c2_i26 += 4;
  }

  _SFD_SCRIPT_CALL(6U, chartInstance->c2_sfEvent, -57);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 40);
  if (CV_EML_IF(0, 1, 2, CV_RELATIONAL_EVAL(4U, 0U, 0, c2_b_xe[0], c2_b_xc[0],
        -1, 4U, c2_b_xe[0] > c2_b_xc[0]))) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 41);
    c2_contact = 1.0;
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 42);
    guard5 = false;
    if (CV_EML_COND(0, 1, 0, CV_RELATIONAL_EVAL(4U, 0U, 1, c2_b_xe[1],
          c2_ywall_r, -1, 4U, c2_b_xe[1] > c2_ywall_r))) {
      if (CV_EML_COND(0, 1, 1, CV_RELATIONAL_EVAL(4U, 0U, 2, c2_b_xe[1],
            c2_ywall_l, -1, 2U, c2_b_xe[1] < c2_ywall_l))) {
        CV_EML_MCDC(0, 1, 0, true);
        CV_EML_IF(0, 1, 3, true);
        _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 43);
        chartInstance->c2_inHole = 1.0;
      } else {
        guard5 = true;
      }
    } else {
      guard5 = true;
    }

    if (guard5 == true) {
      CV_EML_MCDC(0, 1, 0, false);
      CV_EML_IF(0, 1, 3, false);
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 47);
  if (CV_EML_IF(0, 1, 4, CV_RELATIONAL_EVAL(4U, 0U, 3, c2_b_xe[1], c2_b_xc[1],
        -1, 4U, c2_b_xe[1] > c2_b_xc[1]))) {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 48);
    c2_Fext[1] = c2_b_Kp_env[1] * (c2_b_xe[1] - c2_b_xc[1]);
  } else {
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 49);
    if (CV_EML_IF(0, 1, 5, CV_RELATIONAL_EVAL(4U, 0U, 4, c2_b_xe[1], -c2_b_xc[1],
          -1, 2U, c2_b_xe[1] < -c2_b_xc[1]))) {
      _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 50);
      c2_Fext[1] = c2_b_Kp_env[1] * (-c2_b_xe[1] - c2_b_xc[1]);
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 54);
  guard1 = false;
  if (CV_EML_COND(0, 1, 2, c2_contact != 0.0)) {
    if (!CV_EML_COND(0, 1, 3, chartInstance->c2_inHole != 0.0)) {
      CV_EML_MCDC(0, 1, 1, true);
      CV_EML_IF(0, 1, 6, true);
      _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 56);
      c2_Fext[0] = c2_b_Kp_env[0] * (c2_b_xe[0] - c2_b_xc[0]);
      _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 57);
      if (CV_EML_IF(0, 1, 7, CV_RELATIONAL_EVAL(4U, 0U, 5, c2_b_dxe[0], 0.0, -1,
            4U, c2_b_dxe[0] > 0.0))) {
        _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 58);
        c2_Fext[0] += c2_b_Kd_env[0] * c2_b_dxe[0];
      }

      _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 65);
      c2_Fext[1] += 0.0 * c2_b_dxe[1];
      _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 66);
      c2_Fext[2] += 0.0 * c2_b_dxe[2];
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1 == true) {
    CV_EML_MCDC(0, 1, 1, false);
    CV_EML_IF(0, 1, 6, false);
    _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 68);
    guard2 = false;
    if (CV_EML_COND(0, 1, 4, c2_contact != 0.0)) {
      if (CV_EML_COND(0, 1, 5, chartInstance->c2_inHole != 0.0)) {
        CV_EML_MCDC(0, 1, 2, true);
        CV_EML_IF(0, 1, 8, true);
        _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 71);
        if (CV_EML_IF(0, 1, 9, CV_RELATIONAL_EVAL(4U, 0U, 6, c2_b_xe[0],
              c2_b_x_hole[0] + 0.05, -1, 4U, c2_b_xe[0] > c2_b_x_hole[0] + 0.05)))
        {
          _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 72);
          chartInstance->c2_INSERTED = true;
          _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 73);
          c2_Fext[0] += c2_b_Kp_env[0] * (c2_b_xe[0] - (c2_b_x_hole[0] + 0.05));
          _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 74);
          if (CV_EML_IF(0, 1, 10, CV_RELATIONAL_EVAL(4U, 0U, 7, c2_b_dxe[0], 0.0,
                -1, 4U, c2_b_dxe[0] > 0.0))) {
            _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 75);
            c2_Fext[0] += c2_b_Kd_env[0] * c2_b_dxe[0];
          }
        }

        _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 79);
        guard3 = false;
        guard4 = false;
        if (CV_EML_COND(0, 1, 6, CV_RELATIONAL_EVAL(4U, 0U, 8, c2_b_xe[0],
              c2_b_x_hole[0] + 0.016, -1, 4U, c2_b_xe[0] > c2_b_x_hole[0] +
              0.016))) {
          if (CV_EML_COND(0, 1, 7, CV_RELATIONAL_EVAL(4U, 0U, 9, c2_b_xe[0],
                c2_b_x_hole[0] + 0.033, -1, 2U, c2_b_xe[0] < c2_b_x_hole[0] +
                0.033))) {
            if (!CV_EML_COND(0, 1, 8, chartInstance->c2_INSERTED)) {
              CV_EML_MCDC(0, 1, 3, true);
              CV_EML_IF(0, 1, 11, true);
              _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 80);
              c2_Fext[0] += 500.0 * (c2_b_xe[0] - (c2_b_x_hole[0] + 0.016));
            } else {
              guard3 = true;
            }
          } else {
            guard4 = true;
          }
        } else {
          guard4 = true;
        }

        if (guard4 == true) {
          guard3 = true;
        }

        if (guard3 == true) {
          CV_EML_MCDC(0, 1, 3, false);
          CV_EML_IF(0, 1, 11, false);
        }

        _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 84);
        if (CV_EML_IF(0, 1, 12, CV_RELATIONAL_EVAL(4U, 0U, 10, c2_b_xe[1],
              c2_ywall_l, -1, 4U, c2_b_xe[1] > c2_ywall_l))) {
          _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 85);
          c2_Fext[1] += c2_b_Kp_env[1] * (c2_b_xe[1] - c2_ywall_l);
          _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 86);
          if (CV_EML_IF(0, 1, 13, CV_RELATIONAL_EVAL(4U, 0U, 11, c2_b_dxe[1],
                0.0, -1, 4U, c2_b_dxe[1] > 0.0))) {
            _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 87);
            c2_Fext[1] += c2_b_Kd_env[1] * c2_b_dxe[1];
          }

          _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 90);
          c2_Fext[0] += c2_Fext[1] * 0.0 * c2_b_dxe[0];
        } else {
          _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 92);
          if (CV_EML_IF(0, 1, 14, CV_RELATIONAL_EVAL(4U, 0U, 12, c2_b_xe[1],
                c2_ywall_r, -1, 2U, c2_b_xe[1] < c2_ywall_r))) {
            _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 93);
            c2_Fext[1] += c2_b_Kp_env[1] * (c2_b_xe[1] - c2_ywall_r);
            _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 94);
            if (CV_EML_IF(0, 1, 15, CV_RELATIONAL_EVAL(4U, 0U, 13, c2_b_dxe[1],
                  0.0, -1, 2U, c2_b_dxe[1] < 0.0))) {
              _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 95);
              c2_Fext[1] += c2_b_Kd_env[1] * c2_b_dxe[1];
            }

            _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 98);
            c2_Fext[0] += c2_Fext[1] * 0.0 * c2_b_dxe[0];
          }
        }
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2 == true) {
      CV_EML_MCDC(0, 1, 2, false);
      CV_EML_IF(0, 1, 8, false);
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 102);
  for (c2_i28 = 0; c2_i28 < 6; c2_i28++) {
    c2_W0_n[c2_i28] = -c2_Fext[c2_i28];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 103);
  for (c2_i29 = 0; c2_i29 < 16; c2_i29++) {
    c2_b_H0_n[c2_i29] = c2_H0_n[c2_i29];
  }

  for (c2_i30 = 0; c2_i30 < 6; c2_i30++) {
    c2_b_W0_n[c2_i30] = c2_W0_n[c2_i30];
  }

  c2_wtrans(chartInstance, c2_b_H0_n, c2_b_W0_n, c2_dv3);
  for (c2_i31 = 0; c2_i31 < 6; c2_i31++) {
    c2_Wn_n[c2_i31] = c2_dv3[c2_i31];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 104);
  for (c2_i32 = 0; c2_i32 < 6; c2_i32++) {
    c2_b_y[c2_i32] = c2_W0_n[c2_i32];
  }

  for (c2_i33 = 0; c2_i33 < 6; c2_i33++) {
    c2_b_y[c2_i33 + 6] = c2_Wn_n[c2_i33];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, -104);
  _SFD_SYMBOL_SCOPE_POP();
  for (c2_i34 = 0; c2_i34 < 12; c2_i34++) {
    (*chartInstance->c2_y)[c2_i34] = c2_b_y[c2_i34];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
}

static void initSimStructsc2_IPILCO_relativeRPY_S
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c2_my_rpy2r(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
  real_T c2_rpy_in[3], real_T c2_R[9])
{
  uint32_T c2_debug_family_var_map[10];
  c2_s3TYGqYEZfhjEunCONSnN4C c2_opt;
  real_T c2_pitch;
  real_T c2_yaw;
  real_T c2_roll;
  real_T c2_d2r;
  char_T c2_type[3];
  real_T c2_nargin = 2.0;
  real_T c2_nargout = 1.0;
  int32_T c2_i35;
  static char_T c2_cv0[3] = { 'x', 'y', 'z' };

  real_T c2_a[9];
  real_T c2_b[9];
  int32_T c2_i36;
  int32_T c2_i37;
  int32_T c2_i38;
  real_T c2_b_y[9];
  int32_T c2_i39;
  int32_T c2_i40;
  int32_T c2_i41;
  int32_T c2_i42;
  int32_T c2_i43;
  int32_T c2_i44;
  int32_T c2_i45;
  int32_T c2_i46;
  int32_T c2_i47;
  int32_T c2_i48;
  int32_T c2_i49;
  int32_T c2_i50;
  int32_T c2_i51;
  int32_T c2_i52;
  int32_T c2_i53;
  int32_T c2_i54;
  int32_T c2_i55;
  int32_T c2_i56;
  int32_T c2_i57;
  int32_T c2_i58;
  int32_T c2_i59;
  int32_T c2_i60;
  int32_T c2_i61;
  int32_T c2_i62;
  int32_T c2_i63;
  int32_T c2_i64;
  int32_T c2_i65;
  int32_T c2_i66;
  int32_T c2_i67;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 10U, 10U, c2_g_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_opt, 0U, c2_k_sf_marshallOut,
    c2_j_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_pitch, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_yaw, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_roll, 3U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_d2r, 4U, c2_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_type, 5U, c2_j_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 6U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 7U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_rpy_in, 8U, c2_i_sf_marshallOut,
    c2_i_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_R, 9U, c2_h_sf_marshallOut,
    c2_h_sf_marshallIn);
  for (c2_i35 = 0; c2_i35 < 3; c2_i35++) {
    c2_type[c2_i35] = c2_cv0[c2_i35];
  }

  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 50);
  c2_opt.zyx = false;
  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 51);
  c2_opt.deg = false;
  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 54);
  CV_SCRIPT_IF(0, 0, false);
  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 59);
  c2_numcols(chartInstance);
  CV_SCRIPT_IF(0, 1, CV_RELATIONAL_EVAL(14U, 0U, 0, 3.0, 3.0, -1, 0U, 1));
  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 60);
  c2_pitch = c2_rpy_in[1];
  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 61);
  c2_yaw = c2_rpy_in[2];
  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 62);
  c2_roll = c2_rpy_in[0];
  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 68);
  if (CV_SCRIPT_IF(0, 2, c2_opt.deg)) {
    _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 69);
    c2_d2r = 0.017453292519943295;
    _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 70);
    c2_roll *= 0.017453292519943295;
    _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 71);
    c2_pitch *= 0.017453292519943295;
    _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 72);
    c2_yaw *= 0.017453292519943295;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 75);
  if (CV_SCRIPT_IF(0, 3, CV_SCRIPT_MCDC(0, 0, !CV_SCRIPT_COND(0, 0, c2_opt.zyx))))
  {
    _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 77);
    c2_numrows(chartInstance);
    CV_SCRIPT_IF(0, 4, CV_RELATIONAL_EVAL(14U, 0U, 1, 1.0, 1.0, -1, 0U, 1));
    _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 78);
    c2_rotx(chartInstance, c2_roll, c2_a);
    c2_roty(chartInstance, c2_pitch, c2_b);
    c2_eml_scalar_eg(chartInstance);
    c2_eml_scalar_eg(chartInstance);
    c2_threshold(chartInstance);
    for (c2_i36 = 0; c2_i36 < 3; c2_i36++) {
      c2_i37 = 0;
      for (c2_i38 = 0; c2_i38 < 3; c2_i38++) {
        c2_b_y[c2_i37 + c2_i36] = 0.0;
        c2_i39 = 0;
        for (c2_i40 = 0; c2_i40 < 3; c2_i40++) {
          c2_b_y[c2_i37 + c2_i36] += c2_a[c2_i39 + c2_i36] * c2_b[c2_i40 +
            c2_i37];
          c2_i39 += 3;
        }

        c2_i37 += 3;
      }
    }

    c2_rotz(chartInstance, c2_yaw, c2_b);
    c2_eml_scalar_eg(chartInstance);
    c2_eml_scalar_eg(chartInstance);
    for (c2_i41 = 0; c2_i41 < 9; c2_i41++) {
      c2_R[c2_i41] = 0.0;
    }

    for (c2_i42 = 0; c2_i42 < 9; c2_i42++) {
      c2_R[c2_i42] = 0.0;
    }

    for (c2_i43 = 0; c2_i43 < 9; c2_i43++) {
      c2_a[c2_i43] = c2_R[c2_i43];
    }

    for (c2_i44 = 0; c2_i44 < 9; c2_i44++) {
      c2_R[c2_i44] = c2_a[c2_i44];
    }

    c2_threshold(chartInstance);
    for (c2_i45 = 0; c2_i45 < 9; c2_i45++) {
      c2_a[c2_i45] = c2_R[c2_i45];
    }

    for (c2_i46 = 0; c2_i46 < 9; c2_i46++) {
      c2_R[c2_i46] = c2_a[c2_i46];
    }

    for (c2_i47 = 0; c2_i47 < 3; c2_i47++) {
      c2_i48 = 0;
      for (c2_i49 = 0; c2_i49 < 3; c2_i49++) {
        c2_R[c2_i48 + c2_i47] = 0.0;
        c2_i50 = 0;
        for (c2_i51 = 0; c2_i51 < 3; c2_i51++) {
          c2_R[c2_i48 + c2_i47] += c2_b_y[c2_i50 + c2_i47] * c2_b[c2_i51 +
            c2_i48];
          c2_i50 += 3;
        }

        c2_i48 += 3;
      }
    }
  } else {
    _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 87);
    c2_numrows(chartInstance);
    CV_SCRIPT_IF(0, 5, CV_RELATIONAL_EVAL(14U, 0U, 2, 1.0, 1.0, -1, 0U, 1));
    _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 88);
    c2_rotz(chartInstance, c2_roll, c2_a);
    c2_roty(chartInstance, c2_pitch, c2_b);
    c2_eml_scalar_eg(chartInstance);
    c2_eml_scalar_eg(chartInstance);
    c2_threshold(chartInstance);
    for (c2_i52 = 0; c2_i52 < 3; c2_i52++) {
      c2_i53 = 0;
      for (c2_i54 = 0; c2_i54 < 3; c2_i54++) {
        c2_b_y[c2_i53 + c2_i52] = 0.0;
        c2_i55 = 0;
        for (c2_i56 = 0; c2_i56 < 3; c2_i56++) {
          c2_b_y[c2_i53 + c2_i52] += c2_a[c2_i55 + c2_i52] * c2_b[c2_i56 +
            c2_i53];
          c2_i55 += 3;
        }

        c2_i53 += 3;
      }
    }

    c2_rotx(chartInstance, c2_yaw, c2_b);
    c2_eml_scalar_eg(chartInstance);
    c2_eml_scalar_eg(chartInstance);
    for (c2_i57 = 0; c2_i57 < 9; c2_i57++) {
      c2_R[c2_i57] = 0.0;
    }

    for (c2_i58 = 0; c2_i58 < 9; c2_i58++) {
      c2_R[c2_i58] = 0.0;
    }

    for (c2_i59 = 0; c2_i59 < 9; c2_i59++) {
      c2_a[c2_i59] = c2_R[c2_i59];
    }

    for (c2_i60 = 0; c2_i60 < 9; c2_i60++) {
      c2_R[c2_i60] = c2_a[c2_i60];
    }

    c2_threshold(chartInstance);
    for (c2_i61 = 0; c2_i61 < 9; c2_i61++) {
      c2_a[c2_i61] = c2_R[c2_i61];
    }

    for (c2_i62 = 0; c2_i62 < 9; c2_i62++) {
      c2_R[c2_i62] = c2_a[c2_i62];
    }

    for (c2_i63 = 0; c2_i63 < 3; c2_i63++) {
      c2_i64 = 0;
      for (c2_i65 = 0; c2_i65 < 3; c2_i65++) {
        c2_R[c2_i64 + c2_i63] = 0.0;
        c2_i66 = 0;
        for (c2_i67 = 0; c2_i67 < 3; c2_i67++) {
          c2_R[c2_i64 + c2_i63] += c2_b_y[c2_i66 + c2_i63] * c2_b[c2_i67 +
            c2_i64];
          c2_i66 += 3;
        }

        c2_i64 += 3;
      }
    }
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, -92);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c2_numcols(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  uint32_T c2_debug_family_var_map[3];
  real_T c2_c;
  real_T c2_nargin = 1.0;
  real_T c2_nargout = 1.0;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c2_b_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_c, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, 31);
  c2_c = 3.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c2_numrows(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  uint32_T c2_debug_family_var_map[3];
  real_T c2_r;
  real_T c2_nargin = 1.0;
  real_T c2_nargout = 1.0;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c2_c_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_r, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_SCRIPT_FCN(2, 0);
  _SFD_SCRIPT_CALL(2U, chartInstance->c2_sfEvent, 33);
  c2_r = 1.0;
  _SFD_SCRIPT_CALL(2U, chartInstance->c2_sfEvent, -33);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c2_rotx(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                    real_T c2_t, real_T c2_R[9])
{
  uint32_T c2_debug_family_var_map[6];
  real_T c2_ct;
  real_T c2_st;
  real_T c2_nargin = 1.0;
  real_T c2_nargout = 1.0;
  real_T c2_x;
  real_T c2_b_x;
  real_T c2_c_x;
  real_T c2_d_x;
  int32_T c2_i68;
  int32_T c2_i69;
  static real_T c2_dv4[3] = { 1.0, 0.0, 0.0 };

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c2_d_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_ct, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_st, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 3U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_t, 4U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_R, 5U, c2_h_sf_marshallOut,
    c2_h_sf_marshallIn);
  CV_SCRIPT_FCN(3, 0);
  _SFD_SCRIPT_CALL(3U, chartInstance->c2_sfEvent, 33);
  CV_SCRIPT_COND(3, 0, CV_RELATIONAL_EVAL(14U, 3U, 0, 1.0, 1.0, -1, 4U, 0));
  CV_SCRIPT_MCDC(3, 0, false);
  CV_SCRIPT_IF(3, 0, false);
  _SFD_SCRIPT_CALL(3U, chartInstance->c2_sfEvent, 37);
  c2_x = c2_t;
  c2_ct = c2_x;
  c2_b_x = c2_ct;
  c2_ct = c2_b_x;
  c2_ct = muDoubleScalarCos(c2_ct);
  _SFD_SCRIPT_CALL(3U, chartInstance->c2_sfEvent, 38);
  c2_c_x = c2_t;
  c2_st = c2_c_x;
  c2_d_x = c2_st;
  c2_st = c2_d_x;
  c2_st = muDoubleScalarSin(c2_st);
  _SFD_SCRIPT_CALL(3U, chartInstance->c2_sfEvent, 39);
  c2_i68 = 0;
  for (c2_i69 = 0; c2_i69 < 3; c2_i69++) {
    c2_R[c2_i68] = c2_dv4[c2_i69];
    c2_i68 += 3;
  }

  c2_R[1] = 0.0;
  c2_R[4] = c2_ct;
  c2_R[7] = -c2_st;
  c2_R[2] = 0.0;
  c2_R[5] = c2_st;
  c2_R[8] = c2_ct;
  _SFD_SCRIPT_CALL(3U, chartInstance->c2_sfEvent, -39);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c2_roty(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                    real_T c2_t, real_T c2_R[9])
{
  uint32_T c2_debug_family_var_map[6];
  real_T c2_ct;
  real_T c2_st;
  real_T c2_nargin = 1.0;
  real_T c2_nargout = 1.0;
  real_T c2_x;
  real_T c2_b_x;
  real_T c2_c_x;
  real_T c2_d_x;
  int32_T c2_i70;
  int32_T c2_i71;
  static real_T c2_dv5[3] = { 0.0, 1.0, 0.0 };

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c2_e_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_ct, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_st, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 3U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_t, 4U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_R, 5U, c2_h_sf_marshallOut,
    c2_h_sf_marshallIn);
  CV_SCRIPT_FCN(4, 0);
  _SFD_SCRIPT_CALL(4U, chartInstance->c2_sfEvent, 32);
  CV_SCRIPT_COND(4, 0, CV_RELATIONAL_EVAL(14U, 4U, 0, 1.0, 1.0, -1, 4U, 0));
  CV_SCRIPT_MCDC(4, 0, false);
  CV_SCRIPT_IF(4, 0, false);
  _SFD_SCRIPT_CALL(4U, chartInstance->c2_sfEvent, 35);
  c2_x = c2_t;
  c2_ct = c2_x;
  c2_b_x = c2_ct;
  c2_ct = c2_b_x;
  c2_ct = muDoubleScalarCos(c2_ct);
  _SFD_SCRIPT_CALL(4U, chartInstance->c2_sfEvent, 36);
  c2_c_x = c2_t;
  c2_st = c2_c_x;
  c2_d_x = c2_st;
  c2_st = c2_d_x;
  c2_st = muDoubleScalarSin(c2_st);
  _SFD_SCRIPT_CALL(4U, chartInstance->c2_sfEvent, 37);
  c2_R[0] = c2_ct;
  c2_R[3] = 0.0;
  c2_R[6] = c2_st;
  c2_i70 = 0;
  for (c2_i71 = 0; c2_i71 < 3; c2_i71++) {
    c2_R[c2_i70 + 1] = c2_dv5[c2_i71];
    c2_i70 += 3;
  }

  c2_R[2] = -c2_st;
  c2_R[5] = 0.0;
  c2_R[8] = c2_ct;
  _SFD_SCRIPT_CALL(4U, chartInstance->c2_sfEvent, -37);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c2_rotz(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                    real_T c2_t, real_T c2_R[9])
{
  uint32_T c2_debug_family_var_map[6];
  real_T c2_ct;
  real_T c2_st;
  real_T c2_nargin = 1.0;
  real_T c2_nargout = 1.0;
  real_T c2_x;
  real_T c2_b_x;
  real_T c2_c_x;
  real_T c2_d_x;
  int32_T c2_i72;
  int32_T c2_i73;
  static real_T c2_dv6[3] = { 0.0, 0.0, 1.0 };

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c2_f_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_ct, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_st, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 3U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_t, 4U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_R, 5U, c2_h_sf_marshallOut,
    c2_h_sf_marshallIn);
  CV_SCRIPT_FCN(5, 0);
  _SFD_SCRIPT_CALL(5U, chartInstance->c2_sfEvent, 32);
  CV_SCRIPT_COND(5, 0, CV_RELATIONAL_EVAL(14U, 5U, 0, 1.0, 1.0, -1, 4U, 0));
  CV_SCRIPT_MCDC(5, 0, false);
  CV_SCRIPT_IF(5, 0, false);
  _SFD_SCRIPT_CALL(5U, chartInstance->c2_sfEvent, 36);
  c2_x = c2_t;
  c2_ct = c2_x;
  c2_b_x = c2_ct;
  c2_ct = c2_b_x;
  c2_ct = muDoubleScalarCos(c2_ct);
  _SFD_SCRIPT_CALL(5U, chartInstance->c2_sfEvent, 37);
  c2_c_x = c2_t;
  c2_st = c2_c_x;
  c2_d_x = c2_st;
  c2_st = c2_d_x;
  c2_st = muDoubleScalarSin(c2_st);
  _SFD_SCRIPT_CALL(5U, chartInstance->c2_sfEvent, 38);
  c2_R[0] = c2_ct;
  c2_R[3] = -c2_st;
  c2_R[6] = 0.0;
  c2_R[1] = c2_st;
  c2_R[4] = c2_ct;
  c2_R[7] = 0.0;
  c2_i72 = 0;
  for (c2_i73 = 0; c2_i73 < 3; c2_i73++) {
    c2_R[c2_i72 + 2] = c2_dv6[c2_i73];
    c2_i72 += 3;
  }

  _SFD_SCRIPT_CALL(5U, chartInstance->c2_sfEvent, -38);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c2_wtrans(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c2_T[16], real_T c2_W[6], real_T c2_Wt[6])
{
  uint32_T c2_debug_family_var_map[10];
  real_T c2_f[3];
  real_T c2_m[3];
  real_T c2_k[3];
  real_T c2_ft[3];
  real_T c2_mt[3];
  real_T c2_nargin = 2.0;
  real_T c2_nargout = 1.0;
  int32_T c2_i74;
  int32_T c2_i75;
  int32_T c2_i76;
  real_T c2_b_T[16];
  real_T c2_b[3];
  int32_T c2_i77;
  real_T c2_b_f[3];
  int32_T c2_i78;
  real_T c2_b_b[3];
  real_T c2_C[3];
  int32_T c2_i79;
  int32_T c2_i80;
  real_T c2_c_T[16];
  real_T c2_dv7[9];
  int32_T c2_i81;
  int32_T c2_i82;
  int32_T c2_i83;
  int32_T c2_i84;
  real_T c2_a[9];
  int32_T c2_i85;
  int32_T c2_i86;
  int32_T c2_i87;
  int32_T c2_i88;
  int32_T c2_i89;
  int32_T c2_i90;
  int32_T c2_i91;
  int32_T c2_i92;
  int32_T c2_i93;
  int32_T c2_i94;
  int32_T c2_i95;
  real_T c2_d_T[16];
  int32_T c2_i96;
  int32_T c2_i97;
  int32_T c2_i98;
  int32_T c2_i99;
  int32_T c2_i100;
  int32_T c2_i101;
  int32_T c2_i102;
  int32_T c2_i103;
  int32_T c2_i104;
  int32_T c2_i105;
  int32_T c2_i106;
  int32_T c2_i107;
  int32_T c2_i108;
  int32_T c2_i109;
  int32_T c2_i110;
  int32_T c2_i111;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 10U, 10U, c2_p_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_f, 0U, c2_e_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_m, 1U, c2_e_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_k, 2U, c2_e_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_ft, 3U, c2_e_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_mt, 4U, c2_e_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 5U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 6U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_T, 7U, c2_g_sf_marshallOut,
    c2_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_W, 8U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_Wt, 9U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  CV_SCRIPT_FCN(7, 0);
  _SFD_SCRIPT_CALL(7U, chartInstance->c2_sfEvent, 32);
  for (c2_i74 = 0; c2_i74 < 3; c2_i74++) {
    c2_f[c2_i74] = c2_W[c2_i74];
  }

  _SFD_SCRIPT_CALL(7U, chartInstance->c2_sfEvent, 32);
  for (c2_i75 = 0; c2_i75 < 3; c2_i75++) {
    c2_m[c2_i75] = c2_W[c2_i75 + 3];
  }

  _SFD_SCRIPT_CALL(7U, chartInstance->c2_sfEvent, 33);
  for (c2_i76 = 0; c2_i76 < 16; c2_i76++) {
    c2_b_T[c2_i76] = c2_T[c2_i76];
  }

  c2_transl(chartInstance, c2_b_T, c2_b);
  for (c2_i77 = 0; c2_i77 < 3; c2_i77++) {
    c2_b_f[c2_i77] = c2_f[c2_i77];
  }

  for (c2_i78 = 0; c2_i78 < 3; c2_i78++) {
    c2_b_b[c2_i78] = c2_b[c2_i78];
  }

  c2_cross(chartInstance, c2_b_f, c2_b_b, c2_C);
  for (c2_i79 = 0; c2_i79 < 3; c2_i79++) {
    c2_k[c2_i79] = c2_C[c2_i79] + c2_m[c2_i79];
  }

  _SFD_SCRIPT_CALL(7U, chartInstance->c2_sfEvent, 35);
  for (c2_i80 = 0; c2_i80 < 16; c2_i80++) {
    c2_c_T[c2_i80] = c2_T[c2_i80];
  }

  c2_t2r(chartInstance, c2_c_T, c2_dv7);
  c2_i81 = 0;
  for (c2_i82 = 0; c2_i82 < 3; c2_i82++) {
    c2_i83 = 0;
    for (c2_i84 = 0; c2_i84 < 3; c2_i84++) {
      c2_a[c2_i84 + c2_i81] = c2_dv7[c2_i83 + c2_i82];
      c2_i83 += 3;
    }

    c2_i81 += 3;
  }

  for (c2_i85 = 0; c2_i85 < 3; c2_i85++) {
    c2_b[c2_i85] = c2_f[c2_i85];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i86 = 0; c2_i86 < 3; c2_i86++) {
    c2_ft[c2_i86] = 0.0;
  }

  for (c2_i87 = 0; c2_i87 < 3; c2_i87++) {
    c2_ft[c2_i87] = 0.0;
  }

  for (c2_i88 = 0; c2_i88 < 3; c2_i88++) {
    c2_C[c2_i88] = c2_ft[c2_i88];
  }

  for (c2_i89 = 0; c2_i89 < 3; c2_i89++) {
    c2_ft[c2_i89] = c2_C[c2_i89];
  }

  c2_threshold(chartInstance);
  for (c2_i90 = 0; c2_i90 < 3; c2_i90++) {
    c2_C[c2_i90] = c2_ft[c2_i90];
  }

  for (c2_i91 = 0; c2_i91 < 3; c2_i91++) {
    c2_ft[c2_i91] = c2_C[c2_i91];
  }

  for (c2_i92 = 0; c2_i92 < 3; c2_i92++) {
    c2_ft[c2_i92] = 0.0;
    c2_i93 = 0;
    for (c2_i94 = 0; c2_i94 < 3; c2_i94++) {
      c2_ft[c2_i92] += c2_a[c2_i93 + c2_i92] * c2_b[c2_i94];
      c2_i93 += 3;
    }
  }

  _SFD_SCRIPT_CALL(7U, chartInstance->c2_sfEvent, 36);
  for (c2_i95 = 0; c2_i95 < 16; c2_i95++) {
    c2_d_T[c2_i95] = c2_T[c2_i95];
  }

  c2_t2r(chartInstance, c2_d_T, c2_dv7);
  c2_i96 = 0;
  for (c2_i97 = 0; c2_i97 < 3; c2_i97++) {
    c2_i98 = 0;
    for (c2_i99 = 0; c2_i99 < 3; c2_i99++) {
      c2_a[c2_i99 + c2_i96] = c2_dv7[c2_i98 + c2_i97];
      c2_i98 += 3;
    }

    c2_i96 += 3;
  }

  for (c2_i100 = 0; c2_i100 < 3; c2_i100++) {
    c2_b[c2_i100] = c2_k[c2_i100];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i101 = 0; c2_i101 < 3; c2_i101++) {
    c2_mt[c2_i101] = 0.0;
  }

  for (c2_i102 = 0; c2_i102 < 3; c2_i102++) {
    c2_mt[c2_i102] = 0.0;
  }

  for (c2_i103 = 0; c2_i103 < 3; c2_i103++) {
    c2_C[c2_i103] = c2_mt[c2_i103];
  }

  for (c2_i104 = 0; c2_i104 < 3; c2_i104++) {
    c2_mt[c2_i104] = c2_C[c2_i104];
  }

  c2_threshold(chartInstance);
  for (c2_i105 = 0; c2_i105 < 3; c2_i105++) {
    c2_C[c2_i105] = c2_mt[c2_i105];
  }

  for (c2_i106 = 0; c2_i106 < 3; c2_i106++) {
    c2_mt[c2_i106] = c2_C[c2_i106];
  }

  for (c2_i107 = 0; c2_i107 < 3; c2_i107++) {
    c2_mt[c2_i107] = 0.0;
    c2_i108 = 0;
    for (c2_i109 = 0; c2_i109 < 3; c2_i109++) {
      c2_mt[c2_i107] += c2_a[c2_i108 + c2_i107] * c2_b[c2_i109];
      c2_i108 += 3;
    }
  }

  _SFD_SCRIPT_CALL(7U, chartInstance->c2_sfEvent, 38);
  for (c2_i110 = 0; c2_i110 < 3; c2_i110++) {
    c2_Wt[c2_i110] = c2_ft[c2_i110];
  }

  for (c2_i111 = 0; c2_i111 < 3; c2_i111++) {
    c2_Wt[c2_i111 + 3] = c2_mt[c2_i111];
  }

  _SFD_SCRIPT_CALL(7U, chartInstance->c2_sfEvent, -38);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c2_transl(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                      real_T c2_x[16], real_T c2_t1[3])
{
  uint32_T c2_debug_family_var_map[4];
  real_T c2_nargin = 1.0;
  real_T c2_nargout = 1.0;
  boolean_T c2_h;
  real_T c2_d[2];
  real_T c2_b_nargin = 1.0;
  real_T c2_b_nargout = 1.0;
  int32_T c2_i112;
  int32_T c2_i113;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c2_n_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 0U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 1U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_x, 2U, c2_g_sf_marshallOut,
    c2_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_t1, 3U, c2_e_sf_marshallOut,
    c2_g_sf_marshallIn);
  CV_SCRIPT_FCN(8, 0);
  _SFD_SCRIPT_CALL(8U, chartInstance->c2_sfEvent, 53);
  CV_SCRIPT_IF(8, 0, CV_RELATIONAL_EVAL(14U, 8U, 0, 1.0, 1.0, -1, 0U, 1));
  _SFD_SCRIPT_CALL(8U, chartInstance->c2_sfEvent, 54);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c2_m_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_h, 0U, c2_m_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_d, 1U, c2_l_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_nargin, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_nargout, 3U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  CV_SCRIPT_FCN(9, 0);
  _SFD_SCRIPT_CALL(9U, chartInstance->c2_sfEvent, 37);
  for (c2_i112 = 0; c2_i112 < 2; c2_i112++) {
    c2_d[c2_i112] = 4.0;
  }

  _SFD_SCRIPT_CALL(9U, chartInstance->c2_sfEvent, 38);
  CV_SCRIPT_IF(9, 0, CV_RELATIONAL_EVAL(14U, 9U, 0, 2.0, 2.0, -1, 5U, 1));
  _SFD_SCRIPT_CALL(9U, chartInstance->c2_sfEvent, 39);
  c2_h = true;
  _SFD_SCRIPT_CALL(9U, chartInstance->c2_sfEvent, 41);
  CV_SCRIPT_COND(9, 0, c2_h);
  CV_SCRIPT_COND(9, 1, CV_RELATIONAL_EVAL(14U, 9U, 1, 1.0, 1.0, -1, 4U, 0));
  CV_SCRIPT_MCDC(9, 0, false);
  CV_SCRIPT_IF(9, 1, false);
  _SFD_SCRIPT_CALL(9U, chartInstance->c2_sfEvent, -46);
  _SFD_SYMBOL_SCOPE_POP();
  CV_SCRIPT_IF(8, 1, true);
  _SFD_SCRIPT_CALL(8U, chartInstance->c2_sfEvent, 55);
  CV_SCRIPT_IF(8, 2, CV_RELATIONAL_EVAL(14U, 8U, 1, 2.0, 3.0, -1, 0U, 0));
  _SFD_SCRIPT_CALL(8U, chartInstance->c2_sfEvent, 66);
  CV_SCRIPT_COND(8, 0, CV_RELATIONAL_EVAL(14U, 8U, 4, 1.0, 1.0, -1, 0U, 1));
  CV_SCRIPT_MCDC(8, 0, true);
  CV_SCRIPT_IF(8, 5, true);
  _SFD_SCRIPT_CALL(8U, chartInstance->c2_sfEvent, 67);
  for (c2_i113 = 0; c2_i113 < 3; c2_i113++) {
    c2_t1[c2_i113] = c2_x[c2_i113 + 12];
  }

  _SFD_SCRIPT_CALL(8U, chartInstance->c2_sfEvent, -89);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c2_t2r(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                   real_T c2_T[16], real_T c2_R[9])
{
  uint32_T c2_debug_family_var_map[6];
  real_T c2_d[2];
  real_T c2_n;
  real_T c2_nargin = 1.0;
  real_T c2_nargout = 1.0;
  int32_T c2_i114;
  boolean_T c2_b_y;
  int32_T c2_k;
  real_T c2_b_k;
  static boolean_T c2_bv0[2] = { false, true };

  boolean_T c2_b0;
  int32_T c2_i115;
  int32_T c2_i116;
  int32_T c2_i117;
  int32_T c2_i118;
  boolean_T exitg1;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c2_o_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_d, 0U, c2_l_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_n, 1U, c2_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 2U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 3U, c2_d_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_T, 4U, c2_g_sf_marshallOut,
    c2_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_R, 5U, c2_h_sf_marshallOut,
    c2_h_sf_marshallIn);
  CV_SCRIPT_FCN(10, 0);
  _SFD_SCRIPT_CALL(10U, chartInstance->c2_sfEvent, 38);
  for (c2_i114 = 0; c2_i114 < 2; c2_i114++) {
    c2_d[c2_i114] = 4.0;
  }

  _SFD_SCRIPT_CALL(10U, chartInstance->c2_sfEvent, 39);
  CV_SCRIPT_IF(10, 0, CV_RELATIONAL_EVAL(14U, 10U, 0, c2_d[0], c2_d[1], -1, 1U,
    c2_d[0] != c2_d[1]));
  _SFD_SCRIPT_CALL(10U, chartInstance->c2_sfEvent, 42);
  c2_b_y = false;
  c2_k = 0;
  exitg1 = false;
  while ((exitg1 == false) && (c2_k < 2)) {
    c2_b_k = 1.0 + (real_T)c2_k;
    if ((real_T)c2_bv0[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
         _SFD_INTEGER_CHECK("", c2_b_k), 1, 2, 1, 0) - 1] == 0.0) {
      c2_b0 = true;
    } else {
      (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        c2_b_k), 1, 2, 1, 0);
      c2_b0 = false;
    }

    if (!c2_b0) {
      c2_b_y = true;
      exitg1 = true;
    } else {
      c2_k++;
    }
  }

  CV_SCRIPT_IF(10, 1, CV_SCRIPT_MCDC(10, 0, !CV_SCRIPT_COND(10, 0, c2_b_y)));
  _SFD_SCRIPT_CALL(10U, chartInstance->c2_sfEvent, 46);
  c2_n = 4.0;
  _SFD_SCRIPT_CALL(10U, chartInstance->c2_sfEvent, 48);
  CV_SCRIPT_IF(10, 2, CV_RELATIONAL_EVAL(14U, 10U, 1, 2.0, 2.0, -1, 0U, 1));
  _SFD_SCRIPT_CALL(10U, chartInstance->c2_sfEvent, 50);
  c2_i115 = 0;
  c2_i116 = 0;
  for (c2_i117 = 0; c2_i117 < 3; c2_i117++) {
    for (c2_i118 = 0; c2_i118 < 3; c2_i118++) {
      c2_R[c2_i118 + c2_i115] = c2_T[c2_i118 + c2_i116];
    }

    c2_i115 += 3;
    c2_i116 += 4;
  }

  _SFD_SCRIPT_CALL(10U, chartInstance->c2_sfEvent, -55);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber, uint32_T c2_instanceNumber)
{
  (void)c2_machineNumber;
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 0U,
    sf_debug_get_script_id("D:\\victor\\MATLAB\\I-PILCO\\util\\my_rpy2r.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 1U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\common\\numcols.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 2U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\common\\numrows.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 3U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\robot\\rotx.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 4U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\robot\\roty.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 5U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\robot\\rotz.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 6U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\robot\\rt2tr.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 7U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\robot\\wtrans.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 8U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\robot\\transl.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 9U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\common\\ishomog.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, c2_instanceNumber, 10U,
    sf_debug_get_script_id(
    "D:\\victor\\MATLAB\\Toolboxes\\robot-9.10\\rvctools\\robot\\t2r.m"));
}

static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  boolean_T c2_u;
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(boolean_T *)c2_inData;
  c2_b_y = NULL;
  if (!chartInstance->c2_INSERTED_not_empty) {
    sf_mex_assign(&c2_b_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 11, 0U, 0U, 0U, 0), false);
  }

  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static boolean_T c2_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_b_INSERTED, const char_T *c2_identifier)
{
  boolean_T c2_b_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_y = c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_INSERTED),
    &c2_thisId);
  sf_mex_destroy(&c2_b_INSERTED);
  return c2_b_y;
}

static boolean_T c2_b_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  boolean_T c2_b_y;
  boolean_T c2_b1;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_INSERTED_not_empty = false;
  } else {
    chartInstance->c2_INSERTED_not_empty = true;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_b1, 1, 11, 0U, 0, 0U, 0);
    c2_b_y = c2_b1;
  }

  sf_mex_destroy(&c2_u);
  return c2_b_y;
}

static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_INSERTED;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  boolean_T c2_b_y;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_b_INSERTED = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_y = c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_INSERTED),
    &c2_thisId);
  sf_mex_destroy(&c2_b_INSERTED);
  *(boolean_T *)c2_outData = c2_b_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  real_T c2_u;
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(real_T *)c2_inData;
  c2_b_y = NULL;
  if (!chartInstance->c2_inHole_not_empty) {
    sf_mex_assign(&c2_b_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  false);
  } else {
    sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), false);
  }

  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static real_T c2_c_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_b_inHole, const char_T *c2_identifier)
{
  real_T c2_b_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_y = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_inHole),
    &c2_thisId);
  sf_mex_destroy(&c2_b_inHole);
  return c2_b_y;
}

static real_T c2_d_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_b_y;
  real_T c2_d0;
  if (mxIsEmpty(c2_u)) {
    chartInstance->c2_inHole_not_empty = false;
  } else {
    chartInstance->c2_inHole_not_empty = true;
    sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d0, 1, 0, 0U, 0, 0U, 0);
    c2_b_y = c2_d0;
  }

  sf_mex_destroy(&c2_u);
  return c2_b_y;
}

static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_inHole;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_b_y;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_b_inHole = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_y = c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_inHole),
    &c2_thisId);
  sf_mex_destroy(&c2_b_inHole);
  *(real_T *)c2_outData = c2_b_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i119;
  real_T c2_b_inData[12];
  int32_T c2_i120;
  real_T c2_u[12];
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i119 = 0; c2_i119 < 12; c2_i119++) {
    c2_b_inData[c2_i119] = (*(real_T (*)[12])c2_inData)[c2_i119];
  }

  for (c2_i120 = 0; c2_i120 < 12; c2_i120++) {
    c2_u[c2_i120] = c2_b_inData[c2_i120];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 12), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static void c2_e_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_b_y, const char_T *c2_identifier, real_T
  c2_c_y[12])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_y), &c2_thisId, c2_c_y);
  sf_mex_destroy(&c2_b_y);
}

static void c2_f_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[12])
{
  real_T c2_dv8[12];
  int32_T c2_i121;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv8, 1, 0, 0U, 1, 0U, 1, 12);
  for (c2_i121 = 0; c2_i121 < 12; c2_i121++) {
    c2_b_y[c2_i121] = c2_dv8[c2_i121];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_y;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_c_y[12];
  int32_T c2_i122;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_b_y = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_y), &c2_thisId, c2_c_y);
  sf_mex_destroy(&c2_b_y);
  for (c2_i122 = 0; c2_i122 < 12; c2_i122++) {
    (*(real_T (*)[12])c2_outData)[c2_i122] = c2_c_y[c2_i122];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  real_T c2_u;
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(real_T *)c2_inData;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i123;
  real_T c2_b_inData[3];
  int32_T c2_i124;
  real_T c2_u[3];
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i123 = 0; c2_i123 < 3; c2_i123++) {
    c2_b_inData[c2_i123] = (*(real_T (*)[3])c2_inData)[c2_i123];
  }

  for (c2_i124 = 0; c2_i124 < 3; c2_i124++) {
    c2_u[c2_i124] = c2_b_inData[c2_i124];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i125;
  real_T c2_b_inData[6];
  int32_T c2_i126;
  real_T c2_u[6];
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i125 = 0; c2_i125 < 6; c2_i125++) {
    c2_b_inData[c2_i125] = (*(real_T (*)[6])c2_inData)[c2_i125];
  }

  for (c2_i126 = 0; c2_i126 < 6; c2_i126++) {
    c2_u[c2_i126] = c2_b_inData[c2_i126];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 6), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static real_T c2_g_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_b_y;
  real_T c2_d1;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d1, 1, 0, 0U, 0, 0U, 0);
  c2_b_y = c2_d1;
  sf_mex_destroy(&c2_u);
  return c2_b_y;
}

static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_nargout;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_b_y;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_nargout = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_y = c2_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_nargout),
    &c2_thisId);
  sf_mex_destroy(&c2_nargout);
  *(real_T *)c2_outData = c2_b_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static void c2_h_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[6])
{
  real_T c2_dv9[6];
  int32_T c2_i127;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv9, 1, 0, 0U, 1, 0U, 1, 6);
  for (c2_i127 = 0; c2_i127 < 6; c2_i127++) {
    c2_b_y[c2_i127] = c2_dv9[c2_i127];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_Wn_n;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_b_y[6];
  int32_T c2_i128;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_Wn_n = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_Wn_n), &c2_thisId, c2_b_y);
  sf_mex_destroy(&c2_Wn_n);
  for (c2_i128 = 0; c2_i128 < 6; c2_i128++) {
    (*(real_T (*)[6])c2_outData)[c2_i128] = c2_b_y[c2_i128];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_g_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i129;
  int32_T c2_i130;
  int32_T c2_i131;
  real_T c2_b_inData[16];
  int32_T c2_i132;
  int32_T c2_i133;
  int32_T c2_i134;
  real_T c2_u[16];
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_i129 = 0;
  for (c2_i130 = 0; c2_i130 < 4; c2_i130++) {
    for (c2_i131 = 0; c2_i131 < 4; c2_i131++) {
      c2_b_inData[c2_i131 + c2_i129] = (*(real_T (*)[16])c2_inData)[c2_i131 +
        c2_i129];
    }

    c2_i129 += 4;
  }

  c2_i132 = 0;
  for (c2_i133 = 0; c2_i133 < 4; c2_i133++) {
    for (c2_i134 = 0; c2_i134 < 4; c2_i134++) {
      c2_u[c2_i134 + c2_i132] = c2_b_inData[c2_i134 + c2_i132];
    }

    c2_i132 += 4;
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 4, 4), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static void c2_i_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[16])
{
  real_T c2_dv10[16];
  int32_T c2_i135;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv10, 1, 0, 0U, 1, 0U, 2, 4, 4);
  for (c2_i135 = 0; c2_i135 < 16; c2_i135++) {
    c2_b_y[c2_i135] = c2_dv10[c2_i135];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_H0_n;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_b_y[16];
  int32_T c2_i136;
  int32_T c2_i137;
  int32_T c2_i138;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_H0_n = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_H0_n), &c2_thisId, c2_b_y);
  sf_mex_destroy(&c2_H0_n);
  c2_i136 = 0;
  for (c2_i137 = 0; c2_i137 < 4; c2_i137++) {
    for (c2_i138 = 0; c2_i138 < 4; c2_i138++) {
      (*(real_T (*)[16])c2_outData)[c2_i138 + c2_i136] = c2_b_y[c2_i138 +
        c2_i136];
    }

    c2_i136 += 4;
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static void c2_j_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[3])
{
  real_T c2_dv11[3];
  int32_T c2_i139;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv11, 1, 0, 0U, 1, 0U, 1, 3);
  for (c2_i139 = 0; c2_i139 < 3; c2_i139++) {
    c2_b_y[c2_i139] = c2_dv11[c2_i139];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_p0_n;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_b_y[3];
  int32_T c2_i140;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_p0_n = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_p0_n), &c2_thisId, c2_b_y);
  sf_mex_destroy(&c2_p0_n);
  for (c2_i140 = 0; c2_i140 < 3; c2_i140++) {
    (*(real_T (*)[3])c2_outData)[c2_i140] = c2_b_y[c2_i140];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_h_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i141;
  int32_T c2_i142;
  int32_T c2_i143;
  real_T c2_b_inData[9];
  int32_T c2_i144;
  int32_T c2_i145;
  int32_T c2_i146;
  real_T c2_u[9];
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_i141 = 0;
  for (c2_i142 = 0; c2_i142 < 3; c2_i142++) {
    for (c2_i143 = 0; c2_i143 < 3; c2_i143++) {
      c2_b_inData[c2_i143 + c2_i141] = (*(real_T (*)[9])c2_inData)[c2_i143 +
        c2_i141];
    }

    c2_i141 += 3;
  }

  c2_i144 = 0;
  for (c2_i145 = 0; c2_i145 < 3; c2_i145++) {
    for (c2_i146 = 0; c2_i146 < 3; c2_i146++) {
      c2_u[c2_i146 + c2_i144] = c2_b_inData[c2_i146 + c2_i144];
    }

    c2_i144 += 3;
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static void c2_k_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[9])
{
  real_T c2_dv12[9];
  int32_T c2_i147;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv12, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c2_i147 = 0; c2_i147 < 9; c2_i147++) {
    c2_b_y[c2_i147] = c2_dv12[c2_i147];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_R0_n;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_b_y[9];
  int32_T c2_i148;
  int32_T c2_i149;
  int32_T c2_i150;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_R0_n = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_R0_n), &c2_thisId, c2_b_y);
  sf_mex_destroy(&c2_R0_n);
  c2_i148 = 0;
  for (c2_i149 = 0; c2_i149 < 3; c2_i149++) {
    for (c2_i150 = 0; c2_i150 < 3; c2_i150++) {
      (*(real_T (*)[9])c2_outData)[c2_i150 + c2_i148] = c2_b_y[c2_i150 + c2_i148];
    }

    c2_i148 += 3;
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_i_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i151;
  real_T c2_b_inData[3];
  int32_T c2_i152;
  real_T c2_u[3];
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i151 = 0; c2_i151 < 3; c2_i151++) {
    c2_b_inData[c2_i151] = (*(real_T (*)[3])c2_inData)[c2_i151];
  }

  for (c2_i152 = 0; c2_i152 < 3; c2_i152++) {
    c2_u[c2_i152] = c2_b_inData[c2_i152];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 1, 3), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static void c2_l_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId,
  real_T c2_b_y[3])
{
  real_T c2_dv13[3];
  int32_T c2_i153;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv13, 1, 0, 0U, 1, 0U, 2, 1, 3);
  for (c2_i153 = 0; c2_i153 < 3; c2_i153++) {
    c2_b_y[c2_i153] = c2_dv13[c2_i153];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_rpy_in;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_b_y[3];
  int32_T c2_i154;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_rpy_in = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_rpy_in), &c2_thisId, c2_b_y);
  sf_mex_destroy(&c2_rpy_in);
  for (c2_i154 = 0; c2_i154 < 3; c2_i154++) {
    (*(real_T (*)[3])c2_outData)[c2_i154] = c2_b_y[c2_i154];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_j_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i155;
  char_T c2_b_inData[3];
  int32_T c2_i156;
  char_T c2_u[3];
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i155 = 0; c2_i155 < 3; c2_i155++) {
    c2_b_inData[c2_i155] = (*(char_T (*)[3])c2_inData)[c2_i155];
  }

  for (c2_i156 = 0; c2_i156 < 3; c2_i156++) {
    c2_u[c2_i156] = c2_b_inData[c2_i156];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 3),
                false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static const mxArray *c2_k_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  c2_s3TYGqYEZfhjEunCONSnN4C c2_u;
  const mxArray *c2_b_y = NULL;
  boolean_T c2_b_u;
  const mxArray *c2_c_y = NULL;
  boolean_T c2_c_u;
  const mxArray *c2_d_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(c2_s3TYGqYEZfhjEunCONSnN4C *)c2_inData;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_createstruct("structure", 2, 1, 1), false);
  c2_b_u = c2_u.zyx;
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", &c2_b_u, 11, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_b_y, c2_c_y, "zyx", "zyx", 0);
  c2_c_u = c2_u.deg;
  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", &c2_c_u, 11, 0U, 0U, 0U, 0), false);
  sf_mex_addfield(c2_b_y, c2_d_y, "deg", "deg", 0);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static c2_s3TYGqYEZfhjEunCONSnN4C c2_m_emlrt_marshallIn
  (SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance, const mxArray *c2_u,
   const emlrtMsgIdentifier *c2_parentId)
{
  c2_s3TYGqYEZfhjEunCONSnN4C c2_b_y;
  emlrtMsgIdentifier c2_thisId;
  static const char * c2_fieldNames[2] = { "zyx", "deg" };

  c2_thisId.fParent = c2_parentId;
  sf_mex_check_struct(c2_parentId, c2_u, 2, c2_fieldNames, 0U, NULL);
  c2_thisId.fIdentifier = "zyx";
  c2_b_y.zyx = c2_n_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "zyx", "zyx", 0)), &c2_thisId);
  c2_thisId.fIdentifier = "deg";
  c2_b_y.deg = c2_n_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c2_u, "deg", "deg", 0)), &c2_thisId);
  sf_mex_destroy(&c2_u);
  return c2_b_y;
}

static boolean_T c2_n_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  boolean_T c2_b_y;
  boolean_T c2_b2;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_b2, 1, 11, 0U, 0, 0U, 0);
  c2_b_y = c2_b2;
  sf_mex_destroy(&c2_u);
  return c2_b_y;
}

static void c2_j_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_opt;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  c2_s3TYGqYEZfhjEunCONSnN4C c2_b_y;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_opt = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_y = c2_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_opt), &c2_thisId);
  sf_mex_destroy(&c2_opt);
  *(c2_s3TYGqYEZfhjEunCONSnN4C *)c2_outData = c2_b_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_l_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i157;
  real_T c2_b_inData[2];
  int32_T c2_i158;
  real_T c2_u[2];
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i157 = 0; c2_i157 < 2; c2_i157++) {
    c2_b_inData[c2_i157] = (*(real_T (*)[2])c2_inData)[c2_i157];
  }

  for (c2_i158 = 0; c2_i158 < 2; c2_i158++) {
    c2_u[c2_i158] = c2_b_inData[c2_i158];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 1, 2), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static const mxArray *c2_m_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  boolean_T c2_u;
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(boolean_T *)c2_inData;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 11, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

const mxArray *sf_c2_IPILCO_relativeRPY_S_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo = NULL;
  c2_nameCaptureInfo = NULL;
  sf_mex_assign(&c2_nameCaptureInfo, sf_mex_createstruct("structure", 2, 57, 1),
                false);
  c2_info_helper(&c2_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c2_nameCaptureInfo);
  return c2_nameCaptureInfo;
}

static void c2_info_helper(const mxArray **c2_info)
{
  const mxArray *c2_rhs0 = NULL;
  const mxArray *c2_lhs0 = NULL;
  const mxArray *c2_rhs1 = NULL;
  const mxArray *c2_lhs1 = NULL;
  const mxArray *c2_rhs2 = NULL;
  const mxArray *c2_lhs2 = NULL;
  const mxArray *c2_rhs3 = NULL;
  const mxArray *c2_lhs3 = NULL;
  const mxArray *c2_rhs4 = NULL;
  const mxArray *c2_lhs4 = NULL;
  const mxArray *c2_rhs5 = NULL;
  const mxArray *c2_lhs5 = NULL;
  const mxArray *c2_rhs6 = NULL;
  const mxArray *c2_lhs6 = NULL;
  const mxArray *c2_rhs7 = NULL;
  const mxArray *c2_lhs7 = NULL;
  const mxArray *c2_rhs8 = NULL;
  const mxArray *c2_lhs8 = NULL;
  const mxArray *c2_rhs9 = NULL;
  const mxArray *c2_lhs9 = NULL;
  const mxArray *c2_rhs10 = NULL;
  const mxArray *c2_lhs10 = NULL;
  const mxArray *c2_rhs11 = NULL;
  const mxArray *c2_lhs11 = NULL;
  const mxArray *c2_rhs12 = NULL;
  const mxArray *c2_lhs12 = NULL;
  const mxArray *c2_rhs13 = NULL;
  const mxArray *c2_lhs13 = NULL;
  const mxArray *c2_rhs14 = NULL;
  const mxArray *c2_lhs14 = NULL;
  const mxArray *c2_rhs15 = NULL;
  const mxArray *c2_lhs15 = NULL;
  const mxArray *c2_rhs16 = NULL;
  const mxArray *c2_lhs16 = NULL;
  const mxArray *c2_rhs17 = NULL;
  const mxArray *c2_lhs17 = NULL;
  const mxArray *c2_rhs18 = NULL;
  const mxArray *c2_lhs18 = NULL;
  const mxArray *c2_rhs19 = NULL;
  const mxArray *c2_lhs19 = NULL;
  const mxArray *c2_rhs20 = NULL;
  const mxArray *c2_lhs20 = NULL;
  const mxArray *c2_rhs21 = NULL;
  const mxArray *c2_lhs21 = NULL;
  const mxArray *c2_rhs22 = NULL;
  const mxArray *c2_lhs22 = NULL;
  const mxArray *c2_rhs23 = NULL;
  const mxArray *c2_lhs23 = NULL;
  const mxArray *c2_rhs24 = NULL;
  const mxArray *c2_lhs24 = NULL;
  const mxArray *c2_rhs25 = NULL;
  const mxArray *c2_lhs25 = NULL;
  const mxArray *c2_rhs26 = NULL;
  const mxArray *c2_lhs26 = NULL;
  const mxArray *c2_rhs27 = NULL;
  const mxArray *c2_lhs27 = NULL;
  const mxArray *c2_rhs28 = NULL;
  const mxArray *c2_lhs28 = NULL;
  const mxArray *c2_rhs29 = NULL;
  const mxArray *c2_lhs29 = NULL;
  const mxArray *c2_rhs30 = NULL;
  const mxArray *c2_lhs30 = NULL;
  const mxArray *c2_rhs31 = NULL;
  const mxArray *c2_lhs31 = NULL;
  const mxArray *c2_rhs32 = NULL;
  const mxArray *c2_lhs32 = NULL;
  const mxArray *c2_rhs33 = NULL;
  const mxArray *c2_lhs33 = NULL;
  const mxArray *c2_rhs34 = NULL;
  const mxArray *c2_lhs34 = NULL;
  const mxArray *c2_rhs35 = NULL;
  const mxArray *c2_lhs35 = NULL;
  const mxArray *c2_rhs36 = NULL;
  const mxArray *c2_lhs36 = NULL;
  const mxArray *c2_rhs37 = NULL;
  const mxArray *c2_lhs37 = NULL;
  const mxArray *c2_rhs38 = NULL;
  const mxArray *c2_lhs38 = NULL;
  const mxArray *c2_rhs39 = NULL;
  const mxArray *c2_lhs39 = NULL;
  const mxArray *c2_rhs40 = NULL;
  const mxArray *c2_lhs40 = NULL;
  const mxArray *c2_rhs41 = NULL;
  const mxArray *c2_lhs41 = NULL;
  const mxArray *c2_rhs42 = NULL;
  const mxArray *c2_lhs42 = NULL;
  const mxArray *c2_rhs43 = NULL;
  const mxArray *c2_lhs43 = NULL;
  const mxArray *c2_rhs44 = NULL;
  const mxArray *c2_lhs44 = NULL;
  const mxArray *c2_rhs45 = NULL;
  const mxArray *c2_lhs45 = NULL;
  const mxArray *c2_rhs46 = NULL;
  const mxArray *c2_lhs46 = NULL;
  const mxArray *c2_rhs47 = NULL;
  const mxArray *c2_lhs47 = NULL;
  const mxArray *c2_rhs48 = NULL;
  const mxArray *c2_lhs48 = NULL;
  const mxArray *c2_rhs49 = NULL;
  const mxArray *c2_lhs49 = NULL;
  const mxArray *c2_rhs50 = NULL;
  const mxArray *c2_lhs50 = NULL;
  const mxArray *c2_rhs51 = NULL;
  const mxArray *c2_lhs51 = NULL;
  const mxArray *c2_rhs52 = NULL;
  const mxArray *c2_lhs52 = NULL;
  const mxArray *c2_rhs53 = NULL;
  const mxArray *c2_lhs53 = NULL;
  const mxArray *c2_rhs54 = NULL;
  const mxArray *c2_lhs54 = NULL;
  const mxArray *c2_rhs55 = NULL;
  const mxArray *c2_lhs55 = NULL;
  const mxArray *c2_rhs56 = NULL;
  const mxArray *c2_lhs56 = NULL;
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("my_rpy2r"), "name", "name", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_rpy2r.m"), "resolved", "resolved", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464800136U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c2_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_rpy2r.m"), "context", "context", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("numcols"), "name", "name", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/numcols.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793908U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c2_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_rpy2r.m"), "context", "context", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("mrdivide"), "name", "name", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807648U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1370009886U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c2_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c2_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("rdivide"), "name", "name", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363713880U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c2_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c2_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286818796U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c2_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_div"), "name", "name", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1386423952U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c2_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c2_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_rpy2r.m"), "context", "context", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("numrows"), "name", "name", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/numrows.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793908U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c2_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_rpy2r.m"), "context", "context", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("rotx"), "name", "name", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/rotx.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793909U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c2_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/rotx.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("cos"), "name", "name", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c2_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286818722U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c2_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/rotx.m"), "context",
                  "context", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("sin"), "name", "name", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395328504U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c2_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286818736U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c2_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_rpy2r.m"), "context", "context", 15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("roty"), "name", "name", 15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/roty.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793909U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c2_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/roty.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("cos"), "name", "name", 16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c2_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/roty.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("sin"), "name", "name", 17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395328504U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c2_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_rpy2r.m"), "context", "context", 18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c2_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c2_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c2_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c2_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c2_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 23);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  23);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1375980690U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c2_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 24);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 24);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c2_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 25);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c2_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 26);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 26);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c2_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 27);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 27);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c2_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 28);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 28);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c2_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 29);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 29);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c2_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 30);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 30);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c2_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_rpy2r.m"), "context", "context", 31);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("rotz"), "name", "name", 31);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/rotz.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793909U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c2_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/rotz.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("cos"), "name", "name", 32);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c2_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/rotz.m"), "context",
                  "context", 33);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("sin"), "name", "name", 33);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395328504U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c2_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 34);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("rt2tr"), "name", "name", 34);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/rt2tr.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793909U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c2_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/rt2tr.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("numcols"), "name", "name", 35);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/numcols.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793908U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c2_rhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs35, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/rt2tr.m"),
                  "context", "context", 36);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("numrows"), "name", "name", 36);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/numrows.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793908U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c2_rhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs36, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "context", "context", 37);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("wtrans"), "name", "name", 37);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/wtrans.m"),
                  "resolved", "resolved", 37);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793910U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c2_rhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs37, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/wtrans.m"),
                  "context", "context", 38);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("transl"), "name", "name", 38);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/transl.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793910U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c2_rhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs38, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/transl.m"),
                  "context", "context", 39);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("ishomog"), "name", "name", 39);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/ishomog.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793908U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c2_rhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs39, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/ishomog.m"),
                  "context", "context", 40);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("all"), "name", "name", 40);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "resolved",
                  "resolved", 40);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1372582414U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c2_rhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs40, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "context", "context",
                  41);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 41);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c2_rhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs41, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "context", "context",
                  42);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 42);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c2_rhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs42, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "context", "context",
                  43);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.allOrAny"),
                  "name", "name", 43);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1372583158U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c2_rhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs43, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 44);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 44);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c2_rhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs44, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 45);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("isequal"), "name", "name", 45);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 45);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286818758U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c2_rhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs45, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "context",
                  "context", 46);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_isequal_core"), "name",
                  "name", 46);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"),
                  "resolved", "resolved", 46);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286818786U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c2_rhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs46, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs46), "lhs", "lhs",
                  46);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 47);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.constNonSingletonDim"), "name", "name", 47);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/constNonSingletonDim.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 47);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 47);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 47);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 47);
  sf_mex_assign(&c2_rhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs47, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs47), "rhs", "rhs",
                  47);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs47), "lhs", "lhs",
                  47);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/wtrans.m"),
                  "context", "context", 48);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("cross"), "name", "name", 48);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/specfun/cross.m"), "resolved",
                  "resolved", 48);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1286818842U), "fileTimeLo",
                  "fileTimeLo", 48);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 48);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 48);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 48);
  sf_mex_assign(&c2_rhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs48, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs48), "rhs", "rhs",
                  48);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs48), "lhs", "lhs",
                  48);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/wtrans.m"),
                  "context", "context", 49);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("t2r"), "name", "name", 49);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/t2r.m"), "resolved",
                  "resolved", 49);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1464793910U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 49);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 49);
  sf_mex_assign(&c2_rhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs49, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs49), "rhs", "rhs",
                  49);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs49), "lhs", "lhs",
                  49);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/t2r.m"), "context",
                  "context", 50);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("any"), "name", "name", 50);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m"), "resolved",
                  "resolved", 50);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1372582416U), "fileTimeLo",
                  "fileTimeLo", 50);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 50);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 50);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 50);
  sf_mex_assign(&c2_rhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs50, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs50), "rhs", "rhs",
                  50);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs50), "lhs", "lhs",
                  50);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m"), "context", "context",
                  51);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 51);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 51);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 51);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 51);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 51);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 51);
  sf_mex_assign(&c2_rhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs51, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs51), "rhs", "rhs",
                  51);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs51), "lhs", "lhs",
                  51);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m"), "context", "context",
                  52);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 52);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 52);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c2_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs52), "lhs", "lhs",
                  52);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/any.m"), "context", "context",
                  53);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("coder.internal.allOrAny"),
                  "name", "name", 53);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 53);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "resolved", "resolved", 53);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1372583158U), "fileTimeLo",
                  "fileTimeLo", 53);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 53);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 53);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 53);
  sf_mex_assign(&c2_rhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs53, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs53), "rhs", "rhs",
                  53);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs53), "lhs", "lhs",
                  53);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 54);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("isnan"), "name", "name", 54);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 54);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "resolved",
                  "resolved", 54);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1363713858U), "fileTimeLo",
                  "fileTimeLo", 54);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 54);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 54);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 54);
  sf_mex_assign(&c2_rhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs54, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs54), "rhs", "rhs",
                  54);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs54), "lhs", "lhs",
                  54);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m"), "context",
                  "context", 55);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 55);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 55);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 55);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 55);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 55);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 55);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 55);
  sf_mex_assign(&c2_rhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs55, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs55), "rhs", "rhs",
                  55);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs55), "lhs", "lhs",
                  55);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/wtrans.m"),
                  "context", "context", 56);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 56);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 56);
  sf_mex_addfield(*c2_info, c2_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 56);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 56);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 56);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 56);
  sf_mex_addfield(*c2_info, c2_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 56);
  sf_mex_assign(&c2_rhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c2_lhs56, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_rhs56), "rhs", "rhs",
                  56);
  sf_mex_addfield(*c2_info, sf_mex_duplicatearraysafe(&c2_lhs56), "lhs", "lhs",
                  56);
  sf_mex_destroy(&c2_rhs0);
  sf_mex_destroy(&c2_lhs0);
  sf_mex_destroy(&c2_rhs1);
  sf_mex_destroy(&c2_lhs1);
  sf_mex_destroy(&c2_rhs2);
  sf_mex_destroy(&c2_lhs2);
  sf_mex_destroy(&c2_rhs3);
  sf_mex_destroy(&c2_lhs3);
  sf_mex_destroy(&c2_rhs4);
  sf_mex_destroy(&c2_lhs4);
  sf_mex_destroy(&c2_rhs5);
  sf_mex_destroy(&c2_lhs5);
  sf_mex_destroy(&c2_rhs6);
  sf_mex_destroy(&c2_lhs6);
  sf_mex_destroy(&c2_rhs7);
  sf_mex_destroy(&c2_lhs7);
  sf_mex_destroy(&c2_rhs8);
  sf_mex_destroy(&c2_lhs8);
  sf_mex_destroy(&c2_rhs9);
  sf_mex_destroy(&c2_lhs9);
  sf_mex_destroy(&c2_rhs10);
  sf_mex_destroy(&c2_lhs10);
  sf_mex_destroy(&c2_rhs11);
  sf_mex_destroy(&c2_lhs11);
  sf_mex_destroy(&c2_rhs12);
  sf_mex_destroy(&c2_lhs12);
  sf_mex_destroy(&c2_rhs13);
  sf_mex_destroy(&c2_lhs13);
  sf_mex_destroy(&c2_rhs14);
  sf_mex_destroy(&c2_lhs14);
  sf_mex_destroy(&c2_rhs15);
  sf_mex_destroy(&c2_lhs15);
  sf_mex_destroy(&c2_rhs16);
  sf_mex_destroy(&c2_lhs16);
  sf_mex_destroy(&c2_rhs17);
  sf_mex_destroy(&c2_lhs17);
  sf_mex_destroy(&c2_rhs18);
  sf_mex_destroy(&c2_lhs18);
  sf_mex_destroy(&c2_rhs19);
  sf_mex_destroy(&c2_lhs19);
  sf_mex_destroy(&c2_rhs20);
  sf_mex_destroy(&c2_lhs20);
  sf_mex_destroy(&c2_rhs21);
  sf_mex_destroy(&c2_lhs21);
  sf_mex_destroy(&c2_rhs22);
  sf_mex_destroy(&c2_lhs22);
  sf_mex_destroy(&c2_rhs23);
  sf_mex_destroy(&c2_lhs23);
  sf_mex_destroy(&c2_rhs24);
  sf_mex_destroy(&c2_lhs24);
  sf_mex_destroy(&c2_rhs25);
  sf_mex_destroy(&c2_lhs25);
  sf_mex_destroy(&c2_rhs26);
  sf_mex_destroy(&c2_lhs26);
  sf_mex_destroy(&c2_rhs27);
  sf_mex_destroy(&c2_lhs27);
  sf_mex_destroy(&c2_rhs28);
  sf_mex_destroy(&c2_lhs28);
  sf_mex_destroy(&c2_rhs29);
  sf_mex_destroy(&c2_lhs29);
  sf_mex_destroy(&c2_rhs30);
  sf_mex_destroy(&c2_lhs30);
  sf_mex_destroy(&c2_rhs31);
  sf_mex_destroy(&c2_lhs31);
  sf_mex_destroy(&c2_rhs32);
  sf_mex_destroy(&c2_lhs32);
  sf_mex_destroy(&c2_rhs33);
  sf_mex_destroy(&c2_lhs33);
  sf_mex_destroy(&c2_rhs34);
  sf_mex_destroy(&c2_lhs34);
  sf_mex_destroy(&c2_rhs35);
  sf_mex_destroy(&c2_lhs35);
  sf_mex_destroy(&c2_rhs36);
  sf_mex_destroy(&c2_lhs36);
  sf_mex_destroy(&c2_rhs37);
  sf_mex_destroy(&c2_lhs37);
  sf_mex_destroy(&c2_rhs38);
  sf_mex_destroy(&c2_lhs38);
  sf_mex_destroy(&c2_rhs39);
  sf_mex_destroy(&c2_lhs39);
  sf_mex_destroy(&c2_rhs40);
  sf_mex_destroy(&c2_lhs40);
  sf_mex_destroy(&c2_rhs41);
  sf_mex_destroy(&c2_lhs41);
  sf_mex_destroy(&c2_rhs42);
  sf_mex_destroy(&c2_lhs42);
  sf_mex_destroy(&c2_rhs43);
  sf_mex_destroy(&c2_lhs43);
  sf_mex_destroy(&c2_rhs44);
  sf_mex_destroy(&c2_lhs44);
  sf_mex_destroy(&c2_rhs45);
  sf_mex_destroy(&c2_lhs45);
  sf_mex_destroy(&c2_rhs46);
  sf_mex_destroy(&c2_lhs46);
  sf_mex_destroy(&c2_rhs47);
  sf_mex_destroy(&c2_lhs47);
  sf_mex_destroy(&c2_rhs48);
  sf_mex_destroy(&c2_lhs48);
  sf_mex_destroy(&c2_rhs49);
  sf_mex_destroy(&c2_lhs49);
  sf_mex_destroy(&c2_rhs50);
  sf_mex_destroy(&c2_lhs50);
  sf_mex_destroy(&c2_rhs51);
  sf_mex_destroy(&c2_lhs51);
  sf_mex_destroy(&c2_rhs52);
  sf_mex_destroy(&c2_lhs52);
  sf_mex_destroy(&c2_rhs53);
  sf_mex_destroy(&c2_lhs53);
  sf_mex_destroy(&c2_rhs54);
  sf_mex_destroy(&c2_lhs54);
  sf_mex_destroy(&c2_rhs55);
  sf_mex_destroy(&c2_lhs55);
  sf_mex_destroy(&c2_rhs56);
  sf_mex_destroy(&c2_lhs56);
}

static const mxArray *c2_emlrt_marshallOut(const char * c2_u)
{
  const mxArray *c2_b_y = NULL;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c2_u)), false);
  return c2_b_y;
}

static const mxArray *c2_b_emlrt_marshallOut(const uint32_T c2_u)
{
  const mxArray *c2_b_y = NULL;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 7, 0U, 0U, 0U, 0), false);
  return c2_b_y;
}

static void c2_eml_scalar_eg(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c2_threshold(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c2_cross(SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance,
                     real_T c2_a[3], real_T c2_b[3], real_T c2_c[3])
{
  real_T c2_c1;
  real_T c2_c2;
  real_T c2_c3;
  (void)chartInstance;
  c2_c1 = c2_a[1] * c2_b[2] - c2_a[2] * c2_b[1];
  c2_c2 = c2_a[2] * c2_b[0] - c2_a[0] * c2_b[2];
  c2_c3 = c2_a[0] * c2_b[1] - c2_a[1] * c2_b[0];
  c2_c[0] = c2_c1;
  c2_c[1] = c2_c2;
  c2_c[2] = c2_c3;
}

static void c2_b_eml_scalar_eg(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *c2_n_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_u;
  const mxArray *c2_b_y = NULL;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(int32_T *)c2_inData;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c2_mxArrayOutData, c2_b_y, false);
  return c2_mxArrayOutData;
}

static int32_T c2_o_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  int32_T c2_b_y;
  int32_T c2_i159;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_i159, 1, 6, 0U, 0, 0U, 0);
  c2_b_y = c2_i159;
  sf_mex_destroy(&c2_u);
  return c2_b_y;
}

static void c2_k_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_sfEvent;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  int32_T c2_b_y;
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)chartInstanceVoid;
  c2_b_sfEvent = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_y = c2_o_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_sfEvent),
    &c2_thisId);
  sf_mex_destroy(&c2_b_sfEvent);
  *(int32_T *)c2_outData = c2_b_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static uint8_T c2_p_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_b_is_active_c2_IPILCO_relativeRPY_S, const
  char_T *c2_identifier)
{
  uint8_T c2_b_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_y = c2_q_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c2_b_is_active_c2_IPILCO_relativeRPY_S), &c2_thisId);
  sf_mex_destroy(&c2_b_is_active_c2_IPILCO_relativeRPY_S);
  return c2_b_y;
}

static uint8_T c2_q_emlrt_marshallIn(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  uint8_T c2_b_y;
  uint8_T c2_u0;
  (void)chartInstance;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_u0, 1, 3, 0U, 0, 0U, 0);
  c2_b_y = c2_u0;
  sf_mex_destroy(&c2_u);
  return c2_b_y;
}

static void init_dsm_address_info(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc2_IPILCO_relativeRPY_SInstanceStruct
  *chartInstance)
{
  chartInstance->c2_Kp_env = (real_T (*)[6])ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c2_y = (real_T (*)[12])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c2_Kd_env = (real_T (*)[6])ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c2_xe = (real_T (*)[6])ssGetInputPortSignal_wrapper
    (chartInstance->S, 2);
  chartInstance->c2_dxe = (real_T (*)[6])ssGetInputPortSignal_wrapper
    (chartInstance->S, 3);
  chartInstance->c2_xc = (real_T (*)[6])ssGetInputPortSignal_wrapper
    (chartInstance->S, 4);
  chartInstance->c2_x_hole = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 5);
  chartInstance->c2_unused = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 6);
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

void sf_c2_IPILCO_relativeRPY_S_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(958960771U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(120073744U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2035977549U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3325532959U);
}

mxArray* sf_c2_IPILCO_relativeRPY_S_get_post_codegen_info(void);
mxArray *sf_c2_IPILCO_relativeRPY_S_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("I3D4SNwI0xGvYKzuKomrIH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,7,3,dataFields);

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
      pr[0] = (double)(6);
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
      pr[0] = (double)(3);
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
      pr[0] = (double)(1);
      pr[1] = (double)(1);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(12);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxArray* mxPostCodegenInfo =
      sf_c2_IPILCO_relativeRPY_S_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c2_IPILCO_relativeRPY_S_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c2_IPILCO_relativeRPY_S_jit_fallback_info(void)
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

mxArray *sf_c2_IPILCO_relativeRPY_S_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c2_IPILCO_relativeRPY_S_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c2_IPILCO_relativeRPY_S(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x4'type','srcId','name','auxInfo'{{M[1],M[5],T\"y\",},{M[4],M[0],T\"INSERTED\",S'l','i','p'{{M1x2[388 396],M[0],}}},{M[4],M[0],T\"inHole\",S'l','i','p'{{M1x2[381 387],M[0],}}},{M[8],M[0],T\"is_active_c2_IPILCO_relativeRPY_S\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 4, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c2_IPILCO_relativeRPY_S_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _IPILCO_relativeRPY_SMachineNumber_,
           2,
           1,
           1,
           0,
           8,
           0,
           0,
           0,
           0,
           11,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"Kp_env");
          _SFD_SET_DATA_PROPS(1,2,0,1,"y");
          _SFD_SET_DATA_PROPS(2,1,1,0,"Kd_env");
          _SFD_SET_DATA_PROPS(3,1,1,0,"xe");
          _SFD_SET_DATA_PROPS(4,1,1,0,"dxe");
          _SFD_SET_DATA_PROPS(5,1,1,0,"xc");
          _SFD_SET_DATA_PROPS(6,1,1,0,"x_hole");
          _SFD_SET_DATA_PROPS(7,1,1,0,"unused");
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
        _SFD_CV_INIT_EML(0,1,1,16,0,0,0,0,0,9,4);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,3697);
        _SFD_CV_INIT_EML_IF(0,1,0,398,416,-1,436);
        _SFD_CV_INIT_EML_IF(0,1,1,437,457,-1,483);
        _SFD_CV_INIT_EML_IF(0,1,2,1321,1337,-1,1509);
        _SFD_CV_INIT_EML_IF(0,1,3,1359,1396,-1,1505);
        _SFD_CV_INIT_EML_IF(0,1,4,1511,1527,1605,1626);
        _SFD_CV_INIT_EML_IF(0,1,5,1605,1626,-1,1626);
        _SFD_CV_INIT_EML_IF(0,1,6,1692,1713,2250,2274);
        _SFD_CV_INIT_EML_IF(0,1,7,1818,1831,-1,1912);
        _SFD_CV_INIT_EML_IF(0,1,8,2250,2274,-1,2274);
        _SFD_CV_INIT_EML_IF(0,1,9,2326,2353,-1,2605);
        _SFD_CV_INIT_EML_IF(0,1,10,2495,2508,-1,2597);
        _SFD_CV_INIT_EML_IF(0,1,11,2615,2681,-1,2790);
        _SFD_CV_INIT_EML_IF(0,1,12,2822,2840,3213,3235);
        _SFD_CV_INIT_EML_IF(0,1,13,3002,3015,-1,3104);
        _SFD_CV_INIT_EML_IF(0,1,14,3213,3235,-1,3235);
        _SFD_CV_INIT_EML_IF(0,1,15,3362,3375,-1,3464);

        {
          static int condStart[] = { 1362, 1381 };

          static int condEnd[] = { 1377, 1396 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,0,1362,1396,2,0,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 1695, 1707 };

          static int condEnd[] = { 1702, 1713 };

          static int pfixExpr[] = { 0, 1, -1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,1,1695,1713,2,2,&(condStart[0]),&(condEnd[0]),
                                4,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 2257, 2268 };

          static int condEnd[] = { 2264, 2274 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,2,2257,2274,2,4,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          static int condStart[] = { 2618, 2645, 2673 };

          static int condEnd[] = { 2641, 2668, 2681 };

          static int pfixExpr[] = { 0, 1, -3, 2, -1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,3,2618,2681,3,6,&(condStart[0]),&(condEnd[0]),
                                6,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_EML_RELATIONAL(0,1,0,1324,1337,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,1,1362,1377,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,2,1381,1396,-1,2);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,3,1514,1527,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,4,1612,1626,-1,2);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,5,1821,1831,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,6,2329,2353,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,7,2498,2508,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,8,2618,2641,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,9,2645,2668,-1,2);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,10,2825,2840,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,11,3005,3015,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,12,3220,3235,-1,2);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,13,3365,3375,-1,2);
        _SFD_CV_INIT_SCRIPT(0,1,6,0,0,0,2,0,1,1);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"my_rpy2r",2028,-1,3220);
        _SFD_CV_INIT_SCRIPT_IF(0,0,2158,2179,-1,2211);
        _SFD_CV_INIT_SCRIPT_IF(0,1,2248,2271,2342,2405);
        _SFD_CV_INIT_SCRIPT_IF(0,2,2449,2459,-1,2572);
        _SFD_CV_INIT_SCRIPT_IF(0,3,2578,2589,2889,3219);
        _SFD_CV_INIT_SCRIPT_IF(0,4,2618,2639,2702,2884);
        _SFD_CV_INIT_SCRIPT_IF(0,5,2945,2966,3029,3211);
        _SFD_CV_INIT_SCRIPT_FOR(0,0,2761,2783,2872);
        _SFD_CV_INIT_SCRIPT_FOR(0,1,3088,3110,3199);

        {
          static int condStart[] = { 2582 };

          static int condEnd[] = { 2589 };

          static int pfixExpr[] = { 0, -1 };

          _SFD_CV_INIT_SCRIPT_MCDC(0,0,2581,2589,1,0,&(condStart[0]),&(condEnd[0]),
            2,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(0,0,2251,2271,-1,0);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(0,1,2621,2639,-1,0);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(0,2,2948,2966,-1,0);
        _SFD_CV_INIT_SCRIPT(1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(1,0,"numcols",951,-1,991);
        _SFD_CV_INIT_SCRIPT(2,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(2,0,"numrows",946,-1,988);
        _SFD_CV_INIT_SCRIPT(3,1,1,0,0,0,0,0,2,1);
        _SFD_CV_INIT_SCRIPT_FCN(3,0,"rotx",1021,-1,1238);
        _SFD_CV_INIT_SCRIPT_IF(3,0,1052,1087,-1,1118);

        {
          static int condStart[] = { 1055, 1069 };

          static int condEnd[] = { 1065, 1087 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_SCRIPT_MCDC(3,0,1055,1087,2,0,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(3,0,1055,1065,-1,4);
        _SFD_CV_INIT_SCRIPT(4,1,1,0,0,0,0,0,2,1);
        _SFD_CV_INIT_SCRIPT_FCN(4,0,"roty",1021,-1,1228);
        _SFD_CV_INIT_SCRIPT_IF(4,0,1051,1086,-1,1117);

        {
          static int condStart[] = { 1054, 1068 };

          static int condEnd[] = { 1064, 1086 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_SCRIPT_MCDC(4,0,1054,1086,2,0,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(4,0,1054,1064,-1,4);
        _SFD_CV_INIT_SCRIPT(5,1,1,0,0,0,0,0,2,1);
        _SFD_CV_INIT_SCRIPT_FCN(5,0,"rotz",1021,-1,1235);
        _SFD_CV_INIT_SCRIPT_IF(5,0,1051,1086,-1,1117);

        {
          static int condStart[] = { 1054, 1068 };

          static int condEnd[] = { 1064, 1086 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_SCRIPT_MCDC(5,0,1054,1086,2,0,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(5,0,1054,1064,-1,4);
        _SFD_CV_INIT_SCRIPT(6,1,4,0,0,0,1,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(6,0,"rt2tr",1294,-1,1848);
        _SFD_CV_INIT_SCRIPT_IF(6,0,1323,1350,-1,1393);
        _SFD_CV_INIT_SCRIPT_IF(6,1,1398,1425,-1,1493);
        _SFD_CV_INIT_SCRIPT_IF(6,2,1499,1525,-1,1595);
        _SFD_CV_INIT_SCRIPT_IF(6,3,1601,1617,1792,1846);
        _SFD_CV_INIT_SCRIPT_FOR(6,0,1713,1731,1787);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,0,1326,1350,-1,1);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,1,1401,1425,-1,1);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,2,1502,1525,-1,1);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(6,3,1604,1617,-1,4);
        _SFD_CV_INIT_SCRIPT(7,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(7,0,"wtrans",1087,-1,1242);
        _SFD_CV_INIT_SCRIPT(8,1,9,0,0,0,0,0,2,1);
        _SFD_CV_INIT_SCRIPT_FCN(8,0,"transl",1953,-1,3159);
        _SFD_CV_INIT_SCRIPT_IF(8,0,1995,2009,3045,3063);
        _SFD_CV_INIT_SCRIPT_IF(8,1,2018,2031,2716,3036);
        _SFD_CV_INIT_SCRIPT_IF(8,2,2044,2060,2394,2707);
        _SFD_CV_INIT_SCRIPT_IF(8,3,2128,2143,2207,2226);
        _SFD_CV_INIT_SCRIPT_IF(8,4,2207,2226,-1,2226);
        _SFD_CV_INIT_SCRIPT_IF(8,5,2449,2480,2532,2551);
        _SFD_CV_INIT_SCRIPT_IF(8,6,2532,2551,-1,2551);
        _SFD_CV_INIT_SCRIPT_IF(8,7,2716,2737,2872,3036);
        _SFD_CV_INIT_SCRIPT_IF(8,8,3045,3063,-1,3063);

        {
          static int condStart[] = { 2452, 2468 };

          static int condEnd[] = { 2464, 2480 };

          static int pfixExpr[] = { 0, 1, -2 };

          _SFD_CV_INIT_SCRIPT_MCDC(8,0,2452,2480,2,0,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(8,0,1998,2009,-1,0);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(8,1,2047,2060,-1,0);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(8,4,2452,2464,-1,0);
        _SFD_CV_INIT_SCRIPT(9,1,2,0,0,0,0,0,2,1);
        _SFD_CV_INIT_SCRIPT_FCN(9,0,"ishomog",1166,-1,1399);
        _SFD_CV_INIT_SCRIPT_IF(9,0,1220,1237,1367,1398);
        _SFD_CV_INIT_SCRIPT_IF(9,1,1282,1300,-1,1361);

        {
          static int condStart[] = { 1285, 1290 };

          static int condEnd[] = { 1286, 1300 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_SCRIPT_MCDC(9,0,1285,1300,2,0,&(condStart[0]),&(condEnd[0]),
            3,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(9,0,1223,1237,-1,5);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(9,1,1290,1300,-1,4);
        _SFD_CV_INIT_SCRIPT(10,1,3,0,0,0,1,0,1,1);
        _SFD_CV_INIT_SCRIPT_FCN(10,0,"t2r",1170,-1,1743);
        _SFD_CV_INIT_SCRIPT_IF(10,0,1258,1273,-1,1339);
        _SFD_CV_INIT_SCRIPT_IF(10,1,1344,1366,-1,1461);
        _SFD_CV_INIT_SCRIPT_IF(10,2,1517,1533,1595,1742);
        _SFD_CV_INIT_SCRIPT_FOR(10,0,1669,1682,1734);

        {
          static int condStart[] = { 1348 };

          static int condEnd[] = { 1366 };

          static int pfixExpr[] = { 0, -1 };

          _SFD_CV_INIT_SCRIPT_MCDC(10,0,1347,1366,1,0,&(condStart[0]),&(condEnd
            [0]),2,&(pfixExpr[0]));
        }

        _SFD_CV_INIT_SCRIPT_RELATIONAL(10,0,1261,1273,-1,1);
        _SFD_CV_INIT_SCRIPT_RELATIONAL(10,1,1520,1533,-1,0);

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_f_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_c_sf_marshallOut,(MexInFcnForType)
            c2_c_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_f_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_f_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_f_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_f_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_d_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_VALUE_PTR(0U, *chartInstance->c2_Kp_env);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c2_y);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c2_Kd_env);
        _SFD_SET_DATA_VALUE_PTR(3U, *chartInstance->c2_xe);
        _SFD_SET_DATA_VALUE_PTR(4U, *chartInstance->c2_dxe);
        _SFD_SET_DATA_VALUE_PTR(5U, *chartInstance->c2_xc);
        _SFD_SET_DATA_VALUE_PTR(6U, *chartInstance->c2_x_hole);
        _SFD_SET_DATA_VALUE_PTR(7U, chartInstance->c2_unused);
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
  return "lWuscg8unauiBVSYd9Mb8G";
}

static void sf_opaque_initialize_c2_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc2_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c2_IPILCO_relativeRPY_S
    ((SFc2_IPILCO_relativeRPY_SInstanceStruct*) chartInstanceVar);
  initialize_c2_IPILCO_relativeRPY_S((SFc2_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c2_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  enable_c2_IPILCO_relativeRPY_S((SFc2_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_disable_c2_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  disable_c2_IPILCO_relativeRPY_S((SFc2_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c2_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  sf_gateway_c2_IPILCO_relativeRPY_S((SFc2_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c2_IPILCO_relativeRPY_S(SimStruct*
  S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c2_IPILCO_relativeRPY_S
    ((SFc2_IPILCO_relativeRPY_SInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
}

static void sf_opaque_set_sim_state_c2_IPILCO_relativeRPY_S(SimStruct* S, const
  mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c2_IPILCO_relativeRPY_S((SFc2_IPILCO_relativeRPY_SInstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c2_IPILCO_relativeRPY_S(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc2_IPILCO_relativeRPY_SInstanceStruct*) chartInstanceVar
      )->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_IPILCO_relativeRPY_S_optimization_info();
    }

    finalize_c2_IPILCO_relativeRPY_S((SFc2_IPILCO_relativeRPY_SInstanceStruct*)
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
  initSimStructsc2_IPILCO_relativeRPY_S((SFc2_IPILCO_relativeRPY_SInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_IPILCO_relativeRPY_S(SimStruct *S)
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
    initialize_params_c2_IPILCO_relativeRPY_S
      ((SFc2_IPILCO_relativeRPY_SInstanceStruct*)(chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c2_IPILCO_relativeRPY_S(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_IPILCO_relativeRPY_S_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,2,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,2,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,2);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,2,7);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,2,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 7; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,2);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(632146861U));
  ssSetChecksum1(S,(1097249262U));
  ssSetChecksum2(S,(2646857707U));
  ssSetChecksum3(S,(802315861U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c2_IPILCO_relativeRPY_S(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_IPILCO_relativeRPY_S(SimStruct *S)
{
  SFc2_IPILCO_relativeRPY_SInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc2_IPILCO_relativeRPY_SInstanceStruct *)utMalloc(sizeof
    (SFc2_IPILCO_relativeRPY_SInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc2_IPILCO_relativeRPY_SInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.enableChart =
    sf_opaque_enable_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.mdlStart = mdlStart_c2_IPILCO_relativeRPY_S;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c2_IPILCO_relativeRPY_S;
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

void c2_IPILCO_relativeRPY_S_method_dispatcher(SimStruct *S, int_T method, void *
  data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_IPILCO_relativeRPY_S(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_IPILCO_relativeRPY_S(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_IPILCO_relativeRPY_S(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_IPILCO_relativeRPY_S_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
