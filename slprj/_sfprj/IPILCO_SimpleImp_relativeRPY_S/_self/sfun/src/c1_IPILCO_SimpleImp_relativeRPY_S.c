/* Include files */

#include <stddef.h>
#include "blas.h"
#include "IPILCO_SimpleImp_relativeRPY_S_sfun.h"
#include "c1_IPILCO_SimpleImp_relativeRPY_S.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "IPILCO_SimpleImp_relativeRPY_S_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c1_debug_family_names[27] = { "R_e", "R_d", "p_e", "p_d",
  "delta_p", "dp_e", "prevp_d", "ddelta_p", "a_t", "Re_d", "dEUL_de", "a_o",
  "nargin", "nargout", "Kp", "Kd", "He", "xe", "dxe_n", "Hd", "prevHd",
  "prevdp_d", "dt", "prevEUL_de", "a", "dp_d", "EUL_de" };

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

/* Function Declarations */
static void initialize_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void initialize_params_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void enable_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void disable_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void set_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_st);
static void finalize_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void sf_gateway_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void mdl_start_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void c1_chartstep_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void initSimStructsc1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void c1_tr2rt(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
                     *chartInstance, real_T c1_T[16], real_T c1_R[9]);
static void c1_transl(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
                      *chartInstance, real_T c1_x[16], real_T c1_t1[3]);
static void c1_my_tr2rpy(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_m[9], real_T c1_rpy[3]);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber, uint32_T c1_instanceNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static void c1_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_b_EUL_de, const char_T *c1_identifier, real_T c1_y[3]);
static void c1_b_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3]);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_c_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_b_dp_d, const char_T *c1_identifier, real_T c1_y[3]);
static void c1_d_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3]);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_e_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_b_a, const char_T *c1_identifier, real_T c1_y[6]);
static void c1_f_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[6]);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static real_T c1_g_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_g_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_h_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[9]);
static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_i_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[16]);
static void c1_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_h_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_i_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_j_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_j_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3]);
static void c1_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_k_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_l_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static c1_snY5FzvLGhd6UMWwSEwxMX c1_k_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static boolean_T c1_l_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_info_helper(const mxArray **c1_info);
static const mxArray *c1_emlrt_marshallOut(const char * c1_u);
static const mxArray *c1_b_emlrt_marshallOut(const uint32_T c1_u);
static void c1_rdivide(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_x[3], real_T c1_y, real_T c1_z[3]);
static void c1_eml_scalar_eg(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance);
static void c1_threshold(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance);
static void c1_b_eml_scalar_eg(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance);
static void c1_eps(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
                   *chartInstance);
static real_T c1_atan2(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_y, real_T c1_x);
static const mxArray *c1_m_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_m_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_n_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_b_is_active_c1_IPILCO_SimpleImp_relativeRPY_S, const char_T
   *c1_identifier);
static uint8_T c1_o_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void init_dsm_address_info
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);
static void init_simulink_io_address
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c1_is_active_c1_IPILCO_SimpleImp_relativeRPY_S = 0U;
}

static void initialize_params_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c1_update_debugger_state_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  int32_T c1_i0;
  real_T c1_u[3];
  const mxArray *c1_b_y = NULL;
  int32_T c1_i1;
  real_T c1_b_u[6];
  const mxArray *c1_c_y = NULL;
  int32_T c1_i2;
  real_T c1_c_u[3];
  const mxArray *c1_d_y = NULL;
  uint8_T c1_hoistedGlobal;
  uint8_T c1_d_u;
  const mxArray *c1_e_y = NULL;
  c1_st = NULL;
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellmatrix(4, 1), false);
  for (c1_i0 = 0; c1_i0 < 3; c1_i0++) {
    c1_u[c1_i0] = (*chartInstance->c1_EUL_de)[c1_i0];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  for (c1_i1 = 0; c1_i1 < 6; c1_i1++) {
    c1_b_u[c1_i1] = (*chartInstance->c1_a)[c1_i1];
  }

  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", c1_b_u, 0, 0U, 1U, 0U, 2, 6, 1),
                false);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  for (c1_i2 = 0; c1_i2 < 3; c1_i2++) {
    c1_c_u[c1_i2] = (*chartInstance->c1_dp_d)[c1_i2];
  }

  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", c1_c_u, 0, 0U, 1U, 0U, 2, 3, 1),
                false);
  sf_mex_setcell(c1_y, 2, c1_d_y);
  c1_hoistedGlobal =
    chartInstance->c1_is_active_c1_IPILCO_SimpleImp_relativeRPY_S;
  c1_d_u = c1_hoistedGlobal;
  c1_e_y = NULL;
  sf_mex_assign(&c1_e_y, sf_mex_create("y", &c1_d_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c1_y, 3, c1_e_y);
  sf_mex_assign(&c1_st, c1_y, false);
  return c1_st;
}

static void set_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_st)
{
  const mxArray *c1_u;
  real_T c1_dv0[3];
  int32_T c1_i3;
  real_T c1_dv1[6];
  int32_T c1_i4;
  real_T c1_dv2[3];
  int32_T c1_i5;
  chartInstance->c1_doneDoubleBufferReInit = true;
  c1_u = sf_mex_dup(c1_st);
  c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 0)),
                      "EUL_de", c1_dv0);
  for (c1_i3 = 0; c1_i3 < 3; c1_i3++) {
    (*chartInstance->c1_EUL_de)[c1_i3] = c1_dv0[c1_i3];
  }

  c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 1)), "a",
                        c1_dv1);
  for (c1_i4 = 0; c1_i4 < 6; c1_i4++) {
    (*chartInstance->c1_a)[c1_i4] = c1_dv1[c1_i4];
  }

  c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 2)),
                        "dp_d", c1_dv2);
  for (c1_i5 = 0; c1_i5 < 3; c1_i5++) {
    (*chartInstance->c1_dp_d)[c1_i5] = c1_dv2[c1_i5];
  }

  chartInstance->c1_is_active_c1_IPILCO_SimpleImp_relativeRPY_S =
    c1_n_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 3)),
    "is_active_c1_IPILCO_SimpleImp_relativeRPY_S");
  sf_mex_destroy(&c1_u);
  c1_update_debugger_state_c1_IPILCO_SimpleImp_relativeRPY_S(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  int32_T c1_i6;
  int32_T c1_i7;
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
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_IPILCO_SimpleImp_relativeRPY_S(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY
    (_IPILCO_SimpleImp_relativeRPY_SMachineNumber_, chartInstance->chartNumber,
     chartInstance->instanceNumber);
  for (c1_i6 = 0; c1_i6 < 6; c1_i6++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_a)[c1_i6], 0U);
  }

  for (c1_i7 = 0; c1_i7 < 3; c1_i7++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_dp_d)[c1_i7], 1U);
  }

  for (c1_i8 = 0; c1_i8 < 36; c1_i8++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_Kp)[c1_i8], 2U);
  }

  for (c1_i9 = 0; c1_i9 < 36; c1_i9++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_Kd)[c1_i9], 3U);
  }

  for (c1_i10 = 0; c1_i10 < 16; c1_i10++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_He)[c1_i10], 4U);
  }

  for (c1_i11 = 0; c1_i11 < 6; c1_i11++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_xe)[c1_i11], 5U);
  }

  for (c1_i12 = 0; c1_i12 < 6; c1_i12++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_dxe_n)[c1_i12], 6U);
  }

  for (c1_i13 = 0; c1_i13 < 16; c1_i13++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_Hd)[c1_i13], 7U);
  }

  for (c1_i14 = 0; c1_i14 < 16; c1_i14++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_prevHd)[c1_i14], 8U);
  }

  for (c1_i15 = 0; c1_i15 < 3; c1_i15++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_prevdp_d)[c1_i15], 9U);
  }

  _SFD_DATA_RANGE_CHECK(*chartInstance->c1_dt, 10U);
  for (c1_i16 = 0; c1_i16 < 3; c1_i16++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_EUL_de)[c1_i16], 11U);
  }

  for (c1_i17 = 0; c1_i17 < 3; c1_i17++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c1_prevEUL_de)[c1_i17], 12U);
  }
}

static void mdl_start_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_chartstep_c1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  real_T c1_hoistedGlobal;
  int32_T c1_i18;
  real_T c1_b_Kp[36];
  int32_T c1_i19;
  real_T c1_b_Kd[36];
  int32_T c1_i20;
  real_T c1_b_He[16];
  int32_T c1_i21;
  real_T c1_b_xe[6];
  int32_T c1_i22;
  real_T c1_b_dxe_n[6];
  int32_T c1_i23;
  real_T c1_b_Hd[16];
  int32_T c1_i24;
  real_T c1_b_prevHd[16];
  int32_T c1_i25;
  real_T c1_b_prevdp_d[3];
  real_T c1_b_dt;
  int32_T c1_i26;
  real_T c1_b_prevEUL_de[3];
  uint32_T c1_debug_family_var_map[27];
  real_T c1_R_e[9];
  real_T c1_R_d[9];
  real_T c1_p_e[3];
  real_T c1_p_d[3];
  real_T c1_delta_p[3];
  real_T c1_dp_e[3];
  real_T c1_prevp_d[3];
  real_T c1_ddelta_p[3];
  real_T c1_a_t[3];
  real_T c1_Re_d[9];
  real_T c1_dEUL_de[3];
  real_T c1_a_o[3];
  real_T c1_nargin = 10.0;
  real_T c1_nargout = 3.0;
  real_T c1_b_a[6];
  real_T c1_b_dp_d[3];
  real_T c1_b_EUL_de[3];
  int32_T c1_i27;
  real_T c1_c_He[16];
  real_T c1_dv3[9];
  int32_T c1_i28;
  int32_T c1_i29;
  real_T c1_c_Hd[16];
  real_T c1_dv4[9];
  int32_T c1_i30;
  int32_T c1_i31;
  int32_T c1_i32;
  real_T c1_d_Hd[16];
  real_T c1_dv5[3];
  int32_T c1_i33;
  int32_T c1_i34;
  int32_T c1_i35;
  int32_T c1_i36;
  real_T c1_c_prevHd[16];
  real_T c1_dv6[3];
  int32_T c1_i37;
  int32_T c1_i38;
  real_T c1_b_p_d[3];
  real_T c1_b[3];
  int32_T c1_i39;
  int32_T c1_i40;
  int32_T c1_i41;
  int32_T c1_i42;
  int32_T c1_i43;
  int32_T c1_i44;
  real_T c1_c_a[9];
  int32_T c1_i45;
  int32_T c1_i46;
  real_T c1_y[3];
  int32_T c1_i47;
  int32_T c1_i48;
  int32_T c1_i49;
  int32_T c1_i50;
  int32_T c1_i51;
  int32_T c1_i52;
  int32_T c1_i53;
  int32_T c1_i54;
  real_T c1_b_y[3];
  int32_T c1_i55;
  int32_T c1_i56;
  int32_T c1_i57;
  int32_T c1_i58;
  int32_T c1_i59;
  int32_T c1_i60;
  int32_T c1_i61;
  int32_T c1_i62;
  real_T c1_b_b[9];
  int32_T c1_i63;
  int32_T c1_i64;
  int32_T c1_i65;
  real_T c1_C[9];
  int32_T c1_i66;
  int32_T c1_i67;
  int32_T c1_i68;
  int32_T c1_i69;
  int32_T c1_i70;
  int32_T c1_i71;
  int32_T c1_i72;
  int32_T c1_i73;
  int32_T c1_i74;
  real_T c1_b_Re_d[9];
  real_T c1_dv7[3];
  int32_T c1_i75;
  int32_T c1_i76;
  real_T c1_B;
  int32_T c1_i77;
  real_T c1_c_b[3];
  real_T c1_dv8[3];
  int32_T c1_i78;
  int32_T c1_i79;
  int32_T c1_i80;
  int32_T c1_i81;
  int32_T c1_i82;
  int32_T c1_i83;
  int32_T c1_i84;
  int32_T c1_i85;
  int32_T c1_i86;
  int32_T c1_i87;
  int32_T c1_i88;
  int32_T c1_i89;
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
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  c1_hoistedGlobal = *chartInstance->c1_dt;
  for (c1_i18 = 0; c1_i18 < 36; c1_i18++) {
    c1_b_Kp[c1_i18] = (*chartInstance->c1_Kp)[c1_i18];
  }

  for (c1_i19 = 0; c1_i19 < 36; c1_i19++) {
    c1_b_Kd[c1_i19] = (*chartInstance->c1_Kd)[c1_i19];
  }

  for (c1_i20 = 0; c1_i20 < 16; c1_i20++) {
    c1_b_He[c1_i20] = (*chartInstance->c1_He)[c1_i20];
  }

  for (c1_i21 = 0; c1_i21 < 6; c1_i21++) {
    c1_b_xe[c1_i21] = (*chartInstance->c1_xe)[c1_i21];
  }

  for (c1_i22 = 0; c1_i22 < 6; c1_i22++) {
    c1_b_dxe_n[c1_i22] = (*chartInstance->c1_dxe_n)[c1_i22];
  }

  for (c1_i23 = 0; c1_i23 < 16; c1_i23++) {
    c1_b_Hd[c1_i23] = (*chartInstance->c1_Hd)[c1_i23];
  }

  for (c1_i24 = 0; c1_i24 < 16; c1_i24++) {
    c1_b_prevHd[c1_i24] = (*chartInstance->c1_prevHd)[c1_i24];
  }

  for (c1_i25 = 0; c1_i25 < 3; c1_i25++) {
    c1_b_prevdp_d[c1_i25] = (*chartInstance->c1_prevdp_d)[c1_i25];
  }

  c1_b_dt = c1_hoistedGlobal;
  for (c1_i26 = 0; c1_i26 < 3; c1_i26++) {
    c1_b_prevEUL_de[c1_i26] = (*chartInstance->c1_prevEUL_de)[c1_i26];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 27U, 27U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_R_e, 0U, c1_g_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_R_d, 1U, c1_g_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_p_e, 2U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_p_d, 3U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_delta_p, 4U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_dp_e, 5U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_prevp_d, 6U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_ddelta_p, 7U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_a_t, 8U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_Re_d, 9U, c1_g_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_dEUL_de, 10U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_a_o, 11U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 12U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 13U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_Kp, 14U, c1_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_Kd, 15U, c1_f_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_He, 16U, c1_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_xe, 17U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_dxe_n, 18U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_Hd, 19U, c1_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_prevHd, 20U, c1_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_prevdp_d, 21U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_b_dt, 22U, c1_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_b_prevEUL_de, 23U, c1_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_a, 24U, c1_c_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_dp_d, 25U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_b_EUL_de, 26U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 28);
  for (c1_i27 = 0; c1_i27 < 16; c1_i27++) {
    c1_c_He[c1_i27] = c1_b_He[c1_i27];
  }

  c1_tr2rt(chartInstance, c1_c_He, c1_dv3);
  for (c1_i28 = 0; c1_i28 < 9; c1_i28++) {
    c1_R_e[c1_i28] = c1_dv3[c1_i28];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 29);
  for (c1_i29 = 0; c1_i29 < 16; c1_i29++) {
    c1_c_Hd[c1_i29] = c1_b_Hd[c1_i29];
  }

  c1_tr2rt(chartInstance, c1_c_Hd, c1_dv4);
  for (c1_i30 = 0; c1_i30 < 9; c1_i30++) {
    c1_R_d[c1_i30] = c1_dv4[c1_i30];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 32);
  for (c1_i31 = 0; c1_i31 < 3; c1_i31++) {
    c1_p_e[c1_i31] = c1_b_xe[c1_i31];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 33);
  for (c1_i32 = 0; c1_i32 < 16; c1_i32++) {
    c1_d_Hd[c1_i32] = c1_b_Hd[c1_i32];
  }

  c1_transl(chartInstance, c1_d_Hd, c1_dv5);
  for (c1_i33 = 0; c1_i33 < 3; c1_i33++) {
    c1_p_d[c1_i33] = c1_dv5[c1_i33];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 34);
  for (c1_i34 = 0; c1_i34 < 3; c1_i34++) {
    c1_delta_p[c1_i34] = c1_p_d[c1_i34] - c1_p_e[c1_i34];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 36);
  for (c1_i35 = 0; c1_i35 < 3; c1_i35++) {
    c1_dp_e[c1_i35] = c1_b_dxe_n[c1_i35];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 37);
  for (c1_i36 = 0; c1_i36 < 16; c1_i36++) {
    c1_c_prevHd[c1_i36] = c1_b_prevHd[c1_i36];
  }

  c1_transl(chartInstance, c1_c_prevHd, c1_dv6);
  for (c1_i37 = 0; c1_i37 < 3; c1_i37++) {
    c1_prevp_d[c1_i37] = c1_dv6[c1_i37];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 38);
  for (c1_i38 = 0; c1_i38 < 3; c1_i38++) {
    c1_b_p_d[c1_i38] = c1_p_d[c1_i38] - c1_prevp_d[c1_i38];
  }

  c1_rdivide(chartInstance, c1_b_p_d, c1_b_dt, c1_b);
  for (c1_i39 = 0; c1_i39 < 3; c1_i39++) {
    c1_b_dp_d[c1_i39] = c1_b[c1_i39];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 39);
  for (c1_i40 = 0; c1_i40 < 3; c1_i40++) {
    c1_ddelta_p[c1_i40] = c1_b_dp_d[c1_i40] - c1_dp_e[c1_i40];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 43);
  c1_i41 = 0;
  c1_i42 = 0;
  for (c1_i43 = 0; c1_i43 < 3; c1_i43++) {
    for (c1_i44 = 0; c1_i44 < 3; c1_i44++) {
      c1_c_a[c1_i44 + c1_i41] = c1_b_Kd[c1_i44 + c1_i42];
    }

    c1_i41 += 3;
    c1_i42 += 6;
  }

  for (c1_i45 = 0; c1_i45 < 3; c1_i45++) {
    c1_b[c1_i45] = c1_ddelta_p[c1_i45];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i46 = 0; c1_i46 < 3; c1_i46++) {
    c1_y[c1_i46] = 0.0;
    c1_i47 = 0;
    for (c1_i48 = 0; c1_i48 < 3; c1_i48++) {
      c1_y[c1_i46] += c1_c_a[c1_i47 + c1_i46] * c1_b[c1_i48];
      c1_i47 += 3;
    }
  }

  c1_i49 = 0;
  c1_i50 = 0;
  for (c1_i51 = 0; c1_i51 < 3; c1_i51++) {
    for (c1_i52 = 0; c1_i52 < 3; c1_i52++) {
      c1_c_a[c1_i52 + c1_i49] = c1_b_Kp[c1_i52 + c1_i50];
    }

    c1_i49 += 3;
    c1_i50 += 6;
  }

  for (c1_i53 = 0; c1_i53 < 3; c1_i53++) {
    c1_b[c1_i53] = c1_delta_p[c1_i53];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i54 = 0; c1_i54 < 3; c1_i54++) {
    c1_b_y[c1_i54] = 0.0;
    c1_i55 = 0;
    for (c1_i56 = 0; c1_i56 < 3; c1_i56++) {
      c1_b_y[c1_i54] += c1_c_a[c1_i55 + c1_i54] * c1_b[c1_i56];
      c1_i55 += 3;
    }
  }

  for (c1_i57 = 0; c1_i57 < 3; c1_i57++) {
    c1_a_t[c1_i57] = c1_y[c1_i57] + c1_b_y[c1_i57];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 48);
  c1_i58 = 0;
  for (c1_i59 = 0; c1_i59 < 3; c1_i59++) {
    c1_i60 = 0;
    for (c1_i61 = 0; c1_i61 < 3; c1_i61++) {
      c1_c_a[c1_i61 + c1_i58] = c1_R_e[c1_i60 + c1_i59];
      c1_i60 += 3;
    }

    c1_i58 += 3;
  }

  for (c1_i62 = 0; c1_i62 < 9; c1_i62++) {
    c1_b_b[c1_i62] = c1_R_d[c1_i62];
  }

  c1_b_eml_scalar_eg(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  for (c1_i63 = 0; c1_i63 < 9; c1_i63++) {
    c1_Re_d[c1_i63] = 0.0;
  }

  for (c1_i64 = 0; c1_i64 < 9; c1_i64++) {
    c1_Re_d[c1_i64] = 0.0;
  }

  for (c1_i65 = 0; c1_i65 < 9; c1_i65++) {
    c1_C[c1_i65] = c1_Re_d[c1_i65];
  }

  for (c1_i66 = 0; c1_i66 < 9; c1_i66++) {
    c1_Re_d[c1_i66] = c1_C[c1_i66];
  }

  c1_threshold(chartInstance);
  for (c1_i67 = 0; c1_i67 < 9; c1_i67++) {
    c1_C[c1_i67] = c1_Re_d[c1_i67];
  }

  for (c1_i68 = 0; c1_i68 < 9; c1_i68++) {
    c1_Re_d[c1_i68] = c1_C[c1_i68];
  }

  for (c1_i69 = 0; c1_i69 < 3; c1_i69++) {
    c1_i70 = 0;
    for (c1_i71 = 0; c1_i71 < 3; c1_i71++) {
      c1_Re_d[c1_i70 + c1_i69] = 0.0;
      c1_i72 = 0;
      for (c1_i73 = 0; c1_i73 < 3; c1_i73++) {
        c1_Re_d[c1_i70 + c1_i69] += c1_c_a[c1_i72 + c1_i69] * c1_b_b[c1_i73 +
          c1_i70];
        c1_i72 += 3;
      }

      c1_i70 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 49);
  for (c1_i74 = 0; c1_i74 < 9; c1_i74++) {
    c1_b_Re_d[c1_i74] = c1_Re_d[c1_i74];
  }

  c1_my_tr2rpy(chartInstance, c1_b_Re_d, c1_dv7);
  for (c1_i75 = 0; c1_i75 < 3; c1_i75++) {
    c1_b_EUL_de[c1_i75] = c1_dv7[c1_i75];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 50);
  for (c1_i76 = 0; c1_i76 < 3; c1_i76++) {
    c1_b[c1_i76] = c1_b_EUL_de[c1_i76] - c1_b_prevEUL_de[c1_i76];
  }

  c1_B = c1_b_dt;
  for (c1_i77 = 0; c1_i77 < 3; c1_i77++) {
    c1_c_b[c1_i77] = c1_b[c1_i77];
  }

  c1_rdivide(chartInstance, c1_c_b, c1_B, c1_dv8);
  for (c1_i78 = 0; c1_i78 < 3; c1_i78++) {
    c1_dEUL_de[c1_i78] = c1_dv8[c1_i78];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 56);
  c1_i79 = 0;
  c1_i80 = 0;
  for (c1_i81 = 0; c1_i81 < 3; c1_i81++) {
    for (c1_i82 = 0; c1_i82 < 3; c1_i82++) {
      c1_c_a[c1_i82 + c1_i79] = c1_b_Kd[(c1_i82 + c1_i80) + 21];
    }

    c1_i79 += 3;
    c1_i80 += 6;
  }

  for (c1_i83 = 0; c1_i83 < 3; c1_i83++) {
    c1_b[c1_i83] = c1_dEUL_de[c1_i83];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i84 = 0; c1_i84 < 3; c1_i84++) {
    c1_y[c1_i84] = 0.0;
    c1_i85 = 0;
    for (c1_i86 = 0; c1_i86 < 3; c1_i86++) {
      c1_y[c1_i84] += c1_c_a[c1_i85 + c1_i84] * c1_b[c1_i86];
      c1_i85 += 3;
    }
  }

  c1_i87 = 0;
  c1_i88 = 0;
  for (c1_i89 = 0; c1_i89 < 3; c1_i89++) {
    for (c1_i90 = 0; c1_i90 < 3; c1_i90++) {
      c1_c_a[c1_i90 + c1_i87] = c1_b_Kp[(c1_i90 + c1_i88) + 21];
    }

    c1_i87 += 3;
    c1_i88 += 6;
  }

  for (c1_i91 = 0; c1_i91 < 3; c1_i91++) {
    c1_b[c1_i91] = c1_b_EUL_de[c1_i91];
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  c1_threshold(chartInstance);
  for (c1_i92 = 0; c1_i92 < 3; c1_i92++) {
    c1_b_y[c1_i92] = 0.0;
    c1_i93 = 0;
    for (c1_i94 = 0; c1_i94 < 3; c1_i94++) {
      c1_b_y[c1_i92] += c1_c_a[c1_i93 + c1_i92] * c1_b[c1_i94];
      c1_i93 += 3;
    }
  }

  for (c1_i95 = 0; c1_i95 < 3; c1_i95++) {
    c1_a_o[c1_i95] = c1_y[c1_i95] + c1_b_y[c1_i95];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 58);
  for (c1_i96 = 0; c1_i96 < 3; c1_i96++) {
    c1_b_a[c1_i96] = c1_a_t[c1_i96];
  }

  for (c1_i97 = 0; c1_i97 < 3; c1_i97++) {
    c1_b_a[c1_i97 + 3] = c1_a_o[c1_i97];
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -58);
  _SFD_SYMBOL_SCOPE_POP();
  for (c1_i98 = 0; c1_i98 < 6; c1_i98++) {
    (*chartInstance->c1_a)[c1_i98] = c1_b_a[c1_i98];
  }

  for (c1_i99 = 0; c1_i99 < 3; c1_i99++) {
    (*chartInstance->c1_dp_d)[c1_i99] = c1_b_dp_d[c1_i99];
  }

  for (c1_i100 = 0; c1_i100 < 3; c1_i100++) {
    (*chartInstance->c1_EUL_de)[c1_i100] = c1_b_EUL_de[c1_i100];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
}

static void initSimStructsc1_IPILCO_SimpleImp_relativeRPY_S
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void c1_tr2rt(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
                     *chartInstance, real_T c1_T[16], real_T c1_R[9])
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
  int32_T c1_i101;
  int32_T c1_i102;
  int32_T c1_i103;
  int32_T c1_i104;
  int32_T c1_i105;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c1_d_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_n, 0U, c1_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_t, 1U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 2U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 3U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_T, 4U, c1_e_sf_marshallOut,
    c1_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_R, 5U, c1_g_sf_marshallOut,
    c1_e_sf_marshallIn);
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
  c1_i101 = 0;
  c1_i102 = 0;
  for (c1_i103 = 0; c1_i103 < 3; c1_i103++) {
    for (c1_i104 = 0; c1_i104 < 3; c1_i104++) {
      c1_R[c1_i104 + c1_i101] = c1_T[c1_i104 + c1_i102];
    }

    c1_i101 += 3;
    c1_i102 += 4;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, 56);
  for (c1_i105 = 0; c1_i105 < 3; c1_i105++) {
    c1_t[c1_i105] = c1_T[c1_i105 + 12];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c1_sfEvent, -56);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c1_transl(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
                      *chartInstance, real_T c1_x[16], real_T c1_t1[3])
{
  uint32_T c1_debug_family_var_map[4];
  real_T c1_nargin = 1.0;
  real_T c1_nargout = 1.0;
  boolean_T c1_h;
  real_T c1_d[2];
  real_T c1_b_nargin = 1.0;
  real_T c1_b_nargout = 1.0;
  int32_T c1_i106;
  int32_T c1_i107;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c1_f_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 0U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 1U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_x, 2U, c1_e_sf_marshallOut,
    c1_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_t1, 3U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_SCRIPT_FCN(3, 0);
  _SFD_SCRIPT_CALL(3U, chartInstance->c1_sfEvent, 53);
  CV_SCRIPT_IF(3, 0, CV_RELATIONAL_EVAL(14U, 3U, 0, 1.0, 1.0, -1, 0U, 1));
  _SFD_SCRIPT_CALL(3U, chartInstance->c1_sfEvent, 54);
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c1_e_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_h, 0U, c1_i_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_d, 1U, c1_h_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargin, 2U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_b_nargout, 3U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  CV_SCRIPT_FCN(4, 0);
  _SFD_SCRIPT_CALL(4U, chartInstance->c1_sfEvent, 37);
  for (c1_i106 = 0; c1_i106 < 2; c1_i106++) {
    c1_d[c1_i106] = 4.0;
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
  for (c1_i107 = 0; c1_i107 < 3; c1_i107++) {
    c1_t1[c1_i107] = c1_x[c1_i107 + 12];
  }

  _SFD_SCRIPT_CALL(3U, chartInstance->c1_sfEvent, -89);
  _SFD_SYMBOL_SCOPE_POP();
}

static void c1_my_tr2rpy(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_m[9], real_T c1_rpy[3])
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
  int32_T c1_i108;
  static char_T c1_cv0[3] = { 'x', 'y', 'z' };

  int32_T c1_i109;
  int32_T c1_i110;
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
  int32_T c1_i111;
  real_T c1_b_a[3];
  int32_T c1_i112;
  int32_T c1_i113;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 11U, 11U, c1_g_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_opt, 0U, c1_l_sf_marshallOut,
    c1_h_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_s, 1U, c1_h_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_sr, 2U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_cr, 3U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_sp, 4U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_cp, 5U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_type, 6U, c1_k_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 7U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 8U, c1_d_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_m, 9U, c1_g_sf_marshallOut,
    c1_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_rpy, 10U, c1_j_sf_marshallOut,
    c1_g_sf_marshallIn);
  for (c1_i108 = 0; c1_i108 < 3; c1_i108++) {
    c1_type[c1_i108] = c1_cv0[c1_i108];
  }

  CV_SCRIPT_FCN(5, 0);
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 49);
  c1_opt.deg = false;
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 50);
  c1_opt.zyx = false;
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 51);
  CV_SCRIPT_IF(5, 0, false);
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 56);
  for (c1_i109 = 0; c1_i109 < 2; c1_i109++) {
    c1_s[c1_i109] = 3.0;
  }

  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 57);
  CV_SCRIPT_IF(5, 1, CV_RELATIONAL_EVAL(14U, 5U, 0, 2.0, 2.0, -1, 4U, 0));
  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 64);
  for (c1_i110 = 0; c1_i110 < 3; c1_i110++) {
    c1_rpy[c1_i110] = 0.0;
  }

  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 66);
  if (CV_SCRIPT_IF(5, 2, CV_SCRIPT_MCDC(5, 0, !CV_SCRIPT_COND(5, 0, c1_opt.zyx))))
  {
    _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, 68);
    c1_eps(chartInstance);
    c1_x = c1_m[8];
    c1_b_x = c1_x;
    c1_d0 = muDoubleScalarAbs(c1_b_x);
    guard2 = false;
    if (CV_SCRIPT_COND(5, 1, CV_RELATIONAL_EVAL(14U, 5U, 1, c1_d0,
          2.2204460492503131E-16, -1, 2U, c1_d0 < 2.2204460492503131E-16))) {
      c1_eps(chartInstance);
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
    c1_eps(chartInstance);
    c1_i_x = c1_m[0];
    c1_j_x = c1_i_x;
    c1_d2 = muDoubleScalarAbs(c1_j_x);
    guard1 = false;
    if (CV_SCRIPT_COND(5, 3, CV_RELATIONAL_EVAL(14U, 5U, 3, c1_d2,
          2.2204460492503131E-16, -1, 2U, c1_d2 < 2.2204460492503131E-16))) {
      c1_eps(chartInstance);
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
    for (c1_i111 = 0; c1_i111 < 3; c1_i111++) {
      c1_b_a[c1_i111] = c1_rpy[c1_i111];
    }

    for (c1_i112 = 0; c1_i112 < 3; c1_i112++) {
      c1_b_a[c1_i112] *= 180.0;
    }

    for (c1_i113 = 0; c1_i113 < 3; c1_i113++) {
      c1_rpy[c1_i113] = c1_b_a[c1_i113] / 3.1415926535897931;
    }
  }

  _SFD_SCRIPT_CALL(5U, chartInstance->c1_sfEvent, -97);
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
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i114;
  real_T c1_b_inData[3];
  int32_T c1_i115;
  real_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i114 = 0; c1_i114 < 3; c1_i114++) {
    c1_b_inData[c1_i114] = (*(real_T (*)[3])c1_inData)[c1_i114];
  }

  for (c1_i115 = 0; c1_i115 < 3; c1_i115++) {
    c1_u[c1_i115] = c1_b_inData[c1_i115];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_b_EUL_de, const char_T *c1_identifier, real_T c1_y[3])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_EUL_de), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_EUL_de);
}

static void c1_b_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3])
{
  real_T c1_dv9[3];
  int32_T c1_i116;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv9, 1, 0, 0U, 1, 0U, 1, 3);
  for (c1_i116 = 0; c1_i116 < 3; c1_i116++) {
    c1_y[c1_i116] = c1_dv9[c1_i116];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_EUL_de;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[3];
  int32_T c1_i117;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_b_EUL_de = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_EUL_de), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_EUL_de);
  for (c1_i117 = 0; c1_i117 < 3; c1_i117++) {
    (*(real_T (*)[3])c1_outData)[c1_i117] = c1_y[c1_i117];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i118;
  real_T c1_b_inData[3];
  int32_T c1_i119;
  real_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i118 = 0; c1_i118 < 3; c1_i118++) {
    c1_b_inData[c1_i118] = (*(real_T (*)[3])c1_inData)[c1_i118];
  }

  for (c1_i119 = 0; c1_i119 < 3; c1_i119++) {
    c1_u[c1_i119] = c1_b_inData[c1_i119];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 3, 1), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_c_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_b_dp_d, const char_T *c1_identifier, real_T c1_y[3])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_dp_d), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_dp_d);
}

static void c1_d_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3])
{
  real_T c1_dv10[3];
  int32_T c1_i120;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv10, 1, 0, 0U, 1, 0U, 2, 3, 1);
  for (c1_i120 = 0; c1_i120 < 3; c1_i120++) {
    c1_y[c1_i120] = c1_dv10[c1_i120];
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
  int32_T c1_i121;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_b_dp_d = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_dp_d), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_dp_d);
  for (c1_i121 = 0; c1_i121 < 3; c1_i121++) {
    (*(real_T (*)[3])c1_outData)[c1_i121] = c1_y[c1_i121];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i122;
  real_T c1_b_inData[6];
  int32_T c1_i123;
  real_T c1_u[6];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i122 = 0; c1_i122 < 6; c1_i122++) {
    c1_b_inData[c1_i122] = (*(real_T (*)[6])c1_inData)[c1_i122];
  }

  for (c1_i123 = 0; c1_i123 < 6; c1_i123++) {
    c1_u[c1_i123] = c1_b_inData[c1_i123];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 6, 1), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_e_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_b_a, const char_T *c1_identifier, real_T c1_y[6])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_a), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_a);
}

static void c1_f_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[6])
{
  real_T c1_dv11[6];
  int32_T c1_i124;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv11, 1, 0, 0U, 1, 0U, 2, 6, 1);
  for (c1_i124 = 0; c1_i124 < 6; c1_i124++) {
    c1_y[c1_i124] = c1_dv11[c1_i124];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_a;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[6];
  int32_T c1_i125;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_b_a = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_a), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_b_a);
  for (c1_i125 = 0; c1_i125 < 6; c1_i125++) {
    (*(real_T (*)[6])c1_outData)[c1_i125] = c1_y[c1_i125];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  real_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
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
  int32_T c1_i126;
  int32_T c1_i127;
  int32_T c1_i128;
  real_T c1_b_inData[16];
  int32_T c1_i129;
  int32_T c1_i130;
  int32_T c1_i131;
  real_T c1_u[16];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i126 = 0;
  for (c1_i127 = 0; c1_i127 < 4; c1_i127++) {
    for (c1_i128 = 0; c1_i128 < 4; c1_i128++) {
      c1_b_inData[c1_i128 + c1_i126] = (*(real_T (*)[16])c1_inData)[c1_i128 +
        c1_i126];
    }

    c1_i126 += 4;
  }

  c1_i129 = 0;
  for (c1_i130 = 0; c1_i130 < 4; c1_i130++) {
    for (c1_i131 = 0; c1_i131 < 4; c1_i131++) {
      c1_u[c1_i131 + c1_i129] = c1_b_inData[c1_i131 + c1_i129];
    }

    c1_i129 += 4;
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
  int32_T c1_i132;
  int32_T c1_i133;
  int32_T c1_i134;
  real_T c1_b_inData[36];
  int32_T c1_i135;
  int32_T c1_i136;
  int32_T c1_i137;
  real_T c1_u[36];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i132 = 0;
  for (c1_i133 = 0; c1_i133 < 6; c1_i133++) {
    for (c1_i134 = 0; c1_i134 < 6; c1_i134++) {
      c1_b_inData[c1_i134 + c1_i132] = (*(real_T (*)[36])c1_inData)[c1_i134 +
        c1_i132];
    }

    c1_i132 += 6;
  }

  c1_i135 = 0;
  for (c1_i136 = 0; c1_i136 < 6; c1_i136++) {
    for (c1_i137 = 0; c1_i137 < 6; c1_i137++) {
      c1_u[c1_i137 + c1_i135] = c1_b_inData[c1_i137 + c1_i135];
    }

    c1_i135 += 6;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 6, 6), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static real_T c1_g_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d4;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d4, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d4;
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
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_nargout = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_nargout), &c1_thisId);
  sf_mex_destroy(&c1_nargout);
  *(real_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_g_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i138;
  int32_T c1_i139;
  int32_T c1_i140;
  real_T c1_b_inData[9];
  int32_T c1_i141;
  int32_T c1_i142;
  int32_T c1_i143;
  real_T c1_u[9];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i138 = 0;
  for (c1_i139 = 0; c1_i139 < 3; c1_i139++) {
    for (c1_i140 = 0; c1_i140 < 3; c1_i140++) {
      c1_b_inData[c1_i140 + c1_i138] = (*(real_T (*)[9])c1_inData)[c1_i140 +
        c1_i138];
    }

    c1_i138 += 3;
  }

  c1_i141 = 0;
  for (c1_i142 = 0; c1_i142 < 3; c1_i142++) {
    for (c1_i143 = 0; c1_i143 < 3; c1_i143++) {
      c1_u[c1_i143 + c1_i141] = c1_b_inData[c1_i143 + c1_i141];
    }

    c1_i141 += 3;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 3, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_h_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[9])
{
  real_T c1_dv12[9];
  int32_T c1_i144;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv12, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c1_i144 = 0; c1_i144 < 9; c1_i144++) {
    c1_y[c1_i144] = c1_dv12[c1_i144];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_Re_d;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[9];
  int32_T c1_i145;
  int32_T c1_i146;
  int32_T c1_i147;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_Re_d = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_Re_d), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_Re_d);
  c1_i145 = 0;
  for (c1_i146 = 0; c1_i146 < 3; c1_i146++) {
    for (c1_i147 = 0; c1_i147 < 3; c1_i147++) {
      (*(real_T (*)[9])c1_outData)[c1_i147 + c1_i145] = c1_y[c1_i147 + c1_i145];
    }

    c1_i145 += 3;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static void c1_i_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[16])
{
  real_T c1_dv13[16];
  int32_T c1_i148;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv13, 1, 0, 0U, 1, 0U, 2, 4, 4);
  for (c1_i148 = 0; c1_i148 < 16; c1_i148++) {
    c1_y[c1_i148] = c1_dv13[c1_i148];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_T;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[16];
  int32_T c1_i149;
  int32_T c1_i150;
  int32_T c1_i151;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_T = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_T), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_T);
  c1_i149 = 0;
  for (c1_i150 = 0; c1_i150 < 4; c1_i150++) {
    for (c1_i151 = 0; c1_i151 < 4; c1_i151++) {
      (*(real_T (*)[16])c1_outData)[c1_i151 + c1_i149] = c1_y[c1_i151 + c1_i149];
    }

    c1_i149 += 4;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_h_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i152;
  real_T c1_b_inData[2];
  int32_T c1_i153;
  real_T c1_u[2];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i152 = 0; c1_i152 < 2; c1_i152++) {
    c1_b_inData[c1_i152] = (*(real_T (*)[2])c1_inData)[c1_i152];
  }

  for (c1_i153 = 0; c1_i153 < 2; c1_i153++) {
    c1_u[c1_i153] = c1_b_inData[c1_i153];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 1, 2), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_i_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  boolean_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(boolean_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 11, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_j_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i154;
  real_T c1_b_inData[3];
  int32_T c1_i155;
  real_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i154 = 0; c1_i154 < 3; c1_i154++) {
    c1_b_inData[c1_i154] = (*(real_T (*)[3])c1_inData)[c1_i154];
  }

  for (c1_i155 = 0; c1_i155 < 3; c1_i155++) {
    c1_u[c1_i155] = c1_b_inData[c1_i155];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 1, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static void c1_j_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[3])
{
  real_T c1_dv14[3];
  int32_T c1_i156;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv14, 1, 0, 0U, 1, 0U, 2, 1, 3);
  for (c1_i156 = 0; c1_i156 < 3; c1_i156++) {
    c1_y[c1_i156] = c1_dv14[c1_i156];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_rpy;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[3];
  int32_T c1_i157;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_rpy = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_rpy), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_rpy);
  for (c1_i157 = 0; c1_i157 < 3; c1_i157++) {
    (*(real_T (*)[3])c1_outData)[c1_i157] = c1_y[c1_i157];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_k_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i158;
  char_T c1_b_inData[3];
  int32_T c1_i159;
  char_T c1_u[3];
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i158 = 0; c1_i158 < 3; c1_i158++) {
    c1_b_inData[c1_i158] = (*(char_T (*)[3])c1_inData)[c1_i158];
  }

  for (c1_i159 = 0; c1_i159 < 3; c1_i159++) {
    c1_u[c1_i159] = c1_b_inData[c1_i159];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 3), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static const mxArray *c1_l_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  c1_snY5FzvLGhd6UMWwSEwxMX c1_u;
  const mxArray *c1_y = NULL;
  boolean_T c1_b_u;
  const mxArray *c1_b_y = NULL;
  boolean_T c1_c_u;
  const mxArray *c1_c_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
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

static c1_snY5FzvLGhd6UMWwSEwxMX c1_k_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  c1_snY5FzvLGhd6UMWwSEwxMX c1_y;
  emlrtMsgIdentifier c1_thisId;
  static const char * c1_fieldNames[2] = { "deg", "zyx" };

  c1_thisId.fParent = c1_parentId;
  sf_mex_check_struct(c1_parentId, c1_u, 2, c1_fieldNames, 0U, NULL);
  c1_thisId.fIdentifier = "deg";
  c1_y.deg = c1_l_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "deg", "deg", 0)), &c1_thisId);
  c1_thisId.fIdentifier = "zyx";
  c1_y.zyx = c1_l_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getfield
    (c1_u, "zyx", "zyx", 0)), &c1_thisId);
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static boolean_T c1_l_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  boolean_T c1_y;
  boolean_T c1_b0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_b0, 1, 11, 0U, 0, 0U, 0);
  c1_y = c1_b0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_opt;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  c1_snY5FzvLGhd6UMWwSEwxMX c1_y;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_opt = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_opt), &c1_thisId);
  sf_mex_destroy(&c1_opt);
  *(c1_snY5FzvLGhd6UMWwSEwxMX *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray
  *sf_c1_IPILCO_SimpleImp_relativeRPY_S_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  sf_mex_assign(&c1_nameCaptureInfo, sf_mex_createstruct("structure", 2, 53, 1),
                false);
  c1_info_helper(&c1_nameCaptureInfo);
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
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("transl"), "name", "name", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/transl.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464793910U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/robot/transl.m"),
                  "context", "context", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("ishomog"), "name", "name", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/ishomog.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464793908U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/Toolboxes/robot-9.10/rvctools/common/ishomog.m"),
                  "context", "context", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("all"), "name", "name", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372582414U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "context", "context",
                  6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "context", "context",
                  7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 7);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/all.m"), "context", "context",
                  8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.allOrAny"),
                  "name", "name", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583158U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("isequal"), "name", "name", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818758U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_isequal_core"), "name",
                  "name", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818786U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/allOrAny.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.constNonSingletonDim"), "name", "name", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("logical"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/constNonSingletonDim.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
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
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713880U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
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
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818796U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_div"), "name", "name", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1386423952U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
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
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 19);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "context", "context", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980690U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
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
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.xgemm"),
                  "name", "name", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p!below_threshold"),
                  "context", "context", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalarEg"),
                  "name", "name", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalarEg.p"),
                  "resolved", "resolved", 29);
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
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xgemm.p"),
                  "context", "context", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.refblas.xgemm"),
                  "name", "name", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xgemm.p"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
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
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("my_tr2rpy"), "name", "name",
                  31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "resolved", "resolved", 31);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1464800136U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("length"), "name", "name", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1303146206U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("abs"), "name", "name", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 33);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713852U), "fileTimeLo",
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
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818712U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eps"), "name", "name", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "resolved",
                  "resolved", 36);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m"), "context",
                  "context", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_eps"), "name", "name", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "resolved",
                  "resolved", 37);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m"), "context",
                  "context", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("atan2"), "name", "name", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m"), "resolved",
                  "resolved", 39);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m"), "context",
                  "context", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 40);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m"), "context",
                  "context", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 41);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1375980688U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "context", "context", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.scalexpAlloc"),
                  "name", "name", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/scalexpAlloc.p"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m"), "context",
                  "context", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_atan2"), "name",
                  "name", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 43);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818720U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("sin"), "name", "name", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 44);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328504U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818736U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("cos"), "name", "name", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1395328496U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 47);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 47);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1286818722U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 48);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 48);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
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
    "[E]D:/victor/MATLAB/I-PILCO/util/my_tr2rpy.m"), "context", "context", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 49);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807648U), "fileTimeLo",
                  "fileTimeLo", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 49);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370009886U), "mFileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 50);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 50);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
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
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("rdivide"), "name", "name", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 51);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 51);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1363713880U), "fileTimeLo",
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
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(""), "context", "context", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("mrdivide"), "name", "name", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 52);
  sf_mex_addfield(*c1_info, c1_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1410807648U), "fileTimeLo",
                  "fileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(1370009886U), "mFileTimeLo",
                  "mFileTimeLo", 52);
  sf_mex_addfield(*c1_info, c1_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 52);
  sf_mex_assign(&c1_rhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c1_lhs52, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_rhs52), "rhs", "rhs",
                  52);
  sf_mex_addfield(*c1_info, sf_mex_duplicatearraysafe(&c1_lhs52), "lhs", "lhs",
                  52);
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

static void c1_rdivide(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_x[3], real_T c1_y, real_T c1_z[3])
{
  real_T c1_b_y;
  real_T c1_c_y;
  int32_T c1_i160;
  (void)chartInstance;
  c1_b_y = c1_y;
  c1_c_y = c1_b_y;
  for (c1_i160 = 0; c1_i160 < 3; c1_i160++) {
    c1_z[c1_i160] = c1_x[c1_i160] / c1_c_y;
  }
}

static void c1_eml_scalar_eg(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_threshold(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_b_eml_scalar_eg(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c1_eps(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
                   *chartInstance)
{
  (void)chartInstance;
}

static real_T c1_atan2(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
  *chartInstance, real_T c1_y, real_T c1_x)
{
  real_T c1_b_y;
  real_T c1_b_x;
  (void)chartInstance;
  c1_b_y = c1_y;
  c1_b_x = c1_x;
  return muDoubleScalarAtan2(c1_b_y, c1_b_x);
}

static const mxArray *c1_m_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(int32_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, false);
  return c1_mxArrayOutData;
}

static int32_T c1_m_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_y;
  int32_T c1_i161;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_i161, 1, 6, 0U, 0, 0U, 0);
  c1_y = c1_i161;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_i_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_y;
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
    chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_m_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_n_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_b_is_active_c1_IPILCO_SimpleImp_relativeRPY_S, const char_T
   *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_o_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_b_is_active_c1_IPILCO_SimpleImp_relativeRPY_S), &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_IPILCO_SimpleImp_relativeRPY_S);
  return c1_y;
}

static uint8_T c1_o_emlrt_marshallIn
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance, const
   mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  (void)chartInstance;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void init_dsm_address_info
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address
  (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance)
{
  chartInstance->c1_a = (real_T (*)[6])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c1_dp_d = (real_T (*)[3])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 2);
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
  chartInstance->c1_Hd = (real_T (*)[16])ssGetInputPortSignal_wrapper
    (chartInstance->S, 5);
  chartInstance->c1_prevHd = (real_T (*)[16])ssGetInputPortSignal_wrapper
    (chartInstance->S, 6);
  chartInstance->c1_prevdp_d = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 7);
  chartInstance->c1_dt = (real_T *)ssGetInputPortSignal_wrapper(chartInstance->S,
    8);
  chartInstance->c1_EUL_de = (real_T (*)[3])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 3);
  chartInstance->c1_prevEUL_de = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 9);
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

void sf_c1_IPILCO_SimpleImp_relativeRPY_S_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3088004942U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2060234051U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1949604142U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2049975393U);
}

mxArray* sf_c1_IPILCO_SimpleImp_relativeRPY_S_get_post_codegen_info(void);
mxArray *sf_c1_IPILCO_SimpleImp_relativeRPY_S_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("I3pImVMEkVlOgcCXv0N8EG");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,10,3,dataFields);

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
      pr[0] = (double)(4);
      pr[1] = (double)(4);
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
      pr[0] = (double)(3);
      pr[1] = (double)(1);
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
      pr[0] = (double)(1);
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
      pr[1] = (double)(1);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

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
      pr[0] = (double)(3);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxArray* mxPostCodegenInfo =
      sf_c1_IPILCO_SimpleImp_relativeRPY_S_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c1_IPILCO_SimpleImp_relativeRPY_S_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c1_IPILCO_SimpleImp_relativeRPY_S_jit_fallback_info(void)
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

mxArray *sf_c1_IPILCO_SimpleImp_relativeRPY_S_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c1_IPILCO_SimpleImp_relativeRPY_S_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c1_IPILCO_SimpleImp_relativeRPY_S
  (void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x4'type','srcId','name','auxInfo'{{M[1],M[25],T\"EUL_de\",},{M[1],M[16],T\"a\",},{M[1],M[21],T\"dp_d\",},{M[8],M[0],T\"is_active_c1_IPILCO_SimpleImp_relativeRPY_S\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 4, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_IPILCO_SimpleImp_relativeRPY_S_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _IPILCO_SimpleImp_relativeRPY_SMachineNumber_,
           1,
           1,
           1,
           0,
           13,
           0,
           0,
           0,
           0,
           6,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize its own list of scripts */
        init_script_number_translation
          (_IPILCO_SimpleImp_relativeRPY_SMachineNumber_,
           chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,
             _IPILCO_SimpleImp_relativeRPY_SMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _IPILCO_SimpleImp_relativeRPY_SMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,2,0,1,"a");
          _SFD_SET_DATA_PROPS(1,2,0,1,"dp_d");
          _SFD_SET_DATA_PROPS(2,1,1,0,"Kp");
          _SFD_SET_DATA_PROPS(3,1,1,0,"Kd");
          _SFD_SET_DATA_PROPS(4,1,1,0,"He");
          _SFD_SET_DATA_PROPS(5,1,1,0,"xe");
          _SFD_SET_DATA_PROPS(6,1,1,0,"dxe_n");
          _SFD_SET_DATA_PROPS(7,1,1,0,"Hd");
          _SFD_SET_DATA_PROPS(8,1,1,0,"prevHd");
          _SFD_SET_DATA_PROPS(9,1,1,0,"prevdp_d");
          _SFD_SET_DATA_PROPS(10,1,1,0,"dt");
          _SFD_SET_DATA_PROPS(11,2,0,1,"EUL_de");
          _SFD_SET_DATA_PROPS(12,1,1,0,"prevEUL_de");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,2512);
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
          dimVector[0]= 3;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)
            c1_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_f_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_f_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 6;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)NULL);
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
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 4;
          dimVector[1]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_d_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)
            c1_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(12,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_VALUE_PTR(0U, *chartInstance->c1_a);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c1_dp_d);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c1_Kp);
        _SFD_SET_DATA_VALUE_PTR(3U, *chartInstance->c1_Kd);
        _SFD_SET_DATA_VALUE_PTR(4U, *chartInstance->c1_He);
        _SFD_SET_DATA_VALUE_PTR(5U, *chartInstance->c1_xe);
        _SFD_SET_DATA_VALUE_PTR(6U, *chartInstance->c1_dxe_n);
        _SFD_SET_DATA_VALUE_PTR(7U, *chartInstance->c1_Hd);
        _SFD_SET_DATA_VALUE_PTR(8U, *chartInstance->c1_prevHd);
        _SFD_SET_DATA_VALUE_PTR(9U, *chartInstance->c1_prevdp_d);
        _SFD_SET_DATA_VALUE_PTR(10U, chartInstance->c1_dt);
        _SFD_SET_DATA_VALUE_PTR(11U, *chartInstance->c1_EUL_de);
        _SFD_SET_DATA_VALUE_PTR(12U, *chartInstance->c1_prevEUL_de);
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _IPILCO_SimpleImp_relativeRPY_SMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "sBbQaVuKcSnB4A62vlnLUH";
}

static void sf_opaque_initialize_c1_IPILCO_SimpleImp_relativeRPY_S(void
  *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c1_IPILCO_SimpleImp_relativeRPY_S
    ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*) chartInstanceVar);
  initialize_c1_IPILCO_SimpleImp_relativeRPY_S
    ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_IPILCO_SimpleImp_relativeRPY_S(void
  *chartInstanceVar)
{
  enable_c1_IPILCO_SimpleImp_relativeRPY_S
    ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_IPILCO_SimpleImp_relativeRPY_S(void
  *chartInstanceVar)
{
  disable_c1_IPILCO_SimpleImp_relativeRPY_S
    ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_IPILCO_SimpleImp_relativeRPY_S(void
  *chartInstanceVar)
{
  sf_gateway_c1_IPILCO_SimpleImp_relativeRPY_S
    ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*) chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S
  (SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S
    ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*)
     chartInfo->chartInstance);        /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S(SimStruct*
  S, const mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S
    ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*)
     chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c1_IPILCO_SimpleImp_relativeRPY_S(void
  *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*)
                    chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_IPILCO_SimpleImp_relativeRPY_S_optimization_info();
    }

    finalize_c1_IPILCO_SimpleImp_relativeRPY_S
      ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*) chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_IPILCO_SimpleImp_relativeRPY_S
    ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_IPILCO_SimpleImp_relativeRPY_S(SimStruct *S)
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
    initialize_params_c1_IPILCO_SimpleImp_relativeRPY_S
      ((SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*)
       (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c1_IPILCO_SimpleImp_relativeRPY_S(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_IPILCO_SimpleImp_relativeRPY_S_optimization_info();
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
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,10);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,3);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=3; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 10; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(126559395U));
  ssSetChecksum1(S,(1701209393U));
  ssSetChecksum2(S,(4050900960U));
  ssSetChecksum3(S,(1382719332U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c1_IPILCO_SimpleImp_relativeRPY_S(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_IPILCO_SimpleImp_relativeRPY_S(SimStruct *S)
{
  SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct *)utMalloc
    (sizeof(SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct));
  memset(chartInstance, 0, sizeof
         (SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.enableChart =
    sf_opaque_enable_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.disableChart =
    sf_opaque_disable_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_IPILCO_SimpleImp_relativeRPY_S;
  chartInstance->chartInfo.mdlSetWorkWidths =
    mdlSetWorkWidths_c1_IPILCO_SimpleImp_relativeRPY_S;
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

void c1_IPILCO_SimpleImp_relativeRPY_S_method_dispatcher(SimStruct *S, int_T
  method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_IPILCO_SimpleImp_relativeRPY_S(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_IPILCO_SimpleImp_relativeRPY_S(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_IPILCO_SimpleImp_relativeRPY_S(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_IPILCO_SimpleImp_relativeRPY_S_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
