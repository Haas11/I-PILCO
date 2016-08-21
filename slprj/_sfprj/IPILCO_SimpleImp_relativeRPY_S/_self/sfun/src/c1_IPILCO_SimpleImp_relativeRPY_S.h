#ifndef __c1_IPILCO_SimpleImp_relativeRPY_S_h__
#define __c1_IPILCO_SimpleImp_relativeRPY_S_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef struct_snY5FzvLGhd6UMWwSEwxMX
#define struct_snY5FzvLGhd6UMWwSEwxMX

struct snY5FzvLGhd6UMWwSEwxMX
{
  boolean_T deg;
  boolean_T zyx;
};

#endif                                 /*struct_snY5FzvLGhd6UMWwSEwxMX*/

#ifndef typedef_c1_snY5FzvLGhd6UMWwSEwxMX
#define typedef_c1_snY5FzvLGhd6UMWwSEwxMX

typedef struct snY5FzvLGhd6UMWwSEwxMX c1_snY5FzvLGhd6UMWwSEwxMX;

#endif                                 /*typedef_c1_snY5FzvLGhd6UMWwSEwxMX*/

#ifndef typedef_SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
#define typedef_SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_isStable;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_IPILCO_SimpleImp_relativeRPY_S;
  real_T (*c1_a)[6];
  real_T (*c1_dp_d)[3];
  real_T (*c1_Kp)[36];
  real_T (*c1_Kd)[36];
  real_T (*c1_He)[16];
  real_T (*c1_xe)[6];
  real_T (*c1_dxe_n)[6];
  real_T (*c1_Hd)[16];
  real_T (*c1_prevHd)[16];
  real_T (*c1_prevdp_d)[3];
  real_T *c1_dt;
  real_T (*c1_EUL_de)[3];
  real_T (*c1_prevEUL_de)[3];
} SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct;

#endif                                 /*typedef_SFc1_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray
  *sf_c1_IPILCO_SimpleImp_relativeRPY_S_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_IPILCO_SimpleImp_relativeRPY_S_get_check_sum(mxArray *plhs[]);
extern void c1_IPILCO_SimpleImp_relativeRPY_S_method_dispatcher(SimStruct *S,
  int_T method, void *data);

#endif
