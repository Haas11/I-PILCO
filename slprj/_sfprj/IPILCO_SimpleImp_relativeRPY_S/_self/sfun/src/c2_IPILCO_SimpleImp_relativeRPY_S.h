#ifndef __c2_IPILCO_SimpleImp_relativeRPY_S_h__
#define __c2_IPILCO_SimpleImp_relativeRPY_S_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef struct_s3TYGqYEZfhjEunCONSnN4C
#define struct_s3TYGqYEZfhjEunCONSnN4C

struct s3TYGqYEZfhjEunCONSnN4C
{
  boolean_T zyx;
  boolean_T deg;
};

#endif                                 /*struct_s3TYGqYEZfhjEunCONSnN4C*/

#ifndef typedef_c2_s3TYGqYEZfhjEunCONSnN4C
#define typedef_c2_s3TYGqYEZfhjEunCONSnN4C

typedef struct s3TYGqYEZfhjEunCONSnN4C c2_s3TYGqYEZfhjEunCONSnN4C;

#endif                                 /*typedef_c2_s3TYGqYEZfhjEunCONSnN4C*/

#ifndef typedef_SFc2_IPILCO_SimpleImp_relativeRPY_SInstanceStruct
#define typedef_SFc2_IPILCO_SimpleImp_relativeRPY_SInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c2_sfEvent;
  boolean_T c2_isStable;
  boolean_T c2_doneDoubleBufferReInit;
  uint8_T c2_is_active_c2_IPILCO_SimpleImp_relativeRPY_S;
  real_T c2_INHOLE;
  boolean_T c2_INHOLE_not_empty;
  boolean_T c2_INSERTED;
  boolean_T c2_INSERTED_not_empty;
  real_T (*c2_Kp_env)[6];
  real_T (*c2_y)[12];
  real_T (*c2_Kd_env)[6];
  real_T (*c2_xe)[6];
  real_T (*c2_dxe)[6];
  real_T (*c2_xc)[6];
  real_T (*c2_x_hole)[3];
  real_T *c2_unused;
} SFc2_IPILCO_SimpleImp_relativeRPY_SInstanceStruct;

#endif                                 /*typedef_SFc2_IPILCO_SimpleImp_relativeRPY_SInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray
  *sf_c2_IPILCO_SimpleImp_relativeRPY_S_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c2_IPILCO_SimpleImp_relativeRPY_S_get_check_sum(mxArray *plhs[]);
extern void c2_IPILCO_SimpleImp_relativeRPY_S_method_dispatcher(SimStruct *S,
  int_T method, void *data);

#endif
