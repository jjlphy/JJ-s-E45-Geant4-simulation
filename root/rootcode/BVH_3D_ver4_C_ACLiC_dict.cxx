// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIUsersdIjjdIworkdIhyptpcmI11dO0dO2dIrootdIrootcodedIBVH_3D_ver4_C_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/Users/jj/work/hyptpc-11.0.2/root/rootcode/./BVH_3D_ver4.C"

// Header files passed via #pragma extra_include

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_BVH_3D_ver4_C_ACLiC_dict() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_BVH_3D_ver4_C_ACLiC_dict) = RecordReadRules_BVH_3D_ver4_C_ACLiC_dict();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_BVH_3D_ver4_C_ACLiC_dict));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_BVH_3D_ver4_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./BVH_3D_ver4.C",
nullptr
    };
    static const char* includePaths[] = {
"/opt/homebrew/Cellar/root/6.36.00/include/root",
"/opt/homebrew/Cellar/root/6.36.00/etc/root",
"/opt/homebrew/Cellar/root/6.36.00/etc/root/cling",
"/opt/homebrew/Cellar/root/6.36.00/etc/root/cling/plugins/include",
"/opt/homebrew/Cellar/root/6.36.00/include/root",
"/opt/homebrew/Cellar/root/6.36.00/include",
"/opt/homebrew/Cellar/root/6.36.00/include/root",
"/Users/jj/work/hyptpc-11.0.2/root/rootcode/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "BVH_3D_ver4_C_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "BVH_3D_ver4_C_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./BVH_3D_ver4.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"BVH_3D", payloadCode, "@",
"N_BH2", payloadCode, "@",
"N_BVHU", payloadCode, "@",
"N_SCH", payloadCode, "@",
"bind_branch", payloadCode, "@",
"distinguish_dilated", payloadCode, "@",
"draw_dilated_boxes", payloadCode, "@",
"get_unique_hits", payloadCode, "@",
"kExcludedTriplets", payloadCode, "@",
"kIncludedTriplets", payloadCode, "@",
"key3", payloadCode, "@",
"run_analysis", payloadCode, "@",
"sch_dilate_halfwidth", payloadCode, "@",
"u_dilate_halfwidth", payloadCode, "@",
"use_sch_dilation", payloadCode, "@",
"use_u_dilation", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("BVH_3D_ver4_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_BVH_3D_ver4_C_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_BVH_3D_ver4_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_BVH_3D_ver4_C_ACLiC_dict() {
  TriggerDictionaryInitialization_BVH_3D_ver4_C_ACLiC_dict_Impl();
}
