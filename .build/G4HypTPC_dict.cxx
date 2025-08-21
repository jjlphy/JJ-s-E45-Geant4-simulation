// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G4HypTPC_dict
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

// Header files passed as explicit arguments
#include "TParticle.h"
#include "TVector3.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *vectorlEvoidmUgR_Dictionary();
   static void vectorlEvoidmUgR_TClassManip(TClass*);
   static void *new_vectorlEvoidmUgR(void *p = nullptr);
   static void *newArray_vectorlEvoidmUgR(Long_t size, void *p);
   static void delete_vectorlEvoidmUgR(void *p);
   static void deleteArray_vectorlEvoidmUgR(void *p);
   static void destruct_vectorlEvoidmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<void*>*)
   {
      vector<void*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<void*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<void*>", -2, "vector", 389,
                  typeid(vector<void*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvoidmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<void*>) );
      instance.SetNew(&new_vectorlEvoidmUgR);
      instance.SetNewArray(&newArray_vectorlEvoidmUgR);
      instance.SetDelete(&delete_vectorlEvoidmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvoidmUgR);
      instance.SetDestructor(&destruct_vectorlEvoidmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<void*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<void*>","std::vector<void*, std::allocator<void*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<void*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvoidmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<void*>*>(nullptr))->GetClass();
      vectorlEvoidmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvoidmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvoidmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<void*> : new vector<void*>;
   }
   static void *newArray_vectorlEvoidmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<void*>[nElements] : new vector<void*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvoidmUgR(void *p) {
      delete (static_cast<vector<void*>*>(p));
   }
   static void deleteArray_vectorlEvoidmUgR(void *p) {
      delete [] (static_cast<vector<void*>*>(p));
   }
   static void destruct_vectorlEvoidmUgR(void *p) {
      typedef vector<void*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<void*>

namespace ROOT {
   static TClass *vectorlEvectorlEvoidmUgRsPgR_Dictionary();
   static void vectorlEvectorlEvoidmUgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEvoidmUgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEvoidmUgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEvoidmUgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEvoidmUgRsPgR(void *p);
   static void destruct_vectorlEvectorlEvoidmUgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<void*> >*)
   {
      vector<vector<void*> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<void*> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<void*> >", -2, "vector", 389,
                  typeid(vector<vector<void*> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEvoidmUgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<void*> >) );
      instance.SetNew(&new_vectorlEvectorlEvoidmUgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEvoidmUgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEvoidmUgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEvoidmUgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEvoidmUgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<void*> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<void*> >","std::vector<std::vector<void*, std::allocator<void*> >, std::allocator<std::vector<void*, std::allocator<void*> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<void*> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEvoidmUgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<void*> >*>(nullptr))->GetClass();
      vectorlEvectorlEvoidmUgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEvoidmUgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEvoidmUgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<void*> > : new vector<vector<void*> >;
   }
   static void *newArray_vectorlEvectorlEvoidmUgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<void*> >[nElements] : new vector<vector<void*> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEvoidmUgRsPgR(void *p) {
      delete (static_cast<vector<vector<void*> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEvoidmUgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<void*> >*>(p));
   }
   static void destruct_vectorlEvectorlEvoidmUgRsPgR(void *p) {
      typedef vector<vector<void*> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<void*> >

namespace ROOT {
   static TClass *vectorlEvectorlEunsignedsPshortgRsPgR_Dictionary();
   static void vectorlEvectorlEunsignedsPshortgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEunsignedsPshortgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEunsignedsPshortgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEunsignedsPshortgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEunsignedsPshortgRsPgR(void *p);
   static void destruct_vectorlEvectorlEunsignedsPshortgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<unsigned short> >*)
   {
      vector<vector<unsigned short> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<unsigned short> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<unsigned short> >", -2, "vector", 389,
                  typeid(vector<vector<unsigned short> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEunsignedsPshortgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<unsigned short> >) );
      instance.SetNew(&new_vectorlEvectorlEunsignedsPshortgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEunsignedsPshortgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEunsignedsPshortgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEunsignedsPshortgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEunsignedsPshortgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<unsigned short> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<unsigned short> >","std::vector<std::vector<unsigned short, std::allocator<unsigned short> >, std::allocator<std::vector<unsigned short, std::allocator<unsigned short> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned short> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEunsignedsPshortgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned short> >*>(nullptr))->GetClass();
      vectorlEvectorlEunsignedsPshortgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEunsignedsPshortgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEunsignedsPshortgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned short> > : new vector<vector<unsigned short> >;
   }
   static void *newArray_vectorlEvectorlEunsignedsPshortgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned short> >[nElements] : new vector<vector<unsigned short> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEunsignedsPshortgRsPgR(void *p) {
      delete (static_cast<vector<vector<unsigned short> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEunsignedsPshortgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<unsigned short> >*>(p));
   }
   static void destruct_vectorlEvectorlEunsignedsPshortgRsPgR(void *p) {
      typedef vector<vector<unsigned short> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<unsigned short> >

namespace ROOT {
   static TClass *vectorlEvectorlEunsignedsPlonggRsPgR_Dictionary();
   static void vectorlEvectorlEunsignedsPlonggRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEunsignedsPlonggRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEunsignedsPlonggRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEunsignedsPlonggRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEunsignedsPlonggRsPgR(void *p);
   static void destruct_vectorlEvectorlEunsignedsPlonggRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<unsigned long> >*)
   {
      vector<vector<unsigned long> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<unsigned long> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<unsigned long> >", -2, "vector", 389,
                  typeid(vector<vector<unsigned long> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEunsignedsPlonggRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<unsigned long> >) );
      instance.SetNew(&new_vectorlEvectorlEunsignedsPlonggRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEunsignedsPlonggRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEunsignedsPlonggRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEunsignedsPlonggRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEunsignedsPlonggRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<unsigned long> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<unsigned long> >","std::vector<std::vector<unsigned long, std::allocator<unsigned long> >, std::allocator<std::vector<unsigned long, std::allocator<unsigned long> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned long> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEunsignedsPlonggRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned long> >*>(nullptr))->GetClass();
      vectorlEvectorlEunsignedsPlonggRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEunsignedsPlonggRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEunsignedsPlonggRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned long> > : new vector<vector<unsigned long> >;
   }
   static void *newArray_vectorlEvectorlEunsignedsPlonggRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned long> >[nElements] : new vector<vector<unsigned long> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEunsignedsPlonggRsPgR(void *p) {
      delete (static_cast<vector<vector<unsigned long> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEunsignedsPlonggRsPgR(void *p) {
      delete [] (static_cast<vector<vector<unsigned long> >*>(p));
   }
   static void destruct_vectorlEvectorlEunsignedsPlonggRsPgR(void *p) {
      typedef vector<vector<unsigned long> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<unsigned long> >

namespace ROOT {
   static TClass *vectorlEvectorlEunsignedsPintgRsPgR_Dictionary();
   static void vectorlEvectorlEunsignedsPintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEunsignedsPintgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEunsignedsPintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEunsignedsPintgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEunsignedsPintgRsPgR(void *p);
   static void destruct_vectorlEvectorlEunsignedsPintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<unsigned int> >*)
   {
      vector<vector<unsigned int> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<unsigned int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<unsigned int> >", -2, "vector", 389,
                  typeid(vector<vector<unsigned int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEunsignedsPintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<unsigned int> >) );
      instance.SetNew(&new_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<unsigned int> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<unsigned int> >","std::vector<std::vector<unsigned int, std::allocator<unsigned int> >, std::allocator<std::vector<unsigned int, std::allocator<unsigned int> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned int> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEunsignedsPintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned int> >*>(nullptr))->GetClass();
      vectorlEvectorlEunsignedsPintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEunsignedsPintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned int> > : new vector<vector<unsigned int> >;
   }
   static void *newArray_vectorlEvectorlEunsignedsPintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned int> >[nElements] : new vector<vector<unsigned int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      delete (static_cast<vector<vector<unsigned int> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<unsigned int> >*>(p));
   }
   static void destruct_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      typedef vector<vector<unsigned int> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<unsigned int> >

namespace ROOT {
   static TClass *vectorlEvectorlEunsignedsPchargRsPgR_Dictionary();
   static void vectorlEvectorlEunsignedsPchargRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEunsignedsPchargRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEunsignedsPchargRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEunsignedsPchargRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEunsignedsPchargRsPgR(void *p);
   static void destruct_vectorlEvectorlEunsignedsPchargRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<unsigned char> >*)
   {
      vector<vector<unsigned char> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<unsigned char> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<unsigned char> >", -2, "vector", 389,
                  typeid(vector<vector<unsigned char> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEunsignedsPchargRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<unsigned char> >) );
      instance.SetNew(&new_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<unsigned char> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<unsigned char> >","std::vector<std::vector<unsigned char, std::allocator<unsigned char> >, std::allocator<std::vector<unsigned char, std::allocator<unsigned char> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned char> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEunsignedsPchargRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned char> >*>(nullptr))->GetClass();
      vectorlEvectorlEunsignedsPchargRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEunsignedsPchargRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEunsignedsPchargRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned char> > : new vector<vector<unsigned char> >;
   }
   static void *newArray_vectorlEvectorlEunsignedsPchargRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned char> >[nElements] : new vector<vector<unsigned char> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEunsignedsPchargRsPgR(void *p) {
      delete (static_cast<vector<vector<unsigned char> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEunsignedsPchargRsPgR(void *p) {
      delete [] (static_cast<vector<vector<unsigned char> >*>(p));
   }
   static void destruct_vectorlEvectorlEunsignedsPchargRsPgR(void *p) {
      typedef vector<vector<unsigned char> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<unsigned char> >

namespace ROOT {
   static TClass *vectorlEvectorlEshortgRsPgR_Dictionary();
   static void vectorlEvectorlEshortgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEshortgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEshortgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEshortgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEshortgRsPgR(void *p);
   static void destruct_vectorlEvectorlEshortgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<short> >*)
   {
      vector<vector<short> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<short> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<short> >", -2, "vector", 389,
                  typeid(vector<vector<short> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEshortgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<short> >) );
      instance.SetNew(&new_vectorlEvectorlEshortgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEshortgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEshortgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEshortgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEshortgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<short> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<short> >","std::vector<std::vector<short, std::allocator<short> >, std::allocator<std::vector<short, std::allocator<short> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<short> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEshortgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<short> >*>(nullptr))->GetClass();
      vectorlEvectorlEshortgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEshortgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEshortgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<short> > : new vector<vector<short> >;
   }
   static void *newArray_vectorlEvectorlEshortgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<short> >[nElements] : new vector<vector<short> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEshortgRsPgR(void *p) {
      delete (static_cast<vector<vector<short> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEshortgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<short> >*>(p));
   }
   static void destruct_vectorlEvectorlEshortgRsPgR(void *p) {
      typedef vector<vector<short> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<short> >

namespace ROOT {
   static TClass *vectorlEvectorlElonggRsPgR_Dictionary();
   static void vectorlEvectorlElonggRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlElonggRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlElonggRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlElonggRsPgR(void *p);
   static void deleteArray_vectorlEvectorlElonggRsPgR(void *p);
   static void destruct_vectorlEvectorlElonggRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<long> >*)
   {
      vector<vector<long> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<long> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<long> >", -2, "vector", 389,
                  typeid(vector<vector<long> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlElonggRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<long> >) );
      instance.SetNew(&new_vectorlEvectorlElonggRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlElonggRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlElonggRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlElonggRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlElonggRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<long> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<long> >","std::vector<std::vector<long, std::allocator<long> >, std::allocator<std::vector<long, std::allocator<long> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<long> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlElonggRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<long> >*>(nullptr))->GetClass();
      vectorlEvectorlElonggRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlElonggRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlElonggRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<long> > : new vector<vector<long> >;
   }
   static void *newArray_vectorlEvectorlElonggRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<long> >[nElements] : new vector<vector<long> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlElonggRsPgR(void *p) {
      delete (static_cast<vector<vector<long> >*>(p));
   }
   static void deleteArray_vectorlEvectorlElonggRsPgR(void *p) {
      delete [] (static_cast<vector<vector<long> >*>(p));
   }
   static void destruct_vectorlEvectorlElonggRsPgR(void *p) {
      typedef vector<vector<long> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<long> >

namespace ROOT {
   static TClass *vectorlEvectorlEintgRsPgR_Dictionary();
   static void vectorlEvectorlEintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEintgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEintgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEintgRsPgR(void *p);
   static void destruct_vectorlEvectorlEintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<int> >*)
   {
      vector<vector<int> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<int> >", -2, "vector", 389,
                  typeid(vector<vector<int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<int> >) );
      instance.SetNew(&new_vectorlEvectorlEintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEintgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<int> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<int> >","std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<int> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<int> >*>(nullptr))->GetClass();
      vectorlEvectorlEintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEintgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<int> > : new vector<vector<int> >;
   }
   static void *newArray_vectorlEvectorlEintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<int> >[nElements] : new vector<vector<int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEintgRsPgR(void *p) {
      delete (static_cast<vector<vector<int> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEintgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<int> >*>(p));
   }
   static void destruct_vectorlEvectorlEintgRsPgR(void *p) {
      typedef vector<vector<int> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<int> >

namespace ROOT {
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary();
   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p);
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<double> >*)
   {
      vector<vector<double> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<double> >", -2, "vector", 389,
                  typeid(vector<vector<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEdoublegRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<double> >) );
      instance.SetNew(&new_vectorlEvectorlEdoublegRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<double> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<double> >","std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<double> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<double> >*>(nullptr))->GetClass();
      vectorlEvectorlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<double> > : new vector<vector<double> >;
   }
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<double> >[nElements] : new vector<vector<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete (static_cast<vector<vector<double> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete [] (static_cast<vector<vector<double> >*>(p));
   }
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p) {
      typedef vector<vector<double> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<double> >

namespace ROOT {
   static TClass *vectorlEvectorlEconstsPcharmUgRsPgR_Dictionary();
   static void vectorlEvectorlEconstsPcharmUgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEconstsPcharmUgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEconstsPcharmUgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEconstsPcharmUgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEconstsPcharmUgRsPgR(void *p);
   static void destruct_vectorlEvectorlEconstsPcharmUgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<const char*> >*)
   {
      vector<vector<const char*> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<const char*> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<const char*> >", -2, "vector", 389,
                  typeid(vector<vector<const char*> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEconstsPcharmUgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<const char*> >) );
      instance.SetNew(&new_vectorlEvectorlEconstsPcharmUgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEconstsPcharmUgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEconstsPcharmUgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEconstsPcharmUgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEconstsPcharmUgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<const char*> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<const char*> >","std::vector<std::vector<char const*, std::allocator<char const*> >, std::allocator<std::vector<char const*, std::allocator<char const*> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<const char*> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEconstsPcharmUgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<const char*> >*>(nullptr))->GetClass();
      vectorlEvectorlEconstsPcharmUgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEconstsPcharmUgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEconstsPcharmUgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<const char*> > : new vector<vector<const char*> >;
   }
   static void *newArray_vectorlEvectorlEconstsPcharmUgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<const char*> >[nElements] : new vector<vector<const char*> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEconstsPcharmUgRsPgR(void *p) {
      delete (static_cast<vector<vector<const char*> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEconstsPcharmUgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<const char*> >*>(p));
   }
   static void destruct_vectorlEvectorlEconstsPcharmUgRsPgR(void *p) {
      typedef vector<vector<const char*> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<const char*> >

namespace ROOT {
   static TClass *vectorlEvectorlEchargRsPgR_Dictionary();
   static void vectorlEvectorlEchargRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEchargRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEchargRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEchargRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEchargRsPgR(void *p);
   static void destruct_vectorlEvectorlEchargRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<char> >*)
   {
      vector<vector<char> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<char> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<char> >", -2, "vector", 389,
                  typeid(vector<vector<char> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEchargRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<char> >) );
      instance.SetNew(&new_vectorlEvectorlEchargRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEchargRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEchargRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEchargRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEchargRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<char> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<char> >","std::vector<std::vector<char, std::allocator<char> >, std::allocator<std::vector<char, std::allocator<char> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<char> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEchargRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<char> >*>(nullptr))->GetClass();
      vectorlEvectorlEchargRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEchargRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEchargRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<char> > : new vector<vector<char> >;
   }
   static void *newArray_vectorlEvectorlEchargRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<char> >[nElements] : new vector<vector<char> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEchargRsPgR(void *p) {
      delete (static_cast<vector<vector<char> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEchargRsPgR(void *p) {
      delete [] (static_cast<vector<vector<char> >*>(p));
   }
   static void destruct_vectorlEvectorlEchargRsPgR(void *p) {
      typedef vector<vector<char> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<char> >

namespace ROOT {
   static TClass *vectorlEvectorlEcharmUgRsPgR_Dictionary();
   static void vectorlEvectorlEcharmUgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEcharmUgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEcharmUgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEcharmUgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEcharmUgRsPgR(void *p);
   static void destruct_vectorlEvectorlEcharmUgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<char*> >*)
   {
      vector<vector<char*> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<char*> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<char*> >", -2, "vector", 389,
                  typeid(vector<vector<char*> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEcharmUgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<char*> >) );
      instance.SetNew(&new_vectorlEvectorlEcharmUgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEcharmUgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEcharmUgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEcharmUgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEcharmUgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<char*> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<char*> >","std::vector<std::vector<char*, std::allocator<char*> >, std::allocator<std::vector<char*, std::allocator<char*> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<char*> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEcharmUgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<char*> >*>(nullptr))->GetClass();
      vectorlEvectorlEcharmUgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEcharmUgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEcharmUgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<char*> > : new vector<vector<char*> >;
   }
   static void *newArray_vectorlEvectorlEcharmUgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<char*> >[nElements] : new vector<vector<char*> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEcharmUgRsPgR(void *p) {
      delete (static_cast<vector<vector<char*> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEcharmUgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<char*> >*>(p));
   }
   static void destruct_vectorlEvectorlEcharmUgRsPgR(void *p) {
      typedef vector<vector<char*> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<char*> >

namespace ROOT {
   static TClass *vectorlEvectorlEboolgRsPgR_Dictionary();
   static void vectorlEvectorlEboolgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEboolgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEboolgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEboolgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEboolgRsPgR(void *p);
   static void destruct_vectorlEvectorlEboolgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<bool> >*)
   {
      vector<vector<bool> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<bool> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<bool> >", -2, "vector", 389,
                  typeid(vector<vector<bool> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEboolgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<bool> >) );
      instance.SetNew(&new_vectorlEvectorlEboolgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEboolgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEboolgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEboolgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEboolgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<bool> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<bool> >","std::vector<std::vector<bool, std::allocator<bool> >, std::allocator<std::vector<bool, std::allocator<bool> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<bool> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEboolgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<bool> >*>(nullptr))->GetClass();
      vectorlEvectorlEboolgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEboolgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEboolgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<bool> > : new vector<vector<bool> >;
   }
   static void *newArray_vectorlEvectorlEboolgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<bool> >[nElements] : new vector<vector<bool> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEboolgRsPgR(void *p) {
      delete (static_cast<vector<vector<bool> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEboolgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<bool> >*>(p));
   }
   static void destruct_vectorlEvectorlEboolgRsPgR(void *p) {
      typedef vector<vector<bool> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<bool> >

namespace ROOT {
   static TClass *vectorlEvectorlEULong64_tgRsPgR_Dictionary();
   static void vectorlEvectorlEULong64_tgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEULong64_tgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEULong64_tgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEULong64_tgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEULong64_tgRsPgR(void *p);
   static void destruct_vectorlEvectorlEULong64_tgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<ULong64_t> >*)
   {
      vector<vector<ULong64_t> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<ULong64_t> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<ULong64_t> >", -2, "vector", 389,
                  typeid(vector<vector<ULong64_t> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEULong64_tgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<ULong64_t> >) );
      instance.SetNew(&new_vectorlEvectorlEULong64_tgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEULong64_tgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEULong64_tgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEULong64_tgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEULong64_tgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<ULong64_t> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<ULong64_t> >","std::vector<std::vector<unsigned long long, std::allocator<unsigned long long> >, std::allocator<std::vector<unsigned long long, std::allocator<unsigned long long> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<ULong64_t> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEULong64_tgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<ULong64_t> >*>(nullptr))->GetClass();
      vectorlEvectorlEULong64_tgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEULong64_tgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEULong64_tgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<ULong64_t> > : new vector<vector<ULong64_t> >;
   }
   static void *newArray_vectorlEvectorlEULong64_tgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<ULong64_t> >[nElements] : new vector<vector<ULong64_t> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEULong64_tgRsPgR(void *p) {
      delete (static_cast<vector<vector<ULong64_t> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEULong64_tgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<ULong64_t> >*>(p));
   }
   static void destruct_vectorlEvectorlEULong64_tgRsPgR(void *p) {
      typedef vector<vector<ULong64_t> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<ULong64_t> >

namespace ROOT {
   static TClass *vectorlEvectorlELong64_tgRsPgR_Dictionary();
   static void vectorlEvectorlELong64_tgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlELong64_tgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlELong64_tgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlELong64_tgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlELong64_tgRsPgR(void *p);
   static void destruct_vectorlEvectorlELong64_tgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<Long64_t> >*)
   {
      vector<vector<Long64_t> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<Long64_t> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<Long64_t> >", -2, "vector", 389,
                  typeid(vector<vector<Long64_t> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlELong64_tgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<Long64_t> >) );
      instance.SetNew(&new_vectorlEvectorlELong64_tgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlELong64_tgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlELong64_tgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlELong64_tgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlELong64_tgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<Long64_t> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<Long64_t> >","std::vector<std::vector<long long, std::allocator<long long> >, std::allocator<std::vector<long long, std::allocator<long long> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<Long64_t> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlELong64_tgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<Long64_t> >*>(nullptr))->GetClass();
      vectorlEvectorlELong64_tgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlELong64_tgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlELong64_tgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<Long64_t> > : new vector<vector<Long64_t> >;
   }
   static void *newArray_vectorlEvectorlELong64_tgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<Long64_t> >[nElements] : new vector<vector<Long64_t> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlELong64_tgRsPgR(void *p) {
      delete (static_cast<vector<vector<Long64_t> >*>(p));
   }
   static void deleteArray_vectorlEvectorlELong64_tgRsPgR(void *p) {
      delete [] (static_cast<vector<vector<Long64_t> >*>(p));
   }
   static void destruct_vectorlEvectorlELong64_tgRsPgR(void *p) {
      typedef vector<vector<Long64_t> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<Long64_t> >

namespace ROOT {
   static TClass *vectorlEunsignedsPshortgR_Dictionary();
   static void vectorlEunsignedsPshortgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPshortgR(void *p = nullptr);
   static void *newArray_vectorlEunsignedsPshortgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPshortgR(void *p);
   static void deleteArray_vectorlEunsignedsPshortgR(void *p);
   static void destruct_vectorlEunsignedsPshortgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned short>*)
   {
      vector<unsigned short> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned short>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned short>", -2, "vector", 389,
                  typeid(vector<unsigned short>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPshortgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<unsigned short>) );
      instance.SetNew(&new_vectorlEunsignedsPshortgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPshortgR);
      instance.SetDelete(&delete_vectorlEunsignedsPshortgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPshortgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPshortgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned short> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<unsigned short>","std::vector<unsigned short, std::allocator<unsigned short> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<unsigned short>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPshortgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<unsigned short>*>(nullptr))->GetClass();
      vectorlEunsignedsPshortgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPshortgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPshortgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned short> : new vector<unsigned short>;
   }
   static void *newArray_vectorlEunsignedsPshortgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned short>[nElements] : new vector<unsigned short>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPshortgR(void *p) {
      delete (static_cast<vector<unsigned short>*>(p));
   }
   static void deleteArray_vectorlEunsignedsPshortgR(void *p) {
      delete [] (static_cast<vector<unsigned short>*>(p));
   }
   static void destruct_vectorlEunsignedsPshortgR(void *p) {
      typedef vector<unsigned short> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned short>

namespace ROOT {
   static TClass *vectorlEunsignedsPlonggR_Dictionary();
   static void vectorlEunsignedsPlonggR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPlonggR(void *p = nullptr);
   static void *newArray_vectorlEunsignedsPlonggR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPlonggR(void *p);
   static void deleteArray_vectorlEunsignedsPlonggR(void *p);
   static void destruct_vectorlEunsignedsPlonggR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned long>*)
   {
      vector<unsigned long> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned long>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned long>", -2, "vector", 389,
                  typeid(vector<unsigned long>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPlonggR_Dictionary, isa_proxy, 4,
                  sizeof(vector<unsigned long>) );
      instance.SetNew(&new_vectorlEunsignedsPlonggR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPlonggR);
      instance.SetDelete(&delete_vectorlEunsignedsPlonggR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPlonggR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPlonggR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned long> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<unsigned long>","std::vector<unsigned long, std::allocator<unsigned long> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<unsigned long>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPlonggR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<unsigned long>*>(nullptr))->GetClass();
      vectorlEunsignedsPlonggR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPlonggR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPlonggR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned long> : new vector<unsigned long>;
   }
   static void *newArray_vectorlEunsignedsPlonggR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned long>[nElements] : new vector<unsigned long>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPlonggR(void *p) {
      delete (static_cast<vector<unsigned long>*>(p));
   }
   static void deleteArray_vectorlEunsignedsPlonggR(void *p) {
      delete [] (static_cast<vector<unsigned long>*>(p));
   }
   static void destruct_vectorlEunsignedsPlonggR(void *p) {
      typedef vector<unsigned long> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned long>

namespace ROOT {
   static TClass *vectorlEunsignedsPintgR_Dictionary();
   static void vectorlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPintgR(void *p = nullptr);
   static void *newArray_vectorlEunsignedsPintgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPintgR(void *p);
   static void deleteArray_vectorlEunsignedsPintgR(void *p);
   static void destruct_vectorlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned int>*)
   {
      vector<unsigned int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned int>", -2, "vector", 389,
                  typeid(vector<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<unsigned int>) );
      instance.SetNew(&new_vectorlEunsignedsPintgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPintgR);
      instance.SetDelete(&delete_vectorlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPintgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned int> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<unsigned int>","std::vector<unsigned int, std::allocator<unsigned int> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<unsigned int>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<unsigned int>*>(nullptr))->GetClass();
      vectorlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPintgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned int> : new vector<unsigned int>;
   }
   static void *newArray_vectorlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned int>[nElements] : new vector<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPintgR(void *p) {
      delete (static_cast<vector<unsigned int>*>(p));
   }
   static void deleteArray_vectorlEunsignedsPintgR(void *p) {
      delete [] (static_cast<vector<unsigned int>*>(p));
   }
   static void destruct_vectorlEunsignedsPintgR(void *p) {
      typedef vector<unsigned int> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned int>

namespace ROOT {
   static TClass *vectorlEunsignedsPchargR_Dictionary();
   static void vectorlEunsignedsPchargR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPchargR(void *p = nullptr);
   static void *newArray_vectorlEunsignedsPchargR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPchargR(void *p);
   static void deleteArray_vectorlEunsignedsPchargR(void *p);
   static void destruct_vectorlEunsignedsPchargR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned char>*)
   {
      vector<unsigned char> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned char>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned char>", -2, "vector", 389,
                  typeid(vector<unsigned char>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPchargR_Dictionary, isa_proxy, 4,
                  sizeof(vector<unsigned char>) );
      instance.SetNew(&new_vectorlEunsignedsPchargR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPchargR);
      instance.SetDelete(&delete_vectorlEunsignedsPchargR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPchargR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPchargR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned char> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<unsigned char>","std::vector<unsigned char, std::allocator<unsigned char> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<unsigned char>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPchargR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<unsigned char>*>(nullptr))->GetClass();
      vectorlEunsignedsPchargR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPchargR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPchargR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned char> : new vector<unsigned char>;
   }
   static void *newArray_vectorlEunsignedsPchargR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned char>[nElements] : new vector<unsigned char>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPchargR(void *p) {
      delete (static_cast<vector<unsigned char>*>(p));
   }
   static void deleteArray_vectorlEunsignedsPchargR(void *p) {
      delete [] (static_cast<vector<unsigned char>*>(p));
   }
   static void destruct_vectorlEunsignedsPchargR(void *p) {
      typedef vector<unsigned char> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned char>

namespace ROOT {
   static TClass *vectorlEshortgR_Dictionary();
   static void vectorlEshortgR_TClassManip(TClass*);
   static void *new_vectorlEshortgR(void *p = nullptr);
   static void *newArray_vectorlEshortgR(Long_t size, void *p);
   static void delete_vectorlEshortgR(void *p);
   static void deleteArray_vectorlEshortgR(void *p);
   static void destruct_vectorlEshortgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<short>*)
   {
      vector<short> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<short>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<short>", -2, "vector", 389,
                  typeid(vector<short>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEshortgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<short>) );
      instance.SetNew(&new_vectorlEshortgR);
      instance.SetNewArray(&newArray_vectorlEshortgR);
      instance.SetDelete(&delete_vectorlEshortgR);
      instance.SetDeleteArray(&deleteArray_vectorlEshortgR);
      instance.SetDestructor(&destruct_vectorlEshortgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<short> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<short>","std::vector<short, std::allocator<short> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<short>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEshortgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<short>*>(nullptr))->GetClass();
      vectorlEshortgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEshortgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEshortgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<short> : new vector<short>;
   }
   static void *newArray_vectorlEshortgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<short>[nElements] : new vector<short>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEshortgR(void *p) {
      delete (static_cast<vector<short>*>(p));
   }
   static void deleteArray_vectorlEshortgR(void *p) {
      delete [] (static_cast<vector<short>*>(p));
   }
   static void destruct_vectorlEshortgR(void *p) {
      typedef vector<short> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<short>

namespace ROOT {
   static TClass *vectorlElonggR_Dictionary();
   static void vectorlElonggR_TClassManip(TClass*);
   static void *new_vectorlElonggR(void *p = nullptr);
   static void *newArray_vectorlElonggR(Long_t size, void *p);
   static void delete_vectorlElonggR(void *p);
   static void deleteArray_vectorlElonggR(void *p);
   static void destruct_vectorlElonggR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<long>*)
   {
      vector<long> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<long>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<long>", -2, "vector", 389,
                  typeid(vector<long>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElonggR_Dictionary, isa_proxy, 4,
                  sizeof(vector<long>) );
      instance.SetNew(&new_vectorlElonggR);
      instance.SetNewArray(&newArray_vectorlElonggR);
      instance.SetDelete(&delete_vectorlElonggR);
      instance.SetDeleteArray(&deleteArray_vectorlElonggR);
      instance.SetDestructor(&destruct_vectorlElonggR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<long> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<long>","std::vector<long, std::allocator<long> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<long>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElonggR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<long>*>(nullptr))->GetClass();
      vectorlElonggR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElonggR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElonggR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<long> : new vector<long>;
   }
   static void *newArray_vectorlElonggR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<long>[nElements] : new vector<long>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElonggR(void *p) {
      delete (static_cast<vector<long>*>(p));
   }
   static void deleteArray_vectorlElonggR(void *p) {
      delete [] (static_cast<vector<long>*>(p));
   }
   static void destruct_vectorlElonggR(void *p) {
      typedef vector<long> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<long>

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = nullptr);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 389,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr))->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete (static_cast<vector<int>*>(p));
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] (static_cast<vector<int>*>(p));
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 389,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlEconstsPcharmUgR_Dictionary();
   static void vectorlEconstsPcharmUgR_TClassManip(TClass*);
   static void *new_vectorlEconstsPcharmUgR(void *p = nullptr);
   static void *newArray_vectorlEconstsPcharmUgR(Long_t size, void *p);
   static void delete_vectorlEconstsPcharmUgR(void *p);
   static void deleteArray_vectorlEconstsPcharmUgR(void *p);
   static void destruct_vectorlEconstsPcharmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<const char*>*)
   {
      vector<const char*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<const char*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<const char*>", -2, "vector", 389,
                  typeid(vector<const char*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEconstsPcharmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<const char*>) );
      instance.SetNew(&new_vectorlEconstsPcharmUgR);
      instance.SetNewArray(&newArray_vectorlEconstsPcharmUgR);
      instance.SetDelete(&delete_vectorlEconstsPcharmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEconstsPcharmUgR);
      instance.SetDestructor(&destruct_vectorlEconstsPcharmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<const char*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<const char*>","std::vector<char const*, std::allocator<char const*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<const char*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEconstsPcharmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<const char*>*>(nullptr))->GetClass();
      vectorlEconstsPcharmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEconstsPcharmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEconstsPcharmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<const char*> : new vector<const char*>;
   }
   static void *newArray_vectorlEconstsPcharmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<const char*>[nElements] : new vector<const char*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEconstsPcharmUgR(void *p) {
      delete (static_cast<vector<const char*>*>(p));
   }
   static void deleteArray_vectorlEconstsPcharmUgR(void *p) {
      delete [] (static_cast<vector<const char*>*>(p));
   }
   static void destruct_vectorlEconstsPcharmUgR(void *p) {
      typedef vector<const char*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<const char*>

namespace ROOT {
   static TClass *vectorlEchargR_Dictionary();
   static void vectorlEchargR_TClassManip(TClass*);
   static void *new_vectorlEchargR(void *p = nullptr);
   static void *newArray_vectorlEchargR(Long_t size, void *p);
   static void delete_vectorlEchargR(void *p);
   static void deleteArray_vectorlEchargR(void *p);
   static void destruct_vectorlEchargR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<char>*)
   {
      vector<char> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<char>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<char>", -2, "vector", 389,
                  typeid(vector<char>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEchargR_Dictionary, isa_proxy, 4,
                  sizeof(vector<char>) );
      instance.SetNew(&new_vectorlEchargR);
      instance.SetNewArray(&newArray_vectorlEchargR);
      instance.SetDelete(&delete_vectorlEchargR);
      instance.SetDeleteArray(&deleteArray_vectorlEchargR);
      instance.SetDestructor(&destruct_vectorlEchargR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<char> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<char>","std::vector<char, std::allocator<char> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<char>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEchargR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<char>*>(nullptr))->GetClass();
      vectorlEchargR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEchargR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEchargR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<char> : new vector<char>;
   }
   static void *newArray_vectorlEchargR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<char>[nElements] : new vector<char>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEchargR(void *p) {
      delete (static_cast<vector<char>*>(p));
   }
   static void deleteArray_vectorlEchargR(void *p) {
      delete [] (static_cast<vector<char>*>(p));
   }
   static void destruct_vectorlEchargR(void *p) {
      typedef vector<char> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<char>

namespace ROOT {
   static TClass *vectorlEcharmUgR_Dictionary();
   static void vectorlEcharmUgR_TClassManip(TClass*);
   static void *new_vectorlEcharmUgR(void *p = nullptr);
   static void *newArray_vectorlEcharmUgR(Long_t size, void *p);
   static void delete_vectorlEcharmUgR(void *p);
   static void deleteArray_vectorlEcharmUgR(void *p);
   static void destruct_vectorlEcharmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<char*>*)
   {
      vector<char*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<char*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<char*>", -2, "vector", 389,
                  typeid(vector<char*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEcharmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<char*>) );
      instance.SetNew(&new_vectorlEcharmUgR);
      instance.SetNewArray(&newArray_vectorlEcharmUgR);
      instance.SetDelete(&delete_vectorlEcharmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEcharmUgR);
      instance.SetDestructor(&destruct_vectorlEcharmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<char*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<char*>","std::vector<char*, std::allocator<char*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<char*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEcharmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<char*>*>(nullptr))->GetClass();
      vectorlEcharmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEcharmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEcharmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<char*> : new vector<char*>;
   }
   static void *newArray_vectorlEcharmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<char*>[nElements] : new vector<char*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEcharmUgR(void *p) {
      delete (static_cast<vector<char*>*>(p));
   }
   static void deleteArray_vectorlEcharmUgR(void *p) {
      delete [] (static_cast<vector<char*>*>(p));
   }
   static void destruct_vectorlEcharmUgR(void *p) {
      typedef vector<char*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<char*>

namespace ROOT {
   static TClass *vectorlEboolgR_Dictionary();
   static void vectorlEboolgR_TClassManip(TClass*);
   static void *new_vectorlEboolgR(void *p = nullptr);
   static void *newArray_vectorlEboolgR(Long_t size, void *p);
   static void delete_vectorlEboolgR(void *p);
   static void deleteArray_vectorlEboolgR(void *p);
   static void destruct_vectorlEboolgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<bool>*)
   {
      vector<bool> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<bool>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<bool>", -2, "vector", 596,
                  typeid(vector<bool>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEboolgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<bool>) );
      instance.SetNew(&new_vectorlEboolgR);
      instance.SetNewArray(&newArray_vectorlEboolgR);
      instance.SetDelete(&delete_vectorlEboolgR);
      instance.SetDeleteArray(&deleteArray_vectorlEboolgR);
      instance.SetDestructor(&destruct_vectorlEboolgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<bool> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<bool>","std::vector<bool, std::allocator<bool> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<bool>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEboolgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<bool>*>(nullptr))->GetClass();
      vectorlEboolgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEboolgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEboolgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<bool> : new vector<bool>;
   }
   static void *newArray_vectorlEboolgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<bool>[nElements] : new vector<bool>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEboolgR(void *p) {
      delete (static_cast<vector<bool>*>(p));
   }
   static void deleteArray_vectorlEboolgR(void *p) {
      delete [] (static_cast<vector<bool>*>(p));
   }
   static void destruct_vectorlEboolgR(void *p) {
      typedef vector<bool> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<bool>

namespace ROOT {
   static TClass *vectorlEULong64_tgR_Dictionary();
   static void vectorlEULong64_tgR_TClassManip(TClass*);
   static void *new_vectorlEULong64_tgR(void *p = nullptr);
   static void *newArray_vectorlEULong64_tgR(Long_t size, void *p);
   static void delete_vectorlEULong64_tgR(void *p);
   static void deleteArray_vectorlEULong64_tgR(void *p);
   static void destruct_vectorlEULong64_tgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ULong64_t>*)
   {
      vector<ULong64_t> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ULong64_t>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ULong64_t>", -2, "vector", 389,
                  typeid(vector<ULong64_t>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEULong64_tgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ULong64_t>) );
      instance.SetNew(&new_vectorlEULong64_tgR);
      instance.SetNewArray(&newArray_vectorlEULong64_tgR);
      instance.SetDelete(&delete_vectorlEULong64_tgR);
      instance.SetDeleteArray(&deleteArray_vectorlEULong64_tgR);
      instance.SetDestructor(&destruct_vectorlEULong64_tgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ULong64_t> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ULong64_t>","std::vector<unsigned long long, std::allocator<unsigned long long> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ULong64_t>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEULong64_tgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ULong64_t>*>(nullptr))->GetClass();
      vectorlEULong64_tgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEULong64_tgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEULong64_tgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ULong64_t> : new vector<ULong64_t>;
   }
   static void *newArray_vectorlEULong64_tgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ULong64_t>[nElements] : new vector<ULong64_t>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEULong64_tgR(void *p) {
      delete (static_cast<vector<ULong64_t>*>(p));
   }
   static void deleteArray_vectorlEULong64_tgR(void *p) {
      delete [] (static_cast<vector<ULong64_t>*>(p));
   }
   static void destruct_vectorlEULong64_tgR(void *p) {
      typedef vector<ULong64_t> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ULong64_t>

namespace ROOT {
   static TClass *vectorlETVector3gR_Dictionary();
   static void vectorlETVector3gR_TClassManip(TClass*);
   static void *new_vectorlETVector3gR(void *p = nullptr);
   static void *newArray_vectorlETVector3gR(Long_t size, void *p);
   static void delete_vectorlETVector3gR(void *p);
   static void deleteArray_vectorlETVector3gR(void *p);
   static void destruct_vectorlETVector3gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TVector3>*)
   {
      vector<TVector3> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TVector3>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TVector3>", -2, "vector", 389,
                  typeid(vector<TVector3>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETVector3gR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TVector3>) );
      instance.SetNew(&new_vectorlETVector3gR);
      instance.SetNewArray(&newArray_vectorlETVector3gR);
      instance.SetDelete(&delete_vectorlETVector3gR);
      instance.SetDeleteArray(&deleteArray_vectorlETVector3gR);
      instance.SetDestructor(&destruct_vectorlETVector3gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TVector3> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TVector3>","std::vector<TVector3, std::allocator<TVector3> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TVector3>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETVector3gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TVector3>*>(nullptr))->GetClass();
      vectorlETVector3gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETVector3gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETVector3gR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TVector3> : new vector<TVector3>;
   }
   static void *newArray_vectorlETVector3gR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TVector3>[nElements] : new vector<TVector3>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETVector3gR(void *p) {
      delete (static_cast<vector<TVector3>*>(p));
   }
   static void deleteArray_vectorlETVector3gR(void *p) {
      delete [] (static_cast<vector<TVector3>*>(p));
   }
   static void destruct_vectorlETVector3gR(void *p) {
      typedef vector<TVector3> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TVector3>

namespace ROOT {
   static TClass *vectorlETParticlegR_Dictionary();
   static void vectorlETParticlegR_TClassManip(TClass*);
   static void *new_vectorlETParticlegR(void *p = nullptr);
   static void *newArray_vectorlETParticlegR(Long_t size, void *p);
   static void delete_vectorlETParticlegR(void *p);
   static void deleteArray_vectorlETParticlegR(void *p);
   static void destruct_vectorlETParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TParticle>*)
   {
      vector<TParticle> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TParticle>", -2, "vector", 389,
                  typeid(vector<TParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TParticle>) );
      instance.SetNew(&new_vectorlETParticlegR);
      instance.SetNewArray(&newArray_vectorlETParticlegR);
      instance.SetDelete(&delete_vectorlETParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlETParticlegR);
      instance.SetDestructor(&destruct_vectorlETParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TParticle> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TParticle>","std::vector<TParticle, std::allocator<TParticle> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TParticle>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TParticle>*>(nullptr))->GetClass();
      vectorlETParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETParticlegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TParticle> : new vector<TParticle>;
   }
   static void *newArray_vectorlETParticlegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TParticle>[nElements] : new vector<TParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETParticlegR(void *p) {
      delete (static_cast<vector<TParticle>*>(p));
   }
   static void deleteArray_vectorlETParticlegR(void *p) {
      delete [] (static_cast<vector<TParticle>*>(p));
   }
   static void destruct_vectorlETParticlegR(void *p) {
      typedef vector<TParticle> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TParticle>

namespace ROOT {
   static TClass *vectorlELong64_tgR_Dictionary();
   static void vectorlELong64_tgR_TClassManip(TClass*);
   static void *new_vectorlELong64_tgR(void *p = nullptr);
   static void *newArray_vectorlELong64_tgR(Long_t size, void *p);
   static void delete_vectorlELong64_tgR(void *p);
   static void deleteArray_vectorlELong64_tgR(void *p);
   static void destruct_vectorlELong64_tgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Long64_t>*)
   {
      vector<Long64_t> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Long64_t>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Long64_t>", -2, "vector", 389,
                  typeid(vector<Long64_t>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELong64_tgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<Long64_t>) );
      instance.SetNew(&new_vectorlELong64_tgR);
      instance.SetNewArray(&newArray_vectorlELong64_tgR);
      instance.SetDelete(&delete_vectorlELong64_tgR);
      instance.SetDeleteArray(&deleteArray_vectorlELong64_tgR);
      instance.SetDestructor(&destruct_vectorlELong64_tgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Long64_t> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Long64_t>","std::vector<long long, std::allocator<long long> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Long64_t>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELong64_tgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Long64_t>*>(nullptr))->GetClass();
      vectorlELong64_tgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELong64_tgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELong64_tgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Long64_t> : new vector<Long64_t>;
   }
   static void *newArray_vectorlELong64_tgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Long64_t>[nElements] : new vector<Long64_t>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELong64_tgR(void *p) {
      delete (static_cast<vector<Long64_t>*>(p));
   }
   static void deleteArray_vectorlELong64_tgR(void *p) {
      delete [] (static_cast<vector<Long64_t>*>(p));
   }
   static void destruct_vectorlELong64_tgR(void *p) {
      typedef vector<Long64_t> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Long64_t>

namespace {
  void TriggerDictionaryInitialization_libG4HypTPC_Impl() {
    static const char* headers[] = {
"TParticle.h",
"TVector3.h",
nullptr
    };
    static const char* includePaths[] = {
"/gpfs/group/had/sks/Users/jaejin/E72geant4/k18geant4/hyptpc_E45/include",
"/sw/packages/root/6.32.04/include",
"/sw/packages/geant4/11.2.2/include/Geant4",
"/gpfs/group/had/sks/Users/jaejin/E72geant4/k18geant4/hyptpc_E45",
"/gpfs/group/had/sks/Users/jaejin/E72geant4/k18geant4/hyptpc_E45/include",
"/sw/packages/root/6.32.04/include",
"/sw/packages/geant4/11.2.2/include/Geant4",
"/usr/include",
"/usr/include/qt5/",
"/usr/include/qt5/QtGui",
"/usr/include/qt5/QtCore",
"/usr/lib64/qt5//mkspecs/linux-g++",
"/usr/include/qt5/QtOpenGL",
"/usr/include/qt5/QtWidgets",
"/sw/packages/root/6.32.04/include/",
"/gpfs/group/had/sks/Users/jaejin/E72geant4/k18geant4/hyptpc_E45/.build/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libG4HypTPC dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TVector3.h")))  __attribute__((annotate("$clingAutoload$TParticle.h")))  TVector3;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TParticle.h")))  TParticle;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libG4HypTPC dictionary payload"

#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef PTL_BUILD_DLL
  #define PTL_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4VIS_USE_TOOLSSG_QT_GLES
  #define G4VIS_USE_TOOLSSG_QT_GLES 1
#endif
#ifndef G4VIS_USE_TOOLSSG_QT_ZB
  #define G4VIS_USE_TOOLSSG_QT_ZB 1
#endif
#ifndef G4VIS_USE_TOOLSSG_X11_GLES
  #define G4VIS_USE_TOOLSSG_X11_GLES 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4VIS_USE_OPENGLQT
  #define G4VIS_USE_OPENGLQT 1
#endif
#ifndef G4VIS_USE_OPENGLX
  #define G4VIS_USE_OPENGLX 1
#endif
#ifndef QT_GUI_LIB
  #define QT_GUI_LIB 1
#endif
#ifndef QT_CORE_LIB
  #define QT_CORE_LIB 1
#endif
#ifndef QT_NO_DEBUG
  #define QT_NO_DEBUG 1
#endif
#ifndef QT_OPENGL_LIB
  #define QT_OPENGL_LIB 1
#endif
#ifndef QT_WIDGETS_LIB
  #define QT_WIDGETS_LIB 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4UI_USE_QT
  #define G4UI_USE_QT 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif
#ifndef G4LIB_BUILD_DLL
  #define G4LIB_BUILD_DLL 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TParticle.h"
#include "TVector3.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libG4HypTPC",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libG4HypTPC_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libG4HypTPC_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libG4HypTPC() {
  TriggerDictionaryInitialization_libG4HypTPC_Impl();
}
