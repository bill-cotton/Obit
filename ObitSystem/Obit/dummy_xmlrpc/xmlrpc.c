/*  Dummy version of xmlrpc source file xmlrpc.c                           */
/* Includes dummy versions of all needed xmlrpc routines                   */
/* The dummy version of this library is solely to allow compiling and      */
/* linking but not execution of these functions which are stubbed.         */
/* No claims are made about this software except that is is NOT functional */
/* A proper version may be obtained from http://xmlrpc-c.sourceforge.net   */

#include <stdio.h>
#include <stdlib.h>
#include "xmlrpc.h"
#include "xmlrpc-c/client.h"

#define CRASH_AND_BURN  \
{fprintf(stderr,"%s: xmlrpc not implemented\n",routine); exit(1); } 

void xmlrpc_env_init (xmlrpc_env* env)
{
  /* NOP so the program can pretend */
  env->fault_occurred = 0;
}

void xmlrpc_env_clean (xmlrpc_env* env)
{
  /* NOP so the program can pretend */
}

xmlrpc_value * 
xmlrpc_build_value(xmlrpc_env * const env,
                   const char * const format, 
                   ...)
{
  char *routine = "xmlrpc_build_value";
  CRASH_AND_BURN 
  return NULL;
}

void 
xmlrpc_decompose_value(xmlrpc_env *   const envP,
                       xmlrpc_value * const value,
                       const char *   const format, 
                       ...)
{
  char *routine = "xmlrpc_decompose_value";
  CRASH_AND_BURN 
}

extern void xmlrpc_INCREF (xmlrpc_value* value)
{
  char *routine = "xmlrpc_INCREF";
  CRASH_AND_BURN 
}

extern void xmlrpc_DECREF (xmlrpc_value* value)
{
  char *routine = "xmlrpc_DECREF";
  CRASH_AND_BURN 
}

extern xmlrpc_type xmlrpc_value_type (xmlrpc_value* value)
{
  char *routine = " xmlrpc_type";
  CRASH_AND_BURN 
  return XMLRPC_TYPE_DEAD;
}

xmlrpc_value *
xmlrpc_int_new(xmlrpc_env * const envP,
               int          const intValue)
{
  char *routine = "xmlrpc_int_new";
  CRASH_AND_BURN 
  return NULL;
}

void 
xmlrpc_read_int(xmlrpc_env *         const envP,
                const xmlrpc_value * const valueP,
                int *                const intValueP)
{
  char *routine = "xmlrpc_read_int";
  CRASH_AND_BURN 
}

xmlrpc_value *
xmlrpc_bool_new(xmlrpc_env * const envP,
                xmlrpc_bool  const boolValue)
{
  char *routine = "xmlrpc_bool_new";
  CRASH_AND_BURN 
  return NULL;
}

void
xmlrpc_read_bool(xmlrpc_env *         const envP,
                 const xmlrpc_value * const valueP,
                 xmlrpc_bool *        const boolValueP)
{
  char *routine = "xmlrpc_read_bool";
  CRASH_AND_BURN 
}

xmlrpc_value *
xmlrpc_double_new(xmlrpc_env * const envP,
                  double       const doubleValue)
{
  char *routine = "xmlrpc_double_new";
  CRASH_AND_BURN 
  return NULL;
}

void
xmlrpc_read_double(xmlrpc_env *         const envP,
                   const xmlrpc_value * const valueP,
                   xmlrpc_double *      const doubleValueP)
{
  char *routine = "xmlrpc_read_double";
  CRASH_AND_BURN 
}

xmlrpc_value *
xmlrpc_datetime_new_str(xmlrpc_env * const envP,
                        const char * const value)
{
  char *routine = "xmlrpc_datetime_new_str";
  CRASH_AND_BURN 
  return NULL;
}

void
xmlrpc_read_datetime_str(xmlrpc_env *         const envP,
                         const xmlrpc_value * const valueP,
                         const char **        const stringValueP)
{
  char *routine = "xmlrpc_read_datetime_str";
  CRASH_AND_BURN 
}

xmlrpc_value *
xmlrpc_string_new(xmlrpc_env * const envP,
                  const char * const stringValue)
{
  char *routine = "xmlrpc_string_new";
  CRASH_AND_BURN 
  return NULL;
}

xmlrpc_value *
xmlrpc_string_new_lp(xmlrpc_env * const envP, 
                     size_t       const length,
                     const char * const stringValue)
{
  char *routine = "xmlrpc_string_new_lp";
  CRASH_AND_BURN 
  return NULL;
}

void
xmlrpc_read_string(xmlrpc_env *         const envP,
                   const xmlrpc_value * const valueP,
                   const char **        const stringValueP)
{
  char *routine = "xmlrpc_read_string";
  CRASH_AND_BURN 
}


void
xmlrpc_read_string_lp(xmlrpc_env *         const envP,
                      const xmlrpc_value * const valueP,
                      size_t *             const lengthP,
                      const char **        const stringValueP)
{
  char *routine = "xmlrpc_read_string_lp";
  CRASH_AND_BURN 
}

xmlrpc_value *
xmlrpc_base64_new(xmlrpc_env *          const envP,
                  unsigned int          const length,
                  const unsigned char * const bytestringValue)
{
  char *routine = "xmlrpc_base64_new";
  CRASH_AND_BURN 
  return NULL;
}

void
xmlrpc_read_base64(xmlrpc_env *           const envP,
                   const xmlrpc_value *   const valueP,
                   unsigned int *         const lengthP,
                   const unsigned char ** const bytestringValueP)
{
  char *routine = "xmlrpc_read_base64";
  CRASH_AND_BURN 
}

xmlrpc_value *
xmlrpc_array_new(xmlrpc_env * const envP)
{
  char *routine = "xmlrpc_array_new";
  CRASH_AND_BURN 
  return NULL;
}

extern void
xmlrpc_array_append_item(xmlrpc_env*   env,
			 xmlrpc_value* array,
			 xmlrpc_value* value)
{
  char *routine = "xmlrpc_array_append_item";
  CRASH_AND_BURN 
}

xmlrpc_value * 
xmlrpc_array_get_item(xmlrpc_env *         const envP,
                      const xmlrpc_value * const arrayP,
                      int                  const index)
{
  char *routine = "xmlrpc_array_get_item";
  CRASH_AND_BURN 
  return NULL;
}

void
xmlrpc_array_read_item(xmlrpc_env *         const envP,
                       const xmlrpc_value * const arrayP,
                       unsigned int         const index,
                       xmlrpc_value **      const valuePP)
{
  char *routine = "xmlrpc_array_read_item";
  CRASH_AND_BURN 
}

xmlrpc_value *
xmlrpc_struct_new(xmlrpc_env * env)
{
  char *routine = "xmlrpc_struct_new";
  CRASH_AND_BURN 
  return NULL;
}

int
xmlrpc_struct_size (xmlrpc_env   * env, 
                    xmlrpc_value * strct)
{
  char *routine = "xmlrpc_struct_size";
  CRASH_AND_BURN 
  return 1;
}

void 
xmlrpc_struct_set_value(xmlrpc_env *   const env,
                        xmlrpc_value * const strct,
                        const char *   const key,
                        xmlrpc_value * const value)
{
  char *routine = "xmlrpc_struct_set_value";
  CRASH_AND_BURN 
}

int 
xmlrpc_struct_has_key(xmlrpc_env *   const envP,
                      xmlrpc_value * const strctP,
                      const char *   const key)
{
  char *routine = "xmlrpc_struct_has_key";
  CRASH_AND_BURN 
  return 1;
}

void
xmlrpc_struct_find_value(xmlrpc_env *    const envP,
                         xmlrpc_value *  const structP,
                         const char *    const key,
                         xmlrpc_value ** const valuePP)
{
  char *routine = "xmlrpc_struct_find_value";
  CRASH_AND_BURN 
}

void
xmlrpc_struct_read_value(xmlrpc_env *    const envP,
                         xmlrpc_value *  const strctP,
                         const char *    const key,
                         xmlrpc_value ** const valuePP)
{
  char *routine = "xmlrpc_struct_read_value";
  CRASH_AND_BURN 
}

void 
xmlrpc_struct_read_member(xmlrpc_env *    const envP,
                          xmlrpc_value *  const structP,
                          unsigned int    const index,
                          xmlrpc_value ** const keyvalP,
                          xmlrpc_value ** const valueP)

{
  char *routine = "xmlrpc_struct_read_member";
  CRASH_AND_BURN 
}

#define XMLRPC_CLIENT_NO_FLAGS         (0)
#define XMLRPC_CLIENT_SKIP_LIBWWW_INIT (1)

extern void
xmlrpc_client_init(int          const flags,
                   const char * const appname,
                   const char * const appversion)
{
   /* NOP so the program can pretend */
}

extern void
xmlrpc_client_cleanup(void)
{
  char *routine = "xmlrpc_client_cleanup";
  CRASH_AND_BURN 
}

xmlrpc_value * 
xmlrpc_client_call(xmlrpc_env * const envP,
                   const char * const server_url,
                   const char * const method_name,
                   const char * const format,
                   ...)
{
  char *routine = "xmlrpc_client_call";
  CRASH_AND_BURN 
  return NULL;
}

xmlrpc_value * 
xmlrpc_client_call_params(xmlrpc_env *   const envP,
                          const char *   const serverUrl,
                          const char *   const methodName,
                          xmlrpc_value * const paramArrayP)
{
  char *routine = "xmlrpc_client_call_params";
  CRASH_AND_BURN 
  return NULL;
}

/*void 
  xmlrpc_client_call_asynch(const char * const server_url,
  const char * const method_name,
  xmlrpc_response_handler callback,
  void *       const user_data,
  const char * const format,
  ...)
  {
  char *routine = "xmlrpc_client_call_asynch";
  CRASH_AND_BURN 
  return 1;
  }*/

extern void
xmlrpc_client_event_loop_finish_asynch(void)
{
  char *routine = "xmlrpc_client_event_loop_finish_asynch";
  CRASH_AND_BURN 
}

extern void
xmlrpc_client_event_loop_finish_asynch_timeout(unsigned long milliseconds)
{
  char *routine = "xmlrpc_client_event_loop_finish_asynch_timeout";
  CRASH_AND_BURN 
}

xmlrpc_registry *
xmlrpc_registry_new(xmlrpc_env * env)
{
  /* NOP so the program can pretend */
  return NULL;
}

void
xmlrpc_registry_free(xmlrpc_registry * registry)
{
  char *routine = "xmlrpc_registry_free";
  CRASH_AND_BURN 
}

void
xmlrpc_registry_add_method(xmlrpc_env *      env,
                           xmlrpc_registry * registry,
                           const char *      host,
                           const char *      method_name,
                           xmlrpc_method     method,
                           void *            user_data)
{
  /* NOP so the program can pretend */
}

void
xmlrpc_server_abyss(xmlrpc_env *                      const envP,
                    const xmlrpc_server_abyss_parms * const parms,
                    unsigned int                      const parm_size)

{
  /* NOP so the program can pretend */
}

void
xmlrpc_client_setup_global_const(xmlrpc_env * const envP)
{
  /* NOP so the program can pretend */
}

void
xmlrpc_client_teardown_global_const(void)
{
  /* NOP so the program can pretend */
}

void 
xmlrpc_client_create(xmlrpc_env *                      const envP,
                     int                               const flags,
                     const char *                      const appname,
                     const char *                      const appversion,
                     const struct xmlrpc_clientparms * const clientparmsP,
                     unsigned int                      const parmSize,
                     xmlrpc_client **                  const clientPP)
{
  char *routine = "xmlrpc_client_create";
  CRASH_AND_BURN 
}

void 
xmlrpc_client_destroy(xmlrpc_client * const clientP)
{
  char *routine = "xmlrpc_client_destroy";
  CRASH_AND_BURN 
}

void 
xmlrpc_client_call_asynch(const char * const server_url,
                          const char * const method_name,
                          xmlrpc_response_handler callback,
                          void *       const user_data,
                          const char * const format,
                          ...)
{
  char *routine = "xmlrpc_client_call_asynch";
  CRASH_AND_BURN 
}

void
xmlrpc_client_call2f(xmlrpc_env *    const envP,
                     xmlrpc_client * const clientP,
                     const char *    const serverUrl,
                     const char *    const methodName,
                     xmlrpc_value ** const resultPP,
                     const char *    const format,
                     ...)
{
  char *routine = "xmlrpc_client_call2f";
  CRASH_AND_BURN 
}
