/*  Dummy version of xmlrpc include file xmlrpc.h                          */
/* Includes definitions of all needed dummy xmlrpc routines                */
/* The dummy version of this library is solely to allow compiling and      */
/* linking but not execution of these functions which are stubbed.         */
/* No claims are made about this software except that is is NOT functional */
/* A proper version may be obtained from http://xmlrpc-c.sourceforge.net   */

#ifndef _XMLRPC_H
#define _XMLRPC_H

#include <stddef.h>
#include <stdarg.h>

/* &&&&&&&&&&&&&&&&&& type definitions &&&&&&&&&&&&&&&&&&&&&&&&&&& */
typedef signed int xmlrpc_int;  
typedef signed int xmlrpc_int32;
typedef int        xmlrpc_bool;
typedef double     xmlrpc_double;

typedef struct _xmlrpc_env {
    int          fault_occurred;
    xmlrpc_int32 fault_code;
    char*        fault_string;
} xmlrpc_env;

typedef enum {
    XMLRPC_TYPE_INT      = 0,
    XMLRPC_TYPE_BOOL     = 1,
    XMLRPC_TYPE_DOUBLE   = 2,
    XMLRPC_TYPE_DATETIME = 3,
    XMLRPC_TYPE_STRING   = 4,
    XMLRPC_TYPE_BASE64   = 5,
    XMLRPC_TYPE_ARRAY    = 6,
    XMLRPC_TYPE_STRUCT   = 7,
    XMLRPC_TYPE_C_PTR    = 8,
    XMLRPC_TYPE_NIL      = 9,
    XMLRPC_TYPE_DEAD     = 0xDEAD,
} xmlrpc_type;

struct _xmlrpc_value {
    xmlrpc_type _type;
    int _refcount;

  /* stuff deleted */
};

typedef struct _xmlrpc_value xmlrpc_value;


typedef xmlrpc_value *
(*xmlrpc_method)(xmlrpc_env *   env,
                 xmlrpc_value * param_array,
                 void *         user_data);

struct _xmlrpc_registry {
    int _introspection_enabled;
    xmlrpc_value *_methods;
    xmlrpc_value *_default_method;
    xmlrpc_value *_preinvoke_method;
};

typedef struct _xmlrpc_registry xmlrpc_registry;

typedef void ((*runfirstFn)(void *));

typedef struct {
    const char *      config_file_name;
        /* NULL to use preferred proper API-level interface */

    xmlrpc_registry * registryP;

    /* runfirstFn and runfirst_arg are meaningless when 
       config_file_name is NULL
    */
    runfirstFn        runfirst;
    void *            runfirst_arg;

    unsigned int      port_number;
    const char *      log_file_name;
    unsigned int      keepalive_timeout;
    unsigned int      keepalive_max_conn;
    unsigned int      timeout;
    xmlrpc_bool       dont_advertise;

} xmlrpc_server_abyss_parms;

#define XMLRPC_APSIZE(MBRNAME) \
    XMLRPC_STRUCTSIZE(xmlrpc_server_abyss_parms, MBRNAME)

#define _XMLRPC_STRUCT_MEMBER_OFFSET(TYPE, MBRNAME) \
  ((unsigned int)(char*)&((TYPE *)0)->MBRNAME)
#define _XMLRPC_STRUCT_MEMBER_SIZE(TYPE, MBRNAME) \
  sizeof(((TYPE *)0)->MBRNAME)
#define XMLRPC_STRUCTSIZE(TYPE, MBRNAME) \
  (_XMLRPC_STRUCT_MEMBER_OFFSET(TYPE, MBRNAME) + \
  _XMLRPC_STRUCT_MEMBER_SIZE(TYPE, MBRNAME))
#define XMLRPC_FAIL_IF_FAULT(env) \
    do { if ((env)->fault_occurred) goto cleanup; } while (0)

/* &&&&&&&&&&&&&&&&&& function prototypes &&&&&&&&&&&&&&&&&&&&&&&&&&& */
void xmlrpc_env_init (xmlrpc_env* env);
void xmlrpc_env_clean (xmlrpc_env* env);

xmlrpc_value * 
xmlrpc_build_value(xmlrpc_env * const env,
                   const char * const format, 
                   ...);
void 
xmlrpc_decompose_value(xmlrpc_env *   const envP,
                       xmlrpc_value * const value,
                       const char *   const format, 
                       ...);

extern void xmlrpc_INCREF (xmlrpc_value* value);
extern void xmlrpc_DECREF (xmlrpc_value* value);

extern xmlrpc_type xmlrpc_value_type (xmlrpc_value* value);

xmlrpc_value *
xmlrpc_int_new(xmlrpc_env * const envP,
               int          const intValue);

void 
xmlrpc_read_int(xmlrpc_env *         const envP,
                const xmlrpc_value * const valueP,
                int *                const intValueP);

xmlrpc_value *
xmlrpc_bool_new(xmlrpc_env * const envP,
                xmlrpc_bool  const boolValue);

void
xmlrpc_read_bool(xmlrpc_env *         const envP,
                 const xmlrpc_value * const valueP,
                 xmlrpc_bool *        const boolValueP);

xmlrpc_value *
xmlrpc_double_new(xmlrpc_env * const envP,
                  double       const doubleValue);

void
xmlrpc_read_double(xmlrpc_env *         const envP,
                   const xmlrpc_value * const valueP,
                   xmlrpc_double *      const doubleValueP);

xmlrpc_value *
xmlrpc_datetime_new_str(xmlrpc_env * const envP,
                        const char * const value);

void
xmlrpc_read_datetime_str(xmlrpc_env *         const envP,
                         const xmlrpc_value * const valueP,
                         const char **        const stringValueP);

xmlrpc_value *
xmlrpc_string_new(xmlrpc_env * const envP,
                  const char * const stringValue);

xmlrpc_value *
xmlrpc_string_new_lp(xmlrpc_env * const envP, 
                     size_t       const length,
                     const char * const stringValue);

void
xmlrpc_read_string(xmlrpc_env *         const envP,
                   const xmlrpc_value * const valueP,
                   const char **        const stringValueP);


void
xmlrpc_read_string_lp(xmlrpc_env *         const envP,
                      const xmlrpc_value * const valueP,
                      size_t *             const lengthP,
                      const char **        const stringValueP);

xmlrpc_value *
xmlrpc_base64_new(xmlrpc_env *          const envP,
                  unsigned int          const length,
                  const unsigned char * const bytestringValue);

void
xmlrpc_read_base64(xmlrpc_env *           const envP,
                   const xmlrpc_value *   const valueP,
                   unsigned int *         const lengthP,
                   const unsigned char ** const bytestringValueP);

xmlrpc_value *
xmlrpc_array_new(xmlrpc_env * const envP);

extern void
xmlrpc_array_append_item (xmlrpc_env*   env,
                          xmlrpc_value* array,
                          xmlrpc_value* value);

xmlrpc_value * 
xmlrpc_array_get_item(xmlrpc_env *         const envP,
                      const xmlrpc_value * const arrayP,
                      int                  const index);

void
xmlrpc_array_read_item(xmlrpc_env *         const envP,
                       const xmlrpc_value * const arrayP,
                       unsigned int         const index,
                       xmlrpc_value **      const valuePP);

xmlrpc_value *
xmlrpc_struct_new(xmlrpc_env * env);

int
xmlrpc_struct_size (xmlrpc_env   * env, 
                    xmlrpc_value * strct);

void 
xmlrpc_struct_set_value(xmlrpc_env *   const env,
                        xmlrpc_value * const strct,
                        const char *   const key,
                        xmlrpc_value * const value);

int 
xmlrpc_struct_has_key(xmlrpc_env *   const envP,
                      xmlrpc_value * const strctP,
                      const char *   const key);

void
xmlrpc_struct_find_value(xmlrpc_env *    const envP,
                         xmlrpc_value *  const structP,
                         const char *    const key,
                         xmlrpc_value ** const valuePP);

void
xmlrpc_struct_read_value(xmlrpc_env *    const envP,
                         xmlrpc_value *  const strctP,
                         const char *    const key,
                         xmlrpc_value ** const valuePP);

void 
xmlrpc_struct_read_member(xmlrpc_env *    const envP,
                          xmlrpc_value *  const structP,
                          unsigned int    const index,
                          xmlrpc_value ** const keyvalP,
                          xmlrpc_value ** const valueP);

xmlrpc_value * 
xmlrpc_build_value(xmlrpc_env * const env,
                   const char * const format, 
                   ...);

void 
xmlrpc_decompose_value(xmlrpc_env *   const envP,
                       xmlrpc_value * const value,
                       const char *   const format, 
                       ...);


#define XMLRPC_CLIENT_NO_FLAGS         (0)
#define XMLRPC_CLIENT_SKIP_LIBWWW_INIT (1)

extern void
xmlrpc_client_init(int          const flags,
                   const char * const appname,
                   const char * const appversion);

extern void
xmlrpc_client_cleanup(void);

xmlrpc_value * 
xmlrpc_client_call(xmlrpc_env * const envP,
                   const char * const server_url,
                   const char * const method_name,
                   const char * const format,
                   ...);

xmlrpc_value * 
xmlrpc_client_call_params(xmlrpc_env *   const envP,
                          const char *   const serverUrl,
                          const char *   const methodName,
                          xmlrpc_value * const paramArrayP);

/*void 
  xmlrpc_client_call_asynch(const char * const server_url,
  const char * const method_name,
  xmlrpc_response_handler callback,
  void *       const user_data,
  const char * const format,
  ...);*/

extern void
xmlrpc_client_event_loop_finish_asynch(void);

extern void
xmlrpc_client_event_loop_finish_asynch_timeout(unsigned long milliseconds);

xmlrpc_registry *
xmlrpc_registry_new(xmlrpc_env * env);

void
xmlrpc_registry_free(xmlrpc_registry * registry);

void
xmlrpc_registry_add_method(xmlrpc_env *      env,
                           xmlrpc_registry * registry,
                           const char *      host,
                           const char *      method_name,
                           xmlrpc_method     method,
                           void *            user_data);

void
xmlrpc_server_abyss(xmlrpc_env *                      const envP,
                    const xmlrpc_server_abyss_parms * const parms,
                    unsigned int                      const parm_size);

#endif  /*  _XMLRPC_H */

