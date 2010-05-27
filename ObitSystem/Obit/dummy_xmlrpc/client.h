#include <dummy_xmlrpc/base.h>

struct xmlrpc_client;
typedef struct xmlrpc_client xmlrpc_client;

typedef void (*xmlrpc_response_handler) (const char *server_url,
                                         const char *method_name,
                                         xmlrpc_value *param_array,
                                         void *user_data,
                                         xmlrpc_env *fault,
                                         xmlrpc_value *result);

struct xmlrpc_clientparms {
    const char *               transport;
    struct xmlrpc_xportparms * transportparmsP;
        /* Cast a "struct ..._xportparms *" to fit here */
    size_t                     transportparm_size;
};
#include <dummy_xmlrpc/base.h>

struct xmlrpc_client;
typedef struct xmlrpc_client xmlrpc_client;

typedef void (*xmlrpc_response_handler) (const char *server_url,
                                         const char *method_name,
                                         xmlrpc_value *param_array,
                                         void *user_data,
                                         xmlrpc_env *fault,
                                         xmlrpc_value *result);

struct xmlrpc_clientparms {
    const char *               transport;
    struct xmlrpc_xportparms * transportparmsP;
        /* Cast a "struct ..._xportparms *" to fit here */
    size_t                     transportparm_size;
};
