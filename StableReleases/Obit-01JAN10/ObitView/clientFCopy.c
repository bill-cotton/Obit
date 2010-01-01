/* A simple synchronous XML-RPC client written in C. */
/* Test bed for senting files across xmlrpc */
/* Argument 1 = File to copy
   Argument 2 = Server URL
*/

#include <stdlib.h>
#include <stdio.h>

#include <xmlrpc-c/base.h>
#include <xmlrpc-c/client.h>
#include "ObitXML.h"
#include "ObitRPC.h"
#include "ObitRPCUtil.h"
#include "ObitFile.h"
#include "ObitInfoList.h"
#define NAME "XML-RPC C Test Client synch_client"

static void die_if_fault_occurred (xmlrpc_env *env)
{
    if (env->fault_occurred) {
        fprintf(stderr, "XML-RPC Fault: %s (%d)\n",
                env->fault_string, env->fault_code);
        exit(1);
    }
}

/* prototypes */

int 
main(int           const argc, 
     const char ** const argv) {
  
  ObitRPC *client=NULL;
  ObitXML *xml=NULL, *result=NULL;
  ObitErr *err=NULL;
  ObitFile *file=NULL;
  gchar *fileName, *tfileName = "J0854+2006.PCube.fit.gz";
  /*gchar *fileName = "LICENSE";*/
  gchar *serverURL, *tserverURL = "http://localhost:8765/RPC2";
  const char * reply;
  
  /* parse optional arguments, serverURL, filename */
  if (argc>=2) {
    fileName = g_strdup(argv[1]);
  } else fileName = tfileName;

  if (argc>=3) {
    serverURL = g_strdup(argv[2]);
  } else serverURL = tserverURL;

  fprintf (stdout, "Sending %s to %s\n",fileName, serverURL);

  err = newObitErr();
  
  /* start client */
  client = ObitRPCCreateClient ("clientTest", err);

  /*&&&&&&&&&& test ping call &&&&&&&&&&&&&&&&&&&&&&&&*/
  xml = ObitXMLPing2XML(err);
  result = ObitRPCCall (client, serverURL, xml, NULL, NULL, err);
  /* resultP = xmlrpc_client_call(&env, serverURL,
     "ping",
     "(i)", (xmlrpc_int32) 41); */
  ObitErrLog  (err);
  die_if_fault_occurred(&client->envP);
  
  /* print reply. */
  xmlrpc_read_string(&client->envP, result->parmP, &reply);
  die_if_fault_occurred(&client->envP);
  printf("%s\n", reply);
  free((char*)reply);
  /* Dispose of our result value. */
  ObitXMLUnref(xml);
  ObitXMLUnref(result);
  
  /*&&&&&&&&&&&&&&&&&&&&& test copyFile call &&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /* Create File object */
  file = newObitFile("testFile");
  file->fileName = g_strdup(fileName);
  
  /* Send */
  ObitRPCUtilFileSend (file, client, serverURL, err);
  if (err->error) {ObitErrLog (err); return 1;};
  
  ObitFileUnref(file);
  client = ObitRPCUnref(client);
  
  return 0;
}  /* end main */

