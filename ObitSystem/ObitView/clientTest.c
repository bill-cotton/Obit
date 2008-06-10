/* A simple synchronous XML-RPC client written in C. */

#include <stdlib.h>
#include <stdio.h>

/*#include <xmlrpc-c/base.h>
  #include <xmlrpc-c/client.h>*/
#include <xmlrpc.h>
#include <xmlrpc_client.h>
#include "ObitDConCleanWindow.h"
#include "ObitXML.h"
#include "ObitRPC.h"
#include "ObitInfoList.h"
/*#include "xml.h"*/

/* #include "config.h"  information about this build environment */

#define NAME "XML-RPC C Test Client synch_client"
/*#define VERSION "1.0"*/

static void die_if_fault_occurred (xmlrpc_env *env)
{
    if (env->fault_occurred) {
        fprintf(stderr, "XML-RPC Fault: %s (%d)\n",
                env->fault_string, env->fault_code);
        exit(1);
    }
}

/* prototypes */
ObitDConCleanWindow* fakeWindow(ObitErr *err);
static void TellWindow (ObitDConCleanWindow *myWindow, ObitErr *err);

int 
main(int           const argc, 
     const char ** const argv) {

  ObitRPC *client=NULL;
  ObitXML *xml=NULL, *result=NULL;
  ObitErr *err=NULL;
  ObitInfoList *list1=NULL, *list2=NULL;
  gchar *serverURL = "http://localhost:8765/RPC2";
  olong req, field;
  /* debugging XML->infolist */
  /*gint itest; */
  /*gboolean barr[5]={FALSE, TRUE, FALSE, TRUE, FALSE}; */
  /*gfloat ftest, farr[5]={1.234, 5.678, 9.012, 3.456, 7.890}; */
  /*gdouble darr[5]={61.234, 65.678, 69.012, 63.456, 67.890}; */
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType infoType;
  /* xmlrpc_env env;*/
  /*xmlrpc_value * resultP;*/
    const char * reply;
    ObitDConCleanWindow *window;
     /*xmlrpc_value *windowP, *FileInfoP;*/
    /*char dbugBuff[20000];*/

    if (argc-1 > 0) {
        fprintf(stderr, "No arguments");
        exit(0);
    }
    err = newObitErr();

    /* Test InfoList<->XML */
    /*dim[0]=dim[1]=dim[2]=dim[3]=dim[4]=1; */
    /*list1 = newObitInfoList(); */
    /*ftest = 5.3456; */
    /*ObitInfoListPut(list1, "Float", OBIT_float, dim, &ftest, err); */
    /*itest = 9999; */
    /*ObitInfoListPut(list1, "Int", OBIT_long, dim, &itest, err); */
    /*dim[0] = 5; */
    /*ObitInfoListPut(list1, "FloatArr", OBIT_float, dim, farr, err); */
    /*ObitInfoListPut(list1, "DoubleArr", OBIT_double, dim, darr, err); */
    /*ObitInfoListPut(list1, "BoolArr", OBIT_bool, dim, barr, err); */
    /*dim[0] = 14; */
    /*ObitInfoListPut(list1, "strIng", OBIT_string, dim, "Something ugly", err); */
    /*if (list1) ObitInfoListPrint(list1, stdout); */
    /*xml = ObitXMLInfoList2XML (list1, err); */
    /*list2 = ObitXMLXMLInfo2List(xml, err); */
    /*if (err->error) {ObitErrLog (err); return 1;}; */


    /*if (list2) ObitInfoListPrint(list2, stdout); */
    /* Dispose of our result value. */
    /*list1 = ObitInfoListUnref(list1); */
    /*list2 = ObitInfoListUnref(list2); */
   /* xml = ObitXMLUnref(xml); */
    if (err->error) {ObitErrLog (err); return 1;};


    /* start client */
    client = ObitRPCCreateClient ("clientTest", err);
    /* Start up our XML-RPC client library. */
   /*  xmlrpc_client_init(XMLRPC_CLIENT_NO_FLAGS, NAME, VERSION); */

    /* Initialize our error-handling environment. */
    /* xmlrpc_env_init(&env); */

   /*&&&&&&&&&& test ping call &&&&&&&&&&&&&&&&&&&&&&&&*/
    xml = ObitXMLPing2XML(err);
    result = ObitRPCCall (client, serverURL, xml, NULL, NULL, err);
    /* resultP = xmlrpc_client_call(&env, "http://localhost:8765/RPC2",
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

   /*&&&&&&&&&&&&&&&&&&&&& test loadFITS call &&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
   /* resultP = xmlrpc_client_call(&env, "http://localhost:8765/RPC2",
				 "loadFITS",
				 "(s)", "someFile.fits"); */
    /*die_if_fault_occurred(&env); */
    
    /* print reply. */
    /*xmlrpc_read_string(&env, resultP, &reply); */
    /*die_if_fault_occurred(&env); */
    /*printf("%s\n", reply); */
    /*free((char*)reply); */
   
    /* Dispose of our result value. */
    /*xmlrpc_DECREF(resultP); */

    /*&&&&&&&&&& test loadImage call with FITS data &&&&&&&&&&&&&&&&&&&&&&&&*/
    /* FileInfoP = FileInfo2XML (&env, OBIT_IO_FITS, "someFile.fits", 
			      NULL, NULL, 0, 0, err);*/
    /*resultP = xmlrpc_client_call(&env, "http://localhost:8765/RPC2",
      "loadImage",
      "(S)", FileInfoP);*/
    xml = ObitXMLFileInfo2XML (OBIT_IO_FITS, "someFile.fits", 
			       NULL, NULL, 0, 0, 1, 1, err);
    result = ObitRPCCall (client, serverURL, xml, &list1, &list2, err);
    if (err->error) {ObitErrLog (err); return 1;};
    die_if_fault_occurred(&client->envP);

    
    /* print reply. */
    xmlrpc_read_string(&client->envP, result->parmP, &reply);
    die_if_fault_occurred(&client->envP);
    printf("%s\n", reply);
    free((char*)reply);
    if (list1) ObitInfoListPrint(list1, stdout);
    if (list2) ObitInfoListPrint(list2, stdout);
    /* Dispose of our result value. */
    list1 = ObitInfoListUnref(list1);
    list2 = ObitInfoListUnref(list2);
    ObitXMLUnref(xml);
    ObitXMLUnref(result);

    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&& test window editing &&&&&&&&&&&&&&&&&&&&*/
    window = fakeWindow (err);
    xml = ObitXMLWindow2XML (window, 1, err);
    result = ObitRPCCall (client, serverURL, xml, &list1, &list2, err);
    /*windowP = Window2XML (&env, window, 1, err);*/
    ObitErrLog  (err);
   
    /* Have to send argument packaged into an array */
    /*resultP = xmlrpc_client_call(&env, "http://localhost:8765/RPC2",
				 "editWindow",
				 "(S)", windowP);*/

    /* Check Request */
    if (list2) {
      if (ObitInfoListGetTest (list2, "Request", &infoType, dim, &req)) {
	fprintf (stdout, "Request from edit %d\n",req);
	if (ObitInfoListGetTest (list2, "Field", &infoType, dim, &field)) {
	  fprintf (stdout, "Request field %d\n",field);
	}
      }
    }
       

    if (list1) ObitInfoListPrint(list1, stdout);
    if (list2) ObitInfoListPrint(list2, stdout);
    /* Dispose of our result value. */
    list1 = ObitInfoListUnref(list1);
    list2 = ObitInfoListUnref(list2);

    ObitErrLog  (err);
    die_if_fault_occurred(&client->envP);
    window = ObitDConCleanWindowUnref(window);

    /* If this worked then the window will come back edited */
    if (xmlrpc_value_type(result->parmP)==XMLRPC_TYPE_STRUCT) {
      window = ObitXMLXML2Window (result, err);
      /* print reply. */
      TellWindow (window, err);
      
    } else if (xmlrpc_value_type(result->parmP)==XMLRPC_TYPE_STRING) {  
      /* Failed - Print message */
      xmlrpc_read_string(&result->envP, result->parmP, &reply);
      die_if_fault_occurred(&result->envP);
      printf("Failed: %s\n", reply);
      free((char*)reply);
    }
    ObitErrLog  (err);
    die_if_fault_occurred(&result->envP);

    /* Dispose of our result value. */
    ObitXMLUnref(xml);
    ObitXMLUnref(result);
    window = ObitDConCleanWindowUnref(window);

  /*&&&&&&&&&& test loadImage call with AIPS data &&&&&&&&&&&&&&&&&&&&&&&&*/
    /* ObitView    .Test  .   1  user 103  symlinking to AIPSDIR */
    xml = ObitXMLFileInfo2XML (OBIT_IO_AIPS, "ObitView", "Test", 
			       "./AIPSDIR/", 1, 103, 1, 1, err);
    result = ObitRPCCall (client, serverURL, xml, NULL, NULL, err);
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

    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

   /* Clean up our error-handling environment. */
   /* xmlrpc_env_clean(&env); */
    
    /* Shutdown our XML-RPC client library. */
   /* xmlrpc_client_cleanup();. */
    client = ObitRPCUnref(client);

    return 0;
}

/* Fake window for test - match someFile.fits */
ObitDConCleanWindow* fakeWindow(ObitErr *err)
{
  ObitDConCleanWindow *out;
  olong box[4], nax[2], field = 1;
  ObitDConCleanWindowType type;
  
   /* initial structure */
  nax[0] = 61;
  nax[1] = 61;
  out = ObitDConCleanWindowCreate1("CleanWindow", nax, err);

  out->nfield      = 1;
  out->ndim        = 2;
  out->maxId[0]    = 0;
  type = OBIT_DConCleanWindow_round;
  box[0] = 10;
  box[1] = 19;
  box[2] = 46;
  ObitDConCleanWindowAdd (out, field, type, box, err);

  return out;
} /* end  fakeWindow */

/**
 * Print current contents of myWindow
 */
static void TellWindow (ObitDConCleanWindow *myWindow, ObitErr *err)
{
  olong iD, nId;
  olong *win;
  ObitDConCleanWindowType type;

  /* How many potential windows? */
  nId = myWindow->maxId[0];
  /* Extract from Obit Object and draw */
  for (iD=1; iD<=nId; iD++) {
    if (ObitDConCleanWindowInfo(myWindow, 1, iD, &type, &win, err)) {
      /* Which box type */
      if (type==OBIT_DConCleanWindow_rectangle) 
	fprintf (stdout, "Box %d rectangle, blc = [ %d, %d], trc = [ %d, %d]\n", 
		 iD, win[0]+1,win[1]+1, win[2]+1, win[3]+1);
      if (type==OBIT_DConCleanWindow_round)
	fprintf (stdout, "Box %d round, radius =  %d, center = [ %d, %d]\n", 
		 iD, win[0], win[1]+1, win[2]+1);
   }
  } /* end loop over windows */
  /* Any errors? */
  if (err->error) ObitErrLog(err);
} /* end TellWindow */

