/* GObject - GLib Type, Object, Parameter and Signal Library
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <string.h>
#include "doofus.h"

#undef	G_LOG_DOMAIN
#define	G_LOG_DOMAIN "DoofusObject"

/* --- DoofusIface --- */
#define DOOFUS_TYPE_IFACE      (doofus_iface_get_type ())
#define DOOFUS_IFACE(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj),DOOFUS_TYPE_IFACE, DoofusIface))
#define DOOFUS_IS_IFACE(obj)	  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), DOOFUS_TYPE_IFACE))
#define DOOFUS_IFACE_GET_CLASS(obj) (G_TYPE_INSTANCE_GET_INTERFACE ((obj), DOOFUS_TYPE_IFACE, DoofusIfaceClass))
typedef struct _DoofusIface      DoofusIface;
typedef struct _DoofusIfaceClass DoofusIfaceClass;
struct _DoofusIfaceClass
{
  GTypeInterface base_iface;
  void	(*print_string)	(DoofusIface	*tiobj,
			 const gchar	*string);
};
static void	iface_base_init		(DoofusIfaceClass	*iface);
static void	iface_base_finalize	(DoofusIfaceClass	*iface);
static void	print_foo		(DoofusIface	*tiobj,
					 const gchar	*string);
GType
doofus_iface_get_type (void)
{
  static GType doofus_iface_type = 0;

  if (!doofus_iface_type)
    {
      static const GTypeInfo doofus_iface_info =
      {
	sizeof (DoofusIfaceClass),
	(GBaseInitFunc)	iface_base_init,		/* base_init */
	(GBaseFinalizeFunc) iface_base_finalize,	/* base_finalize */
      };

      doofus_iface_type = g_type_register_static (G_TYPE_INTERFACE, "DoofusIface", &doofus_iface_info, 0);
      g_type_interface_add_prerequisite (doofus_iface_type, G_TYPE_OBJECT);
    }

  return doofus_iface_type;
}
static guint iface_base_init_count = 0;
static void
iface_base_init (DoofusIfaceClass *iface)
{
  iface_base_init_count++;
  if (iface_base_init_count == 1)
    {
      /* add signals here */
    }
}
static void
iface_base_finalize (DoofusIfaceClass *iface)
{
  iface_base_init_count--;
  if (iface_base_init_count == 0)
    {
      /* destroy signals here */
    }
}
static void
print_foo (DoofusIface   *tiobj,
	   const gchar *string)
{
  if (!string)
    string = "<NULL>";
  g_print ("Iface-FOO: \"%s\" from %p\n", string, tiobj);
}
static void
doofus_object_doofus_iface_init (gpointer giface,
			     gpointer iface_data)
{
  DoofusIfaceClass *iface = giface;

  g_assert (iface_data == GUINT_TO_POINTER (42));

  g_assert (G_TYPE_FROM_INTERFACE (iface) == DOOFUS_TYPE_IFACE);

  /* assert iface_base_init() was already called */
  g_assert (iface_base_init_count > 0);

  /* initialize stuff */
  iface->print_string = print_foo;
}
void
iface_print_string (DoofusIface   *tiobj,
		    const gchar *string)
{
  DoofusIfaceClass *iface;

  g_return_if_fail (DOOFUS_IS_IFACE (tiobj));
  g_return_if_fail (G_IS_OBJECT (tiobj)); /* ensured through prerequisite */

  iface = DOOFUS_IFACE_GET_CLASS (tiobj);
  g_object_ref (tiobj);
  iface->print_string (tiobj, string);
  g_object_unref (tiobj);
}


/* --- DoofusObject --- */
typedef struct _DoofusObjectClass DoofusObjectClass;
struct _DoofusObjectClass
{
  GObjectClass parent_class;

  gchar* (*doofus_signal) (DoofusObject *tobject,
			 DoofusIface  *iface_object,
			 gpointer    tdata);
};
static void	doofus_object_class_init	(DoofusObjectClass	*class);
static void	doofus_object_init	(DoofusObject		*tobject);
static gboolean	doofus_signal_accumulator	(GSignalInvocationHint	*ihint,
					 GValue            	*return_accu,
					 const GValue       	*handler_return,
					 gpointer                data);
static gchar*	doofus_object_test_signal	(DoofusObject		*tobject,
					 DoofusIface		*iface_object,
					 gpointer		 tdata);
GType
doofus_object_get_type (void)
{
  static GType doofus_object_type = 0;

  if (!doofus_object_type)
    {
      static const GTypeInfo doofus_object_info =
      {
	sizeof (DoofusObjectClass),
	NULL,           /* base_init */
	NULL,           /* base_finalize */
	(GClassInitFunc) doofus_object_class_init,
	NULL,           /* class_finalize */
	NULL,           /* class_data */
	sizeof (DoofusObject),
	5,              /* n_preallocs */
	(GInstanceInitFunc) doofus_object_init,
      };
      GInterfaceInfo iface_info = { doofus_object_test_iface_init, NULL, GUINT_TO_POINTER (42) };

      doofus_object_type = g_type_register_static (G_TYPE_OBJECT, "DoofusObject", &doofus_object_info, 0);
      g_type_add_interface_static (doofus_object_type, DOOFUS_TYPE_IFACE, &iface_info);
    }

  return doofus_object_type;
}
static void
doofus_object_class_init (DoofusObjectClass *class)
{
  /*  GObjectClass *gobject_class = G_OBJECT_CLASS (class); */

  class->doofus_signal = doofus_object_test_signal;

  g_signal_new ("test-signal",
		G_OBJECT_CLASS_TYPE (class),
		G_SIGNAL_RUN_FIRST | G_SIGNAL_RUN_LAST | G_SIGNAL_RUN_CLEANUP,
		G_STRUCT_OFFSET (DoofusObjectClass, test_signal),
		test_signal_accumulator, NULL,
		g_cclosure_marshal_STRING__OBJECT_POINTER,
		G_TYPE_STRING, 2, DOOFUS_TYPE_IFACE, G_TYPE_POINTER);
}
static void
doofus_object_init (DoofusObject *tobject)
{
}
static gboolean
doofus_signal_accumulator (GSignalInvocationHint *ihint,
			 GValue                *return_accu,
			 const GValue          *handler_return,
			 gpointer               data)
{
  gchar *accu_string = g_value_get_string (return_accu);
  gchar *new_string = g_value_get_string (handler_return);
  gchar *result_string;

  if (accu_string)
    result_string = g_strconcat (accu_string, new_string, NULL);
  else if (new_string)
    result_string = g_strdup (new_string);
  else
    result_string = NULL;

  g_value_set_string_take_ownership (return_accu, result_string);

  return TRUE;
}
static gchar*
doofus_object_test_signal (DoofusObject *tobject,
			 DoofusIface  *iface_object,
			 gpointer    tdata)
{
  g_message ("::test_signal default_handler called");

  g_return_val_if_fail (DOOFUS_IS_IFACE (iface_object), NULL);
  
  return g_strdup ("<default_handler>");
}


/* --- DoofusIface for DerivedObject --- */
static void
print_bar (DoofusIface   *tiobj,
	   const gchar *string)
{
  DoofusIfaceClass *parent_iface;

  g_return_if_fail (DOOFUS_IS_IFACE (tiobj));

  if (!string)
    string = "<NULL>";
  g_print ("Iface-BAR: \"%s\" from %p\n", string, tiobj);

  g_print ("chaining: ");
  parent_iface = g_type_interface_peek_parent (DOOFUS_IFACE_GET_CLASS (tiobj));
  parent_iface->print_string (tiobj, string);

  g_assert (g_type_interface_peek_parent (parent_iface) == NULL);
}

static void
derived_object_doofus_iface_init (gpointer giface,
				gpointer iface_data)
{
  DoofusIfaceClass *iface = giface;

  g_assert (iface_data == GUINT_TO_POINTER (87));

  g_assert (G_TYPE_FROM_INTERFACE (iface) == DOOFUS_TYPE_IFACE);

  /* assert doofus_object_test_iface_init() was already called */
  g_assert (iface->print_string == print_foo);

  /* override stuff */
  iface->print_string = print_bar;
}


/* --- DerivedObject --- */
#define DERIVED_TYPE_OBJECT            (derived_object_get_type ())
#define DERIVED_OBJECT(object)         (G_TYPE_CHECK_INSTANCE_CAST ((object), DERIVED_TYPE_OBJECT, DerivedObject))
#define DERIVED_OBJECT_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST ((klass), DERIVED_TYPE_OBJECT, DerivedObjectClass))
#define DERIVED_IS_OBJECT(object)      (G_TYPE_CHECK_INSTANCE_TYPE ((object), DERIVED_TYPE_OBJECT))
#define DERIVED_IS_OBJECT_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), DERIVED_TYPE_OBJECT))
#define DERIVED_OBJECT_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj), DERIVED_TYPE_OBJECT, DerivedObjectClass))
typedef struct _DoofusObject      DerivedObject;
typedef struct _DoofusObjectClass DerivedObjectClass;
GType
derived_object_get_type (void)
{
  static GType derived_object_type = 0;

  if (!derived_object_type)
    {
      static const GTypeInfo derived_object_info =
      {
	sizeof (DerivedObjectClass),
	NULL,           /* base_init */
	NULL,           /* base_finalize */
	NULL,		/* class_init */
	NULL,           /* class_finalize */
	NULL,           /* class_data */
	sizeof (DerivedObject),
	5,              /* n_preallocs */
	NULL,		/* instance_init */
      };
      GInterfaceInfo iface_info = { derived_object_doofus_iface_init, NULL, GUINT_TO_POINTER (87) };

      derived_object_type = g_type_register_static (DOOFUS_TYPE_OBJECT, "DerivedObject", &derived_object_info, 0);
      g_type_add_interface_static (derived_object_type, DOOFUS_TYPE_IFACE, &iface_info);
    }

  return derived_object_type;
}


/* --- main --- */
int
xxxmain (int   argc,
      char *argv[])
{
  GTypeInfo info = { 0, };
  GTypeFundamentalInfo finfo = { 0, };
  GType type;
  DoofusObject *sigarg;
  DerivedObject *dobject;
  gchar *string = NULL;

  g_log_set_always_fatal (g_log_set_always_fatal (G_LOG_FATAL_MASK) |
			  G_LOG_LEVEL_WARNING |
			  G_LOG_LEVEL_CRITICAL);
  g_type_init_with_debug_flags (G_TYPE_DEBUG_OBJECTS | G_TYPE_DEBUG_SIGNALS);

  /* test new fundamentals */
  g_assert (G_TYPE_MAKE_FUNDAMENTAL (G_TYPE_RESERVED_USER_FIRST) == g_type_fundamental_next ());
  type = g_type_register_fundamental (g_type_fundamental_next (), "FooShadow1", &info, &finfo, 0);
  g_assert (G_TYPE_MAKE_FUNDAMENTAL (G_TYPE_RESERVED_USER_FIRST + 1) == g_type_fundamental_next ());
  type = g_type_register_fundamental (g_type_fundamental_next (), "FooShadow2", &info, &finfo, 0);
  g_assert (G_TYPE_MAKE_FUNDAMENTAL (G_TYPE_RESERVED_USER_FIRST + 2) == g_type_fundamental_next ());
  g_assert (g_type_from_name ("FooShadow1") == G_TYPE_MAKE_FUNDAMENTAL (G_TYPE_RESERVED_USER_FIRST));
  g_assert (g_type_from_name ("FooShadow2") == G_TYPE_MAKE_FUNDAMENTAL (G_TYPE_RESERVED_USER_FIRST + 1));

  /* to test past class initialization interface setups, create the class here */
  g_type_class_ref (DOOFUS_TYPE_OBJECT);

  dobject = g_object_new (DERIVED_TYPE_OBJECT, NULL);
  sigarg = g_object_new (DOOFUS_TYPE_OBJECT, NULL);

  g_print ("MAIN: emit test-signal:\n");
  g_signal_emit_by_name (dobject, "test-signal", sigarg, NULL, &string);
  g_message ("signal return: \"%s\"", string);
  g_assert (strcmp (string, "<default_handler><default_handler>") == 0);
  g_free (string);

  g_print ("MAIN: call iface print-string on test and derived object:\n");
  iface_print_string (DOOFUS_IFACE (sigarg), "iface-string-from-test-type");
  iface_print_string (DOOFUS_IFACE (dobject), "iface-string-from-derived-type");
  
  g_object_unref (sigarg);
  g_object_unref (dobject);

  g_message ("%s done", argv[0]);

  return 0;
}
