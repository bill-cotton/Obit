/* header class for demo of generating gobject/gtype class */
#include	<glib-object.h>
#ifndef DOOFUS_H
#define DOOFUS_H
struct _DoofusObject
{
  GObject parent_instance;
};
typedef struct _DoofusObject      Doofus;

/* prototypes */
GType doofus_object_get_type (void);

/* macroes */
#define DOOFUS_TYPE_OBJECT            (doofus_object_get_type ())
#define DOOFUS_OBJECT(object)         (G_TYPE_CHECK_INSTANCE_CAST ((object), DOOFUS_TYPE_OBJECT, DoofusObject))
#define DOOFUS_OBJECT_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST ((klass), DOOFUS_TYPE_OBJECT, DoofusObjectClass))
#define DOOFUS_IS_OBJECT(object)      (G_TYPE_CHECK_INSTANCE_TYPE ((object), DOOFUS_TYPE_OBJECT))
#define DOOFUS_IS_OBJECT_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), DOOFUS_TYPE_OBJECT))
#define DOOFUS_OBJECT_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj), DOOFUS_TYPE_OBJECT, DoofusObjectClass))
#endif /* DOOFUS_H */
