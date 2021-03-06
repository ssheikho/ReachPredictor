#ifndef SW_WRAPPER_H
#define SW_WRAPPER_H

#include "ExternalUpdate.h"

#include <gtk/gtk.h>

/*	On update()
	0 - Do nothing
	1 - requestFullChildWidth();
*/

class SWWrapper : public ExternalUpdate {
	public:
		SWWrapper();
		~SWWrapper();

		void update();

		void setChild(GtkWidget *child);
		GtkWidget *getWidget();

		void requestFullChildWidth();

	protected:
		GtkWidget *_scrolledWindow, *_child;
		int _onUpdate;
};

#endif
