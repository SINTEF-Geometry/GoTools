
#ifndef _GVACTIONSHEET_H
#define _GVACTIONSHEET_H

#include <QObject>
#include "GoTools/utils/config.h"

class QWIDGET;
class gvObserver;

/** gvActionSheet:
 */

class gvActionSheet : public QObject
{
  Q_OBJECT;


public:
  virtual ~gvActionSheet()
    {;}
  virtual void showDialog(gvObserver* observer) = 0;
  
 protected:
  shared_ptr<QWidget> form_;
};

#endif // _GVACTIONSHEET_H

