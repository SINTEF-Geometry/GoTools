
#ifndef _GVACTIONSHEET_H
#define _GVACTIONSHEET_H

#include <QObject>
#include <memory>

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
  std::shared_ptr<QWidget> form_;
};

#endif // _GVACTIONSHEET_H

