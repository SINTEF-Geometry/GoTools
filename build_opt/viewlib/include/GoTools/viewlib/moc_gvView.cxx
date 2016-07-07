/****************************************************************************
** Meta object code from reading C++ file 'gvView.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../../viewlib/include/GoTools/viewlib/gvView.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'gvView.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_gvView[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      13,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: signature, parameters, type, tag, flags
      13,    8,    7,    7, 0x05,
      47,   32,    7,    7, 0x05,
      82,   72,    7,    7, 0x05,

 // slots: signature, parameters, type, tag, flags
     105,  100,    7,    7, 0x0a,
     124,  100,    7,    7, 0x0a,
     147,  100,    7,    7, 0x0a,
     161,  100,    7,    7, 0x0a,
     179,  100,    7,    7, 0x0a,
     197,  100,    7,    7, 0x0a,
     218,  100,    7,    7, 0x0a,
     240,  100,    7,    7, 0x0a,
     262,  100,    7,    7, 0x0a,
     288,  284,    7,    7, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_gvView[] = {
    "gvView\0\0name\0objectPicked(uint)\0"
    "names,numnames\0objectsPicked(uint*,int)\0"
    "xrel,yrel\0feedback(int,int)\0mode\0"
    "setWireframe(bool)\0setSelectionmode(bool)\0"
    "setAxis(bool)\0setBackCull(bool)\0"
    "setSpecular(bool)\0setPerspective(bool)\0"
    "setFeedbackmode(bool)\0setBlendingmode(bool)\0"
    "setGetClickmode(bool)\0x,y\0setCenter(int,int)\0"
};

void gvView::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        gvView *_t = static_cast<gvView *>(_o);
        switch (_id) {
        case 0: _t->objectPicked((*reinterpret_cast< uint(*)>(_a[1]))); break;
        case 1: _t->objectsPicked((*reinterpret_cast< uint*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 2: _t->feedback((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 3: _t->setWireframe((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 4: _t->setSelectionmode((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: _t->setAxis((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 6: _t->setBackCull((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 7: _t->setSpecular((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 8: _t->setPerspective((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 9: _t->setFeedbackmode((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 10: _t->setBlendingmode((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: _t->setGetClickmode((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: _t->setCenter((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData gvView::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject gvView::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_gvView,
      qt_meta_data_gvView, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &gvView::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *gvView::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *gvView::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_gvView))
        return static_cast<void*>(const_cast< gvView*>(this));
    if (!strcmp(_clname, "gvObserver"))
        return static_cast< gvObserver*>(const_cast< gvView*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int gvView::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 13)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 13;
    }
    return _id;
}

// SIGNAL 0
void gvView::objectPicked(unsigned int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void gvView::objectsPicked(unsigned int * _t1, int _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void gvView::feedback(int _t1, int _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}
QT_END_MOC_NAMESPACE
