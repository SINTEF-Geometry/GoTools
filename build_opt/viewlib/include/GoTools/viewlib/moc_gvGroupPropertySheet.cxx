/****************************************************************************
** Meta object code from reading C++ file 'gvGroupPropertySheet.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../../viewlib/include/GoTools/viewlib/gvGroupPropertySheet.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'gvGroupPropertySheet.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_gvGroupPropertySheet[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      35,   22,   21,   21, 0x05,

 // slots: signature, parameters, type, tag, flags
      76,   21,   21,   21, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_gvGroupPropertySheet[] = {
    "gvGroupPropertySheet\0\0members,name\0"
    "value_changed(std::vector<int>&,QString)\0"
    "accept()\0"
};

void gvGroupPropertySheet::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        gvGroupPropertySheet *_t = static_cast<gvGroupPropertySheet *>(_o);
        switch (_id) {
        case 0: _t->value_changed((*reinterpret_cast< std::vector<int>(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 1: _t->accept(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData gvGroupPropertySheet::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject gvGroupPropertySheet::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_gvGroupPropertySheet,
      qt_meta_data_gvGroupPropertySheet, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &gvGroupPropertySheet::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *gvGroupPropertySheet::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *gvGroupPropertySheet::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_gvGroupPropertySheet))
        return static_cast<void*>(const_cast< gvGroupPropertySheet*>(this));
    if (!strcmp(_clname, "gvPropertySheet"))
        return static_cast< gvPropertySheet*>(const_cast< gvGroupPropertySheet*>(this));
    return QObject::qt_metacast(_clname);
}

int gvGroupPropertySheet::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void gvGroupPropertySheet::value_changed(std::vector<int> & _t1, QString _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
