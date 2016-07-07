/****************************************************************************
** Meta object code from reading C++ file 'gvApplicationVolAndLR.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../../../viewlib/include/GoTools/viewlib/vol_and_lr/gvApplicationVolAndLR.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'gvApplicationVolAndLR.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_gvApplicationVolAndLR[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      23,   22,   22,   22, 0x0a,
      36,   22,   22,   22, 0x0a,
      58,   22,   22,   22, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_gvApplicationVolAndLR[] = {
    "gvApplicationVolAndLR\0\0view_reset()\0"
    "translate_to_origin()\0move_vertices_to_origin()\0"
};

void gvApplicationVolAndLR::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        gvApplicationVolAndLR *_t = static_cast<gvApplicationVolAndLR *>(_o);
        switch (_id) {
        case 0: _t->view_reset(); break;
        case 1: _t->translate_to_origin(); break;
        case 2: _t->move_vertices_to_origin(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObjectExtraData gvApplicationVolAndLR::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject gvApplicationVolAndLR::staticMetaObject = {
    { &gvApplication::staticMetaObject, qt_meta_stringdata_gvApplicationVolAndLR,
      qt_meta_data_gvApplicationVolAndLR, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &gvApplicationVolAndLR::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *gvApplicationVolAndLR::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *gvApplicationVolAndLR::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_gvApplicationVolAndLR))
        return static_cast<void*>(const_cast< gvApplicationVolAndLR*>(this));
    return gvApplication::qt_metacast(_clname);
}

int gvApplicationVolAndLR::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = gvApplication::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
