/****************************************************************************
** Meta object code from reading C++ file 'gvApplication.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../../viewlib/include/GoTools/viewlib/gvApplication.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'gvApplication.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_gvApplication[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      40,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x0a,
      22,   14,   14,   14, 0x0a,
      48,   14,   14,   14, 0x0a,
      68,   14,   14,   14, 0x0a,
      85,   14,   14,   14, 0x0a,
      92,   14,   14,   14, 0x0a,
     100,   14,   14,   14, 0x0a,
     113,   14,   14,   14, 0x0a,
     134,   14,   14,   14, 0x0a,
     151,   14,   14,   14, 0x0a,
     163,   14,   14,   14, 0x0a,
     175,   14,   14,   14, 0x0a,
     191,   14,   14,   14, 0x0a,
     211,   14,   14,   14, 0x0a,
     234,   14,   14,   14, 0x0a,
     257,  253,   14,   14, 0x0a,
     286,   14,   14,   14, 0x0a,
     314,   14,   14,   14, 0x0a,
     331,   14,   14,   14, 0x0a,
     355,   14,   14,   14, 0x0a,
     381,   14,   14,   14, 0x0a,
     398,   14,   14,   14, 0x0a,
     416,   14,   14,   14, 0x0a,
     432,   14,   14,   14, 0x0a,
     445,   14,   14,   14, 0x0a,
     459,   14,   14,   14, 0x0a,
     476,   14,   14,   14, 0x0a,
     498,   14,   14,   14, 0x0a,
     518,   14,   14,   14, 0x0a,
     539,   14,   14,   14, 0x0a,
     568,  563,   14,   14, 0x0a,
     610,  595,   14,   14, 0x0a,
     647,   14,   14,   14, 0x0a,
     664,   14,   14,   14, 0x0a,
     685,   14,   14,   14, 0x0a,
     705,   14,   14,   14, 0x0a,
     732,  724,   14,   14, 0x09,
     780,  760,   14,   14, 0x09,
     827,  814,   14,   14, 0x09,
     884,  864,   14,   14, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_gvApplication[] = {
    "gvApplication\0\0open()\0reload_last_opened_file()\0"
    "save_selection_as()\0close_document()\0"
    "quit()\0about()\0view_reset()\0"
    "view_reset_visible()\0view_wireframe()\0"
    "view_axis()\0view_cull()\0view_specular()\0"
    "view_orthographic()\0toggle_blending_mode()\0"
    "view_focus_point()\0x,y\0"
    "view_focus_point_cb(int,int)\0"
    "display_object_properties()\0"
    "assign_texture()\0set_curve_resolutions()\0"
    "set_surface_resolutions()\0enable_objects()\0"
    "disable_objects()\0toggle_enable()\0"
    "select_all()\0select_none()\0select_inverse()\0"
    "select_all_surfaces()\0select_all_curves()\0"
    "select_all_visible()\0toggle_selection_mode()\0"
    "name\0toggle_select_object(uint)\0"
    "names,numnames\0toggle_multiselect_object(uint*,int)\0"
    "group_selected()\0dismiss_selections()\0"
    "show_control_nets()\0set_random_color()\0"
    "new_res\0changeCurveResolutions(int)\0"
    "new_u_res,new_v_res\0"
    "changeSurfaceResolutions(int,int)\0"
    "members,name\0add_group(std::vector<int>&,QString)\0"
    "new_objs,new_colors\0"
    "add_objects(ObjContainer&,ColContainer&)\0"
};

void gvApplication::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        gvApplication *_t = static_cast<gvApplication *>(_o);
        switch (_id) {
        case 0: _t->open(); break;
        case 1: _t->reload_last_opened_file(); break;
        case 2: _t->save_selection_as(); break;
        case 3: _t->close_document(); break;
        case 4: _t->quit(); break;
        case 5: _t->about(); break;
        case 6: _t->view_reset(); break;
        case 7: _t->view_reset_visible(); break;
        case 8: _t->view_wireframe(); break;
        case 9: _t->view_axis(); break;
        case 10: _t->view_cull(); break;
        case 11: _t->view_specular(); break;
        case 12: _t->view_orthographic(); break;
        case 13: _t->toggle_blending_mode(); break;
        case 14: _t->view_focus_point(); break;
        case 15: _t->view_focus_point_cb((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 16: _t->display_object_properties(); break;
        case 17: _t->assign_texture(); break;
        case 18: _t->set_curve_resolutions(); break;
        case 19: _t->set_surface_resolutions(); break;
        case 20: _t->enable_objects(); break;
        case 21: _t->disable_objects(); break;
        case 22: _t->toggle_enable(); break;
        case 23: _t->select_all(); break;
        case 24: _t->select_none(); break;
        case 25: _t->select_inverse(); break;
        case 26: _t->select_all_surfaces(); break;
        case 27: _t->select_all_curves(); break;
        case 28: _t->select_all_visible(); break;
        case 29: _t->toggle_selection_mode(); break;
        case 30: _t->toggle_select_object((*reinterpret_cast< uint(*)>(_a[1]))); break;
        case 31: _t->toggle_multiselect_object((*reinterpret_cast< uint*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 32: _t->group_selected(); break;
        case 33: _t->dismiss_selections(); break;
        case 34: _t->show_control_nets(); break;
        case 35: _t->set_random_color(); break;
        case 36: _t->changeCurveResolutions((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 37: _t->changeSurfaceResolutions((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 38: _t->add_group((*reinterpret_cast< std::vector<int>(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 39: _t->add_objects((*reinterpret_cast< ObjContainer(*)>(_a[1])),(*reinterpret_cast< ColContainer(*)>(_a[2]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData gvApplication::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject gvApplication::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_gvApplication,
      qt_meta_data_gvApplication, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &gvApplication::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *gvApplication::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *gvApplication::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_gvApplication))
        return static_cast<void*>(const_cast< gvApplication*>(this));
    return QWidget::qt_metacast(_clname);
}

int gvApplication::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 40)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 40;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
