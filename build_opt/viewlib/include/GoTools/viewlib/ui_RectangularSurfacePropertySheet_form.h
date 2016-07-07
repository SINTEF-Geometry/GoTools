/********************************************************************************
** Form generated from reading UI file 'RectangularSurfacePropertySheet_form.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_RECTANGULARSURFACEPROPERTYSHEET_FORM_H
#define UI_RECTANGULARSURFACEPROPERTYSHEET_FORM_H

#include <Qt3Support/Q3GroupBox>
#include <Qt3Support/Q3MimeSourceFactory>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLCDNumber>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_RectangularSurfacePropertySheet_form
{
public:
    Q3GroupBox *box;
    QLabel *UresLabel;
    QSlider *UresSlider;
    QLCDNumber *UresLCDNumber;
    QSlider *VresSlider;
    QLCDNumber *VresLCDNumber;
    QLabel *VresLabel;
    QWidget *widget;
    QHBoxLayout *hboxLayout;
    QSpacerItem *spacerItem;
    QPushButton *ApplyButton;
    QPushButton *CloseButton;
    QSpacerItem *spacerItem1;
    QWidget *widget1;
    QHBoxLayout *hboxLayout1;
    QCheckBox *VisibleCheck;
    QCheckBox *TurnOrientationCheck;
    QSpacerItem *spacerItem2;

    void setupUi(QWidget *RectangularSurfacePropertySheet_form)
    {
        if (RectangularSurfacePropertySheet_form->objectName().isEmpty())
            RectangularSurfacePropertySheet_form->setObjectName(QString::fromUtf8("RectangularSurfacePropertySheet_form"));
        RectangularSurfacePropertySheet_form->resize(544, 496);
        box = new Q3GroupBox(RectangularSurfacePropertySheet_form);
        box->setObjectName(QString::fromUtf8("box"));
        box->setGeometry(QRect(0, 0, 420, 220));
        box->setOrientation(Qt::Vertical);
        UresLabel = new QLabel(box);
        UresLabel->setObjectName(QString::fromUtf8("UresLabel"));
        UresLabel->setGeometry(QRect(335, 77, 63, 23));
        UresLabel->setWordWrap(false);
        UresSlider = new QSlider(box);
        UresSlider->setObjectName(QString::fromUtf8("UresSlider"));
        UresSlider->setGeometry(QRect(22, 77, 237, 23));
        UresSlider->setMinimum(2);
        UresSlider->setMaximum(120);
        UresSlider->setValue(10);
        UresSlider->setOrientation(Qt::Horizontal);
        UresLCDNumber = new QLCDNumber(box);
        UresLCDNumber->setObjectName(QString::fromUtf8("UresLCDNumber"));
        UresLCDNumber->setGeometry(QRect(265, 77, 64, 23));
        UresLCDNumber->setFrameShape(QFrame::Box);
        UresLCDNumber->setFrameShadow(QFrame::Raised);
        UresLCDNumber->setProperty("intValue", QVariant(10));
        VresSlider = new QSlider(box);
        VresSlider->setObjectName(QString::fromUtf8("VresSlider"));
        VresSlider->setGeometry(QRect(22, 116, 238, 23));
        VresSlider->setMinimum(2);
        VresSlider->setMaximum(120);
        VresSlider->setValue(10);
        VresSlider->setOrientation(Qt::Horizontal);
        VresLCDNumber = new QLCDNumber(box);
        VresLCDNumber->setObjectName(QString::fromUtf8("VresLCDNumber"));
        VresLCDNumber->setGeometry(QRect(266, 116, 64, 23));
        VresLCDNumber->setProperty("intValue", QVariant(10));
        VresLabel = new QLabel(box);
        VresLabel->setObjectName(QString::fromUtf8("VresLabel"));
        VresLabel->setGeometry(QRect(336, 116, 62, 23));
        VresLabel->setWordWrap(false);
        widget = new QWidget(box);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(22, 155, 381, 30));
        hboxLayout = new QHBoxLayout(widget);
#ifndef Q_OS_MAC
        hboxLayout->setSpacing(6);
#endif
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        spacerItem = new QSpacerItem(20, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout->addItem(spacerItem);

        ApplyButton = new QPushButton(widget);
        ApplyButton->setObjectName(QString::fromUtf8("ApplyButton"));

        hboxLayout->addWidget(ApplyButton);

        CloseButton = new QPushButton(widget);
        CloseButton->setObjectName(QString::fromUtf8("CloseButton"));

        hboxLayout->addWidget(CloseButton);

        spacerItem1 = new QSpacerItem(20, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout->addItem(spacerItem1);

        widget1 = new QWidget(box);
        widget1->setObjectName(QString::fromUtf8("widget1"));
        widget1->setGeometry(QRect(21, 41, 381, 30));
        hboxLayout1 = new QHBoxLayout(widget1);
#ifndef Q_OS_MAC
        hboxLayout1->setSpacing(6);
#endif
        hboxLayout1->setContentsMargins(0, 0, 0, 0);
        hboxLayout1->setObjectName(QString::fromUtf8("hboxLayout1"));
        VisibleCheck = new QCheckBox(widget1);
        VisibleCheck->setObjectName(QString::fromUtf8("VisibleCheck"));

        hboxLayout1->addWidget(VisibleCheck);

        TurnOrientationCheck = new QCheckBox(widget1);
        TurnOrientationCheck->setObjectName(QString::fromUtf8("TurnOrientationCheck"));

        hboxLayout1->addWidget(TurnOrientationCheck);

        spacerItem2 = new QSpacerItem(20, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout1->addItem(spacerItem2);


        retranslateUi(RectangularSurfacePropertySheet_form);
        QObject::connect(UresSlider, SIGNAL(valueChanged(int)), UresLCDNumber, SLOT(display(int)));
        QObject::connect(VresSlider, SIGNAL(valueChanged(int)), VresLCDNumber, SLOT(display(int)));

        QMetaObject::connectSlotsByName(RectangularSurfacePropertySheet_form);
    } // setupUi

    void retranslateUi(QWidget *RectangularSurfacePropertySheet_form)
    {
        RectangularSurfacePropertySheet_form->setWindowTitle(QApplication::translate("RectangularSurfacePropertySheet_form", "Rectangular surface properties", 0, QApplication::UnicodeUTF8));
        box->setTitle(QApplication::translate("RectangularSurfacePropertySheet_form", "Spline surface properties", 0, QApplication::UnicodeUTF8));
        UresLabel->setText(QApplication::translate("RectangularSurfacePropertySheet_form", "U resolution", 0, QApplication::UnicodeUTF8));
        VresLabel->setText(QApplication::translate("RectangularSurfacePropertySheet_form", "V resolution", 0, QApplication::UnicodeUTF8));
        ApplyButton->setText(QApplication::translate("RectangularSurfacePropertySheet_form", "Apply", 0, QApplication::UnicodeUTF8));
        CloseButton->setText(QApplication::translate("RectangularSurfacePropertySheet_form", "Close", 0, QApplication::UnicodeUTF8));
        VisibleCheck->setText(QApplication::translate("RectangularSurfacePropertySheet_form", "Visible", 0, QApplication::UnicodeUTF8));
        TurnOrientationCheck->setText(QApplication::translate("RectangularSurfacePropertySheet_form", "Turn orientation", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class RectangularSurfacePropertySheet_form: public Ui_RectangularSurfacePropertySheet_form {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RECTANGULARSURFACEPROPERTYSHEET_FORM_H
