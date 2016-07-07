/********************************************************************************
** Form generated from reading UI file 'SurfaceResolutionSheet_form.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SURFACERESOLUTIONSHEET_FORM_H
#define UI_SURFACERESOLUTIONSHEET_FORM_H

#include <Qt3Support/Q3GroupBox>
#include <Qt3Support/Q3MimeSourceFactory>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLCDNumber>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSlider>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SurfaceResolutionSheet_form
{
public:
    Q3GroupBox *box;
    QLabel *TextLabel1_2;
    QLabel *TextLabel1;
    QSlider *UresSlider;
    QSlider *VresSlider;
    QWidget *Layout1;
    QHBoxLayout *hboxLayout;
    QPushButton *OkButton;
    QPushButton *CancelButton;
    QPushButton *button200x200;
    QLCDNumber *LCDNumberU;
    QLCDNumber *LCDNumberV;

    void setupUi(QDialog *SurfaceResolutionSheet_form)
    {
        if (SurfaceResolutionSheet_form->objectName().isEmpty())
            SurfaceResolutionSheet_form->setObjectName(QString::fromUtf8("SurfaceResolutionSheet_form"));
        SurfaceResolutionSheet_form->resize(500, 167);
        SurfaceResolutionSheet_form->setSizeGripEnabled(true);
        box = new Q3GroupBox(SurfaceResolutionSheet_form);
        box->setObjectName(QString::fromUtf8("box"));
        box->setGeometry(QRect(10, 0, 491, 161));
        TextLabel1_2 = new QLabel(box);
        TextLabel1_2->setObjectName(QString::fromUtf8("TextLabel1_2"));
        TextLabel1_2->setGeometry(QRect(400, 60, 80, 20));
        TextLabel1_2->setWordWrap(false);
        TextLabel1 = new QLabel(box);
        TextLabel1->setObjectName(QString::fromUtf8("TextLabel1"));
        TextLabel1->setGeometry(QRect(400, 20, 80, 20));
        TextLabel1->setWordWrap(false);
        UresSlider = new QSlider(box);
        UresSlider->setObjectName(QString::fromUtf8("UresSlider"));
        UresSlider->setGeometry(QRect(10, 20, 300, 20));
        UresSlider->setMinimum(2);
        UresSlider->setMaximum(1000);
        UresSlider->setValue(20);
        UresSlider->setOrientation(Qt::Horizontal);
        VresSlider = new QSlider(box);
        VresSlider->setObjectName(QString::fromUtf8("VresSlider"));
        VresSlider->setGeometry(QRect(10, 60, 300, 20));
        VresSlider->setMinimum(2);
        VresSlider->setMaximum(1000);
        VresSlider->setValue(20);
        VresSlider->setOrientation(Qt::Horizontal);
        Layout1 = new QWidget(box);
        Layout1->setObjectName(QString::fromUtf8("Layout1"));
        Layout1->setGeometry(QRect(60, 100, 210, 50));
        hboxLayout = new QHBoxLayout(Layout1);
        hboxLayout->setSpacing(6);
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        OkButton = new QPushButton(Layout1);
        OkButton->setObjectName(QString::fromUtf8("OkButton"));
        OkButton->setAutoDefault(true);
        OkButton->setDefault(true);

        hboxLayout->addWidget(OkButton);

        CancelButton = new QPushButton(Layout1);
        CancelButton->setObjectName(QString::fromUtf8("CancelButton"));
        CancelButton->setAutoDefault(true);

        hboxLayout->addWidget(CancelButton);

        button200x200 = new QPushButton(box);
        button200x200->setObjectName(QString::fromUtf8("button200x200"));
        button200x200->setGeometry(QRect(330, 110, 129, 34));
        LCDNumberU = new QLCDNumber(box);
        LCDNumberU->setObjectName(QString::fromUtf8("LCDNumberU"));
        LCDNumberU->setGeometry(QRect(320, 20, 64, 23));
        LCDNumberU->setProperty("value", QVariant(20));
        LCDNumberU->setProperty("intValue", QVariant(20));
        LCDNumberV = new QLCDNumber(box);
        LCDNumberV->setObjectName(QString::fromUtf8("LCDNumberV"));
        LCDNumberV->setGeometry(QRect(320, 60, 64, 23));
        LCDNumberV->setProperty("value", QVariant(20));
        LCDNumberV->setProperty("intValue", QVariant(20));

        retranslateUi(SurfaceResolutionSheet_form);
        QObject::connect(OkButton, SIGNAL(clicked()), SurfaceResolutionSheet_form, SLOT(accept()));
        QObject::connect(CancelButton, SIGNAL(clicked()), SurfaceResolutionSheet_form, SLOT(reject()));
        QObject::connect(UresSlider, SIGNAL(valueChanged(int)), LCDNumberU, SLOT(display(int)));
        QObject::connect(VresSlider, SIGNAL(valueChanged(int)), LCDNumberV, SLOT(display(int)));

        QMetaObject::connectSlotsByName(SurfaceResolutionSheet_form);
    } // setupUi

    void retranslateUi(QDialog *SurfaceResolutionSheet_form)
    {
        SurfaceResolutionSheet_form->setWindowTitle(QApplication::translate("SurfaceResolutionSheet_form", "Surface Resolution", 0, QApplication::UnicodeUTF8));
        box->setTitle(QApplication::translate("SurfaceResolutionSheet_form", "GroupBox2", 0, QApplication::UnicodeUTF8));
        TextLabel1_2->setText(QApplication::translate("SurfaceResolutionSheet_form", "V resolution", 0, QApplication::UnicodeUTF8));
        TextLabel1->setText(QApplication::translate("SurfaceResolutionSheet_form", "U resolution", 0, QApplication::UnicodeUTF8));
        OkButton->setText(QApplication::translate("SurfaceResolutionSheet_form", "&OK", 0, QApplication::UnicodeUTF8));
        CancelButton->setText(QApplication::translate("SurfaceResolutionSheet_form", "&Cancel", 0, QApplication::UnicodeUTF8));
        CancelButton->setShortcut(QApplication::translate("SurfaceResolutionSheet_form", "Alt+C", 0, QApplication::UnicodeUTF8));
        button200x200->setText(QApplication::translate("SurfaceResolutionSheet_form", "200x200", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class SurfaceResolutionSheet_form: public Ui_SurfaceResolutionSheet_form {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SURFACERESOLUTIONSHEET_FORM_H
