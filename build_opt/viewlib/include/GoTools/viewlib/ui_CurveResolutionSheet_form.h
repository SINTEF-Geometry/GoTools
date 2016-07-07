/********************************************************************************
** Form generated from reading UI file 'CurveResolutionSheet_form.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CURVERESOLUTIONSHEET_FORM_H
#define UI_CURVERESOLUTIONSHEET_FORM_H

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

class Ui_CurveResolutionSheet_form
{
public:
    Q3GroupBox *box;
    QLabel *TextLabel1;
    QSlider *ResSlider;
    QPushButton *button5000;
    QLCDNumber *LCDNumber1;
    QWidget *Layout1;
    QHBoxLayout *hboxLayout;
    QPushButton *OkButton;
    QPushButton *CancelButton;

    void setupUi(QDialog *CurveResolutionSheet_form)
    {
        if (CurveResolutionSheet_form->objectName().isEmpty())
            CurveResolutionSheet_form->setObjectName(QString::fromUtf8("CurveResolutionSheet_form"));
        CurveResolutionSheet_form->resize(510, 159);
        CurveResolutionSheet_form->setSizeGripEnabled(true);
        box = new Q3GroupBox(CurveResolutionSheet_form);
        box->setObjectName(QString::fromUtf8("box"));
        box->setGeometry(QRect(0, 10, 501, 131));
        box->setOrientation(Qt::Vertical);
        TextLabel1 = new QLabel(box);
        TextLabel1->setObjectName(QString::fromUtf8("TextLabel1"));
        TextLabel1->setGeometry(QRect(410, 30, 70, 20));
        TextLabel1->setWordWrap(false);
        ResSlider = new QSlider(box);
        ResSlider->setObjectName(QString::fromUtf8("ResSlider"));
        ResSlider->setGeometry(QRect(20, 30, 300, 20));
        ResSlider->setMinimum(2);
        ResSlider->setMaximum(10000);
        ResSlider->setValue(100);
        ResSlider->setOrientation(Qt::Horizontal);
        button5000 = new QPushButton(box);
        button5000->setObjectName(QString::fromUtf8("button5000"));
        button5000->setGeometry(QRect(340, 80, 129, 34));
        LCDNumber1 = new QLCDNumber(box);
        LCDNumber1->setObjectName(QString::fromUtf8("LCDNumber1"));
        LCDNumber1->setGeometry(QRect(330, 30, 64, 23));
        LCDNumber1->setProperty("value", QVariant(100));
        LCDNumber1->setProperty("intValue", QVariant(100));
        Layout1 = new QWidget(box);
        Layout1->setObjectName(QString::fromUtf8("Layout1"));
        Layout1->setGeometry(QRect(70, 70, 210, 50));
        hboxLayout = new QHBoxLayout(Layout1);
        hboxLayout->setSpacing(6);
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        OkButton = new QPushButton(Layout1);
        OkButton->setObjectName(QString::fromUtf8("OkButton"));
        OkButton->setAutoDefault(true);
        OkButton->setDefault(true);

        hboxLayout->addWidget(OkButton);

        CancelButton = new QPushButton(Layout1);
        CancelButton->setObjectName(QString::fromUtf8("CancelButton"));
        CancelButton->setAutoDefault(true);

        hboxLayout->addWidget(CancelButton);


        retranslateUi(CurveResolutionSheet_form);
        QObject::connect(OkButton, SIGNAL(clicked()), CurveResolutionSheet_form, SLOT(accept()));
        QObject::connect(CancelButton, SIGNAL(clicked()), CurveResolutionSheet_form, SLOT(reject()));
        QObject::connect(ResSlider, SIGNAL(valueChanged(int)), LCDNumber1, SLOT(display(int)));

        QMetaObject::connectSlotsByName(CurveResolutionSheet_form);
    } // setupUi

    void retranslateUi(QDialog *CurveResolutionSheet_form)
    {
        CurveResolutionSheet_form->setWindowTitle(QApplication::translate("CurveResolutionSheet_form", "Curve Resolution", 0, QApplication::UnicodeUTF8));
        box->setTitle(QApplication::translate("CurveResolutionSheet_form", "GroupBox3", 0, QApplication::UnicodeUTF8));
        TextLabel1->setText(QApplication::translate("CurveResolutionSheet_form", "Resolution", 0, QApplication::UnicodeUTF8));
        button5000->setText(QApplication::translate("CurveResolutionSheet_form", "5000", 0, QApplication::UnicodeUTF8));
        OkButton->setText(QApplication::translate("CurveResolutionSheet_form", "&OK", 0, QApplication::UnicodeUTF8));
        CancelButton->setText(QApplication::translate("CurveResolutionSheet_form", "&Cancel", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class CurveResolutionSheet_form: public Ui_CurveResolutionSheet_form {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CURVERESOLUTIONSHEET_FORM_H
