/********************************************************************************
** Form generated from reading UI file 'RectangularVolumePropertySheet_form.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_RECTANGULARVOLUMEPROPERTYSHEET_FORM_H
#define UI_RECTANGULARVOLUMEPROPERTYSHEET_FORM_H

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

class Ui_RectangularVolumePropertySheet_form
{
public:
    Q3GroupBox *box;
    QLabel *TextLabel1;
    QSlider *ResSlider;
    QLCDNumber *LCDNumber1;
    QWidget *Layout1;
    QHBoxLayout *hboxLayout;
    QPushButton *ApplyButton;
    QPushButton *CloseButton;

    void setupUi(QDialog *RectangularVolumePropertySheet_form)
    {
        if (RectangularVolumePropertySheet_form->objectName().isEmpty())
            RectangularVolumePropertySheet_form->setObjectName(QString::fromUtf8("RectangularVolumePropertySheet_form"));
        RectangularVolumePropertySheet_form->resize(510, 159);
        RectangularVolumePropertySheet_form->setSizeGripEnabled(true);
        box = new Q3GroupBox(RectangularVolumePropertySheet_form);
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
        ResSlider->setMaximum(1000);
        ResSlider->setValue(3);
        ResSlider->setOrientation(Qt::Horizontal);
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
        ApplyButton = new QPushButton(Layout1);
        ApplyButton->setObjectName(QString::fromUtf8("ApplyButton"));
        ApplyButton->setAutoDefault(true);
        ApplyButton->setDefault(true);

        hboxLayout->addWidget(ApplyButton);

        CloseButton = new QPushButton(Layout1);
        CloseButton->setObjectName(QString::fromUtf8("CloseButton"));
        CloseButton->setAutoDefault(true);

        hboxLayout->addWidget(CloseButton);


        retranslateUi(RectangularVolumePropertySheet_form);
        QObject::connect(ResSlider, SIGNAL(valueChanged(int)), LCDNumber1, SLOT(display(int)));

        QMetaObject::connectSlotsByName(RectangularVolumePropertySheet_form);
    } // setupUi

    void retranslateUi(QDialog *RectangularVolumePropertySheet_form)
    {
        RectangularVolumePropertySheet_form->setWindowTitle(QApplication::translate("RectangularVolumePropertySheet_form", "Volume Resolution", 0, QApplication::UnicodeUTF8));
        box->setTitle(QApplication::translate("RectangularVolumePropertySheet_form", "GroupBox3", 0, QApplication::UnicodeUTF8));
        TextLabel1->setText(QApplication::translate("RectangularVolumePropertySheet_form", "Resolution", 0, QApplication::UnicodeUTF8));
        ApplyButton->setText(QApplication::translate("RectangularVolumePropertySheet_form", "&Apply", 0, QApplication::UnicodeUTF8));
        CloseButton->setText(QApplication::translate("RectangularVolumePropertySheet_form", "&Close", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class RectangularVolumePropertySheet_form: public Ui_RectangularVolumePropertySheet_form {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RECTANGULARVOLUMEPROPERTYSHEET_FORM_H
