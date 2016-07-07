/********************************************************************************
** Form generated from reading UI file 'SplineCurvePropertySheet_form.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SPLINECURVEPROPERTYSHEET_FORM_H
#define UI_SPLINECURVEPROPERTYSHEET_FORM_H

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
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SplineCurvePropertySheet_form
{
public:
    Q3GroupBox *box;
    QWidget *layout13;
    QVBoxLayout *vboxLayout;
    QHBoxLayout *hboxLayout;
    QCheckBox *VisibleCheck;
    QSpacerItem *Spacer2_5;
    QHBoxLayout *hboxLayout1;
    QSlider *ResSlider;
    QLCDNumber *ResLCDNumber;
    QLabel *ResLabel;
    QHBoxLayout *hboxLayout2;
    QLabel *ColorLabel;
    QSpacerItem *spacer4;
    QLabel *redLabel;
    QSpinBox *redSpinBox;
    QLabel *greenLabel;
    QSpinBox *greenSpinBox;
    QLabel *blueLabel;
    QSpinBox *blueSpinBox;
    QLabel *alphaLabel;
    QSpinBox *alphaSpinBox;
    QHBoxLayout *hboxLayout3;
    QSpacerItem *Spacer2;
    QPushButton *ApplyButton;
    QPushButton *CloseButton;
    QSpacerItem *Spacer2_2;

    void setupUi(QWidget *SplineCurvePropertySheet_form)
    {
        if (SplineCurvePropertySheet_form->objectName().isEmpty())
            SplineCurvePropertySheet_form->setObjectName(QString::fromUtf8("SplineCurvePropertySheet_form"));
        SplineCurvePropertySheet_form->resize(588, 480);
        box = new Q3GroupBox(SplineCurvePropertySheet_form);
        box->setObjectName(QString::fromUtf8("box"));
        box->setGeometry(QRect(30, 20, 420, 210));
        layout13 = new QWidget(box);
        layout13->setObjectName(QString::fromUtf8("layout13"));
        layout13->setGeometry(QRect(0, 30, 410, 160));
        vboxLayout = new QVBoxLayout(layout13);
        vboxLayout->setSpacing(6);
        vboxLayout->setContentsMargins(0, 0, 0, 0);
        vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
        vboxLayout->setContentsMargins(0, 0, 0, 0);
        hboxLayout = new QHBoxLayout();
        hboxLayout->setSpacing(6);
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        VisibleCheck = new QCheckBox(layout13);
        VisibleCheck->setObjectName(QString::fromUtf8("VisibleCheck"));

        hboxLayout->addWidget(VisibleCheck);

        Spacer2_5 = new QSpacerItem(20, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout->addItem(Spacer2_5);


        vboxLayout->addLayout(hboxLayout);

        hboxLayout1 = new QHBoxLayout();
        hboxLayout1->setSpacing(6);
        hboxLayout1->setContentsMargins(0, 0, 0, 0);
        hboxLayout1->setObjectName(QString::fromUtf8("hboxLayout1"));
        ResSlider = new QSlider(layout13);
        ResSlider->setObjectName(QString::fromUtf8("ResSlider"));
        ResSlider->setMinimum(2);
        ResSlider->setMaximum(500);
        ResSlider->setValue(10);
        ResSlider->setOrientation(Qt::Horizontal);

        hboxLayout1->addWidget(ResSlider);

        ResLCDNumber = new QLCDNumber(layout13);
        ResLCDNumber->setObjectName(QString::fromUtf8("ResLCDNumber"));
        ResLCDNumber->setProperty("intValue", QVariant(10));

        hboxLayout1->addWidget(ResLCDNumber);

        ResLabel = new QLabel(layout13);
        ResLabel->setObjectName(QString::fromUtf8("ResLabel"));
        ResLabel->setWordWrap(false);

        hboxLayout1->addWidget(ResLabel);


        vboxLayout->addLayout(hboxLayout1);

        hboxLayout2 = new QHBoxLayout();
        hboxLayout2->setSpacing(6);
        hboxLayout2->setContentsMargins(0, 0, 0, 0);
        hboxLayout2->setObjectName(QString::fromUtf8("hboxLayout2"));
        ColorLabel = new QLabel(layout13);
        ColorLabel->setObjectName(QString::fromUtf8("ColorLabel"));
        ColorLabel->setWordWrap(false);

        hboxLayout2->addWidget(ColorLabel);

        spacer4 = new QSpacerItem(20, 30, QSizePolicy::Minimum, QSizePolicy::Expanding);

        hboxLayout2->addItem(spacer4);

        redLabel = new QLabel(layout13);
        redLabel->setObjectName(QString::fromUtf8("redLabel"));
        redLabel->setWordWrap(false);

        hboxLayout2->addWidget(redLabel);

        redSpinBox = new QSpinBox(layout13);
        redSpinBox->setObjectName(QString::fromUtf8("redSpinBox"));
        redSpinBox->setMaximum(255);

        hboxLayout2->addWidget(redSpinBox);

        greenLabel = new QLabel(layout13);
        greenLabel->setObjectName(QString::fromUtf8("greenLabel"));
        greenLabel->setWordWrap(false);

        hboxLayout2->addWidget(greenLabel);

        greenSpinBox = new QSpinBox(layout13);
        greenSpinBox->setObjectName(QString::fromUtf8("greenSpinBox"));
        greenSpinBox->setMaximum(255);

        hboxLayout2->addWidget(greenSpinBox);

        blueLabel = new QLabel(layout13);
        blueLabel->setObjectName(QString::fromUtf8("blueLabel"));
        blueLabel->setWordWrap(false);

        hboxLayout2->addWidget(blueLabel);

        blueSpinBox = new QSpinBox(layout13);
        blueSpinBox->setObjectName(QString::fromUtf8("blueSpinBox"));
        blueSpinBox->setMaximum(255);

        hboxLayout2->addWidget(blueSpinBox);

        alphaLabel = new QLabel(layout13);
        alphaLabel->setObjectName(QString::fromUtf8("alphaLabel"));
        alphaLabel->setWordWrap(false);

        hboxLayout2->addWidget(alphaLabel);

        alphaSpinBox = new QSpinBox(layout13);
        alphaSpinBox->setObjectName(QString::fromUtf8("alphaSpinBox"));
        alphaSpinBox->setMaximum(255);

        hboxLayout2->addWidget(alphaSpinBox);


        vboxLayout->addLayout(hboxLayout2);

        hboxLayout3 = new QHBoxLayout();
        hboxLayout3->setSpacing(6);
        hboxLayout3->setContentsMargins(0, 0, 0, 0);
        hboxLayout3->setObjectName(QString::fromUtf8("hboxLayout3"));
        Spacer2 = new QSpacerItem(20, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout3->addItem(Spacer2);

        ApplyButton = new QPushButton(layout13);
        ApplyButton->setObjectName(QString::fromUtf8("ApplyButton"));

        hboxLayout3->addWidget(ApplyButton);

        CloseButton = new QPushButton(layout13);
        CloseButton->setObjectName(QString::fromUtf8("CloseButton"));

        hboxLayout3->addWidget(CloseButton);

        Spacer2_2 = new QSpacerItem(20, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout3->addItem(Spacer2_2);


        vboxLayout->addLayout(hboxLayout3);


        retranslateUi(SplineCurvePropertySheet_form);
        QObject::connect(ResSlider, SIGNAL(valueChanged(int)), ResLCDNumber, SLOT(display(int)));

        QMetaObject::connectSlotsByName(SplineCurvePropertySheet_form);
    } // setupUi

    void retranslateUi(QWidget *SplineCurvePropertySheet_form)
    {
        SplineCurvePropertySheet_form->setWindowTitle(QApplication::translate("SplineCurvePropertySheet_form", "Splinecurve properties", 0, QApplication::UnicodeUTF8));
        box->setTitle(QApplication::translate("SplineCurvePropertySheet_form", "Spline curve properties", 0, QApplication::UnicodeUTF8));
        VisibleCheck->setText(QApplication::translate("SplineCurvePropertySheet_form", "Visible", 0, QApplication::UnicodeUTF8));
        ResLabel->setText(QApplication::translate("SplineCurvePropertySheet_form", "Resolution", 0, QApplication::UnicodeUTF8));
        ColorLabel->setText(QApplication::translate("SplineCurvePropertySheet_form", "Color", 0, QApplication::UnicodeUTF8));
        redLabel->setText(QApplication::translate("SplineCurvePropertySheet_form", "R", 0, QApplication::UnicodeUTF8));
        greenLabel->setText(QApplication::translate("SplineCurvePropertySheet_form", "G", 0, QApplication::UnicodeUTF8));
        blueLabel->setText(QApplication::translate("SplineCurvePropertySheet_form", "B", 0, QApplication::UnicodeUTF8));
        alphaLabel->setText(QApplication::translate("SplineCurvePropertySheet_form", "A", 0, QApplication::UnicodeUTF8));
        ApplyButton->setText(QApplication::translate("SplineCurvePropertySheet_form", "Apply", 0, QApplication::UnicodeUTF8));
        CloseButton->setText(QApplication::translate("SplineCurvePropertySheet_form", "Close", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class SplineCurvePropertySheet_form: public Ui_SplineCurvePropertySheet_form {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SPLINECURVEPROPERTYSHEET_FORM_H
