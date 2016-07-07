/********************************************************************************
** Form generated from reading UI file 'PointCloudPropertySheet_form.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_POINTCLOUDPROPERTYSHEET_FORM_H
#define UI_POINTCLOUDPROPERTYSHEET_FORM_H

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

class Ui_PointCloudPropertySheet_form
{
public:
    Q3GroupBox *box;
    QWidget *layoutWidget1;
    QHBoxLayout *hboxLayout;
    QSlider *RenderSlider;
    QLCDNumber *RenderLCDNumber;
    QLabel *RenderLabel;
    QWidget *layoutWidget2;
    QHBoxLayout *hboxLayout1;
    QCheckBox *VisibleCheck;
    QCheckBox *IdCheck;
    QSpacerItem *spacerItem;
    QWidget *layoutWidget3;
    QHBoxLayout *hboxLayout2;
    QSpacerItem *spacerItem1;
    QPushButton *ApplyButton;
    QPushButton *CloseButton;
    QSpacerItem *spacerItem2;
    QWidget *layoutWidget;
    QHBoxLayout *hboxLayout3;
    QSlider *PointsizeSlider;
    QLCDNumber *PointsizeLCDNumber;
    QLabel *PointsizeLabel;

    void setupUi(QWidget *PointCloudPropertySheet_form)
    {
        if (PointCloudPropertySheet_form->objectName().isEmpty())
            PointCloudPropertySheet_form->setObjectName(QString::fromUtf8("PointCloudPropertySheet_form"));
        PointCloudPropertySheet_form->resize(544, 496);
        QSizePolicy sizePolicy(static_cast<QSizePolicy::Policy>(7), static_cast<QSizePolicy::Policy>(7));
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(PointCloudPropertySheet_form->sizePolicy().hasHeightForWidth());
        PointCloudPropertySheet_form->setSizePolicy(sizePolicy);
        box = new Q3GroupBox(PointCloudPropertySheet_form);
        box->setObjectName(QString::fromUtf8("box"));
        box->setGeometry(QRect(10, 10, 284, 143));
        box->setAlignment(Qt::AlignLeading);
        box->setOrientation(Qt::Vertical);
        layoutWidget1 = new QWidget(box);
        layoutWidget1->setObjectName(QString::fromUtf8("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(20, 50, 251, 30));
        hboxLayout = new QHBoxLayout(layoutWidget1);
#ifndef Q_OS_MAC
        hboxLayout->setSpacing(6);
#endif
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        RenderSlider = new QSlider(layoutWidget1);
        RenderSlider->setObjectName(QString::fromUtf8("RenderSlider"));
        RenderSlider->setMinimum(0);
        RenderSlider->setMaximum(100);
        RenderSlider->setValue(100);
        RenderSlider->setOrientation(Qt::Horizontal);

        hboxLayout->addWidget(RenderSlider);

        RenderLCDNumber = new QLCDNumber(layoutWidget1);
        RenderLCDNumber->setObjectName(QString::fromUtf8("RenderLCDNumber"));
        RenderLCDNumber->setProperty("intValue", QVariant(100));

        hboxLayout->addWidget(RenderLCDNumber);

        RenderLabel = new QLabel(layoutWidget1);
        RenderLabel->setObjectName(QString::fromUtf8("RenderLabel"));
        RenderLabel->setWordWrap(false);

        hboxLayout->addWidget(RenderLabel);

        layoutWidget2 = new QWidget(box);
        layoutWidget2->setObjectName(QString::fromUtf8("layoutWidget2"));
        layoutWidget2->setGeometry(QRect(20, 90, 154, 30));
        hboxLayout1 = new QHBoxLayout(layoutWidget2);
#ifndef Q_OS_MAC
        hboxLayout1->setSpacing(6);
#endif
        hboxLayout1->setContentsMargins(0, 0, 0, 0);
        hboxLayout1->setObjectName(QString::fromUtf8("hboxLayout1"));
        VisibleCheck = new QCheckBox(layoutWidget2);
        VisibleCheck->setObjectName(QString::fromUtf8("VisibleCheck"));

        hboxLayout1->addWidget(VisibleCheck);

        IdCheck = new QCheckBox(layoutWidget2);
        IdCheck->setObjectName(QString::fromUtf8("IdCheck"));

        hboxLayout1->addWidget(IdCheck);

        spacerItem = new QSpacerItem(20, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout1->addItem(spacerItem);

        layoutWidget3 = new QWidget(box);
        layoutWidget3->setObjectName(QString::fromUtf8("layoutWidget3"));
        layoutWidget3->setGeometry(QRect(180, 100, 100, 30));
        hboxLayout2 = new QHBoxLayout(layoutWidget3);
#ifndef Q_OS_MAC
        hboxLayout2->setSpacing(6);
#endif
        hboxLayout2->setContentsMargins(0, 0, 0, 0);
        hboxLayout2->setObjectName(QString::fromUtf8("hboxLayout2"));
        spacerItem1 = new QSpacerItem(20, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout2->addItem(spacerItem1);

        ApplyButton = new QPushButton(layoutWidget3);
        ApplyButton->setObjectName(QString::fromUtf8("ApplyButton"));

        hboxLayout2->addWidget(ApplyButton);

        CloseButton = new QPushButton(layoutWidget3);
        CloseButton->setObjectName(QString::fromUtf8("CloseButton"));

        hboxLayout2->addWidget(CloseButton);

        spacerItem2 = new QSpacerItem(20, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout2->addItem(spacerItem2);

        layoutWidget = new QWidget(box);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(24, 20, 241, 30));
        hboxLayout3 = new QHBoxLayout(layoutWidget);
#ifndef Q_OS_MAC
        hboxLayout3->setSpacing(6);
#endif
        hboxLayout3->setContentsMargins(0, 0, 0, 0);
        hboxLayout3->setObjectName(QString::fromUtf8("hboxLayout3"));
        PointsizeSlider = new QSlider(layoutWidget);
        PointsizeSlider->setObjectName(QString::fromUtf8("PointsizeSlider"));
        PointsizeSlider->setMinimum(0);
        PointsizeSlider->setMaximum(100);
        PointsizeSlider->setValue(100);
        PointsizeSlider->setOrientation(Qt::Horizontal);

        hboxLayout3->addWidget(PointsizeSlider);

        PointsizeLCDNumber = new QLCDNumber(layoutWidget);
        PointsizeLCDNumber->setObjectName(QString::fromUtf8("PointsizeLCDNumber"));
        PointsizeLCDNumber->setProperty("intValue", QVariant(100));

        hboxLayout3->addWidget(PointsizeLCDNumber);

        PointsizeLabel = new QLabel(layoutWidget);
        PointsizeLabel->setObjectName(QString::fromUtf8("PointsizeLabel"));
        PointsizeLabel->setWordWrap(false);

        hboxLayout3->addWidget(PointsizeLabel);


        retranslateUi(PointCloudPropertySheet_form);
        QObject::connect(RenderSlider, SIGNAL(valueChanged(int)), RenderLCDNumber, SLOT(display(int)));
        QObject::connect(PointsizeSlider, SIGNAL(valueChanged(int)), PointsizeLCDNumber, SLOT(display(int)));

        QMetaObject::connectSlotsByName(PointCloudPropertySheet_form);
    } // setupUi

    void retranslateUi(QWidget *PointCloudPropertySheet_form)
    {
        PointCloudPropertySheet_form->setWindowTitle(QApplication::translate("PointCloudPropertySheet_form", "Point cloud properties", 0, QApplication::UnicodeUTF8));
        box->setTitle(QApplication::translate("PointCloudPropertySheet_form", "Point cloud properties", 0, QApplication::UnicodeUTF8));
        RenderLabel->setText(QApplication::translate("PointCloudPropertySheet_form", "Render percentage", 0, QApplication::UnicodeUTF8));
        VisibleCheck->setText(QApplication::translate("PointCloudPropertySheet_form", "Visible", 0, QApplication::UnicodeUTF8));
        IdCheck->setText(QApplication::translate("PointCloudPropertySheet_form", "Show id", 0, QApplication::UnicodeUTF8));
        ApplyButton->setText(QApplication::translate("PointCloudPropertySheet_form", "Apply", 0, QApplication::UnicodeUTF8));
        CloseButton->setText(QApplication::translate("PointCloudPropertySheet_form", "Close", 0, QApplication::UnicodeUTF8));
        PointsizeLabel->setText(QApplication::translate("PointCloudPropertySheet_form", "Point size", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class PointCloudPropertySheet_form: public Ui_PointCloudPropertySheet_form {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_POINTCLOUDPROPERTYSHEET_FORM_H
