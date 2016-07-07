/********************************************************************************
** Form generated from reading UI file 'gvGroupPropertySheet_form.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GVGROUPPROPERTYSHEET_FORM_H
#define UI_GVGROUPPROPERTYSHEET_FORM_H

#include <Qt3Support/Q3Frame>
#include <Qt3Support/Q3GroupBox>
#include <Qt3Support/Q3MimeSourceFactory>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_gvGroupPropertySheet_form
{
public:
    Q3GroupBox *box;
    Q3Frame *Frame14;
    Q3Frame *TolGroup;
    QLineEdit *GroupName;
    QLabel *Name;
    QPushButton *ApplyButton;
    QPushButton *CloseButton;

    void setupUi(QWidget *gvGroupPropertySheet_form)
    {
        if (gvGroupPropertySheet_form->objectName().isEmpty())
            gvGroupPropertySheet_form->setObjectName(QString::fromUtf8("gvGroupPropertySheet_form"));
        gvGroupPropertySheet_form->resize(289, 119);
        box = new Q3GroupBox(gvGroupPropertySheet_form);
        box->setObjectName(QString::fromUtf8("box"));
        box->setGeometry(QRect(10, 0, 280, 110));
        Frame14 = new Q3Frame(box);
        Frame14->setObjectName(QString::fromUtf8("Frame14"));
        Frame14->setGeometry(QRect(10, 330, 260, 50));
        Frame14->setFrameShape(Q3Frame::NoFrame);
        Frame14->setFrameShadow(Q3Frame::Raised);
        TolGroup = new Q3Frame(box);
        TolGroup->setObjectName(QString::fromUtf8("TolGroup"));
        TolGroup->setGeometry(QRect(10, 30, 260, 80));
        TolGroup->setFrameShape(Q3Frame::NoFrame);
        TolGroup->setFrameShadow(Q3Frame::Raised);
        GroupName = new QLineEdit(TolGroup);
        GroupName->setObjectName(QString::fromUtf8("GroupName"));
        GroupName->setGeometry(QRect(90, 0, 120, 25));
        Name = new QLabel(TolGroup);
        Name->setObjectName(QString::fromUtf8("Name"));
        Name->setGeometry(QRect(30, 0, 50, 20));
        Name->setWordWrap(false);
        ApplyButton = new QPushButton(TolGroup);
        ApplyButton->setObjectName(QString::fromUtf8("ApplyButton"));
        ApplyButton->setGeometry(QRect(0, 40, 119, 35));
        CloseButton = new QPushButton(TolGroup);
        CloseButton->setObjectName(QString::fromUtf8("CloseButton"));
        CloseButton->setGeometry(QRect(130, 40, 119, 35));

        retranslateUi(gvGroupPropertySheet_form);

        QMetaObject::connectSlotsByName(gvGroupPropertySheet_form);
    } // setupUi

    void retranslateUi(QWidget *gvGroupPropertySheet_form)
    {
        gvGroupPropertySheet_form->setWindowTitle(QApplication::translate("gvGroupPropertySheet_form", "gvGroup Properties", 0, QApplication::UnicodeUTF8));
        box->setTitle(QApplication::translate("gvGroupPropertySheet_form", "gvGroup Properties", 0, QApplication::UnicodeUTF8));
        Name->setText(QApplication::translate("gvGroupPropertySheet_form", "Name", 0, QApplication::UnicodeUTF8));
        ApplyButton->setText(QApplication::translate("gvGroupPropertySheet_form", "Apply", 0, QApplication::UnicodeUTF8));
        CloseButton->setText(QApplication::translate("gvGroupPropertySheet_form", "Close", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class gvGroupPropertySheet_form: public Ui_gvGroupPropertySheet_form {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GVGROUPPROPERTYSHEET_FORM_H
