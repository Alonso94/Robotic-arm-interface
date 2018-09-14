#include "solutions.h"
#include "ui_solutions.h"

#include <QDebug>
#include <QtWidgets>

Solutions::Solutions(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Solutions)
{
    ui->setupUi(this);
}

Solutions::~Solutions()
{
    delete ui;
}


void Solutions::on_pushButton_clicked()
{
    QButtonGroup group;
    QList<QRadioButton*> all=this->findChildren<QRadioButton*>();
    for(int i = 0; i < all.size(); ++i)
    {
        group.addButton(all[i],i);
    }
    id=group.checkedId();
    this->close();
}
