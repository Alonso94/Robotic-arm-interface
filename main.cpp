#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.setFixedSize(760,350);
    w.setWindowTitle("Robotic Arm");
    w.show();

    return a.exec();
}
