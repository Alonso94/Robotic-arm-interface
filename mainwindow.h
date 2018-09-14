#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSerialPort>
#include "solutions.h"
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_sr1_valueChanged(int value);
    void on_sr2_valueChanged(int value);
    void on_sr3_valueChanged(int value);
    void on_sr4_valueChanged(int value);
    void on_sr5_valueChanged(int value);
    void on_sr6_valueChanged(int value);
    void update_servos(QString);
    void on_calc_clicked();

    void on_go_clicked();

    void on_update_clicked();

private:
    Ui::MainWindow *ui;
    QSerialPort *arduino;
    static const quint16 vid=9025;
    static const quint16 pid=66;
    QString arduino_port;
    bool arduino_av;
    Solutions sol;
};

#endif // MAINWINDOW_H
