#ifndef SOLUTIONS_H
#define SOLUTIONS_H

#include <QWidget>

namespace Ui {
class Solutions;
}

class Solutions : public QWidget
{
    Q_OBJECT

public:
    explicit Solutions(QWidget *parent = 0);
    ~Solutions();
    Ui::Solutions *ui;
    int id=-1;
private slots:
    void on_pushButton_clicked();
};

#endif // SOLUTIONS_H
