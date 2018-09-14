#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>
#include <QtWidgets>
#include <QSerialPort>
#include <QSerialPortInfo>
#include <QVector>

typedef long double ld;
typedef QVector<ld> vld;
ld EPS=1e-9,pi=acos(-1);

void zero(ld t[4][4]){
    for(int i=0;i<4;++i)
        for(int j=0;j<4;++j)if(fabs(t[i][j])<EPS)t[i][j]=0;
}
struct point{
    ld x,y,z;
    ld theta,psi,phi;
    ld beta,alpha,gamma;
    bool euler,RPY;
    point(){}
    point(ld x,ld y,ld z):x(x),y(y),z(z){}
    void calc(ld e[4][4]){
        x=e[0][3];
        y=e[1][3];
        z=e[2][3];
        euler=RPY=0;
        if(fabs(e[2][2]-1)>EPS&&fabs(e[2][2]+1)>EPS){
            euler=1;
            theta=atan2(e[2][2],hypot(e[0][2],e[1][2]));
            psi=atan2(e[0][2],-e[1][2]);
            phi=atan2(e[2][0],e[2][1]);
        }
        if(fabs(e[2][0]-1)>EPS&&fabs(e[2][0]+1)>EPS){
            RPY=1;
            beta=atan2(-e[2][0],hypot(e[0][0],e[1][0]));
            alpha=atan2(e[1][0],e[0][0]);
            gamma=atan2(e[2][1],e[2][2]);
        }
    }
};
struct joint{
    ld d,alpha,a,theta;
    joint(){}
    joint(ld d,ld a,ld alpha):d(d),a(a),alpha(alpha){}
    ld T[4][4];
    void calcT(){
        T[0][0]=cos(theta);
        T[0][1]=-sin(theta)*cos(alpha);
        T[0][2]=sin(theta)*sin(alpha);
        T[0][3]=a*cos(theta);
        T[1][0]=sin(theta);
        T[1][1]=cos(theta)*cos(alpha);
        T[1][2]=-cos(theta)*sin(alpha);
        T[1][3]=a*sin(theta);
        T[2][0]=0;
        T[2][1]=sin(alpha);
        T[2][2]=cos(alpha);
        T[2][3]=d;
        T[3][0]=T[3][1]=T[3][2]=0;
        T[3][3]=1;
        zero(T);
    }
}DH[8];
void mul(ld t1[4][4],ld t2[4][4],ld res[4][4]){
    for(int i=0;i<4;++i){
        for(int j=0;j<4;++j){
            res[i][j]=0;
            for(int k=0;k<4;++k)res[i][j]+=t1[i][k]*t2[k][j];
        }
    }
    zero(res);
}
void DKM(vld &v,point &p,ld endeff[4][4]){
    v[3]=v[7]=0;
    for(int i=0;i<8;++i){
        DH[i].theta=v[i];
        DH[i].calcT();
    }
    ld temp[4][4];
    mul(DH[0].T,DH[1].T,temp);
    for(int i=2;i<8;++i){
        mul(temp,DH[i].T,endeff);
        for(int k=0;k<4;++k)
            for(int j=0;j<4;++j)temp[k][j]=endeff[k][j];
    }
    p.calc(endeff);
}
void IKM(QVector<QVector<ld> >&sols,ld t08[4][4])
{
    QVector< QVector<ld> >ret(6);
    bool notacceptedsols[8]={0};
    ld t07[4][4],t87[4][4];
    DH[7].theta=0;
    DH[7].calcT();
    for(int i=0;i<4;++i)
        for(int j=0;j<4;++j)
            t87[i][j]=DH[7].T[i][j];
    t87[2][3]*=-1;
    mul(t08,t87,t07);
    ld dx=t07[0][3],dy=t07[1][3],dz=t07[2][3];
    ret[0].push_back(atan2(dy,dx));
    ret[0].push_back(atan2(-dy,-dx));
    ld l1=DH[0].d,l2=DH[1].a,l3=DH[3].d,l4=DH[7].d;
    if(fabs(ret[0][0])>pi/2+EPS)
        notacceptedsols[0]=notacceptedsols[1]=notacceptedsols[2]=notacceptedsols[3]=1;
    if(fabs(ret[0][1])>pi/2+EPS)
        notacceptedsols[4]=notacceptedsols[5]=notacceptedsols[6]=notacceptedsols[7]=1;
    for(int i=0;i<2;++i){
        ld a=dx*cos(ret[0][i])+dy*sin(ret[0][i]);
        ld b=dz-l1;
        ld c=(l2*l2+a*a+b*b-l3*l3)/(2*l2);
        if(fabs(a)<EPS)a=0;
        if(fabs(b)<EPS)b=0;
        if(fabs(c)<EPS)c=0;
        ret[1].push_back(atan2(a,b)-atan2(c,sqrt(a*a+b*b-c*c)));
        ret[2].push_back(atan2(a-l2*cos(ret[1].back()),l2*sin(ret[1].back())+b)-ret[1].back());
        if(fabs(ret[1].back())>pi/2+EPS||fabs(ret[2].back())>pi/2+EPS)
            notacceptedsols[2*i]=notacceptedsols[2*i+1]=1;
        ret[1].push_back(atan2(a,b)-atan2(c,-sqrt(a*a+b*b-c*c)));
        ret[2].push_back(atan2(a-l2*cos(ret[1].back()),l2*sin(ret[1].back())+b)-ret[1].back());
        if(fabs(ret[1].back())>pi/2+EPS||fabs(ret[2].back())>pi/2+EPS)
            notacceptedsols[2*i+2]=notacceptedsols[2*i+3]=1;
    }

    DH[3].theta=0;
    DH[3].calcT();
    for(int i=0;i<4;++i){
        DH[0].theta=ret[0][i/2];
        DH[1].theta=ret[1][i];
        DH[2].theta=ret[2][i];
        for(int j=0;j<3;++j)
            DH[j].calcT();
        ld t04[4][4],temp[4][4];
        mul(DH[0].T,DH[1].T,t04);
        mul(t04,DH[2].T,temp);
        mul(temp,DH[3].T,t04);
        ld t40[4][4];
        for(int i=0;i<3;++i)
            for(int j=0;j<3;++j)
                t40[i][j]=t04[j][i];
        for(int i=0;i<4;++i)
            t40[3][i]=t04[3][i];
        for(int i=0;i<3;++i){
            t40[i][3]=0;
            for(int j=0;j<3;++j)
                t40[i][3]-=t40[i][j]*t04[j][3];
        }
        zero(t40);
        ld k[4][4];
        mul(t40,t07,k);
        ret[4].push_back(atan2(hypot(k[2][0],k[2][1]),k[2][2]));
        if(fabs(ret[4].back())>EPS){
            ret[3].push_back(atan2(-k[1][2]/sin(ret[4].back()),-k[0][2]/sin(ret[4].back())));
            ret[5].push_back(atan2(-k[2][1]/sin(ret[4].back()),k[2][0]/sin(ret[4].back())));
        }
        else {
            ld theta=atan2(k[1][0],k[0][0]);
            ret[3].push_back(theta/2);
            ret[5].push_back(theta/2);
        }
        if(fabs(ret[3].back())>pi/2+EPS||fabs(ret[4].back())>pi/2+EPS||fabs(ret[5].back())>pi/2+EPS)
            notacceptedsols[2*i]=1;
        ret[4].push_back(atan2(-hypot(k[2][0],k[2][1]),k[2][2]));
        if(fabs(ret[4].back())>EPS){
            ret[3].push_back(atan2(-k[1][2]/sin(ret[4].back()),-k[0][2])/sin(ret[4].back()));
            ret[5].push_back(atan2(-k[2][1]/sin(ret[4].back()),k[2][0]/sin(ret[4].back())));
        }
        else {
            ld theta=atan2(k[1][0],k[0][0]);
            ret[3].push_back(theta/2);
            ret[5].push_back(theta/2);
        }
        if(fabs(ret[3].back())>pi/2+EPS||fabs(ret[4].back())>pi/2+EPS||fabs(ret[5].back())>pi/2+EPS)
            notacceptedsols[2*i+1]=1;
    }
    int ind=0;
    for(int i=0;i<8;++i){
        if(notacceptedsols[i])continue;
        sols.push_back(QVector<ld>(8,0.0));
        for(int j=0;j<6;++j){
            sols[ind][j+(j>=3)]=ret[j][i*ret[j].size()/8];
            if(sols[ind][j+(j>=3)]<-pi+EPS)
                sols[ind][j+(j>=3)]+=2*pi;
            if(sols[ind][j+(j>=3)]>pi+EPS)
                sols[ind][j+(j>=3)]-=2*pi;
        }
        ++ind;
    }
}


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    arduino_av=false;
    arduino_port="";
    arduino=new QSerialPort;
    qDebug() << "The number of available ports: " << QSerialPortInfo::availablePorts().length()<<endl;
    foreach(const QSerialPortInfo &serialPortInfo,QSerialPortInfo::availablePorts())
    {
        if(serialPortInfo.hasVendorIdentifier() && serialPortInfo.hasProductIdentifier())
            if(serialPortInfo.vendorIdentifier()==vid && serialPortInfo.productIdentifier()==pid)
            {
                arduino_port=serialPortInfo.portName();
                arduino_av=true;
            }
    }
    if(arduino_av)
    {
        arduino->setPortName(arduino_port);
        arduino->open(QSerialPort::WriteOnly);
        arduino->setBaudRate(QSerialPort::Baud9600);
        arduino->setDataBits(QSerialPort::Data8);
        arduino->setParity(QSerialPort::NoParity);
        arduino->setStopBits(QSerialPort::OneStop);
        arduino->setFlowControl(QSerialPort::NoFlowControl);
    }
    else
    {
        QMessageBox::warning(this,"Port error","Couldn't find arduino");
    }
    DH[0]=joint(107,0,-pi/2.0);
    DH[1]=joint(0,350,0);
    DH[2]=joint(0,0,pi/2.0);
    DH[3]=joint(270,0,0);
    DH[4]=joint(0,0,pi/2.0);
    DH[5]=joint(0,0,-pi/2.0);
    DH[6]=joint(0,0,0);
    DH[7]=joint(104,0,0);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_sr1_valueChanged(int value)
{
    ui->s1->setText(QString::number((double)(value-1350)*180/1800,'f',3)+QString(" deg"));
    MainWindow::update_servos(QString("a%1").arg(value));
}
void MainWindow::on_sr2_valueChanged(int value)
{
    ui->s2->setText(QString::number((double)(value-1450)*180/2100,'f',3)+QString(" deg"));
    MainWindow::update_servos(QString("b%1").arg(value));
}
void MainWindow::on_sr3_valueChanged(int value)
{
    ui->s3->setText(QString::number((double)(value-1500)*270/2000,'f',3)+QString(" deg"));
    MainWindow::update_servos(QString("c%1").arg(value));
}
void MainWindow::on_sr4_valueChanged(int value)
{
    ui->s4->setText(QString::number((double)(value-1500)*270/2000,'f',3)+QString(" deg"));
    MainWindow::update_servos(QString("d%1").arg(value));
}
void MainWindow::on_sr5_valueChanged(int value)
{
    ui->s5->setText(QString::number((double)(value-1811)*180/1600,'f',3)+QString(" deg"));
    MainWindow::update_servos(QString("e%1").arg(value));
}
void MainWindow::on_sr6_valueChanged(int value)
{
    ui->s6->setText(QString::number((double)(value-1500)*180/1600,'f',3)+QString(" deg"));
    MainWindow::update_servos(QString("f%1").arg(value));
}
void MainWindow::update_servos(QString command)
{
    qDebug()<<command<<endl;
    if(arduino->isWritable())
    {

        arduino->write(command.toStdString().c_str());
    }
    else
    {
       qDebug() <<"Couldn't write to serial";
    }
}

void MainWindow::on_calc_clicked()
{
    if(ui->s1->isModified()) ui->sr1->setValue(1350+1800*ui->s1->text().toDouble()/180);
    if(ui->s2->isModified()) ui->sr2->setValue(1450+2100*ui->s2->text().toDouble()/180);
    if(ui->s3->isModified()) ui->sr3->setValue(1500+2000*ui->s4->text().toDouble()/270);
    if(ui->s4->isModified()) ui->sr4->setValue(1500+2000*ui->s4->text().toDouble()/270);
    if(ui->s5->isModified()) ui->sr5->setValue(1811+1600*ui->s5->text().toDouble()/180);
    if(ui->s6->isModified()) ui->sr6->setValue(1500+1600*ui->s6->text().toDouble()/180);
    vld v(8,0.0);
    v[0]=(ui->sr1->value()-1350)*pi/1800;
    v[1]=(ui->sr2->value()-1450)*pi/2100;
    v[2]=(ui->sr3->value()-1500)*3*pi/(2*2000);
    v[4]=(ui->sr4->value()-1500)*3*pi/(2*2000);
    v[5]=(ui->sr5->value()-1811)*pi/1600;
    v[6]=(ui->sr6->value()-1500)*pi/1600;
    point p;
    ld t[4][4];
    DKM(v,p,t);
    ui->table->setItem(1,0,new QTableWidgetItem(QString::number(p.x,'f',3)));
    ui->table->setItem(1,1,new QTableWidgetItem(QString::number(p.y,'f',3)));
    ui->table->setItem(1,2,new QTableWidgetItem(QString::number(p.z,'f',3)));
    if(p.euler)
    {
        ui->table->setItem(2,0,new QTableWidgetItem("Theta"));
        ui->table->setItem(2,1,new QTableWidgetItem("Psi"));
        ui->table->setItem(2,2,new QTableWidgetItem("Phi"));
        ui->table->setItem(3,0,new QTableWidgetItem(QString::number((180/pi)*p.theta,'f',3)));
        ui->table->setItem(3,1,new QTableWidgetItem(QString::number((180/pi)*p.psi,'f',3)));
        ui->table->setItem(3,2,new QTableWidgetItem(QString::number((180/pi)*p.phi,'f',3)));
    }
    else if(p.RPY)
    {
        ui->table->setItem(2,0,new QTableWidgetItem("Beta"));
        ui->table->setItem(2,1,new QTableWidgetItem("Alpha"));
        ui->table->setItem(2,2,new QTableWidgetItem("Gamma"));
        ui->table->setItem(3,0,new QTableWidgetItem(QString::number((180/pi)*p.beta,'f',3)));
        ui->table->setItem(3,1,new QTableWidgetItem(QString::number((180/pi)*p.alpha,'f',3)));
        ui->table->setItem(3,2,new QTableWidgetItem(QString::number((180/pi)*p.gamma,'f',3)));
    }
}

void MainWindow::on_go_clicked()
{
    ld t[4][4];
    QVector<QVector<ld>> sols;
    t[0][3]=ui->x->text().toDouble();
    t[1][3]=ui->y->text().toDouble();
    t[2][3]=ui->z->text().toDouble();
    t[3][0]=t[3][1]=t[3][2]=0;
    t[3][3]=1;
    if(!ui->rpy->isChecked())
    {
        ld th=(pi/180)*ui->t->text().toDouble();
        ld p=(pi/180)*ui->p->text().toDouble();
        ld h=(pi/180)*ui->h->text().toDouble();
        t[0][0]=cos(p)*cos(h)-sin(p)*cos(th)*sin(h);
        t[0][1]=-cos(p)*sin(h)-sin(p)*cos(th)*cos(h);
        t[0][2]=sin(p)*sin(th);
        t[1][0]=sin(p)*cos(h)+cos(p)*cos(th)*sin(h);
        t[1][1]=-sin(p)*sin(h)+cos(p)*cos(th)*cos(h);
        t[1][2]=-cos(p)*sin(th);
        t[2][0]=sin(th)*sin(h);
        t[2][1]=sin(th)*cos(h);
        t[2][2]=cos(th);
    }
    else
    {
        ld b=ui->b->text().toDouble();
        ld al=ui->a->text().toDouble();
        ld g=ui->g->text().toDouble();
        t[0][0]=cos(al)*cos(b);
        t[0][1]=-sin(al)*cos(g)+cos(al)*sin(b)*sin(g);
        t[0][2]=sin(al)*sin(g)+cos(al)*sin(b)*cos(g);
        t[1][0]=sin(al)*cos(b);
        t[1][1]=cos(al)*cos(g)+sin(al)*sin(b)*sin(g);
        t[1][2]=-cos(al)*sin(g)+sin(al)*sin(b)*cos(g);
        t[2][0]=-sin(b);
        t[2][1]=cos(b)*sin(g);
        t[2][2]=cos(b)*cos(g);
    }
    zero(t);
    IKM(sols,t);
    if(sols.size()==0) {QMessageBox::warning(this,"Solutions","No Solutions");return;}
    QRadioButton *radio;
    sol.setWindowTitle("Solutions");
    sol.show();
    sol.findChild<QTextBrowser*>("tt")->clear();
    foreach (QRadioButton * w, sol.findChildren<QRadioButton *>() ) delete w;
    for(int i=0;i<sols.size();++i)
    {
       sol.findChild<QTextBrowser*>("tt")->append(QString("Solution# %1 : ").arg(i+1));
       radio=new QRadioButton(tr("Solution #%1").arg(i+1));
       sol.findChild<QVBoxLayout*>("vl")->addWidget(radio);
       for(int j=0;j<7;++j)
           if(j!=3)
               sol.findChild<QTextBrowser*>("tt")->append(QString::number((180/pi)*sols[i][j],'f',3));
       sol.findChild<QTextBrowser*>("tt")->append("\n");
    }
    QEventLoop loop;
    QObject::connect(sol.findChild<QPushButton*>("pushButton"), SIGNAL(clicked()), &loop, SLOT(quit()));
    loop.exec();
    int id=sol.id;
    if(id==-1) {QMessageBox::warning(this,"IKM","No Solution checked");return;}
    ui->sr1->setValue(1350+1800*sols[id][0]/pi);
    ui->sr2->setValue(1450+2100*sols[id][1]/pi);
    ui->sr3->setValue(1500+4000*sols[id][2]/(3*pi));
    ui->sr4->setValue(1500+4000*sols[id][4]/(3*pi));
    ui->sr5->setValue(1811+1600*sols[id][5]/pi);
    ui->sr6->setValue(1500+1600*sols[id][6]/pi);
}

void MainWindow::on_update_clicked()
{
    ld d=ui->DH->item(1,1)->text().toDouble();
    ld a=ui->DH->item(1,2)->text().toDouble();
    ld alpha=(pi/180)*ui->DH->item(1,3)->text().toDouble();
    DH[0]=joint(d,a,alpha);
    d=ui->DH->item(2,1)->text().toDouble();
    a=ui->DH->item(2,2)->text().toDouble();
    alpha=(pi/180)*ui->DH->item(2,3)->text().toDouble();
    DH[1]=joint(d,a,alpha);
    d=ui->DH->item(3,1)->text().toDouble();
    a=ui->DH->item(3,2)->text().toDouble();
    alpha=(pi/180)*ui->DH->item(3,3)->text().toDouble();
    DH[2]=joint(d,a,alpha);
    d=ui->DH->item(4,1)->text().toDouble();
    a=ui->DH->item(4,2)->text().toDouble();
    alpha=(pi/180)*ui->DH->item(4,3)->text().toDouble();
    DH[3]=joint(d,a,alpha);
    d=ui->DH->item(5,1)->text().toDouble();
    a=ui->DH->item(5,2)->text().toDouble();
    alpha=(pi/180)*ui->DH->item(5,3)->text().toDouble();
    DH[4]=joint(d,a,alpha);
    d=ui->DH->item(6,1)->text().toDouble();
    a=ui->DH->item(6,2)->text().toDouble();
    alpha=(pi/180)*ui->DH->item(6,3)->text().toDouble();
    DH[5]=joint(d,a,alpha);
    d=ui->DH->item(7,1)->text().toDouble();
    a=ui->DH->item(7,2)->text().toDouble();
    alpha=(pi/180)*ui->DH->item(7,3)->text().toDouble();
    DH[6]=joint(d,a,alpha);
    d=ui->DH->item(8,1)->text().toDouble();
    a=ui->DH->item(8,2)->text().toDouble();
    alpha=(pi/180)*ui->DH->item(8,3)->text().toDouble();
    DH[7]=joint(d,a,alpha);
}
