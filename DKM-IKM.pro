#-------------------------------------------------
#
# Project created by QtCreator 2016-07-02T01:06:29
#
#-------------------------------------------------

QT       += core gui serialport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = DKM-IKM
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    solutions.cpp

HEADERS  += mainwindow.h \
    solutions.h

FORMS    += mainwindow.ui \
    solutions.ui
