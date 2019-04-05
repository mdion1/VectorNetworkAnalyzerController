QT += core serialport
QT -= gui

CONFIG += c++11

TARGET = NetworkVectorAnalyzerController
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    gpib_comms.cpp

HEADERS += \
    gpib_comms.h
