MOC_DIR = ./moc
OBJECTS_DIR = ./obj

QT+= gui declarative

TARGET = image_abstraction

INCLUDEPATH += src/flst  \
    src/synth \
    src/mw3 \
    src/kdtree \
    src/kdtree/ann/include
DEPENDPATH += src/flst \
    src/synth \
    src/mw3 \
    src/kdtree \
    src/kdtree/ann/include

HEADERS       =  src/flst/*.h \
    src/mw3/*.h \
    src/synth/*.h \
    src/kdtree/*.h \
    src/kdtree/ann/include/ANN/*.h \
    TreeOfShapes.h
SOURCES       = main.cpp \
    src/flst/*.c \
    src/mw3/*.c \
    src/synth/*.c \
    src/kdtree/ann/src/*.cpp \
    TreeOfShapes.cpp
LIBS = -L/usr/lib/x86_64-linux-gnu/

QMAKE_CFLAGS += -std=c99
QMAKE_CXXFLAGS_RELEASE = -std=c++0x
QMAKE_CXXFLAGS_RELEASE += -Ofast

QMAKE_CXXFLAGS_DEBUG = -std=c++0x
QMAKE_CXXFLAGS_DEBUG += -Ofast
